// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double my_abs ( double x ) {
  if ( 0.0 <= x ) {
    return x;
  }
  else {
    return ( -x );
  }
}

// Numerically stable log(exp(a) + exp(b)) (audit0.10 §1.11). R_NegInf comes
// from Rmath.h. Both -Inf is the only case that needs guarding directly (it
// would otherwise produce -Inf - -Inf = NaN below).
double log_sum_exp2(double a, double b) {
  double m = std::max(a, b);
  if (m == R_NegInf) {
    return R_NegInf;
  }
  return m + log(exp(a - m) + exp(b - m));
}

// Stable 3-category softmax for one observation's state probabilities, with
// an implicit zero logit for the count (baseline) category (audit0.10 §1.11):
// subtracting the max logit before exponentiating keeps every exp() argument
// <= 0, avoiding overflow for a large linear predictor (coef_limit allows up
// to 50). Mathematically identical to the un-shifted softmax.
void fill_props_row(arma::mat &props, int i, double eta_z, double eta_pl) {
  double m = std::max(0.0, std::max(eta_z, eta_pl));
  double base = exp(-m);
  double d_z = exp(eta_z - m);
  double d_pl = exp(eta_pl - m);
  double denom = base + d_z + d_pl;
  props(i,0) = d_z / denom;
  props(i,1) = base / denom;
  props(i,2) = d_pl / denom;
}

//[[Rcpp::export]]
double ell_nb_i_fun(arma::vec beta_nb, double alpha_nb, arma::vec x_nb_ext_i, int y_i, double offset_nb_i = 0.0){
  // int n_beta_nb = beta_nb.size();

  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;

  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a + offset_nb_i);

  // audit0.10 §1.11: closed form of the O(y) sums below (a rising-factorial /
  // Pochhammer identity), sum_{j=0}^{y-1} log(j + 1/alpha) - sum_{j=1}^{y}
  // log(j) = lgamma(y + 1/alpha) - lgamma(1/alpha) - lgamma(y + 1); valid at
  // y = 0 too (gives 0, matching the previous loops which never ran there).
  // round8 0.9: for a large y with a small 1/alpha the three lgamma() terms
  // above partially cancel (max relative error ~1.5e-10 at y=1e6, alpha=5,
  // against a Kahan-summed reference loop -- the naive-summation loop drifts
  // by ~1e-8 there on its own and is not a reliable reference at that scale).
  // -log(y + 1/alpha) - lbeta(y + 1, 1/alpha) is algebraically identical (via
  // Gamma(y+1+r) = (y+r)*Gamma(y+r)) but R::lbeta() avoids the cancellation,
  // and was at least as accurate as the lgamma form in every case checked.
  double r = 1 / alpha_nb;
  double ell_nb_i = -log(y_i + r) - R::lbeta(y_i + 1, r);

  ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y_i*log(1 + alpha_nb*mu_i) + y_i*log(alpha_nb) + y_i*log(mu_i);

  return(ell_nb_i);
}

// round9 E.1 (audit §5.5): Poisson count-state log-pmf, the alpha_nb -> 0
// limit of ell_nb_i_fun() above -- no dispersion parameter, no digamma.
//[[Rcpp::export]]
double ell_pois_i_fun(arma::vec beta_nb, arma::vec x_nb_ext_i, int y_i, double offset_nb_i = 0.0){
  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;
  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a + offset_nb_i);
  return y_i*log(mu_i) - mu_i - R::lgammafn(y_i + 1);
}

// round9 E.2 (audit §5.5): derivatives of G(mu, alpha) = -log(1 - f0), where
// f0 = P(Y=0) under the untruncated count distribution -- the term a hurdle
// zero process adds to the count state's log-density for every y > 0 row
// (the count state is zero-truncated: its emission density is
// f_count(y)/(1-f0)). f0 = (1+alpha*mu)^(-1/alpha) (NB) or exp(-mu)
// (Poisson, the alpha -> 0 limit); with logf0 = log(f0):
//   d logf0/dmu       = -1/(1+alpha*mu)                                (NB)
//                      = -1                                            (Poisson)
//   d logf0/dalpha    = log(1+alpha*mu)/alpha^2 - mu/(alpha*(1+alpha*mu))  (NB only)
//   d2 logf0/dmu2     = alpha/(1+alpha*mu)^2                           (NB)
//                      = 0                                             (Poisson)
//   d2 logf0/dmu dalpha = mu/(1+alpha*mu)^2                            (NB only)
//   d2 logf0/dalpha2  = mu*(2+3*alpha*mu)/(alpha^2*(1+alpha*mu)^2)
//                        - 2*log(1+alpha*mu)/alpha^3                   (NB only)
// G = h(logf0) with h(L) = -log(1-exp(L)); h'(L) = f0/(1-f0),
// h''(L) = f0/(1-f0)^2. Chain rule: dG/dx = h'*Lx,
// d2G/dxdy = h''*Lx*Ly + h'*Lxy for x, y in {mu, alpha}. All of the above
// (and the formulas below) were verified against central finite differences
// over a grid of mu in {0.5,...,50} x alpha in {0.05,...,3} (max relative
// error < 1e-5) before being committed -- see tests/testthat/test-hurdle.R.
struct hurdle_trunc_derivs {
  double f0, G, dGdmu, dGdalpha, d2Gdmu2, d2Gdmudalpha, d2Gdalpha2;
};
hurdle_trunc_derivs hurdle_trunc_derivs_fun(double mu, double alpha_nb, int family_count) {
  double logf0, q_mu, q_alpha = 0.0, Q_mumu, Q_mualpha = 0.0, Q_alphaalpha = 0.0;
  if (family_count == 1) {
    logf0 = -mu;
    q_mu = -1.0;
    Q_mumu = 0.0;
  } else {
    double oam = 1 + alpha_nb*mu;
    logf0 = -(1/alpha_nb) * log(oam);
    q_mu = -1.0/oam;
    Q_mumu = alpha_nb/(oam*oam);
    q_alpha = log(oam)/(alpha_nb*alpha_nb) - mu/(alpha_nb*oam);
    Q_mualpha = mu/(oam*oam);
    Q_alphaalpha = mu*(2+3*alpha_nb*mu)/(alpha_nb*alpha_nb*oam*oam) - 2*log(oam)/(alpha_nb*alpha_nb*alpha_nb);
  }
  double f0 = exp(logf0);
  // log(1-f0) via -expm1(logf0), the same overflow/cancellation-free idiom
  // used elsewhere in this file (ell_pl_i_fun()).
  double log_1m_f0 = log(-expm1(logf0));
  double G = -log_1m_f0;
  double p = f0/(1-f0);      // h'(L)
  double hpp = p*(1+p);      // h''(L) = f0/(1-f0)^2
  double dGdmu = p * q_mu;
  double dGdalpha = p * q_alpha;
  double d2Gdmu2 = hpp*q_mu*q_mu + p*Q_mumu;
  double d2Gdmudalpha = hpp*q_mu*q_alpha + p*Q_mualpha;
  double d2Gdalpha2 = hpp*q_alpha*q_alpha + p*Q_alphaalpha;
  hurdle_trunc_derivs out = {f0, G, dGdmu, dGdalpha, d2Gdmu2, d2Gdmudalpha, d2Gdalpha2};
  return out;
}

//[[Rcpp::export]]
arma::vec delldtheta_nb_i_fun(arma::vec beta_nb, double alpha_nb, arma::vec x_nb_ext_i, int y_i, double offset_nb_i = 0.0){
  int n_beta_nb = beta_nb.size();

  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;

  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a + offset_nb_i);

  double delldalpha_nb_i = 0;

  if(y_i==0){
    delldalpha_nb_i = log(1 + alpha_nb*mu_i)/(alpha_nb*alpha_nb) - mu_i/(alpha_nb*(1+alpha_nb*mu_i));
  }else{

    // round8 A.2 (audit §5.4): sum_{j=0}^{y-1} 1/(j + 1/alpha) is the
    // digamma difference digamma(y + r) - digamma(r) with r = 1/alpha
    // (verified numerically against the loop over y in {0,1,5,100,1e4,1e5}
    // x alpha in {1e-3,0.01,0.5,1,5,50}, max relative error < 1e-9).
    double r_nb = 1 / alpha_nb;
    double sum_inv_j_plus_r = R::digamma(y_i + r_nb) - R::digamma(r_nb);

    delldalpha_nb_i = log(1 + alpha_nb*mu_i) - sum_inv_j_plus_r;
    delldalpha_nb_i = delldalpha_nb_i/(alpha_nb*alpha_nb);

    delldalpha_nb_i = delldalpha_nb_i + (y_i-mu_i)/(alpha_nb*(1+alpha_nb*mu_i));
  }

  double b = (y_i-mu_i)/(1+alpha_nb*mu_i);
  arma::vec delldbeta_nb_i = x_nb_ext_i*b;
  arma::vec delldtheta_nb_i = zeros<vec>(n_beta_nb+1);
  for(int k=0; k<n_beta_nb; k++){
    delldtheta_nb_i(k) = delldbeta_nb_i(k);
  }
  delldtheta_nb_i(n_beta_nb) = delldalpha_nb_i;

  return(delldtheta_nb_i);
}

//[[Rcpp::export]]
arma::mat d2elldtheta2_nb_i_fun(arma::vec beta_nb, double alpha_nb, arma::vec x_nb_ext_i, int y_i, double offset_nb_i = 0.0){
  int n_beta_nb = beta_nb.size();

  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;

  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a + offset_nb_i);

  mat hessian_nb_i = zeros<mat>(n_beta_nb+1,n_beta_nb+1);
  double d2elldalpha2_nb_i = 0;

  double b = -1.0*mu_i*(1+alpha_nb*y_i)/((1+alpha_nb*mu_i)*(1+alpha_nb*mu_i));

  arma::mat d2elldbeta2_nb_i = x_nb_ext_i*trans(x_nb_ext_i);// as.numeric(-1.0*mu.i*(1+alpha.nb.old*y[i])/(1+alpha.nb.old*mu.i)^2)*

  d2elldbeta2_nb_i = b*d2elldbeta2_nb_i;

  double b2 = -1.0*mu_i*(y_i-mu_i)/((1+alpha_nb*mu_i)*(1+alpha_nb*mu_i));
  arma::vec d2elldalphadbeta_nb_i = b2*x_nb_ext_i;

  if(y_i>0){
    // round8 A.2: this loop is sum_{j=0}^{y-1} (j / (1 + alpha*j))^2, which
    // is NOT sum 1/(j+r)^2 (r = 1/alpha). Writing j/(1+alpha*j) =
    // 1/alpha - 1/(alpha^2*(j+r)) (check: alpha*(1+alpha*j) = alpha^2*(j+r))
    // and squaring gives, after summing over j,
    //   sum_j (j/(1+alpha*j))^2 = y/alpha^2
    //     - (2/alpha^3) * (digamma(y+r) - digamma(r))
    //     + (1/alpha^4) * (trigamma(r) - trigamma(y+r))
    // (the 2/alpha^3 and 1/alpha^4 factors are the chain rule through
    // r = 1/alpha, dr/dalpha = -1/alpha^2, applied twice). Verified
    // numerically against the loop, but this closed form cancels
    // catastrophically for a small alpha*y (not just a small y) -- e.g.
    // y=5, alpha=1e-3 (alpha*y=5e-3) loses ~9 significant digits, and,
    // originally missed, y=31, alpha=1e-6 (alpha*y=3.1e-5) is off by 18%
    // even though y alone is well above the old y > 30 threshold -- because
    // the true value is then a small difference of O(1/alpha^4) terms; the
    // cancellation is governed by the product alpha*y, not y alone (round9
    // 0.4, review §4). Switched the guard accordingly: a fresh verification
    // grid down to alpha = 1e-6 (dev/round9 report) shows the closed form's
    // relative error crossing below 1e-9 once alpha*y exceeds roughly 0.02;
    // alpha*y > 0.1 keeps a wide margin (worst observed relative error
    // ~1e-11 there, several orders of magnitude under the ~1e-9/1e-10 bar
    // used elsewhere in this file) while still using the closed form for
    // every case that matters for cost (a large y with an alpha not
    // vanishingly small). A large y together with an extremely small alpha
    // (alpha*y <= 0.1, e.g. the near-Poisson regime Part E.1 targets) still
    // takes the O(y) loop -- accepted per the round9 brief rather than
    // adding a series expansion, since the loop is only slow, not wrong.
    if (alpha_nb * y_i <= 0.1) {
      for (int j = 0; j < y_i; j++) {
        d2elldalpha2_nb_i = d2elldalpha2_nb_i - (j/(1+alpha_nb*j))*(j/(1+alpha_nb*j));
      }
    } else {
      double r_nb2 = 1 / alpha_nb;
      double loop_term =
        y_i / (alpha_nb*alpha_nb)
        - (2 / (alpha_nb*alpha_nb*alpha_nb)) * (R::digamma(y_i + r_nb2) - R::digamma(r_nb2))
        + (1 / (alpha_nb*alpha_nb*alpha_nb*alpha_nb)) * (R::trigamma(r_nb2) - R::trigamma(y_i + r_nb2));
      d2elldalpha2_nb_i = d2elldalpha2_nb_i - loop_term;
    }
  }

  d2elldalpha2_nb_i = d2elldalpha2_nb_i - 2/(alpha_nb*alpha_nb*alpha_nb)*log(1 + alpha_nb*mu_i) + (2*(1/(alpha_nb*alpha_nb))*mu_i)/(1+alpha_nb*mu_i) + (y_i+1/alpha_nb)*mu_i*mu_i/((1+alpha_nb*mu_i)*(1+alpha_nb*mu_i));

  hessian_nb_i.submat(0,0,n_beta_nb-1,n_beta_nb-1) = d2elldbeta2_nb_i;
  hessian_nb_i.submat(n_beta_nb,0,n_beta_nb,n_beta_nb-1) = trans(d2elldalphadbeta_nb_i);
  hessian_nb_i.submat(0,n_beta_nb,n_beta_nb-1,n_beta_nb) = d2elldalphadbeta_nb_i;
  hessian_nb_i.submat(n_beta_nb,n_beta_nb,n_beta_nb,n_beta_nb) = d2elldalpha2_nb_i;

  return(hessian_nb_i);
}

//[[Rcpp::export]]
double ell_pl_i_fun(arma::vec beta_pl,double c_pl, arma::vec x_pl_ext_i, double y_i){
  arma::mat xtb_pl_i = trans(x_pl_ext_i)*beta_pl;
  double a = xtb_pl_i.eval()(0,0);
  double exp_xtb_pl_i = exp(a);
  //OBS!! The restriction that y.i>c.pl will be taken care of outside
  // audit0.10 §1.11: log((C/y)^a - (C/(y+1))^a) in a form that never computes
  // (C/y)^a directly (which overflows/underflows for large y or extreme a).
  // Algebraically identical: (C/y)^a - (C/(y+1))^a = (C/y)^a * (1 - (y/(y+1))^a),
  // so log(diff) = a*log(C/y) + log(1 - exp(a*log(y/(y+1)))), and
  // 1 - exp(u) = -expm1(u) is accurate for u near 0 (y/(y+1) near 1).
  double ell_pl_i = exp_xtb_pl_i * log(c_pl / y_i) +
    log(-expm1(exp_xtb_pl_i * log(y_i / (y_i + 1))));
  return(ell_pl_i);
}

//[[Rcpp::export]]
arma::vec delldbeta_pl_i_fun_approx(arma::vec beta_pl,double c_pl, arma::vec x_pl_ext_i, double y_i){
  arma::mat xtb_pl_i = trans(x_pl_ext_i)*beta_pl;
  double a = xtb_pl_i.eval()(0,0);
  double exp_xtb_pl_i = exp(a);
  arma::vec delldbeta_pl_i = x_pl_ext_i*(1 + log(c_pl)*exp_xtb_pl_i - log(y_i)*exp_xtb_pl_i);
  return(delldbeta_pl_i);
}

//[[Rcpp::export]]
arma::mat d2elldbeta2_pl_i_fun_approx(arma::vec beta_pl,double c_pl, arma::vec x_pl_ext_i, double y_i){
  arma::mat xtb_pl_i = trans(x_pl_ext_i)*beta_pl;
  double a = xtb_pl_i.eval()(0,0);
  double exp_xtb_pl_i = exp(a);
  arma::mat hessian_pl_i = x_pl_ext_i*trans(x_pl_ext_i)*(log(c_pl)*exp_xtb_pl_i - log(y_i)*exp_xtb_pl_i);
  return(hessian_pl_i);
}

// Shared derivative pieces for the exact discretised-Pareto gradient/Hessian
// (audit0.10 §1.8): with alpha = exp(x'b), u = (C/y)^alpha, v = (C/(y+1))^alpha,
// L1 = log(C/y), L2 = log(C/(y+1)), l = log(u - v):
//   dl/dalpha    = (u*L1 - v*L2) / (u - v)
//   d2l/dalpha2  = (u*L1^2 - v*L2^2) / (u - v) - (dl/dalpha)^2
// round8 0.2 (review §3): forming u and v separately here (even though the
// denominator already used the cancellation-free product below) underflowed
// both to exactly 0 for a sharply peaked block -- large alpha * |L1| -- giving
// 0/0 = NaN despite the log-likelihood itself (ell_pl_i_fun()) staying finite.
// Factor u out of the numerator the same way: with r = v/u = exp(alpha*(L2-L1))
// in (0, 1), (u*L1 - v*L2)/(u - v) = (L1 - r*L2)/(1 - r), and 1 - r is
// -expm1(alpha*(L2-L1)) for the same reason ell_pl_i_fun() uses it.
struct pareto_exact_derivs {
  double alpha_i, dl_dalpha, d2l_dalpha2;
};
pareto_exact_derivs pareto_exact_derivs_fun(arma::vec beta_pl, double c_pl,
                                            arma::vec x_pl_ext_i, double y_i) {
  double lp = (trans(x_pl_ext_i)*beta_pl).eval()(0,0);
  double alpha_i = exp(lp);
  double L1 = log(c_pl / y_i);
  double L2 = log(c_pl / (y_i + 1));
  double r = exp(alpha_i * (L2 - L1));
  double one_minus_r = -expm1(alpha_i * (L2 - L1));
  double dl_dalpha = (L1 - r * L2) / one_minus_r;
  double d2l_dalpha2 = (L1 * L1 - r * L2 * L2) / one_minus_r - dl_dalpha * dl_dalpha;
  pareto_exact_derivs out = {alpha_i, dl_dalpha, d2l_dalpha2};
  return out;
}

//[[Rcpp::export]]
arma::vec delldbeta_pl_i_fun_exact(arma::vec beta_pl,double c_pl, arma::vec x_pl_ext_i, double y_i){
  pareto_exact_derivs d = pareto_exact_derivs_fun(beta_pl, c_pl, x_pl_ext_i, y_i);
  arma::vec delldbeta_pl_i = x_pl_ext_i * (d.alpha_i * d.dl_dalpha);
  return(delldbeta_pl_i);
}

//[[Rcpp::export]]
arma::mat d2elldbeta2_pl_i_fun_exact(arma::vec beta_pl,double c_pl, arma::vec x_pl_ext_i, double y_i){
  pareto_exact_derivs d = pareto_exact_derivs_fun(beta_pl, c_pl, x_pl_ext_i, y_i);
  arma::mat hessian_pl_i = x_pl_ext_i*trans(x_pl_ext_i) *
    (d.alpha_i * d.dl_dalpha + d.alpha_i * d.alpha_i * d.d2l_dalpha2);
  return(hessian_pl_i);
}

//[[Rcpp::export]]
double log_lik_fun(arma::vec gamma_z, arma::vec gamma_pl,arma::vec beta_nb, double alpha_nb, arma::vec beta_pl, double c_pl,arma::mat x_mult_z_ext,arma::mat x_mult_pl_ext,arma::mat x_nb_ext, arma::mat x_pl_ext, arma::vec y, arma::vec offset_nb, arma::vec offset_zc, arma::vec offset_pl_mult, arma::vec w, int family_count = 0, int family_zero = 0, bool has_weights = true){

  int n = x_mult_z_ext.n_rows;
  arma::mat props = zeros<mat>(n,3) ;

  // round8 A.3 (audit §5.4): one matrix-vector product per linear predictor
  // instead of trans(X.submat(i,...)) allocated fresh every row. This
  // function is called from every stats::optimise() line search in
  // R/em_step.R (dozens of evaluations per em_step() call across the four
  // blocks), so it's the single highest-leverage target in this file.
  // round9 D.1: offset_zc / offset_pl_mult (zero-inflation / EVI-inflation
  // component offsets) enter additively into the pre-softmax linear
  // predictor, exactly like offset_nb enters mu_i below -- an offset shifts
  // eta, not d eta / d beta, so every gradient/Hessian formula in this file
  // keeps its existing form.
  arma::vec eta_z_vec = x_mult_z_ext * gamma_z + offset_zc;
  arma::vec eta_pl_mult_vec = x_mult_pl_ext * gamma_pl + offset_pl_mult;
  for(int i=0; i<n; i++){
    fill_props_row(props, i, eta_z_vec(i), eta_pl_mult_vec(i));
  }

  arma::vec eta_nb_vec = x_nb_ext * beta_nb;
  arma::vec eta_pl_vec = x_pl_ext * beta_pl;
  // round9 E.1: r_nb is only meaningful for family_count == 0 (nbinom); the
  // Poisson branch below never uses it.
  double r_nb = (family_count == 0) ? 1 / alpha_nb : 0.0;

  double func_val = 0;

  for(int i=0; i<n; i++){
    double mu_i = exp(eta_nb_vec(i) + offset_nb(i));
    // round9 E.1: Poisson count state is the alpha_nb -> 0 limit of the NB
    // one below -- see ell_pois_i_fun().
    double ell_nb_i;
    if (family_count == 1) {
      ell_nb_i = y(i)*log(mu_i) - mu_i - R::lgammafn(y(i) + 1);
    } else {
      // audit0.10 §1.11 / round8 0.9: closed form, see ell_nb_i_fun().
      ell_nb_i = -log(y(i) + r_nb) - R::lbeta(y(i) + 1, r_nb);
      ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y(i)*log(1 + alpha_nb*mu_i) + y(i)*log(alpha_nb) + y(i)*log(mu_i);
    }

    // audit0.10 §1.11: combine the mixture terms with a log-sum-exp instead of
    // log(p1*exp(.) + p2*exp(.)), which underflows exp(.) to exactly 0 (losing
    // all information) well before the sum itself would over/underflow.
    // round9 D.2: w(i) is a frequency/analytic weight multiplying this row's
    // contribution to the total log-likelihood.
    // round10 0.5 (review §5): has_weights guards every w(i) * / w % site in
    // this file. w(i) == 1 makes the multiplication itself an IEEE754 no-op,
    // but the 0 < y < c branch below used to also regroup a chained
    // `func_val + a + b` into `func_val + w(i)*(a + b)` -- floating-point
    // addition is not associative, so that alone was a real bit-level change
    // that moved default (unweighted) fits off their pre-round9-D.2 optimum
    // by up to 8.6e-8 on a flat likelihood ridge. has_weights = false now
    // takes the exact pre-D.2 expression at every site (see that branch's
    // own comment for why it needs an if/else rather than this file's usual
    // term-then-conditional-multiply pattern). This restores true
    // bit-identity on every case checked directly against pre-D.2 fits; on
    // some other inputs the compiler's own vectorisation of this loop is
    // independently sensitive to the branch's presence and can still differ
    // by up to a couple of ULPs, which no source-level fix here reaches.
    if (family_zero == 1) {
      // round9 E.2: hurdle zero process -- the zero state owns y=0 entirely
      // (no separate count-state density factor there), and the count
      // state is zero-truncated for y>0: -log(1-f0) (see
      // hurdle_trunc_derivs_fun()) is added to its log-density in both the
      // 0<y<c and y>=c branches.
      if (y(i) == 0) {
        func_val = func_val + (has_weights ? w(i) * log(props(i,0)) : log(props(i,0)));
      } else {
        double ell_nb_trunc_i = ell_nb_i + hurdle_trunc_derivs_fun(mu_i, alpha_nb, family_count).G;
        if (y(i) < c_pl) {
          double term = log(props(i,1)) + ell_nb_trunc_i;
          func_val = func_val + (has_weights ? w(i) * term : term);
        } else {
          double exp_xtb_pl_i = exp(eta_pl_vec(i));
          double ell_pl_i = exp_xtb_pl_i * log(c_pl / y(i)) +
            log(-expm1(exp_xtb_pl_i * log(y(i) / (y(i) + 1))));
          double term = log_sum_exp2(log(props(i,1)) + ell_nb_trunc_i, log(props(i,2)) + ell_pl_i);
          func_val = func_val + (has_weights ? w(i) * term : term);
        }
      }
    } else if(y(i)==0){
      double term = log_sum_exp2(log(props(i,0)), log(props(i,1)) + ell_nb_i);
      func_val = func_val + (has_weights ? w(i) * term : term);
    }else if(y(i)>0 && y(i)<c_pl){
      // round10 0.5: NOT the term-then-conditional-multiply pattern used
      // everywhere else in this function -- the pre-D.2 expression here was
      // `func_val + log(props(i,1)) + ell_nb_i` (two chained additions, left
      // to right, no parentheses), not `func_val + (log(props(i,1)) +
      // ell_nb_i)`. Floating-point addition is not associative, so
      // regrouping it that way (as introducing a shared `term` variable
      // would) is itself a real bit-level change, independent of whether w
      // is 1 -- this was the actual source of round9 D.2's identity drift,
      // not the multiplication by w. An explicit if/else keeps the
      // has_weights = false branch textually identical to the pre-D.2 code.
      if (has_weights) {
        func_val = func_val + w(i) * (log(props(i,1)) + ell_nb_i);
      } else {
        func_val = func_val + log(props(i,1)) + ell_nb_i;
      }
    }else{
      double exp_xtb_pl_i = exp(eta_pl_vec(i));
      // audit0.10 §1.11: stable Pareto log-pmf, see ell_pl_i_fun().
      double ell_pl_i = exp_xtb_pl_i * log(c_pl / y(i)) +
        log(-expm1(exp_xtb_pl_i * log(y(i) / (y(i) + 1))));

      double term = log_sum_exp2(log(props(i,1)) + ell_nb_i, log(props(i,2)) + ell_pl_i);
      func_val = func_val + (has_weights ? w(i) * term : term);
    }
  }

  return(func_val);
}

// round8 A.1 (audit §5.4): em_profile_c() used to call log_lik_fun() once per
// C_EV candidate, and every call redid the O(n) work above (state
// probabilities, the whole NB log-likelihood) even though neither depends on
// C. Here that work happens once; only the branch selection and, for
// y >= c, the Pareto term are evaluated per candidate. No monotone-tail
// incremental accumulation (sorting y once and updating as c rises) -- the
// benchmark (A.4) didn't call for it; if it ever does, that's its own commit.
//[[Rcpp::export]]
arma::vec log_lik_profile_fun(arma::vec gamma_z, arma::vec gamma_pl,
                              arma::vec beta_nb, double alpha_nb, arma::vec beta_pl,
                              arma::vec c_candidates,
                              arma::mat x_mult_z_ext, arma::mat x_mult_pl_ext,
                              arma::mat x_nb_ext, arma::mat x_pl_ext,
                              arma::vec y, arma::vec offset_nb,
                              arma::vec offset_zc, arma::vec offset_pl_mult,
                              arma::vec w, int family_count = 0, int family_zero = 0,
                              bool has_weights = true){

  int n = x_mult_z_ext.n_rows;
  int n_mult_z = x_mult_z_ext.n_cols;
  int n_mult_pl = x_mult_pl_ext.n_cols;
  int n_nb = x_nb_ext.n_cols;
  int n_pl = x_pl_ext.n_cols;
  int n_c = c_candidates.n_elem;

  arma::mat props = zeros<mat>(n,3);
  arma::vec ell_nb = zeros<vec>(n);
  arma::vec alpha_pl_vec = zeros<vec>(n);

  double r_nb = (family_count == 0) ? 1 / alpha_nb : 0.0;

  // Steps 1-3: state probabilities, the count log-likelihood (ell_nb_i, see
  // ell_nb_i_fun()/ell_pois_i_fun()) and the Pareto shape (alpha_pl_vec) --
  // none of these depend on C, so each is computed once per observation.
  for (int i = 0; i < n; i++) {
    double eta_z = (trans(gamma_z)*trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))).eval()(0,0) + offset_zc(i);
    double eta_pl = (trans(gamma_pl)*trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))).eval()(0,0) + offset_pl_mult(i);
    fill_props_row(props, i, eta_z, eta_pl);

    arma::mat x_nb_ext_i = trans(x_nb_ext.submat(i,0,i,n_nb-1));
    double xtb_nb_i = (trans(x_nb_ext_i)*beta_nb).eval()(0,0);
    double mu_i = exp(xtb_nb_i + offset_nb(i));
    double ell_nb_i;
    if (family_count == 1) {
      ell_nb_i = y(i)*log(mu_i) - mu_i - R::lgammafn(y(i) + 1);
    } else {
      ell_nb_i = -log(y(i) + r_nb) - R::lbeta(y(i) + 1, r_nb);
      ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y(i)*log(1 + alpha_nb*mu_i) + y(i)*log(alpha_nb) + y(i)*log(mu_i);
    }
    // round9 E.2: precompute the zero-truncation correction once per
    // observation too (it doesn't depend on C), same as ell_nb_i itself.
    if (family_zero == 1 && y(i) > 0) {
      ell_nb_i = ell_nb_i + hurdle_trunc_derivs_fun(mu_i, alpha_nb, family_count).G;
    }
    ell_nb(i) = ell_nb_i;

    arma::mat x_pl_ext_i = trans(x_pl_ext.submat(i,0,i,n_pl-1));
    double xtb_pl_i = (trans(x_pl_ext_i)*beta_pl).eval()(0,0);
    alpha_pl_vec(i) = exp(xtb_pl_i);
  }

  arma::vec loglik = zeros<vec>(n_c);

  // Step 4-5: per candidate, branch on y_i == 0 / y_i < c / y_i >= c exactly
  // as log_lik_fun() does, reusing props/ell_nb/alpha_pl_vec from above; only
  // the Pareto term (y_i >= c) is evaluated fresh per candidate, combined via
  // the same log_sum_exp2() as log_lik_fun().
  for (int g = 0; g < n_c; g++) {
    double c_pl = c_candidates(g);
    double func_val = 0;
    for (int i = 0; i < n; i++) {
      // round9 E.2: under a hurdle zero process, y=0 rows are entirely the
      // zero state's (no C-dependence there), and ell_nb(i) already carries
      // the zero-truncation correction for y>0 rows (precomputed above).
      // round10 0.5: has_weights guards the multiply, see log_lik_fun() above.
      if (family_zero == 1) {
        if (y(i) == 0) {
          func_val += has_weights ? w(i) * log(props(i,0)) : log(props(i,0));
        } else if (y(i) < c_pl) {
          double term = log(props(i,1)) + ell_nb(i);
          func_val += has_weights ? w(i) * term : term;
        } else {
          double a = alpha_pl_vec(i);
          double ell_pl_i = a * log(c_pl / y(i)) +
            log(-expm1(a * log(y(i) / (y(i) + 1))));
          double term = log_sum_exp2(log(props(i,1)) + ell_nb(i), log(props(i,2)) + ell_pl_i);
          func_val += has_weights ? w(i) * term : term;
        }
      } else if (y(i) == 0) {
        double term = log_sum_exp2(log(props(i,0)), log(props(i,1)) + ell_nb(i));
        func_val += has_weights ? w(i) * term : term;
      } else if (y(i) > 0 && y(i) < c_pl) {
        double term = log(props(i,1)) + ell_nb(i);
        func_val += has_weights ? w(i) * term : term;
      } else {
        double a = alpha_pl_vec(i);
        double ell_pl_i = a * log(c_pl / y(i)) +
          log(-expm1(a * log(y(i) / (y(i) + 1))));
        double term = log_sum_exp2(log(props(i,1)) + ell_nb(i), log(props(i,2)) + ell_pl_i);
        func_val += has_weights ? w(i) * term : term;
      }
    }
    loglik(g) = func_val;
  }

  return loglik;
}

// Bounded Newton step -H^{-1}*g, robust to a singular Hessian (audit0.10
// §1.3). H is the Hessian of the (expected complete-data log-likelihood)
// objective being MAXIMISED, so it is negative (semi-)definite near the
// optimum; -inv(H)*g is therefore an ascent step. arma::solve(...,
// solve_opts::no_approx) returns false on failure instead of throwing, so a
// singular H never raises a C++ exception here. On failure, retry once with
// a small ridge that pushes H further into negative-definite territory
// (subtracting a positive multiple of the identity, consistent with H's
// existing sign convention) rather than flipping it. If that also fails,
// return a vector of NA_REAL of the same length as g; every caller below
// then propagates that into the corresponding "_old" parameter via ordinary
// (NaN-propagating) vector addition, which is exactly what em_step() (R/
// em_step.R) already checks for to keep that block's old value.
arma::mat safe_newton_step(const arma::mat &H, const arma::mat &g) {
  arma::mat step;
  bool ok = arma::solve(step, H, -g, arma::solve_opts::no_approx);
  if (!ok) {
    double scale = std::max(1.0, arma::abs(H).max());
    arma::mat H_ridge = H - 1e-8 * scale * arma::eye(H.n_rows, H.n_cols);
    ok = arma::solve(step, H_ridge, -g, arma::solve_opts::no_approx);
  }
  if (!ok) {
    step = arma::mat(g.n_rows, g.n_cols);
    step.fill(NA_REAL);
  }
  return step;
}

//[[Rcpp::export]]
List update_bfgs_fun(arma::vec gamma_z_in, arma::vec gamma_pl_in,arma::vec beta_nb_in, double alpha_nb_in, arma::vec beta_pl_in, double c_pl,arma::mat x_mult_z_ext,arma::mat x_mult_pl_ext,arma::mat x_nb_ext, arma::mat x_pl_ext, arma::vec y, double max_upd_par, int no_m_bfgs_steps, arma::vec offset_nb, arma::vec offset_zc, arma::vec offset_pl_mult, arma::vec w, int family_count = 0, bool exact_pl = false, int family_zero = 0, bool has_weights = true){

  int n = x_mult_z_ext.n_rows;
  int n_mult_z = x_mult_z_ext.n_cols;
  int n_mult_pl = x_mult_pl_ext.n_cols;
  int n_nb = x_nb_ext.n_cols;
  int n_pl = x_pl_ext.n_cols;
  arma::mat props = zeros<mat>(n,3) ;

  arma::vec gamma_z_old = gamma_z_in;
  arma::vec gamma_z_after_bfgs = gamma_z_in;
  arma::vec gamma_pl_old = gamma_pl_in;
  arma::vec beta_nb_old = beta_nb_in;
  arma::vec beta_nb_after_bfgs = beta_nb_in;
  double alpha_nb_old = alpha_nb_in;
// double alpha_nb_after_bfgs = alpha_nb_in;
  arma::vec beta_pl_old = beta_pl_in;
  arma::vec beta_pl_after_bfgs = beta_pl_in;


  // round8 A.3: precomputed linear predictor instead of trans(X.submat(i,...))
  // per row (same pattern as log_lik_fun()). round9 D.1: offsets enter
  // additively, same as log_lik_fun()/log_lik_profile_fun() above.
  arma::vec eta_z_vec0 = x_mult_z_ext * gamma_z_old + offset_zc;
  arma::vec eta_pl_mult_vec0 = x_mult_pl_ext * gamma_pl_old + offset_pl_mult;
  for(int i=0; i<n; i++){
    fill_props_row(props, i, eta_z_vec0(i), eta_pl_mult_vec0(i));
  }

  arma::mat resp = zeros<mat>(n,3);
  double marginal_yx_i = 0;

  for (int i=0; i<n; i++){
  // Marginal probability mass function of y and x (sum over components) - marginal._x_i
    if(y(i)==0){
      // round9 E.2: under a hurdle zero process the count state cannot
      // produce y=0 at all (it's zero-truncated), so the posterior is
      // degenerate -- forced to the zero state, not derived from a
      // props/density ratio.
      if (family_zero == 1) {
        resp(i,0) = 1;
        resp(i,1) = 0;
        resp(i,2) = 0;
      } else {
        double d_nb = (family_count == 1)
          ? exp(ell_pois_i_fun(beta_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i)))
          : exp(ell_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i)));
        marginal_yx_i = props(i,0) + props(i,1)*d_nb;
        resp(i,0) = props(i,0)/marginal_yx_i;
        resp(i,1) = props(i,1)*d_nb/marginal_yx_i;
        resp(i,2) = 0;
      }
    }else if(y(i)>0 && y(i)<c_pl){
      resp(i,0) = 0;
      resp(i,1) = 1;
      resp(i,2) = 0;
    }else if(y(i)>=c_pl){
      double d_nb = (family_count == 1)
        ? exp(ell_pois_i_fun(beta_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i)))
        : exp(ell_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i)));
      if (family_zero == 1) {
        // round9 E.2: the count state's emission density is zero-truncated
        // (f_count(y)/(1-f0)); exp(G) = 1/(1-f0).
        double mu_i = exp((trans(beta_nb_old)*trans(x_nb_ext.submat(i,0,i,n_nb-1))).eval()(0,0) + offset_nb(i));
        d_nb = d_nb * exp(hurdle_trunc_derivs_fun(mu_i, alpha_nb_old, family_count).G);
      }
      double d_pl = exp(ell_pl_i_fun(beta_pl_old,c_pl,trans(x_pl_ext.submat(i,0,i,n_pl-1)),y(i)));
      marginal_yx_i = props(i,1)*d_nb  + props(i,2)*d_pl;
      resp(i,0) = 0;
      resp(i,1) = props(i,1)*d_nb/marginal_yx_i;
      resp(i,2) = props(i,2)*d_pl/marginal_yx_i;
    }
  }

  //Calculate function value before the algorithm starts
  double func_val_before_bfgs = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb,offset_zc,offset_pl_mult,w,family_count,family_zero,has_weights);

  arma::mat d2Qdtheta2_nb = zeros<mat>(n_nb+1,n_nb+1);
  arma::mat dQdtheta_nb = zeros<mat>(n_nb+1,1);
  double maxabschange = 0;
  arma::mat change_nb_bfgs = zeros<mat>(n_nb+1,1);
  arma::mat d2Qdbeta2_pl = zeros<mat>(n_pl,n_pl);
  arma::mat dQdbeta_pl = zeros<mat>(n_pl,1);
  arma::mat change_pl_bfgs = zeros<mat>(n_pl,1);
  arma::mat dQdgamma_z = zeros<mat>(n_mult_z,1);
  arma::mat d2Qdgamma2_z = zeros<mat>(n_mult_z,n_mult_z);
  arma::mat change_mult_z_bfgs = zeros<mat>(n_mult_z,1);
  arma::mat dQdgamma_pl = zeros<mat>(n_mult_pl,1);
  arma::mat d2Qdgamma2_pl = zeros<mat>(n_mult_pl,n_mult_pl);
  arma::mat change_mult_pl_bfgs = zeros<mat>(n_mult_pl,1);



  //Update beta_nb and alpha_nb with BFGS
  // round8 A.3: the per-row loop rebuilt mu_i (via a fresh submat/trans) and
  // allocated a new gradient/Hessian block on every row, every sub-iteration,
  // by calling delldtheta_nb_i_fun()/d2elldtheta2_nb_i_fun(). Precompute the
  // linear predictor and the per-row weights in one pass, then accumulate the
  // beta block as X'w (gradient) and X' diag(w) X (Hessian) -- the digamma/
  // trigamma alpha terms (A.2, including its y <= 30 loop fallback) are
  // scalar and still computed per row, but with no matrix allocation.
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    arma::vec eta_nb_i = x_nb_ext * beta_nb_old;
    arma::vec mu_vec = exp(eta_nb_i + offset_nb);

    arma::vec grad_beta, H_cross;
    arma::mat H_beta;
    double grad_alpha = 0.0, H_alpha = -1.0;

    if (family_count == 1) {
      // round9 E.1: Poisson count state -- no dispersion parameter. Score
      // X'(w*resp*(y-mu)), Hessian -X' diag(w*resp*mu) X: simpler and
      // better conditioned than the NB block below (no digamma, no
      // dispersion row/col). grad_alpha/H_cross stay 0 and H_alpha stays a
      // decoupled -1, so the Newton system is block-diagonal and the
      // alpha coordinate's step is exactly 0 every time -- alpha_nb_old
      // never moves, and R drops it from par.all entirely for this family.
      arma::vec b_grad = y - mu_vec;
      arma::vec b_hess = -mu_vec;
      if (family_zero == 1) {
        // round9 E.2: hurdle correction -- the count state's emission
        // density is zero-truncated for y>0; resp.col(1) is already
        // exactly 0 for y=0 rows (forced in the E-step above), so this
        // addition only affects rows that actually contribute to this
        // block.
        for (int i = 0; i < n; i++) {
          hurdle_trunc_derivs hd = hurdle_trunc_derivs_fun(mu_vec(i), alpha_nb_old, family_count);
          b_grad(i) += hd.dGdmu * mu_vec(i);
          b_hess(i) += hd.d2Gdmu2 * mu_vec(i) * mu_vec(i) + hd.dGdmu * mu_vec(i);
        }
      }
      // round10 0.5: has_weights guards the multiply, see log_lik_fun() above.
      // Armadillo's %/eGlue expression templates give the two ternary
      // branches different C++ types, so this is an if/else into a
      // concrete arma::vec rather than a ?: expression.
      arma::vec w_beta, w_hess_beta;
      if (has_weights) {
        w_beta = w % resp.col(1) % b_grad;
        w_hess_beta = w % resp.col(1) % b_hess;
      } else {
        w_beta = resp.col(1) % b_grad;
        w_hess_beta = resp.col(1) % b_hess;
      }
      grad_beta = trans(x_nb_ext) * w_beta;
      H_beta = trans(x_nb_ext) * (x_nb_ext.each_col() % w_hess_beta);
      H_cross = zeros<vec>(n_nb);
    } else {
      arma::vec one_plus_am = 1 + alpha_nb_old*mu_vec;

      arma::vec b_grad = (y - mu_vec) / one_plus_am;
      arma::vec b_hess = -1.0*mu_vec % (1 + alpha_nb_old*y) / arma::square(one_plus_am);
      arma::vec b2_hess = -1.0*mu_vec % (y - mu_vec) / arma::square(one_plus_am);

      double r_nb = 1 / alpha_nb_old;
      arma::vec delldalpha_vec(n), d2elldalpha2_vec(n);
      for (int i=0; i<n; i++){
        double mu_i = mu_vec(i);
        double oam_i = one_plus_am(i);
        int yi = (int) y(i);

        // audit0.10 §1.11 / round8 A.2: digamma(y+r)-digamma(r) is exactly 0 at
        // y=0 (same argument on both sides), so this one formula covers both
        // branches of the original delldtheta_nb_i_fun().
        double sum_inv = R::digamma(y(i) + r_nb) - R::digamma(r_nb);
        delldalpha_vec(i) = (log(oam_i) - sum_inv)/(alpha_nb_old*alpha_nb_old) + (y(i)-mu_i)/(alpha_nb_old*oam_i);

        // round9 0.4 (review §4): the cancellation this closed form suffers
        // from is governed by alpha*y, not y alone -- see the comment at the
        // matching guard in d2elldtheta2_nb_i_fun() above for the
        // verification grid and the alpha_nb_old * y(i) <= 0.1 cutoff.
        double loop_term = 0;
        if (alpha_nb_old * y(i) > 0.1) {
          loop_term = y(i)/(alpha_nb_old*alpha_nb_old)
            - (2/(alpha_nb_old*alpha_nb_old*alpha_nb_old))*sum_inv
            + (1/(alpha_nb_old*alpha_nb_old*alpha_nb_old*alpha_nb_old))*(R::trigamma(r_nb) - R::trigamma(y(i)+r_nb));
        } else if (yi > 0) {
          for (int j=0; j<yi; j++) {
            loop_term += (j/(1+alpha_nb_old*j))*(j/(1+alpha_nb_old*j));
          }
        }
        d2elldalpha2_vec(i) = -loop_term - 2/(alpha_nb_old*alpha_nb_old*alpha_nb_old)*log(oam_i)
          + (2/(alpha_nb_old*alpha_nb_old))*mu_i/oam_i + (y(i)+r_nb)*mu_i*mu_i/(oam_i*oam_i);

        if (family_zero == 1) {
          // round9 E.2: hurdle correction, same rationale as the Poisson
          // block above -- resp.col(1) is 0 for y=0 rows, so this addition
          // is inert there regardless.
          hurdle_trunc_derivs hd = hurdle_trunc_derivs_fun(mu_i, alpha_nb_old, family_count);
          b_grad(i) += hd.dGdmu * mu_i;
          b_hess(i) += hd.d2Gdmu2 * mu_i * mu_i + hd.dGdmu * mu_i;
          b2_hess(i) += hd.d2Gdmudalpha * mu_i;
          delldalpha_vec(i) += hd.dGdalpha;
          d2elldalpha2_vec(i) += hd.d2Gdalpha2;
        }
      }

      // round9 D.2: w (the frequency/analytic weight, distinct from the
      // w_* local names below which predate it) multiplies resp the same
      // way throughout -- score = X'(w*resp*u), Hessian = -X' diag(w*resp*v) X.
      // round10 0.5: has_weights guards the multiply, see log_lik_fun() above.
      // (if/else into concrete arma::vec, not ?: -- see the Poisson block's
      // comment above on why.)
      arma::vec w_beta, w_hess_beta, w_hess_cross;
      if (has_weights) {
        w_beta = w % resp.col(1) % b_grad;
        grad_alpha = arma::sum(w % resp.col(1) % delldalpha_vec);
        w_hess_beta = w % resp.col(1) % b_hess;
        w_hess_cross = w % resp.col(1) % b2_hess;
        H_alpha = arma::sum(w % resp.col(1) % d2elldalpha2_vec);
      } else {
        w_beta = resp.col(1) % b_grad;
        grad_alpha = arma::sum(resp.col(1) % delldalpha_vec);
        w_hess_beta = resp.col(1) % b_hess;
        w_hess_cross = resp.col(1) % b2_hess;
        H_alpha = arma::sum(resp.col(1) % d2elldalpha2_vec);
      }
      grad_beta = trans(x_nb_ext) * w_beta;
      H_beta = trans(x_nb_ext) * (x_nb_ext.each_col() % w_hess_beta);
      H_cross = trans(x_nb_ext) * w_hess_cross;
    }

    dQdtheta_nb = zeros<mat>(n_nb+1,1);
    dQdtheta_nb.submat(0,0,n_nb-1,0) = grad_beta;
    dQdtheta_nb(n_nb,0) = grad_alpha;

    d2Qdtheta2_nb = zeros<mat>(n_nb+1,n_nb+1);
    d2Qdtheta2_nb.submat(0,0,n_nb-1,n_nb-1) = H_beta;
    d2Qdtheta2_nb.submat(n_nb,0,n_nb,n_nb-1) = trans(H_cross);
    d2Qdtheta2_nb.submat(0,n_nb,n_nb-1,n_nb) = H_cross;
    d2Qdtheta2_nb(n_nb,n_nb) = H_alpha;

    change_nb_bfgs = safe_newton_step(d2Qdtheta2_nb, dQdtheta_nb);
    if (change_nb_bfgs.has_nan()) {
      // A permanently singular Hessian: propagate NA into beta_nb_old /
      // alpha_nb_old (so em_step() falls back to the pre-step values) and
      // stop iterating this block -- further steps from NaN parameters are
      // wasted work.
      beta_nb_old = beta_nb_old + change_nb_bfgs.submat(0,0,n_nb-1,0);
      alpha_nb_old = alpha_nb_old + change_nb_bfgs.submat(n_nb,0,n_nb,0).eval()(0,0);
      break;
    }

    maxabschange = max(abs(change_nb_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_nb_bfgs = max_upd_par/maxabschange*change_nb_bfgs;
    }

    beta_nb_old = beta_nb_old + change_nb_bfgs.submat(0,0,n_nb-1,0);
    alpha_nb_old = alpha_nb_old + change_nb_bfgs.submat(n_nb,0,n_nb,0).eval()(0,0);

  }

  double func_val_after_nb = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_old,alpha_nb_old,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb,offset_zc,offset_pl_mult,w,family_count,family_zero,has_weights);

  //Update beta_pl with BFGS
  // round8 A.3: restrict to the y >= c_pl rows once (arma::find()), then
  // accumulate as X'w / X' diag(w) X on that subset instead of a per-row
  // loop that allocated a fresh gradient/Hessian block (and re-evaluated
  // x_i'beta_pl) every row, every sub-iteration. The nonlinear per-row terms
  // (dl_dalpha etc.) are still scalar, but only over the active subset,
  // which review's own framing puts at "typically a few percent of n".
  arma::uvec pl_idx = arma::find(y >= c_pl);
  int n_pl_active = pl_idx.n_elem;
  arma::mat x_pl_sub;
  arma::vec y_pl_sub, resp2_sub, w_pl_sub;
  if (n_pl_active > 0) {
    x_pl_sub = x_pl_ext.rows(pl_idx);
    y_pl_sub = y.elem(pl_idx);
    arma::vec resp2_full = resp.col(2);
    resp2_sub = resp2_full.elem(pl_idx);
    w_pl_sub = w.elem(pl_idx);
  }

  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    d2Qdbeta2_pl = zeros<mat>(n_pl,n_pl);
    dQdbeta_pl = zeros<mat>(n_pl,1);

    if (n_pl_active > 0) {
      arma::vec eta_pl_sub = x_pl_sub * beta_pl_old;
      arma::vec alpha_sub = exp(eta_pl_sub);
      arma::vec grad_w(n_pl_active), hess_w(n_pl_active);

      for (int k = 0; k < n_pl_active; k++) {
        double y_k = y_pl_sub(k);
        double alpha_k = alpha_sub(k);
        if (exact_pl) {
          double L1 = log(c_pl / y_k);
          double L2 = log(c_pl / (y_k + 1));
          double r = exp(alpha_k * (L2 - L1));
          double one_minus_r = -expm1(alpha_k * (L2 - L1));
          double dl_dalpha = (L1 - r*L2) / one_minus_r;
          double d2l_dalpha2 = (L1*L1 - r*L2*L2)/one_minus_r - dl_dalpha*dl_dalpha;
          grad_w(k) = alpha_k * dl_dalpha;
          hess_w(k) = alpha_k*dl_dalpha + alpha_k*alpha_k*d2l_dalpha2;
        } else {
          double term = log(c_pl)*alpha_k - log(y_k)*alpha_k;
          grad_w(k) = 1 + term;
          hess_w(k) = term;
        }
      }

      // round10 0.5: has_weights guards the multiply, see log_lik_fun() above.
      arma::vec w_grad, w_hess;
      if (has_weights) {
        w_grad = w_pl_sub % resp2_sub % grad_w;
        w_hess = w_pl_sub % resp2_sub % hess_w;
      } else {
        w_grad = resp2_sub % grad_w;
        w_hess = resp2_sub % hess_w;
      }
      dQdbeta_pl = x_pl_sub.t() * w_grad;
      d2Qdbeta2_pl = x_pl_sub.t() * (x_pl_sub.each_col() % w_hess);
    }

    change_pl_bfgs = safe_newton_step(d2Qdbeta2_pl, dQdbeta_pl);
    if (change_pl_bfgs.has_nan()) {
      beta_pl_old = beta_pl_old + change_pl_bfgs;
      break;
    }

    maxabschange = max(abs(change_pl_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_pl_bfgs = max_upd_par/maxabschange*change_pl_bfgs;
    }
    beta_pl_old = beta_pl_old + change_pl_bfgs;

  }

  double func_val_after_pl = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_old,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb,offset_zc,offset_pl_mult,w,family_count,family_zero,has_weights);

  //Update gamma_z with BFGS
  // round8 A.3: X'w / X' diag(w) X instead of a per-row loop building a
  // fresh outer product x_i x_i' (and re-evaluating both linear predictors)
  // every row, every sub-iteration. gamma_pl_old is fixed throughout this
  // block (the gamma_pl block runs after), so its linear predictor is
  // computed once, outside the loop.
  arma::vec eta_pl_mult_fixed = x_mult_pl_ext * gamma_pl_old + offset_pl_mult;
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    arma::vec eta_z_i = x_mult_z_ext * gamma_z_old + offset_zc;
    arma::vec denom_vec = 1 + exp(eta_z_i) + exp(eta_pl_mult_fixed);
    // round10 0.5: has_weights guards the multiply, see log_lik_fun() above.
    arma::vec z_grad_base = resp.col(0) - exp(eta_z_i)/denom_vec;
    arma::vec w_grad_z = has_weights ? (w % z_grad_base) : z_grad_base;
    dQdgamma_z = x_mult_z_ext.t() * w_grad_z;
    arma::vec z_hess_base = 1.0/denom_vec;
    arma::vec w_hess_z = has_weights ? (w % z_hess_base) : z_hess_base;
    d2Qdgamma2_z = -1.0 * (x_mult_z_ext.t() * (x_mult_z_ext.each_col() % w_hess_z));

    change_mult_z_bfgs = safe_newton_step(d2Qdgamma2_z, dQdgamma_z);
    if (change_mult_z_bfgs.has_nan()) {
      gamma_z_old = gamma_z_old + change_mult_z_bfgs;
      break;
    }

    maxabschange = max(abs(change_mult_z_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_mult_z_bfgs = max_upd_par/maxabschange*change_mult_z_bfgs;
    }
    gamma_z_old = gamma_z_old + change_mult_z_bfgs;

  }

  double func_val_after_mult_z = log_lik_fun(gamma_z_old,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb,offset_zc,offset_pl_mult,w,family_count,family_zero,has_weights);

  //Update gamma_pl with BFGS
  // round8 A.3: same pattern as the gamma_z block; gamma_z_old is now fixed
  // (already updated above), so its linear predictor is computed once.
  arma::vec eta_z_fixed = x_mult_z_ext * gamma_z_old + offset_zc;
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    arma::vec eta_pl_mult_i = x_mult_pl_ext * gamma_pl_old + offset_pl_mult;
    arma::vec denom_vec2 = 1 + exp(eta_z_fixed) + exp(eta_pl_mult_i);
    // round10 0.5: has_weights guards the multiply, see log_lik_fun() above.
    arma::vec pl_grad_base = resp.col(2) - exp(eta_pl_mult_i)/denom_vec2;
    arma::vec w_grad_pl = has_weights ? (w % pl_grad_base) : pl_grad_base;
    dQdgamma_pl = x_mult_pl_ext.t() * w_grad_pl;
    arma::vec pl_hess_base = 1.0/denom_vec2;
    arma::vec w_hess_pl = has_weights ? (w % pl_hess_base) : pl_hess_base;
    d2Qdgamma2_pl = -1.0 * (x_mult_pl_ext.t() * (x_mult_pl_ext.each_col() % w_hess_pl));

    change_mult_pl_bfgs = safe_newton_step(d2Qdgamma2_pl, dQdgamma_pl);
    if (change_mult_pl_bfgs.has_nan()) {
      gamma_pl_old = gamma_pl_old + change_mult_pl_bfgs;
      break;
    }

    maxabschange = max(abs(change_mult_pl_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_mult_pl_bfgs = max_upd_par/maxabschange*change_mult_pl_bfgs;
    }
    gamma_pl_old = gamma_pl_old + change_mult_pl_bfgs;

  }

  double func_val_after_mult_pl = log_lik_fun(gamma_z_in,gamma_pl_old,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb,offset_zc,offset_pl_mult,w,family_count,family_zero,has_weights);

  return Rcpp::List::create(
    Rcpp::Named("props") = props,
    Rcpp::Named("resp") = resp,
    Rcpp::Named("beta_nb_old") = beta_nb_old,
    Rcpp::Named("alpha_nb_old") = alpha_nb_old,
    Rcpp::Named("beta_pl_old") = beta_pl_old,
    Rcpp::Named("gamma_z_old") = gamma_z_old,
    Rcpp::Named("gamma_pl_old") = gamma_pl_old,
    Rcpp::Named("change_nb_bfgs") = change_nb_bfgs,
    Rcpp::Named("change_pl_bfgs") = change_pl_bfgs,
    Rcpp::Named("change_mult_z_bfgs") = change_mult_z_bfgs,
    Rcpp::Named("change_mult_pl_bfgs") = change_mult_pl_bfgs,
    Rcpp::Named("func_val_before_bfgs") = func_val_before_bfgs,
    Rcpp::Named("func_val_after_nb") = func_val_after_nb,
    Rcpp::Named("func_val_after_pl") = func_val_after_pl,
    Rcpp::Named("func_val_after_mult_z") = func_val_after_mult_z,
    Rcpp::Named("func_val_after_mult_pl") = func_val_after_mult_pl
  ) ;

}


