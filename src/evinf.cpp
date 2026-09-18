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
  double ell_nb_i = lgamma(y_i + 1/alpha_nb) - lgamma(1/alpha_nb) - lgamma(y_i + 1);

  ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y_i*log(1 + alpha_nb*mu_i) + y_i*log(alpha_nb) + y_i*log(mu_i);

  return(ell_nb_i);
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

    delldalpha_nb_i = log(1 + alpha_nb*mu_i);
    for(int j=0; j<y_i; j++){
      delldalpha_nb_i = delldalpha_nb_i - 1/(j+1/alpha_nb);
    }
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
    for(int j=0; j<y_i; j++){
      d2elldalpha2_nb_i = d2elldalpha2_nb_i - (j/(1+alpha_nb*j))*(j/(1+alpha_nb*j));
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
double log_lik_fun(arma::vec gamma_z, arma::vec gamma_pl,arma::vec beta_nb, double alpha_nb, arma::vec beta_pl, double c_pl,arma::mat x_mult_z_ext,arma::mat x_mult_pl_ext,arma::mat x_nb_ext, arma::mat x_pl_ext, arma::vec y, arma::vec offset_nb){

  int n = x_mult_z_ext.n_rows;
  int n_mult_z = x_mult_z_ext.n_cols;
  int n_mult_pl = x_mult_pl_ext.n_cols;
  int n_nb = x_nb_ext.n_cols;
  int n_pl = x_pl_ext.n_cols;
  arma::mat props = zeros<mat>(n,3) ;

  for(int i=0; i<n; i++){
    double eta_z = (trans(gamma_z)*trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))).eval()(0,0);
    double eta_pl = (trans(gamma_pl)*trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))).eval()(0,0);
    fill_props_row(props, i, eta_z, eta_pl);
  }

  double func_val = 0;

  for(int i=0; i<n; i++){
    arma::mat x_nb_ext_i = trans(x_nb_ext.submat(i,0,i,n_nb-1));
    double xtb_nb_i = (trans(x_nb_ext_i)*beta_nb).eval()(0,0);
    double mu_i = exp(xtb_nb_i + offset_nb(i));
    // audit0.10 §1.11: closed form, see ell_nb_i_fun().
    double ell_nb_i = lgamma(y(i) + 1/alpha_nb) - lgamma(1/alpha_nb) - lgamma(y(i) + 1);
    ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y(i)*log(1 + alpha_nb*mu_i) + y(i)*log(alpha_nb) + y(i)*log(mu_i);

    // audit0.10 §1.11: combine the mixture terms with a log-sum-exp instead of
    // log(p1*exp(.) + p2*exp(.)), which underflows exp(.) to exactly 0 (losing
    // all information) well before the sum itself would over/underflow.
    if(y(i)==0){
      func_val = func_val + log_sum_exp2(log(props(i,0)), log(props(i,1)) + ell_nb_i);
    }else if(y(i)>0 && y(i)<c_pl){
      func_val = func_val + log(props(i,1)) + ell_nb_i;
    }else{
      arma::mat x_pl_ext_i = trans(x_pl_ext.submat(i,0,i,n_pl-1));
      double xtb_pl_i = (trans(x_pl_ext_i)*beta_pl).eval()(0,0);
      double exp_xtb_pl_i = exp(xtb_pl_i);
      // audit0.10 §1.11: stable Pareto log-pmf, see ell_pl_i_fun().
      double ell_pl_i = exp_xtb_pl_i * log(c_pl / y(i)) +
        log(-expm1(exp_xtb_pl_i * log(y(i) / (y(i) + 1))));

      func_val = func_val + log_sum_exp2(log(props(i,1)) + ell_nb_i, log(props(i,2)) + ell_pl_i);
    }
  }

  return(func_val);
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
List update_bfgs_fun(arma::vec gamma_z_in, arma::vec gamma_pl_in,arma::vec beta_nb_in, double alpha_nb_in, arma::vec beta_pl_in, double c_pl,arma::mat x_mult_z_ext,arma::mat x_mult_pl_ext,arma::mat x_nb_ext, arma::mat x_pl_ext, arma::vec y, double max_upd_par, int no_m_bfgs_steps, arma::vec offset_nb, bool exact_pl = false){

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


  double denominator = 0;
  for(int i=0; i<n; i++){
    double eta_z = (trans(gamma_z_old)*trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))).eval()(0,0);
    double eta_pl = (trans(gamma_pl_old)*trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))).eval()(0,0);
    fill_props_row(props, i, eta_z, eta_pl);
  }

  arma::mat resp = zeros<mat>(n,3);
  double marginal_yx_i = 0;

  for (int i=0; i<n; i++){
  // Marginal probability mass function of y and x (sum over components) - marginal._x_i
    if(y(i)==0){
      double d_nb = exp(ell_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i)));
      marginal_yx_i = props(i,0) + props(i,1)*d_nb;
      resp(i,0) = props(i,0)/marginal_yx_i;
      resp(i,1) = props(i,1)*d_nb/marginal_yx_i;
      resp(i,2) = 0;
    }else if(y(i)>0 && y(i)<c_pl){
      resp(i,0) = 0;
      resp(i,1) = 1;
      resp(i,2) = 0;
    }else if(y(i)>=c_pl){
      double d_nb = exp(ell_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i)));
      double d_pl = exp(ell_pl_i_fun(beta_pl_old,c_pl,trans(x_pl_ext.submat(i,0,i,n_pl-1)),y(i)));
      marginal_yx_i = props(i,1)*d_nb  + props(i,2)*d_pl;
      resp(i,0) = 0;
      resp(i,1) = props(i,1)*d_nb/marginal_yx_i;
      resp(i,2) = props(i,2)*d_pl/marginal_yx_i;
    }
  }

  //Calculate function value before the algorithm starts
  double func_val_before_bfgs = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb);

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
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    d2Qdtheta2_nb = zeros<mat>(n_nb+1,n_nb+1);
    dQdtheta_nb = zeros<mat>(n_nb+1,1);

    for(int i=0; i<n; i++){

      dQdtheta_nb = dQdtheta_nb + delldtheta_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i))*resp(i,1);
      d2Qdtheta2_nb = d2Qdtheta2_nb + d2elldtheta2_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i),offset_nb(i))*resp(i,1);
    }

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

  double func_val_after_nb = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_old,alpha_nb_old,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb);

  //Update beta_pl with BFGS
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    d2Qdbeta2_pl = zeros<mat>(n_pl,n_pl);
    dQdbeta_pl = zeros<mat>(n_pl,1);

    for(int i=0; i<n; i++){
      if(y(i)>=c_pl){
        arma::vec x_pl_ext_i = trans(x_pl_ext.submat(i,0,i,n_pl-1));
        if (exact_pl) {
          dQdbeta_pl = dQdbeta_pl + delldbeta_pl_i_fun_exact(beta_pl_old,c_pl,x_pl_ext_i,y(i))*resp(i,2);
          d2Qdbeta2_pl = d2Qdbeta2_pl + d2elldbeta2_pl_i_fun_exact(beta_pl_old,c_pl,x_pl_ext_i,y(i))*resp(i,2);
        } else {
          dQdbeta_pl = dQdbeta_pl + delldbeta_pl_i_fun_approx(beta_pl_old,c_pl,x_pl_ext_i,y(i))*resp(i,2);
          d2Qdbeta2_pl = d2Qdbeta2_pl + d2elldbeta2_pl_i_fun_approx(beta_pl_old,c_pl,x_pl_ext_i,y(i))*resp(i,2);
        }
      }
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

  double func_val_after_pl = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_old,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb);

  //Update gamma_z with BFGS
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    dQdgamma_z = zeros<mat>(n_mult_z,1);
    d2Qdgamma2_z = zeros<mat>(n_mult_z,n_mult_z);

    for(int i=0; i<n; i++){
      double xtgamma_z_i = (x_mult_z_ext.submat(i,0,i,n_mult_z-1)*gamma_z_old).eval()(0,0);
      double xtgamma_pl_i = (x_mult_pl_ext.submat(i,0,i,n_mult_pl-1)*gamma_pl_old).eval()(0,0);
      arma::mat xtx_gamma_z_i = trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))*x_mult_z_ext.submat(i,0,i,n_mult_z-1);
      denominator = (1 + exp(xtgamma_z_i) + exp(xtgamma_pl_i));
      arma::mat dQdgamma_z_i = trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))*(resp(i,0) - exp(xtgamma_z_i)/denominator);
      arma::mat d2Qdgamma2_z_i = -1.0*xtx_gamma_z_i/denominator;

      dQdgamma_z = dQdgamma_z + dQdgamma_z_i;
      d2Qdgamma2_z = d2Qdgamma2_z + d2Qdgamma2_z_i;
    }

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

  double func_val_after_mult_z = log_lik_fun(gamma_z_old,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb);

  //Update gamma_pl with BFGS
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    dQdgamma_pl = zeros<mat>(n_mult_pl,1);
    d2Qdgamma2_pl = zeros<mat>(n_mult_pl,n_mult_pl);

    for(int i=0; i<n; i++){
      double xtgamma_z_i = (x_mult_z_ext.submat(i,0,i,n_mult_z-1)*gamma_z_old).eval()(0,0);
      double xtgamma_pl_i = (x_mult_pl_ext.submat(i,0,i,n_mult_pl-1)*gamma_pl_old).eval()(0,0);
      arma::mat xtx_gamma_pl_i = trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))*x_mult_pl_ext.submat(i,0,i,n_mult_pl-1);
      denominator = (1 + exp(xtgamma_z_i) + exp(xtgamma_pl_i));
      arma::mat dQdgamma_pl_i = trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))*(resp(i,2) - exp(xtgamma_pl_i)/denominator);
      arma::mat d2Qdgamma2_pl_i = -1.0*xtx_gamma_pl_i/denominator;

      dQdgamma_pl = dQdgamma_pl + dQdgamma_pl_i;
      d2Qdgamma2_pl = d2Qdgamma2_pl + d2Qdgamma2_pl_i;
    }

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

  double func_val_after_mult_pl = log_lik_fun(gamma_z_in,gamma_pl_old,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y,offset_nb);

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


