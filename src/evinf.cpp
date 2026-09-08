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

//[[Rcpp::export]]
double ell_nb_i_fun(arma::vec beta_nb, double alpha_nb, arma::vec x_nb_ext_i, int y_i){
  // int n_beta_nb = beta_nb.size();

  double ell_nb_i = 0;

  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;

  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a);

  if(y_i>0){
    for(int j=0; j<y_i; j++){
      ell_nb_i = ell_nb_i + log(j + 1/alpha_nb);
    }
  }

  if(y_i>0){
    for(int j=1; j<=y_i; j++){
      ell_nb_i = ell_nb_i - log(j);
    }
  }

  ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y_i*log(1 + alpha_nb*mu_i) + y_i*log(alpha_nb) + y_i*log(mu_i);

  return(ell_nb_i);
}

//[[Rcpp::export]]
arma::vec delldtheta_nb_i_fun(arma::vec beta_nb, double alpha_nb, arma::vec x_nb_ext_i, int y_i){
  int n_beta_nb = beta_nb.size();

  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;

  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a);

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
arma::mat d2elldtheta2_nb_i_fun(arma::vec beta_nb, double alpha_nb, arma::vec x_nb_ext_i, int y_i){
  int n_beta_nb = beta_nb.size();

  arma::mat xtb_nb_i = trans(x_nb_ext_i)*beta_nb;

  double a = xtb_nb_i.eval()(0,0);
  double mu_i = exp(a);

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
  double cdivy = c_pl/y_i;
  double cdivyp1 = c_pl/(y_i+1);
  //OBS!! The restriction that y.i>c.pl will be taken care of outside
  double ell_pl_i = log(pow(cdivy,exp_xtb_pl_i)-pow(cdivyp1,exp_xtb_pl_i));
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

//[[Rcpp::export]]
double log_lik_fun(arma::vec gamma_z, arma::vec gamma_pl,arma::vec beta_nb, double alpha_nb, arma::vec beta_pl, double c_pl,arma::mat x_mult_z_ext,arma::mat x_mult_pl_ext,arma::mat x_nb_ext, arma::mat x_pl_ext, arma::vec y){

  int n = x_mult_z_ext.n_rows;
  int n_mult_z = x_mult_z_ext.n_cols;
  int n_mult_pl = x_mult_pl_ext.n_cols;
  int n_nb = x_nb_ext.n_cols;
  int n_pl = x_pl_ext.n_cols;
  arma::mat props = zeros<mat>(n,3) ;

  double denominator = 0;
  for(int i=0; i<n; i++){

    double d_z = exp(trans(gamma_z)*trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))).eval()(0,0);
    double d_pl = exp(trans(gamma_pl)*trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))).eval()(0,0);

    denominator = 1 + d_z + d_pl;
    props(i,0) = d_z/denominator;
    props(i,1) = 1/denominator;
    props(i,2) = d_pl/denominator;
  }

  double func_val = 0;

  for(int i=0; i<n; i++){
    arma::mat x_nb_ext_i = trans(x_nb_ext.submat(i,0,i,n_nb-1));
    double xtb_nb_i = (trans(x_nb_ext_i)*beta_nb).eval()(0,0);
    double mu_i = exp(xtb_nb_i);
    double ell_nb_i = 0;
    if(y(i)>0){
      for(int j=0; j<y(i); j++){
        ell_nb_i = ell_nb_i + log(j + 1/alpha_nb);
      }
    }

    if(y(i)>0){
      for(int j=1; j<=y(i); j++){
        ell_nb_i = ell_nb_i - log(j);
      }
    }

    ell_nb_i = ell_nb_i - (1/alpha_nb)*log(1 + alpha_nb*mu_i) - y(i)*log(1 + alpha_nb*mu_i) + y(i)*log(alpha_nb) + y(i)*log(mu_i);

    if(y(i)==0){
      func_val = func_val + log(props(i,0) + props(i,1)*exp(ell_nb_i));
    }else if(y(i)>0 && y(i)<c_pl){
      func_val = func_val + log(props(i,1)*exp(ell_nb_i));
    }else{
       arma::mat x_pl_ext_i = trans(x_pl_ext.submat(i,0,i,n_pl-1));
       double xtb_pl_i = (trans(x_pl_ext_i)*beta_pl).eval()(0,0);

      double exp_xtb_pl_i = exp(xtb_pl_i);
      double cdivy = c_pl/y(i);
      double cdivyp1 = c_pl/(y(i)+1);
      double ell_pl_i = log(pow(cdivy,exp_xtb_pl_i) - pow(cdivyp1,exp_xtb_pl_i));

      func_val = func_val + log(props(i,1)*exp(ell_nb_i) + props(i,2)*exp(ell_pl_i));
    }
  }

  return(func_val);
}

//[[Rcpp::export]]
List update_bfgs_fun(arma::vec gamma_z_in, arma::vec gamma_pl_in,arma::vec beta_nb_in, double alpha_nb_in, arma::vec beta_pl_in, double c_pl,arma::mat x_mult_z_ext,arma::mat x_mult_pl_ext,arma::mat x_nb_ext, arma::mat x_pl_ext, arma::vec y, double max_upd_par, int no_m_bfgs_steps){

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

    double d_z = exp(trans(gamma_z_old)*trans(x_mult_z_ext.submat(i,0,i,n_mult_z-1))).eval()(0,0);
    double d_pl = exp(trans(gamma_pl_old)*trans(x_mult_pl_ext.submat(i,0,i,n_mult_pl-1))).eval()(0,0);

    denominator = 1 + d_z + d_pl;
    props(i,0) = d_z/denominator;
    props(i,1) = 1/denominator;
    props(i,2) = d_pl/denominator;
  }

  arma::mat resp = zeros<mat>(n,3);
  double marginal_yx_i = 0;

  for (int i=0; i<n; i++){
  // Marginal probability mass function of y and x (sum over components) - marginal._x_i
    if(y(i)==0){
      double d_nb = exp(ell_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i)));
      marginal_yx_i = props(i,0) + props(i,1)*d_nb;
      resp(i,0) = props(i,0)/marginal_yx_i;
      resp(i,1) = props(i,1)*d_nb/marginal_yx_i;
      resp(i,2) = 0;
    }else if(y(i)>0 && y(i)<c_pl){
      resp(i,0) = 0;
      resp(i,1) = 1;
      resp(i,2) = 0;
    }else if(y(i)>=c_pl){
      double d_nb = exp(ell_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i)));
      double d_pl = exp(ell_pl_i_fun(beta_pl_old,c_pl,trans(x_pl_ext.submat(i,0,i,n_pl-1)),y(i)));
      marginal_yx_i = props(i,1)*d_nb  + props(i,2)*d_pl;
      resp(i,0) = 0;
      resp(i,1) = props(i,1)*d_nb/marginal_yx_i;
      resp(i,2) = props(i,2)*d_pl/marginal_yx_i;
    }
  }

  //Calculate function value before the algorithm starts
  double func_val_before_bfgs = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y);

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

      dQdtheta_nb = dQdtheta_nb + delldtheta_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i))*resp(i,1);
      d2Qdtheta2_nb = d2Qdtheta2_nb + d2elldtheta2_nb_i_fun(beta_nb_old,alpha_nb_old,trans(x_nb_ext.submat(i,0,i,n_nb-1)),y(i))*resp(i,1);
    }

    change_nb_bfgs = -inv(d2Qdtheta2_nb)*dQdtheta_nb;

    maxabschange = max(abs(change_nb_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_nb_bfgs = max_upd_par/maxabschange*change_nb_bfgs;
    }

    beta_nb_old = beta_nb_old + change_nb_bfgs.submat(0,0,n_nb-1,0);
    alpha_nb_old = alpha_nb_old + change_nb_bfgs.submat(n_nb,0,n_nb,0).eval()(0,0);

  }

  double func_val_after_nb = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_old,alpha_nb_old,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y);

  //Update beta_pl with BFGS
  for(int i_bfgs=1; i_bfgs<=no_m_bfgs_steps; i_bfgs++){
    d2Qdbeta2_pl = zeros<mat>(n_pl,n_pl);
    dQdbeta_pl = zeros<mat>(n_pl,1);

    for(int i=0; i<n; i++){
      if(y(i)>=c_pl){
        dQdbeta_pl = dQdbeta_pl + delldbeta_pl_i_fun_approx(beta_pl_old,c_pl,trans(x_pl_ext.submat(i,0,i,n_pl-1)),y(i))*resp(i,2);
        d2Qdbeta2_pl = d2Qdbeta2_pl + d2elldbeta2_pl_i_fun_approx(beta_pl_old,c_pl,trans(x_pl_ext.submat(i,0,i,n_pl-1)),y(i))*resp(i,2);
      }
    }

    change_pl_bfgs = -inv(d2Qdbeta2_pl)*dQdbeta_pl;

    maxabschange = max(abs(change_pl_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_pl_bfgs = max_upd_par/maxabschange*change_pl_bfgs;
    }
    beta_pl_old = beta_pl_old + change_pl_bfgs;

  }

  double func_val_after_pl = log_lik_fun(gamma_z_in,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_old,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y);

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

    change_mult_z_bfgs = -inv(d2Qdgamma2_z)*dQdgamma_z;


    maxabschange = max(abs(change_mult_z_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_mult_z_bfgs = max_upd_par/maxabschange*change_mult_z_bfgs;
    }
    gamma_z_old = gamma_z_old + change_mult_z_bfgs;

  }

  double func_val_after_mult_z = log_lik_fun(gamma_z_old,gamma_pl_in,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y);

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

    change_mult_pl_bfgs = -inv(d2Qdgamma2_pl)*dQdgamma_pl;

    maxabschange = max(abs(change_mult_pl_bfgs)).eval()(0,0);
    if(maxabschange>max_upd_par){
      change_mult_pl_bfgs = max_upd_par/maxabschange*change_mult_pl_bfgs;
    }
    gamma_pl_old = gamma_pl_old + change_mult_pl_bfgs;

  }

  double func_val_after_mult_pl = log_lik_fun(gamma_z_in,gamma_pl_old,beta_nb_in,alpha_nb_in,beta_pl_in,c_pl,x_mult_z_ext,x_mult_pl_ext,x_nb_ext,x_pl_ext,y);

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


