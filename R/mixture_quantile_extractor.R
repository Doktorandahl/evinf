
# Pareto(scale = C, shape = a) CDF, vectorised over the shape (mistr::ppareto is
# not: its scalar `shape <= 0` guard errors on a vector).
ppareto_vec <- function(q, C, a) {
  out <- 1 - (C / q)^a
  out[q < C] <- 0
  out[!is.finite(out)] <- 0
  out
}

mixture_p <- function(x,pl_alphas,C,nb_mu,nb_alpha,probabilities){
  p <- probabilities[,1] + probabilities[,2] * pnbinom(x,mu=nb_mu,size=1/nb_alpha) + probabilities[,3] * ppareto_vec(x,C,pl_alphas)
  return(p)
}
