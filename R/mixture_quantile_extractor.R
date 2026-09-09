
# CDF of the *discretised* Pareto used by the extreme-value component of the
# mixture (audit N3). The EV pmf is (C/y)^a - (C/(y+1))^a for integer y >= C, so
# F(y) = 1 - (C/(y+1))^a for y >= C and 0 below. Vectorised over the shape `a`
# (mistr::ppareto is not: its scalar `shape <= 0` guard errors on a vector).
ppareto_vec <- function(q, C, a) {
  out <- 1 - (C / (q + 1))^a
  out[q < C] <- 0
  out[!is.finite(out)] <- 0
  out
}

mixture_p <- function(x,pl_alphas,C,nb_mu,nb_alpha,probabilities){
  p <- probabilities[,1] + probabilities[,2] * pnbinom(x,mu=nb_mu,size=1/nb_alpha) + probabilities[,3] * ppareto_vec(x,C,pl_alphas)
  return(p)
}
