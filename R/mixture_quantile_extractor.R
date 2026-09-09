mixture_quantiles_2 <- function(quantile=0.975,pl_alphas,C,nb_mu,nb_alpha,probabilities,theoretical_max=1e20){
  q <- round(mistr::qpareto(probabilities[,3]*quantile,C,pl_alphas))
  q[which(q>theoretical_max)] <- theoretical_max
  ch <- q
  n <- 0
  p <- rep(0,length(q))
  repeat{
    old_p <- p
    p <- mixture_p(q,pl_alphas,C,nb_mu,nb_alpha,probabilities)
    ch[which(p<quantile & abs(p-old_p)<0.01 & quantile-p>0.01)] <- ch[which(p<quantile & abs(p-old_p)<0.01 & quantile-p>0.01)]*2
    ch[which(round(p,8)>=quantile & ch>0)] <- -floor(ch[which(round(p,8)>=quantile & ch>0)]/2)
    ch[which(round(p,8)<=quantile & ch<0)] <- -round(ch[which(round(p,8)<=quantile & ch<0)]/2)
    ch[which(ch<0)] <- 0
    ch[which(ch>theoretical_max)]<-0
    q[which(q<0)] <- 0
    q <- q + ch
    q[which(q>theoretical_max)] <- theoretical_max
    q[which(ch==0 & p<quantile)] <- q[which(ch==0 & p<quantile)] + 1
    n<-n+1
    if(all(round(p,8)>=quantile | q>=theoretical_max) & all(ch>=0)){
      break
    }
  }
  return(q)
}

mixture_quantiles <- function(quantile=0.975,pl_alphas,C,nb_mu,nb_alpha,probabilities,theoretical_max=1e20){
  q <- round(mistr::qpareto(probabilities[,3]*quantile,C,pl_alphas))
  q[which(q>theoretical_max)] <- theoretical_max
  ch <- q
  n <- 0
  p <- rep(0,length(q))
  repeat{
    old_p <- p
    p <- mixture_p(q,pl_alphas,C,nb_mu,nb_alpha,probabilities)
    ch[which(p<quantile & abs(p-old_p)<0.01 & quantile-p>0.01)] <- ch[which(p<quantile & abs(p-old_p)<0.01 & quantile-p>0.01)]*2
    ch[which(round(p,8)>=quantile & ch>0)] <- -floor(ch[which(round(p,8)>=quantile & ch>0)]/2)
    ch[which(round(p,8)<=quantile & ch<0)] <- -round(ch[which(round(p,8)<=quantile & ch<0)]/2)
    ch[which(q<0)] <- 0
    ch[which(ch>theoretical_max)]<-0
    q[which(q<0)] <- 0
    q <- q + ch
    q[which(q>theoretical_max)] <- theoretical_max
    q[which(ch==0 & p<quantile)] <- q[which(ch==0 & p<quantile)] + 1
    n<-n+1
    if(all(round(p,8)>=quantile | q>=theoretical_max) & all(ch>=0)){
      break
    }
  }
  return(q)
}


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
