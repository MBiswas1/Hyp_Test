#' This function gives us test statistic and corresponding p-value for testing if a scalar predictor is correlated to a functional response.
#' @param data has n data points comprising X, Y, Z, and t
#' @param X is a matrix. Each column of X corresponds to q many covariates of an individual
#' @param Z is a matrix. Each column of Z corresponds to r many nuisance covariates of an individual.
#' @param t is the list of n time points at which the data is available
#' @param Y is a list of n responses at time points corresponding to t
#' @param ker = 1 corresponds to linear kernel for X
#'     = 2 corresponds to quadratic kernel for X
#'     = 0 will correspond to gaussian kernel in the future
#' @param sig is indicator for whether we are using sigma in the test statistic
#' @details A functional response variable Y is regressed on scalar variables X and Z to determine if X is a significant contributor to Y.
#' @return p.val is the p-value for the given data.
#' @return test.stat is the test statistic obtained.
#' @export p_val
#' @author Mityl Biswas

p_val <- function(data, kernel = 0, use.sig = 1 )
{
  require(mvtnorm)
  require(mgcv)
  require(face)
  X <- data$X
  Z <- data$Z
  t.star <- data$t
  Y <- data$Y
  n <- ncol(X)
  m.i <- lengths(t.star)
  seq <- seq(1:n)
  # Z1 <- sapply(seq,function(i) rep(Z[2,i], m.i[i]))
  # Z2 <- sapply(seq,function(i) rep(Z[3,i], m.i[i]))
  # Z1 <- c(unlist(Z1))
  # Z2 <- c(unlist(Z2))
  q <- nrow(Z)
  Z0 <- matrix(NA,q,sum(m.i))
  for(i in 1:q)
    Z0[i,] <- unlist(sapply(seq,function(j) rep(Z[i,j], m.i[j])))
  sim <- 10000
  T1 <- rep(NA, sim)
  t.2 <- unlist(t.star)
  index <- as.numeric(ordered(t.2))
  index <- tolist(index, m.i)
  Y0 <- c(unlist(Y))
  form <- "Y0 ~ s(t.2)"
  for(i in 1:q)
    form <- paste(form, "+s(t.2, by = Z0[",i,",])", sep = "")
  t <- sort(unique(t.2))
  m <- length(t)
  E <- exp(-outer(t, t, function(a, b){abs(a-b)}))
  fit <- gam(as.formula(form))
  Y.center <- fit$residuals

  # seq2 <- seq(1:sum(m.i))
  # seq.t <- cbind(t.star, seq2)
  # t.seq <- seq.t[order(t.star),]

  if(kernel == 0)
  { K <- ker_Xg(X, E)}else
  if(kernel == 1)
 { K <- ker_X(X, E)}else
   if(kernel == 2)
   {K = ker_X2(X, E)}
  index.star <- unlist(sapply(seq, function(i) unlist(index[i])+m*(i-1)))
  K.0 <- kern_miss(index.star, K)

  if(use.sig == 0)
   { Y.sigi <- Y.center}  else
     {
       # Assume est.sigma works
       est.sigma <- est_sig(Y, t, t.star, index)
       sig <- est.sigma$sig
       sigi <- est.sigma$sigi
       Y.sigi <- sigi %*% unlist(Y.center)
       }
       T0 <- test_stat(Y.sigi, K.0)
       m.i.cum <- cumsum(m.i)
       m.i.cum.low <- c(0, m.i.cum[-n])+1


#Get new K
for(i in 1: sim)
{
  sam <- sample(n)
  m.i.cum.new <- m.i.cum[sam]
  m.i.cum.low.new <- m.i.cum.low[sam]
  index.new <- index[sam]
  index.star2 <- unlist(sapply(seq, function(i) unlist(index.new[i])+m*(i-1)))
  index.star3 <- unlist(sapply(seq, function(i) seq(m.i.cum.low.new[i], m.i.cum.new[i])))
  K1 <- kern_miss(index.star2, K)
  Y.new <- Y.center[index.star3]
  # Get test statistic
  T1[i] <- as.vector(test_stat(Y.new, K1))
}
       p.val<- sum(T1>as.numeric(T0))/sim
       return(list(p.val = p.val, test.stat = T0))
}
