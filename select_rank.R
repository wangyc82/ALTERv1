select_rank <-
function(Y,X,constant=4){
    #require(MASS)
  Dim <- dim(Y)
  n <- Dim[1]
  q <- Dim[2]
  p <- dim(X)[2]
  
  sqX <- t(X)%*%X
  ytilde <- t(X)%*%Y
  sqX_inv <- ginv(sqX)
  if(q<=p){
  P <- t(ytilde)%*%sqX_inv%*%ytilde
  P <- (P+t(P))/2
  Eigen <- eigen(P,symmetric=TRUE,only.values=TRUE)
  Eigen$values[Eigen$values<=0] = 0
  d <- Eigen$values
  }
  if(q>p){
      Eigen <- eigen(X%*%sqX_inv%*%t(X))
      Eigen$values[Eigen$values<=0] <- 0
      temp <- Eigen$vectors%*%diag(sqrt(Eigen$values))%*%t(Eigen$vectors)
      temp <- t(Y)%*%temp
      temp <- svd(temp,nu=0)
      d <- temp$d^2
      }
  sigmahat <- mad(as.vector(Y))
  cut.off <- constant*max(q,n)*sigmahat^2
  return(sum(d>cut.off))
}
