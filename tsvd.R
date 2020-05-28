tsvd <-
function(Y, X, rank,ini, BICtype = 2, thr = NULL, control=list(thr.lower = -10,thr.upper=-1,thr.by=0.5)){
## ------------------------------------------------------------------------
## MAIN FUNCTION
## ------------------------------------------------------------------------
##
## DESCRIPTION:
##         This is the main algorithm function
##       
## 
## INPUT ARGUMENTS:
## X            The design matrix
## Y            The response matrix
## rank         The desired rank 
## thr		The two input thresholding paramters eg: thr <- c(0.1,0.1)
##			Default: NUll            
## 			If not specified and p<100, the parameter space is thr1=thr2=seq(0.04,0.24,by=0.02)
##			If not specified and p>=100, the parameter space is thr1=thr2=seq(0.01,0.14,by=0.01)
##
## OUTPUT ARGUMENTS:
## u            The estimated p*rank real unitary matrix (matrix of left-singular vectors)
## v  		The estimated q*rank real unitary matrix (matrix of right-singular vectors)	
## d         	The estimated square matrix with singular values on the diagonal (dimention=rank*rank)
## C		The estimated coefficient matrix(beta = u%*%d%*%v)
## opt          The best paramter pair choosed with BIC
##
##
##
##          
## AUTHORS:
##   
##
## REFERENCE:
##   
##   
##   
##
##
## DEVELOPMENT: 
##              
##
  
 


  Dim <- dim(X);
  n <- Dim[1];
  p <- Dim[2];
  q <- dim(Y)[2];
  ytilde <- t(X)%*%Y;          							# storing values that only need to compute once X'X and X'Y
  sqX <- t(X)%*%X;
  SVD <-svd(X,nu = min(n,p),nv = min(n,p))  #singular value decomposition with n left and right singular vectors


  thr.lower <- control$thr.lower
  thr.upper <- control$thr.upper
  thr.by <- control$thr.by
################## search for optimal thresholding parameters
################## --------------------------------------------
#source("sira.r");						# source the iteration algorithm
## a function to calculate df
df <- function(u){
    
    if(is.vector(u)) u <- matrix(u,nrow=1)
    ind <- colSums(u!=0)
    const <- dim(u)[2]-sum(ind==0)
    return(max(const*(const-1)/2,0))
}



	if(!is.null(thr))							# use the specified pair of parameters
		{
		if(is.matrix(thr)) {Thr1<- thr[,1]; Thr2<-thr[,2];}
		else {Thr1<-thr[1]; Thr2<-thr[2];}
		}else {									# search through the parameter space
		  if(length(thr.lower)==1) {thr.lower <- rep(thr.lower,2)}
          if(length(thr.upper)==1) {thr.upper <- rep(thr.upper,2)}
          if(length(thr.by)==1)    {thr.by <- rep(thr.by,2)}
		  Thr1 <- exp(seq(thr.lower[1],thr.upper[1],by=thr.by[1]))
          Thr2 <- exp(seq(thr.lower[2],thr.upper[2],by=thr.by[2]))
		}
	

	BIC <- matrix(100000,length(Thr1),length(Thr2))
	for(i in 1:length(Thr1))
		{
		  for(j in 1:length(Thr2))
			{
			thr <- c(Thr1[i],Thr2[j]);

################
 #ERROR HANDLING
  # ##http://stackoverflow.com/questions/8093914/skip-to-next-value-of-loop-upon-error-in-r-trycatch
      
  possibleError <- tryCatch({
    res <- sira(ytilde,sqX,SVD$d[1],ini, thr, iter.max = 200, rel.diff=0.001);
  }, error=function(e) e) 
  
  if(inherits(possibleError, "error")) {print(paste("iteration stops for tuning parameters =",signif(thr[1],3),signif(thr[2],3)));  next}
     ans <- res$est

      Df <- sum((ans$u!=0)) + sum((ans$v!=0)) - df(ans$u) - df(ans$v) - ans$rank;
      if(rank == 1) C <- ans$d*ans$u%*%t(ans$v);
      if(rank >1)   C <- ans$u%*%diag(ans$d)%*%t(ans$v);
      SSE <- sum((Y[,ans$v.index]-X[,ans$u.index]%*%C)^2)
      if(length(ans$v.index)<q) SSE <- SSE + sum(Y[,setdiff(1:q,ans$v.index)]^2)
      BIC[i,j] <- log(SSE) + log(q*n)/(q*n)*Df

      
      if(BICtype==2){
          
          dfu <- sum((ans$u!=0)) - df(ans$u)
          dfv <- sum((ans$v!=0)) - df(ans$v)
          BIC[i,j] <- log(SSE) + log(ans$rank*n)/(ans$rank*n)*dfu  + log(q*n)/(q*n)*dfv
      }
             
			}## for j
		}## for i
  
	  loc <- which.min(BIC);						# parameter choosing using BIC
	  BICbest <- BIC[loc];
          loc1 <- loc%%length(Thr1);
	  if(loc1==0) loc1 <- length(Thr1);
	  loc2 <- (loc-loc1)/length(Thr1) + 1;
################## --------------------------------------------		  
################## end search for optimal thresholding parameters
  
	  

	  
################## obtain the estimates
################## --------------------------------------------	
	  thr.opt <- c(Thr1[loc1],Thr2[loc2]);
    
	  res <- sira(ytilde,sqX,SVD$d[1],ini,thr.opt,iter.max = 1000,rel.diff=0.0001);
	  if(!res$convergence) {cat("For ", thr.opt[1]," ",thr.opt[2]," ,estimate fails to converge during iteration!\n",sep="");}
################## --------------------------------------------		  
################## algorithm finished
      ans = res$est
      C <- matrix(0,p,q)
      if(rank == 1) C[ans$u.index,ans$v.index] <- ans$d*ans$u%*%t(ans$v);
      if(rank >1)   C[ans$u.index,ans$v.index] <- ans$u%*%diag(ans$d)%*%t(ans$v);
   
   result <- list(est = ans,C = C,thr = thr.opt, BICtype = BICtype,BIC = BICbest);
   return(result);
	
 }
