sira <-
function(ytilde,sqX,X_norm,ini,thr,iter.max = 1000,rel.diff=0.001){
## ------------------------------------------------------------------------
## Sparse Iteration Regression Algorithm(SIRA)
## ------------------------------------------------------------------------
##
## DESCRIPTION:
##         This function will peform Sparse Iteration Regression Algorithm outlined as in the paper
##       
##
## INPUT ARGUMENTS:
## u           	 	The p*rank real unitary matrix 
## v            	The q*rank real unitary matrix
## d            	The square matrix with singular values on the diagonal
## rank 		The pre-specified designed number for rank
## thr          	The pre-specified pair of thresholding parameter
## iter.max     	The maximum number of iteration
## rel.diff     	The difference between result from time t and (t-1)
##
##
## OUTPUT ARGUMENTS:
## beta  		The estimated coefficient matrix resulted from (SIRA)
## u			The estimated p*rank real unitary matrix resulted from (SIRA)
## v			The estimated q*rank real unitary matrix resulted from (SIRA)
## d			The estimated square matrix with singular values on the diagonal resulted from (SIRA)
## convergence		The indicator of the convergence(TRUE/FALSE) for pre-specified paramter "thr"
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

p <- dim(ytilde)[1];
q <- dim(ytilde)[2];
if(length(thr)==1) {thr <- rep(thr,2)}
rank <- ini$rank
u <- ini$u
v <- ini$v
d <- ini$d
u.index <- ini$u.index
v.index <- ini$v.index
convergence <- FALSE;
#source("soda.r");
################## iteration starts here
################## --------------------------------------------
u.temp <- matrix(0,p,rank);
v.temp <- matrix(0,q,rank);
u.temp[u.index,] <- u; v.temp[v.index,] <- v; d.temp <- d;

iter <- 1;
while(iter<=iter.max&&convergence==FALSE){
  if(rank==1) {ud.ini<- d.temp*u.temp}
  if(rank>1)  {ud.ini<- u.temp%*%diag(d.temp)}
  ytilde.u <- ytilde[,v.index]%*%v;
  
################## Step in (4)	    
  u <- (diag(p)-sqX/X_norm^2)%*%ud.ini + ytilde.u/X_norm^2;
################## Step in (5)	 
  for(r in 1:rank)  {u[abs(u[,r])<= thr[1]*sqrt(sum(u[,r]^2)),r] <- 0 }
################## Step in (6)	
  U <- soda(u);
  u <- U$Q;
  u.index <- U$index;
################## Step in (1) 
  A <- t(u)%*%sqX[u.index,u.index]%*%u
  tuning <- 0*min(diag(A)[diag(A)>0]) ##no change
  v <- t(solve(A+tuning*diag(rank))%*%t(u)%*%ytilde[u.index,]);

################## Step in (2)
  for(r in 1:rank)  {v[abs(v[,r])<=thr[2]*sqrt(sum(v[,r]^2)),r] <- 0}
################## Step in (3)		  
  V2 <- soda(v);
  v2 <- V2$Q;
  v.index <- V2$index;
  d <- abs(diag(t(v[v.index,])%*%v2));
  v <- v2;
  Ord <- order(d,decreasing=TRUE)
  if(sum((1:rank-Ord)^2)>0){
      u <- u[,Ord]
      v <- v[,Ord]
      d <- d[Ord]
    }

		iter <- iter + 1;
		ufull <- matrix(0,p,rank);
		vfull <- matrix(0,q,rank);
        u[abs(u)<=0.000001] <- 0
        v[abs(v)<=0.000001] <- 0
		ufull[u.index,] <- u;
		vfull[v.index,] <- v;
		diff <- sum((ufull-u.temp)^2)/sum(u^2) + sum((vfull-v.temp)^2)/sum(v^2) + sum((d-d.temp)^2)/sum(d^2);
		ind <- sum((ufull!=0&&u.temp==0)) + sum((ufull==0&&u.temp!=0)) + sum((vfull!=0&&v.temp==0)) + sum((vfull==0&&v.temp!=0));

		u.temp <- ufull; v.temp <- vfull; d.temp <- d;
		if(diff <= rel.diff&&ind==0) convergence <- TRUE;
  }
################## iteration ends here
################## --------------------------------------------
ans <- list(rank = rank,
            u= u,
            v= v,
            d = d,
            u.index = u.index,
            v.index = v.index);
class(ans)="sparseSVD";
return(list(est=ans,convergence=convergence));
}
