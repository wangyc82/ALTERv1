tsvd_ini <-
function(Y, X, rank, ini = NULL, method=1){

    Dim <- dim(X);
    n <- Dim[1];
    p <- Dim[2];
    q <- dim(Y)[2];
    ytilde <- t(X)%*%Y;
    sqX <- t(X)%*%X;
    SVD <-svd(X,nu = min(n,p),nv = min(n,p));
 	
    if(method==1){## ridge regression estimate ##
    if(!is.null(ini)) cat("ridge regression estimate!\n")
	tuning <-mean((SVD$d[1:min(n,p)])^2)/p*0.1
	est <- SVD$v%*%diag(1/(SVD$d+tuning))%*%t(SVD$v)%*%ytilde
	rowSS <- apply(est^2,1,sum)
	u.index <- (1:p)[rowSS!=0]
	colSS <- apply(est^2,2,sum)
	v.index <- (1:q)[colSS!=0]
	est0 <- est[u.index,v.index]
	
	ans.ini <- svd(est0,nu=rank,nv=rank);
	ini <- list(    rank = rank,
                    u=ans.ini$u[,1:rank],
	                v=ans.ini$v[,1:rank],
	                d=ans.ini$d[1:rank],
	                u.index= u.index,
	                v.index= v.index)
  
	}
    
    if(method==2) {## with V, use lasso for U and D ##
        #require(glmnet)
        v <- matrix(0,q,rank)
        v[ini$v.index,] <- ini$v
        Yu <- Y%*%v
        u0 <- matrix(NA,p,rank)
        for(i in 1:rank){
            Lambda = exp(seq(-15,10,by=0.5))
            Beta = glmnet(X,Yu[,i],standardize=FALSE,lambda=Lambda)$beta
            BIC = log(apply((X%*%Beta-Yu[,i])^2,2,sum)) + log(n)/n*apply(Beta,2,function(x){sum(x!=0)})
            k0 = which.min(BIC)
            u0[,i] = Beta[,k0]
        }
        est <- u0%*%t(v)
        
        colSS <- apply(v^2,1,sum)
        v.index <- (1:q)[colSS!=0]
        
        std <- function(x){
            if(sum(x^2)>0)
            return(x/sqrt(sum(x^2)));
            return(x);
        }
        
        u <- apply(u0,2,std)
        rowSS <- apply(est^2,1,sum)
        u.index <- (1:p)[rowSS!=0]
        
        d <- sqrt(apply(u0^2,2,sum))
        for(i in 1:rank) d[i] = max(d[i],min(d[d>0]))
        
        ini <-list(rank = rank,
        u=u[u.index,1:rank],
        v=v[v.index,1:rank],
        d=d[1:rank],
        u.index= u.index,
        v.index= v.index)
        
    }
class(ini)="sparseSVD";
return(ini)
}
