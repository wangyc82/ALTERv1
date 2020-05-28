tsvd_wrapper <- 
function(Y,X,rank=NULL,BICtype=2,thr=NULL, control = list(thr.lower = -10,
         thr.upper = -1, thr.by = 0.5)){
         	
            if(is.null(rank)) rank <- select_rank(Y,X)
         	ini1 <- tsvd_ini(Y,X,rank,method=1)
         	ans1 <- tsvd(Y,X,rank,ini=ini1,BICtype,thr,control)
            ini2 <- tsvd_ini(Y,X,rank,ini=ans1$est,method=2)
            ans2 <- tsvd(Y,X,rank,ini=ini2,BICtype,thr,control)
            return(ans2)
         	
         }