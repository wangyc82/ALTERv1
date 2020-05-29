Correlations<-cTSVD(MouseRPKMs,rise_lncRNA_coding,TF.gene.R1) 
{
  # this procedure is to get the lncRNA-TF associations based on developemtnal RNA-seq data (take human data, MouseRPKMs as example)
  # rise_lncRNA_coding is RISE data
  # TF.gene.R1 is TTRUEST TF target genes
  
  mouse_gene_exp<-data.matrix(MouseRPKMs[,-1])
  # convert ensemble ID to gene symbol for mouse
  
  library("EnsDb.Mmusculus.v79")
  ensembleIDs<-MouseRPKMs$X1
  gene_symbol<-select(EnsDb.Mmusculus.v79, keys=ensembleIDs, columns=c("GENEID","SYMBOL","GENENAME","ENTREZID"),keytype="GENEID")
  lab<-match(ensembleIDs,gene_symbol$GENEID)
  mouse_gene_symbol<-gene_symbol$SYMBOL[lab]
  mouse_gene_exp<-mouse_gene_exp[-which(mouse_gene_symbol %in% "NA"),]
  
  # the prcedure for 3step TSVD
  `%notin%` <- Negate(`%in%`)
  
  tissue_type<-sapply(1:length(strsplit(colnames(mouse_gene_exp),"[.]")),function(x) strsplit(colnames(mouse_gene_exp),"[.]")[[x]][1])
  TT<-names(table(as.factor(tissue_type)))[which(table(as.factor(tissue_type))>40)] # only for organs with bio-samples larger than 40.
  source('~/Documents/tsvd/R/select_rank.R')
  source('~/Documents/tsvd/R/tsvd.R')
  source('~/Documents/tsvd/R/tsvd_ini.R')
  source('~/Documents/tsvd/R/sira.R')
  source('~/Documents/tsvd/R/soda.R')
  
  dev.stage<-list()
  Cgl<-list()
  Cgt<-list()
  Clt<-list()
  for (m in 1:length(TT)) {
    mouse_gene_exp_T<-mouse_gene_exp[,which(tissue_type %in% TT[m])]
    dev.stage[[m]]<-sapply(1:length(strsplit(colnames(mouse_gene_exp_T),"[.]")),function(x) strsplit(colnames(mouse_gene_exp_T),"[.]")[[x]][2])
    mouse_gene_exp_T1<-mouse_gene_exp_T[which(apply(mouse_gene_exp_T,1,sd)>1.5),]
    TF.gene.R1<-TF.gene.R[which(TF.gene.R$TF %in% intersect(rownames(mouse_gene_exp_T1),levels(as.factor(TF.gene.R$TF)))),]
    mouse_gene_exp_T2<-mouse_gene_exp_T1[which(rownames(mouse_gene_exp_T1) %notin% levels(as.factor(TF.gene.R1$TF))),]
    mouse_TF_exp_T<-mouse_gene_exp_T1[intersect(rownames(mouse_gene_exp_T1),levels(as.factor(TF.gene.R$TF))),]
    mouse_lncRNA_exp_known_T<-mouse_lncRNA_exp_known[,which(tissue_type %in% TT[m])]
    l<-sapply(1:nrow(mouse_lncRNA_exp_known_T),function(x) length(which(mouse_lncRNA_exp_known_T[x,]==0)))
    mouse_lncRNA_exp_known_T1<-mouse_lncRNA_exp_known_T[which(l<ncol(mouse_lncRNA_exp_known_T)),]
    ini0<-matrix(0,nrow(mouse_lncRNA_exp_known_T1),nrow(mouse_gene_exp_T2))
    rownames(ini0)<-rownames(mouse_lncRNA_exp_known_T1)
    colnames(ini0)<-rownames(mouse_gene_exp_T2)
    for (k in 1:nrow(rise_lncRNA_coding)) {ini0[which(rownames(ini0) %in% rise_lncRNA_coding$gene_name1[k]),which(colnames(ini0) %in% rise_lncRNA_coding$gene_name[k])]<-1}
    rank<-select_rank(t(mouse_gene_exp_T2),t(mouse_lncRNA_exp_known_T1))
    est0<-svd(ini0,nu=rank,nv=rank)
    ini1 <- tsvd_ini(t(mouse_gene_exp_T2),t(mouse_lncRNA_exp_known_T1),rank,method=1)
    ini1$u<-est0$u[,1:rank]
    ini1$v<-est0$v[,1:rank]
    ini1$d<-est0$d[1:rank]
    ans1 <- tsvd(t(mouse_gene_exp_T2),t(mouse_lncRNA_exp_known_T1),rank,ini=ini1)
    Cgl[[m]]<-ans1
    ini01<-matrix(0,nrow(mouse_TF_exp_T),nrow(mouse_gene_exp_T2))
    rownames(ini01)<-rownames(mouse_TF_exp_T)
    colnames(ini01)<-rownames(mouse_gene_exp_T2)
    for (k in 1:nrow(TF.gene.R1)) {ini01[which(rownames(ini01) %in% TF.gene.R1$TF[k]),which(colnames(ini01) %in% TF.gene.R1$Target[k])]<-1}
    rank1<-select_rank(t(mouse_gene_exp_T2),t(mouse_TF_exp_T))
    est01<-svd(ini01,nu=rank1,nv=rank1)
    ini2 <- tsvd_ini(t(mouse_gene_exp_T2),t(mouse_TF_exp_T),rank1,method=1)
    ini2$u<-est01$u[,1:rank1]
    ini2$v<-est01$v[,1:rank1]
    ini2$d<-est01$d[1:rank1]
    ans2 <- tsvd(t(mouse_gene_exp_T2),t(mouse_TF_exp_T),rank1,ini=ini2)
    Cgt[[m]]<-ans2
    ini02<-cor(t(ans2$C),t(ans1$C))
    ini02[is.na(ini02)]<-0
    rank2<-select_rank(t(mouse_lncRNA_exp_known_T1),t(mouse_TF_exp_T))
    est02<-svd(ini02,nu=rank2,nv=rank2)
    ini3 <- tsvd_ini(t(mouse_lncRNA_exp_known_T1),t(mouse_TF_exp_T),rank2,method=1)
    ini3$u<-est02$u[,1:rank2]
    ini3$v<-est02$v[,1:rank2]
    ini3$d<-est02$d[1:rank2]
    ans3 <- tsvd(t(mouse_lncRNA_exp_known_T1),t(mouse_TF_exp_T),rank2,ini=ini3)
    Clt[[m]]<-ans3
    rm(mouse_gene_exp_T,mouse_gene_exp_T1,TF.gene.R1,mouse_gene_exp_T2,mouse_TF_exp_T,mouse_lncRNA_exp_known_T,l,mouse_lncRNA_exp_known_T1)
    rm(ini1,ini2,ini3,est0,est01,est02,ini0,ini01,ini02,rank,rank1,rank2,ans1,ans2,ans3)
    }
  return(list(corGandL=Cgl, corGandT=Cgt, corLandT=Clt))
  }

