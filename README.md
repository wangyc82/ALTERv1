ALTERv1

this procedure is to get the expanded lncRNA target genes, expanded TF target genes, and the lncRNA-TF associations

the main function is cTSVD, rise_lncRNA_coding is initial lncRNA target genes from RISE (mouse-RISE-data.RData), and TF.gene.R1 is TF target genes from TTRUEST (mouse-TTRUST-data.RData)

Usage

correlations<-cTSVD(RNA.data,rise_lncRNA_coding,TF.gene.R1)

It will generate lncRNA-gene relations (corGandL), TF-gene relations (corGandT), and lncRNA-TF associations (corLandT)

It also includes the related R function in TSVD package, and mouse RNA-seq data (MouseRPKMs) as an example.

The RNA-seq data across human and mouse development stages can be downloaded from supllementary data in Nature paper, SuppData2_RPKMs.

The Nature paper is Cardoso-Moreira, M. et al. (2019) Gene expression across mammalian organ development. Nature, 571(7766), 505-509.
