#Script by: Hitesh Kore
#This script computes the diferentially expressed genes/transcripts across two different conditions i.e. normal and treatment using DESeq2
#It is based on a blog on the URL: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#Rscript DeSeqDifferentialExpression.R <Salmon_outdir> <Design_matrix> <geneTranscriptmap> <outdir> <quantification_type>
#quantification- gene/transcript
#Note: tximport library reads the salmon quantification output in each sample directory and compiles a sample-wise quantification matrix for DESeq2 calculations  

args = commandArgs(trailingOnly=TRUE)

#DeSeq2 using Salmon
DESeqSalmon<- function(sdir,smat,trmap,atype) {
 rm(list = ls()) #clear workspace
#Libraries required
library(tximport)
library(DESeq2)
tx2gene=read.delim(trmap,sep="\t",header=FALSE)
colnames(tx2gene) <-c("TXNAME","GENEID")
##sample matrix
sample_matrix=read.table(smat,header=T)

##path for quant.sf files
files=c()
for(i in 1:length(sample_matrix$samples)){
 files[i]=paste(sdir,sample_matrix$samples[i],"_1\\quant.sf", sep="")
 
 }


file.exists(files)

#labelling files with sample names

names(files)<-sample_matrix$samples
#files


	if	(atype=="transcript"){
		tximp=tximport(files, type="salmon",txOut=TRUE)
		#print(tximp$counts[1:5,1:5])
	
	}

	if	(atype=="gene"){
		##gene level summarization
		tximp <- tximport(files, type = "salmon", tx2gene = tx2gene)
		print(tximp$counts[1:5,1:5])
	
	}
	
	ddsTxi=DESeqDataSetFromTximport(tximp, colData=sample_matrix, design=~condition)
	ddsTxi.filtered=ddsTxi[apply(counts(ddsTxi,normalized=TRUE),1,median) > 1, ] ##Removing lowely expressed genes
	#dim(counts(ddsTxi.filtered))
	diff=DESeq(ddsTxi.filtered)
	res=results(diff)
	summary(res)
	#class(res)
	res.lna=na.omit(res)
	sig.res.lna=res.lna[res.lna$padj < 0.05, ]
	write.table(sig.res.lna,file=paste(outdir,"Sig_diff_LFC2.csv",sep=""),sep="\t")

}




if (length(args)==0) {
  stop("Invalid Arguments", call.=FALSE)
} else if (length(args)==5) {
  # default output file
 
  sdir=args[1]
  smat=args[2]
  trmap=args[3]
  outdir=args[4]
  atype=args[5]  
  DESeqSalmon(sdir,smat,trmap,atype)
  
}

