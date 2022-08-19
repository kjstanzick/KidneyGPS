######################
###generate genelist with gene correct names

setwd("/stk05236/clean_GPS")

genes = fread( "/stk05236/inputs/genes/glist-hg19.tw.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)


hgnc_genes = read.table("/stk05236/inputs/hgnc/hgnc_custom_220818.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

hgnc_genes = hgnc_genes[which(hgnc_genes$Status!=""),]

for(i in 1:nrow(hgnc_genes)){
  
  if(hgnc_genes$Previous.symbols[i]!=""){
    prev_genes = unlist(strsplit(hgnc_genes$Previous.symbols[i], ", "))
    if(i==1) {
      gene = c(hgnc_genes$Approved.symbol[i],prev_genes)
      real_genes = hgnc_genes[i,c(1,2,7:10)]
      for(x in 1:length(prev_genes)){
        real_genes = rbind(real_genes,hgnc_genes[i,c(1,2,7:10)])
      }
    } else {
      gene = c(gene,hgnc_genes$Approved.symbol[i],prev_genes)
      for(x in 1:(length(prev_genes)+1)){
        real_genes = rbind(real_genes,hgnc_genes[i,c(1,2,7:10)])
      }
    }
  } else {
    if(i==1) {
      gene = c(hgnc_genes$Approved.symbol[i])
      real_genes = hgnc_genes[i,c(1,2,7:10)]
      
    } else {
      gene = c(gene,hgnc_genes$Approved.symbol[i])
      real_genes = rbind(real_genes,hgnc_genes[i,c(1,2,7:10)])
    }
  }
}

real_genes['genes']=gene

real_genes = real_genes[which(!is.na(real_genes$Approved.symbol)),]

real_genes = real_genes[which(!is.na(real_genes$genes)),]

for(i in 1:nrow(genes)){
  
  if(!any(genes$Gene[i]==real_genes$Approved.symbol) & any(genes$Gene[i]==real_genes$genes)){
    
    genes$Gene[i]=real_genes$Approved.symbol[which(real_genes$genes==genes$Gene[i])[1]]
  }
}

genes_sorted <- genes[order(genes$Chr, genes$Pos1), ]
write.table(genes_sorted, file="./inputs/glist.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(real_genes, file="./inputs/genes_an_ids.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
