setwd("/stk05236/")

cred_var_99 = read.table("./lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

neptune_files = list.files("./inputs/nephQTL")

##############################################
# add cpaid:
library(data.table)

ref = fread( "./lmm/03_eval/european/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

ref_cpaid = ref[,c(1,15)]

cred_var_99 = merge (cred_var_99,ref_cpaid, by.x="rsid", by.y="RSID",all.x=TRUE,all.y=FALSE)


	



#########################################################
###############
#ENTREZ_ID CHROM POS SNP REF ALT BETA SE T_STAT P cpaid FDR.bh Approved.symbol
#
#benötigt wird 5,6,7,8,10,11,12
#
#peerEQTL_G.2017-09-25.reformat.bgzip.cis.cpaid_hugo_fdr.txt  -> glomerular   
#peerEQTL_T.2017-09-25.reformat.bgzip.cis.cpaid_hugo_fdr.txt  -> tubulointerstitial 
#
#
#########################################################
library(data.table)

neptune_spalten = c(5,6,7,8,10,11,12)

neptune_table = data.frame()

input_tissue = c("glomerular","tubulointerstitial")


for (y in 1:length(neptune_files)){
	
	input_file = paste("./inputs/nephQTL/",neptune_files[y],sep="")
		
	input = fread(input_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
		
	input = input[which(input$FDR.bh<0.05),]		## nominale Signifikanz False Discovery Rate
	
	idx = which (input$cpaid %in% cred_var_99$MarkerName)
	
	neptune_table = input[idx,]
	
	
	#for (i in 1:nrow(cred_var_99)){
	
		
		#if (any(input$cpaid == cred_var_99$cpaid[i])){
		
			#idx = which(input$cpaid == cred_var_99$cpaid[i])
		
			#if(nrow(neptune_table)==0){
				
				#neptune_table = input[idx,neptune_spalten]
					
			#}else{
				
				#neptune_table = rbind(neptune_table,input[idx,neptune_spalten], make.row.names = FALSE, stringsAsFactors = FALSE)
					
			#}
				 
				
	tissue = rep(input_tissue[y],length(idx))
	
	neptune_table['tissue'] <- tissue
	
	neptune_table = merge(neptune_table,cred_var_99, by.x = "cpaid", by.y="MarkerName", all.x=TRUE, all.y=FALSE)
	
	output=paste("./lmm/all_loci/neptune/cred_var_with_nephQTL_",tissue[y],".txt",sep="")
	
	write.table(neptune_table, file=output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
	#write.table(neptune_table, file="./eQTL/neptune/all_cred_var_with_nephQTL.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, append=TRUE)
}
			


