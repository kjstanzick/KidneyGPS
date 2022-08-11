setwd("/stk05236/")

cred_var_99 = read.table("./lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gtex_files = list.files("./inputs/gtex")

gtex_files = gtex_files[-which(grepl("10rows",gtex_files))]

##########noch .rep und .out datei rausschmeißen, stehen ganz oben:

gtex_files = gtex_files[-c(1,2)]

#########################################################
###############
#gene_id	variant_id	tss_distance	ma_samples	ma_count	maf	pval_nominal	slope	slope_se	pval_nominal_threshold	min_pval_nominal	pval_beta	ref	alt	num_alt_per_site	rs_id_dbSNP151_GRCh38p7	gene
#
#benötigt wird 6,7,8,9,13,14,15,16
#
#
#
#
#########################################################
library(data.table)

#gtex_spalten = c(6,7,8,9,13,14,15,16)

gtex_table = data.frame()

for (y in 1:length(gtex_files)){

	input_file = paste("./inputs/gtex/",gtex_files[y],sep="")
		
	input = fread(input_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
		
	input = input[-which(is.na(input$rs_id_dbSNP151_GRCh38p7)),]
	
	idx = which (input$rs_id_dbSNP151_GRCh38p7 %in% cred_var_99$rsid)
	
	gtex_table = input[idx,]
	
	#for (i in 1:nrow(cred_var_99)){
		
		#idx = which(input$rs_id_dbSNP151_GRCh38p7 == cred_var_99$rsid[i])
		#if (any(input$rs_id_dbSNP151_GRCh38p7 == cred_var_99$rsid[i])){
				
			#print(paste("i = ", i))
			
			# stop("Hallo")
			
			#if(nrow(gtex_table)==0){
				
				#gtex_table = input[idx,gtex_spalten]
					
			#}else{
				
				#gtex_table = rbind(gtex_table,input[idx,gtex_spalten], make.row.names = FALSE, stringsAsFactors = FALSE)
					
			#}
				
			
			#region = c(region,rep(cred_var_99$region_id[i],length(idx)))
			
			#signal = c(signal,rep(cred_var_99$signal_id[i],length(idx)))
			
			#ea = c(ea,rep(toupper(cred_var_99$ea[i]),length(idx)))
			
		#}
		
		
	#}
	
	
	input_tissue = unlist(strsplit(gtex_files[y],".",fixed=TRUE))
				
	tissue = rep(input_tissue[1],nrow(gtex_table))
	
	gtex_table['tissue'] <- tissue
	
	gtex_table = merge(gtex_table, cred_var_99, by.x = "rs_id_dbSNP151_GRCh38p7", by.y= "rsid")
	
	output = paste("./lmm/all_loci/gtex_eqtl/cred_var_wakefield_with_gtex_",input_tissue[1],".txt",sep="")
	
	write.table(gtex_table, file=output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
	#write.table(gtex_table, file="./eQTL/gtex/all_cred_var_zscore_with_gtex.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, append = TRUE)
	#print(y)
}
			
			
				
		
			



