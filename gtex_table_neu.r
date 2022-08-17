setwd("/stk05236/clean_GPS")

cred_var_99 = read.table("./05_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gtex_files = list.files("/stk05236/inputs/gtex")

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

dir.create("./08_eQTLs/gtex")
#gtex_spalten = c(6,7,8,9,13,14,15,16)

gtex_table = data.frame()

for (y in 1:length(gtex_files)){

	input_file = paste("/stk05236/inputs/gtex/",gtex_files[y],sep="")
		
	input = fread(input_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
		
	input = input[-which(is.na(input$rs_id_dbSNP151_GRCh38p7)),]
	
	idx = which (input$rs_id_dbSNP151_GRCh38p7 %in% cred_var_99$rsid)
	
	gtex_table = input[idx,]	
	
	input_tissue = unlist(strsplit(gtex_files[y],".",fixed=TRUE))
				
	tissue = rep(input_tissue[1],nrow(gtex_table))
	
	gtex_table['tissue'] <- tissue
	
	gtex_table = merge(gtex_table, cred_var_99, by.x = "rs_id_dbSNP151_GRCh38p7", by.y= "rsid")
	
	output = paste("./08_eQTLs/gtex_/cred_var_with_gtex_",input_tissue[1],".txt",sep="")
	
	write.table(gtex_table, file=output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
}
			
gtex_files = list.files("./08_eQTLs/gtex")

library(data.table)

for (i in 1:length(gtex_files)){

	input_file = paste("./08_eQTLs/gtex/",gtex_files[i],sep="")
	
	input = fread(input_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
	
	if (i==1) { gtex_table = input}
	
	else{gtex_table = rbind(gtex_table, input, make.row.names = FALSE, stringsAsFactors = FALSE)}
	
	
}

write.table(gtex_table, file="./08_eQTLs/gtex/all_cred_var_with_gtex.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)		
				
		
			



