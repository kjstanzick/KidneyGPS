setwd("/stk05236/")

gtex_files = list.files("./lmm/all_loci/gtex_sqtl")

library(data.table)

for (i in 1:length(gtex_files)){

	input_file = paste("./lmm/all_loci/gtex_sqtl/",gtex_files[i],sep="")
	
	input = fread(input_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
	
	if (i==1) { gtex_table = input}
	
	else{gtex_table = rbind(gtex_table, input, make.row.names = FALSE, stringsAsFactors = FALSE)}
	
	
}

write.table(gtex_table, file="./lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)