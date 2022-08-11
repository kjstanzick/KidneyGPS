setwd("/stk05236/")

library(data.table)

variant_99 = read.table("./lmm/new_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#rsid [1]

gtex = fread("./lmm/new_loci/gtex/all_cred_var_wakefield_with_gtex.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#rsid [1]

sqtl = fread("./lmm/new_loci/gtex/all_cred_var_wakefield_with_gtex_sqtl.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#rsid[1]

#glo = fread("./eQTL/neptune/cred_var_with_nephQTL_glomerular.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#rsid[15]

#intest = fread("./eQTL/neptune/cred_var_with_nephQTL_tubulointerstitial.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#rsid[15]

#effect = fread("./variant_effect_predictor/results_99_2020_27_02/variants_with_effect.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#rsid[1]

gtex_kidney = gtex[which(gtex$tissue == "Kidney_Cortex"),]

gtex_without_kidney = gtex[which(gtex$tissue != "Kidney_Cortex"),]

sqtl_kidney = sqtl[which(sqtl$tissue == "Kidney_Cortex"),]

sqtl_without_kidney = sqtl[-which(sqtl$tissue == "Kidney_Cortex"),]

gtex_kidney = gtex_kidney[,c(1,8,9,10,14,15,17,18,19,20,21,23,31,32)]
	
gtex_without_kidney = gtex_without_kidney[,c(1,8,9,10,14,15,17,18,19,20,21,23,31,32)]	

sqtl_kidney = sqtl_kidney[,c(1,9,10,11,15,16,18,19,20,21,22,24,32,33)]	

sqtl_without_kidney = sqtl_without_kidney[,c(1,9,10,11,15,16,18,19,20,21,22,24,32,33)]

gtex_kidney$ea = toupper(gtex_kidney$ea)

gtex_without_kidney$ea = toupper(gtex_without_kidney$ea)

sqtl_kidney$ea = toupper(sqtl_kidney$ea)

sqtl_without_kidney$ea = toupper(sqtl_without_kidney$ea)

write.table(gtex_kidney, file="./lmm/new_loci/gtex/gtex_kidney_eqtl.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(gtex_without_kidney, file="./lmm/new_loci/gtex/gtex_without_kidney_eqtl.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(sqtl_kidney, file="./lmm/new_loci/gtex/gtex_kidney_sqtl.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(sqtl_without_kidney, file="./lmm/new_loci/gtex/gtex_without_kidney_sqtl.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
