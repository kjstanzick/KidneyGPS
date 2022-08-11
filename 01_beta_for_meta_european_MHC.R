library(data.table)

## path to main directory
setwd("/stk05236/clean_GPS/")

## eGFRcrea lead variants
meta = fread( "/stk05236/lmm/03_eval/metal_eGFR_meta1.TBL.Indep.500k.5e8.indepX_edited_MHC-tw.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

## eGFRcrea in Europeans
meta_european = fread( "/stk05236/lmm/02_gwama/european/metal_eGFR_meta_ea1.TBL",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)



meta_new <- meta

beta_meta <- meta_new$Effect

beta_meta_european <- vector(length=nrow(meta_new))

for (i in 1:nrow(meta_new)){
	
	if(any(meta_new$MarkerName[i]==meta_european$MarkerName)){
		beta_meta_european[i] <-meta_european$Effect[which(meta_new$MarkerName[i]==meta_european$MarkerName)]
		if (meta_european$Allele1[which(meta_new$MarkerName[i]==meta_european$MarkerName)] != meta_new$Allele1[i]){
			beta_meta_european[i] <- 0-(beta_meta_european[i])
		}
	}else{beta_meta_european[i] <-NA}
	
	print(i)
}

meta_new['beta_european'] <- beta_meta_european


write.table(meta_new, file="./01_EUR/nea_meta_MHC_plus_european_beta.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)