library(data.table)

setwd("/stk05236/clean_GPS/")

meta = fread( "/stk05236/lmm/all_loci/01_EUR/nea_meta_MHC_plus_european_beta.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

cys = fread( "/stk05236/lmm/GWAMA_cys/02_gwama/metal_eGFRcys_meta1.TBL",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#bun = fread( "./lmm/GWAMA_bun/02_gwama/metal_bun_meta1.TBL",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
bun = fread( "/stk05236/lmm/GWAMA_bun/02_gwama/metal_bun_meta_all1.TBL",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
names(cys) [names(cys)=="P-value"] <- "P.value"
names(bun) [names(bun)=="P-value"] <- "P.value"

p_cys <- vector(length=nrow(meta))

p_bun <- vector(length=nrow(meta))


for (i in 1:nrow(meta)){

	if(any(meta$MarkerName[i]==cys$MarkerName)){
	
		p_cys[i] <- as.numeric(cys$P.value[which(meta$MarkerName[i]==cys$MarkerName)])
	}else{p_cys[i]<-NA}
	
	if(any(meta$MarkerName[i]==bun$MarkerName)){
		p_bun[i] <-as.numeric(bun$P.value[which(meta$MarkerName[i]==bun$MarkerName)])
	}else{p_meta[i] <-NA}
		
	print(i)
}

beta_cys <- vector(length=nrow(meta))
beta_bun <- vector(length=nrow(meta))

for (i in 1:nrow(meta)){
	
	if(any(meta$MarkerName[i]==cys$MarkerName)){
		beta_cys[i] <-as.numeric(cys$Effect[which(meta$MarkerName[i]==cys$MarkerName)])
		if (cys$Allele1[which(meta$MarkerName[i]==cys$MarkerName)] != meta$Allele1[i]){
			beta_cys[i] <- 0-(beta_cys[i])
		}
	}else{beta_cys[i] <-NA}
	
	if(any(meta$MarkerName[i]==bun$MarkerName)){
		beta_bun[i] <-as.numeric(bun$Effect[which(meta$MarkerName[i]==bun$MarkerName)])
		if (bun$Allele1[which(meta$MarkerName[i]==bun$MarkerName)] != meta$Allele1[i]){
			beta_bun[i] <- 0-(beta_bun[i])
		}
	}else{beta_bun[i] <-NA}
	
	print(i)
}

effect_direction_cys <- vector(length=nrow(meta))
effect_direction_bun <- vector(length=nrow(meta))

for (i in 1:nrow(meta)){

	
	if (beta_cys[i]*meta$Effect[i] > 0){
		effect_direction_cys[i] <- "same"
	}
	if (beta_cys[i]*meta$Effect[i] < 0){
		effect_direction_cys[i] <- "opposite"
	}
	if (beta_cys[i]*meta$Effect[i] == 0){
		effect_direction_cys[i] <- "no_effect"
	}
	
	
	
	
	if (beta_bun[i]*meta$Effect[i] > 0){
		effect_direction_bun[i] <- "same"
	}
	if (beta_bun[i]*meta$Effect[i] < 0){
		effect_direction_bun[i] <- "opposite"
	}
	if (beta_bun[i]*meta$Effect[i] == 0){
		effect_direction_bun[i] <- "no_effect"
	}
	
}

meta['p_cys'] <-p_cys
meta['beta_cys'] <-beta_cys
meta['effect_direction_cys'] <-effect_direction_cys
meta['p_bun'] <-p_bun
meta['beta_bun'] <-beta_bun
meta['effect_direction_bun'] <-effect_direction_bun

##nominal Signifikant: 0,05

nominal_significant_cys <-  vector(length=nrow(meta))
nominal_significant_bun <-  vector(length=nrow(meta))
nominal_significant_both <-  vector(length=nrow(meta))

for (i in 1:nrow(meta)){

	if(p_cys[i] < (0.05)){
		nominal_significant_cys[i] <- "yes"
	}else{
		nominal_significant_cys[i] <- "no"
	}
	
	if(p_bun[i] < (0.05)){
		nominal_significant_bun[i] <- "yes"
	}else{
		nominal_significant_bun[i] <- "no"
	}
	if(p_bun[i] < (0.05)){
		if(p_cys[i] < (0.05)){
		nominal_significant_both[i] <- "yes"
		}else{
			nominal_significant_both[i] <- "no"
		}
	}else{
		nominal_significant_both[i] <- "no"
	}
}

meta['nominal_significant_cys'] <-nominal_significant_cys
meta['nominal_significant_bun'] <-nominal_significant_bun
meta['nominal_significant_both'] <-nominal_significant_both



write.table(meta, file="./02_cys_bun/new_meta_MHC_cys_bun_validation.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)