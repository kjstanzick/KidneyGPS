library("data.table")
all_sig=fread("/stk05236/lmm/03_eval/metal_eGFR_meta1.TBL.Indep.500k.5e8.indep.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
cred_var=read.table("/stk05236/lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
library(data.table)
region_table=read.table("/stk05236/lmm/all_loci/05_region_table/region_table_MHC_cred_var_cys_bun.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
CADD=read.table("/stk05236/lmm/all_loci/cadd_lmm/results/cadd_in_cred_var_no_dup.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NEPTUNE_glo=fread("/stk05236/lmm/all_loci/neptune/cred_var_with_nephQTL_glomerular.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
NEPTUNE_tub=fread("/stk05236/lmm/all_loci/neptune/cred_var_with_nephQTL_tubulointerstitial.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
NEPTUNE_tub=NEPTUNE_tub[which(!(is.na(NEPTUNE_tub$Approved.symbol))),]
GTEx_eQTL=fread("/stk05236/lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup_kidney.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
GTEx_sQTL=fread("/stk05236/lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup_kidney.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
MGI=read.table( "/stk05236/lmm/all_loci/stats_and_results/candidate_genes_with_MGI_right_name.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill=TRUE, quote="",comment="")
OMIM=read.table("/stk05236/lmm/all_loci/stats_and_results/candidate_genes_with_OMIM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GPS=read.table("/stk05236/lmm/all_loci/06_GPS_table/GPS_signal_based.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE)

ref=fread("/stk05236/lmm/03_eval/european/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

cred_var_m=merge(cred_var, ref[,c(2,3,6,7,8,15)], by.x="rsid", by.y="RSID", all.x=TRUE, all.y=F, sort=FALSE)
cred_var_m$ea=toupper(cred_var_m$ea)
cred_var_m$Allele1=toupper(cred_var_m$Allele1)
cred_var_m$Allele2=toupper(cred_var_m$Allele2)
#ea und beta_cond auf eGFR lowering drehen
not_lowering=which(cred_var_m$beta>0)

for(i in not_lowering){
	cred_var_m$beta[i]=(-1)*cred_var_m$beta[i]
	cred_var_m$eaf[i]=1-cred_var_m$eaf[i]
	temp_allele=c(cred_var_m$Allele1[i],cred_var_m$Allele2[i])
	cred_var_m$ea[i]=temp_allele[which(temp_allele!=cred_var_m$ea[i])]
}

#Allele1, Allele2 und Effect uncond auf ea anpassen

for(i in 1:nrow(cred_var_m)){
	if(cred_var_m$Allele1[i]!=cred_var_m$ea[i]){
		cred_var_m$Effect[i]=(-1)*cred_var_m$Effect[i]  
		cred_var_m$Allele1[i] <- paste(cred_var_m$Allele1[i],cred_var_m$Allele2[i], sep = "")
		cred_var_m$Allele2[i] <- substr(cred_var_m$Allele1[i],0,nchar(cred_var_m$Allele1[i]) - nchar(cred_var_m$Allele2[i]))
		cred_var_m$Allele1[i] <- substr(cred_var_m$Allele1[i],nchar(cred_var_m$Allele2[i]) + 1, nchar(cred_var_m$Allele1[i]))
	}
}

NEPTUNE_glo_m=merge(NEPTUNE_glo[,c(1:15,24)],cred_var_m,by="rsid",all.x=T, all.y=F, sort=FALSE)
NEPTUNE_glo_m$ALT=toupper(NEPTUNE_glo_m$ALT)
NEPTUNE_glo_m$REF=toupper(NEPTUNE_glo_m$REF)
for (i in 1:nrow(NEPTUNE_glo_m)){
	if(NEPTUNE_glo_m$ALT[i]!=NEPTUNE_glo_m$ea[i]){
		NEPTUNE_glo_m$BETA[i] = (-1)*NEPTUNE_glo_m$BETA[i]
		NEPTUNE_glo_m$ALT[i] <- paste(NEPTUNE_glo_m$ALT[i],NEPTUNE_glo_m$REF[i], sep = "")
		NEPTUNE_glo_m$REF[i] <- substr(NEPTUNE_glo_m$ALT[i],0,nchar(NEPTUNE_glo_m$ALT[i]) - nchar(NEPTUNE_glo_m$REF[i]))
		NEPTUNE_glo_m$ALT[i] <- substr(NEPTUNE_glo_m$ALT[i],nchar(NEPTUNE_glo_m$REF[i]) + 1, nchar(NEPTUNE_glo_m$ALT[i]))
	}
}

NEPTUNE_tub_m=merge(NEPTUNE_tub[,c(1:15,24)],cred_var_m,by="rsid",all.x=T, all.y=F, sort=FALSE)
NEPTUNE_tub_m$ALT=toupper(NEPTUNE_tub_m$ALT)
NEPTUNE_tub_m$REF=toupper(NEPTUNE_tub_m$REF)
for (i in 1:nrow(NEPTUNE_tub_m)){
	if(NEPTUNE_tub_m$ALT[i]!=NEPTUNE_tub_m$ea[i]){
		NEPTUNE_tub_m$BETA[i] = (-1)*NEPTUNE_tub_m$BETA[i]
		NEPTUNE_tub_m$ALT[i] <- paste(NEPTUNE_tub_m$ALT[i],NEPTUNE_tub_m$REF[i], sep = "")
		NEPTUNE_tub_m$REF[i] <- substr(NEPTUNE_tub_m$ALT[i],0,nchar(NEPTUNE_tub_m$ALT[i]) - nchar(NEPTUNE_tub_m$REF[i]))
		NEPTUNE_tub_m$ALT[i] <- substr(NEPTUNE_tub_m$ALT[i],nchar(NEPTUNE_tub_m$REF[i]) + 1, nchar(NEPTUNE_tub_m$ALT[i]))
	}
}

GTEx_eQTL_m=merge(GTEx_eQTL[,c(1:18,27)],cred_var_m,by.x="rs_id_dbSNP151_GRCh38p7",by.y="rsid",all.x=T, all.y=F, sort=FALSE)
GTEx_eQTL_m$alt=toupper(GTEx_eQTL_m$alt)
GTEx_eQTL_m$ref=toupper(GTEx_eQTL_m$ref)

for (i in 1:nrow(GTEx_eQTL_m)){
	if(GTEx_eQTL_m$alt[i]!=GTEx_eQTL_m$ea[i]){
		GTEx_eQTL_m$slope[i] = (-1)*GTEx_eQTL_m$slope[i]
		GTEx_eQTL_m$alt[i] <- paste(GTEx_eQTL_m$alt[i],GTEx_eQTL_m$ref[i], sep = "")
		GTEx_eQTL_m$ref[i] <- substr(GTEx_eQTL_m$alt[i],0,nchar(GTEx_eQTL_m$alt[i]) - nchar(GTEx_eQTL_m$ref[i]))
		GTEx_eQTL_m$alt[i] <- substr(GTEx_eQTL_m$alt[i],nchar(GTEx_eQTL_m$ref[i]) + 1, nchar(GTEx_eQTL_m$alt[i]))
	}
}

GTEx_sQTL_m=merge(GTEx_sQTL[,c(1:19,28)],cred_var_m,by.x="rs_id_dbSNP151_GRCh38p7",by.y="rsid",all.x=T, all.y=F, sort=FALSE)
GTEx_sQTL_m$alt=toupper(GTEx_sQTL_m$alt)
GTEx_sQTL_m$ref=toupper(GTEx_sQTL_m$ref)
for (i in 1:nrow(GTEx_sQTL_m)){
	if(GTEx_sQTL_m$alt[i]!=GTEx_sQTL_m$ea[i]){
		GTEx_sQTL_m$slope[i] = (-1)*GTEx_sQTL_m$slope[i]
		GTEx_sQTL_m$alt[i] <- paste(GTEx_sQTL_m$alt[i],GTEx_sQTL_m$ref[i], sep = "")
		GTEx_sQTL_m$ref[i] <- substr(GTEx_sQTL_m$alt[i],0,nchar(GTEx_sQTL_m$alt[i]) - nchar(GTEx_sQTL_m$ref[i]))
		GTEx_sQTL_m$alt[i] <- substr(GTEx_sQTL_m$alt[i],nchar(GTEx_sQTL_m$ref[i]) + 1, nchar(GTEx_sQTL_m$alt[i]))
	}
}

not_lowering_all=which(all_sig$Effect>0)
all_sig$Allele1=toupper(all_sig$Allele1)
all_sig$Allele2=toupper(all_sig$Allele2)
for(i in not_lowering_all){
	all_sig$Effect[i]=(-1)*all_sig$Effect[i]
	all_sig$Freq1[i]=1-all_sig$Freq1[i]
	temp_allele=c(all_sig$Allele1[i],all_sig$Allele2[i])
	all_sig$Allele1[i]=temp_allele[which(temp_allele!=all_sig$Allele1[i])]
	all_sig$Allele2[i]=temp_allele[which(temp_allele!=all_sig$Allele2[i])]
}

write.table(cred_var_m, file="/stk05236/lmm/all_loci/07_cred_var/all_cred_var_99_egfr_lowering.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(NEPTUNE_glo_m, file="/stk05236/lmm/all_loci/neptune/cred_var_with_nephQTL_glomerular_egfr_lowering.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(NEPTUNE_tub_m, file="/stk05236/lmm/all_loci/neptune/cred_var_with_nephQTL_tubulointerstitial_egfr_lowering.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(GTEx_eQTL_m, file="/stk05236/lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup_kidney_egfr_lowering.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(GTEx_sQTL_m, file="/stk05236/lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup_kidney_egfr_lowering.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(all_sig, file="/stk05236/lmm/03_eval/metal_eGFR_meta1.TBL.Indep.500k.5e8.indep_egfr_lowering.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)