setwd("/stk05236/")

library(data.table)

GPS_table = read.table( "./lmm/all_loci/06_GPS_table/table_for_GPS.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#id = rep("k",nrow(GPS_table))

#id = paste(id,GPS_table$locus_id,sep="")

#GPS_table$locus_id = id

Gene_Score = rep(NA,nrow(GPS_table))

GPS_table['Gene_Score'] = Gene_Score

#
#	Gen-Namen anpassen
#

real_genes = read.table("./inputs/synonym_table.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

real_genes = real_genes[which(!is.na(real_genes$approved_symbol)),]

real_genes = real_genes[which(!is.na(real_genes$genes)),]

for(i in 1:nrow(GPS_table)){

	if(!any(GPS_table$Gene[i]==real_genes$approved_symbol) & any(GPS_table$Gene[i]==real_genes$genes)){
	
		GPS_table$Gene[i]=real_genes$approved_symbol[which(real_genes$genes==GPS_table$Gene[i])[1]]
	}
}



print("CADD")

######### CADD
cadd_15 = read.table( "./lmm/all_loci/cadd_lmm/results/cadd_in_cred_var.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment="")

cred_var_99 = read.table("./lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rsids_dup=cadd_15$RSID[which(duplicated(cadd_15$RSID))]

all_dups = which(cadd_15$RSID%in%rsids_dup)

more_cons = cadd_15[all_dups,]

low_cons_score = c()

for (i in 1:nrow(more_cons)){

	x = more_cons[which(more_cons$RSID[i]==more_cons$RSID),]
	
	if(any(x$ConsScore> more_cons$ConsScore[i])){
	
		low_cons_score[i] = TRUE
	}else{
		low_cons_score[i] = FALSE
	}
}

cadd_15 = cadd_15[-which(cadd_15$rows2%in%more_cons$rows2[which(low_cons_score)]),]

cadd_15 = cadd_15[-which(duplicated(cadd_15[,c(1:8,9:11,14,16,20)])),]

category_1 = c()

category_2 =c()

category_3 = c()

for (i in 1:nrow(GPS_table)){
	
	x = cadd_15[which(cadd_15$GeneName==GPS_table$Gene[i]),]
	
	a=cred_var_99[which(cred_var_99$region_id==GPS_table$locus_id[i]),]
	
	
	b=c(rep(0,max(a$signal_id)))
	
	d=merge(x,a,by.x="RSID",by.y="rsid",all.x=TRUE,all.y=TRUE)
	
	d=d[which(!is.na(d$ConsScore)),]
	
	interesting = d[which(d$ConsScore==8 | d$ConsScore==7),]
	
	medium = d[which(d$ConsScore==6 | d$ConsScore==5),]
	
	low = d[which(d$ConsScore<=4 ),]
	
	if(nrow(interesting)>0){
		for(y in 1:nrow(interesting)){
	
			b[as.numeric(interesting$signal_id[y])]=b[as.numeric(interesting$signal_id[y])]+1
		}
		
		category_1[i]=paste(b,collapse="|")
	}else{	
		category_1[i]=0
	}
	
	b=c(rep(0,max(a$signal_id)))
	
	if(nrow(medium)>0){
		for(y in 1:nrow(medium)){
	
			b[as.numeric(medium$signal_id[y])]=b[as.numeric(medium$signal_id[y])]+1
		}
		category_2[i]=paste(b,collapse="|")
	}else{	
		category_2[i]=0
	}
	
	b=c(rep(0,max(a$signal_id)))
	
	if(nrow(low)>0){
		for(y in 1:nrow(low)){
	
			b[as.numeric(low$signal_id[y])]=b[as.numeric(low$signal_id[y])]+1
		}
		category_3[i]=paste(b,collapse="|")	
	}else{	
	category_3[i]=0
	}
}

GPS_table['stop-gained/stop-lost/non-synonymus'] = category_1

GPS_table['canonical-splice/noncoding-change/synonymous/splice-site']=category_2

GPS_table['Cadd15_other']=category_3


entries15 = length(which(category_1!=0 | category_2!=0 |category_3!=0))

print(paste("CADD15_entries: ", entries15))

cadd_15_neu=merge(cadd_15, cred_var_99, by.x="RSID",by.y="rsid",all.x=TRUE,all.y=TRUE)
	
cadd_15_neu=cadd_15_neu[which(!is.na(d$ConsScore)),]

write.table(cadd_15_neu, file= "./lmm/all_loci/cadd_lmm/results/cadd_in_cred_var_no_dup.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)







#######hier weiter bearbeiten, Variable für eQTLs einführen, als Spalte zu GPS hinzufügen, Statistik für jede Variable als Variable erstellen und als Tabelle speichern



print("Neptune")
###
########## neptune

#glo = fread("./lmm/known_loci/NEPTUNE/cred_var_with_nephQTL_glomerular.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#intest = fread("./lmm/known_loci/NEPTUNE/cred_var_with_nephQTL_tubulointerstitial.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

glo = fread("./lmm/all_loci/neptune/cred_var_with_nephQTL_glomerular.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

intest = fread("./lmm/all_loci/neptune/cred_var_with_nephQTL_tubulointerstitial.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

names(glo)[13] = "Gene"

glo = glo[which(!is.na(glo$Gene)),]

names(intest)[13] = "Gene"

intest = intest[which(!is.na(intest$Gene)),]

for(i in 1:nrow(glo)){

	if(!any(glo$Gene[i]==real_genes$approved_symbol) & any(glo$Gene[i]==real_genes$genes)){
	
		glo$Gene[i]=real_genes$approved_symbol[which(real_genes$genes==glo$Gene[i])[1]]
	}
}
for(i in 1:nrow(intest)){

	if(!any(intest$Gene[i]==real_genes$approved_symbol) & any(intest$Gene[i]==real_genes$genes)){
	
		intest$Gene[i]=real_genes$approved_symbol[which(real_genes$genes==intest$Gene[i])[1]]
	}
}

#write.table(glo, file= "./lmm/all_loci/neptune/cred_var_with_nephQTL_glomerular_right_gene_name.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(intest, file= "./lmm/all_loci/neptune/cred_var_with_nephQTL_tubulointerstitial_right_gene_name.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

NEPTUNE_glomerulus = c()

NEPTUNE_tubulointerstitium = c()

for (i in 1:nrow(GPS_table)){

	if(any(GPS_table$Gene[i] == glo$Gene)) {
	
		a = rep(0,max(cred_var_99$signal_id[which(cred_var_99$region_id==GPS_table$locus_id[i])]))
		
		b=glo[which(glo$Gene==GPS_table$Gene[i]&glo$region==GPS_table$locus_id[i]),]
		
		for(y in 1:nrow(b)){
		
			a[b$signal_id[y]]=a[b$signal_id[y]]+1
		}
	
		NEPTUNE_glomerulus[i] =paste(a,collapse="|")
		
	}else{
	
		NEPTUNE_glomerulus[i] =0
		
	}
	
	if(any(GPS_table$Gene[i] == intest$Gene)) {
	
		a = rep(0,max(cred_var_99$signal_id[which(cred_var_99$region_id==GPS_table$locus_id[i])]))
		
		b=intest[which(intest$Gene==GPS_table$Gene[i] & intest$region_id==GPS_table$locus_id[i]),]
		
		for(y in 1:nrow(b)){
		
			a[b$signal_id[y]]=a[b$signal_id[y]]+1
		}
	
		NEPTUNE_tubulointerstitium[i] =paste(a,collapse="|")
		
	}else{
	
		NEPTUNE_tubulointerstitium[i] =0
		
	}
	
}

GPS_table['NEPTUNE_glomerulus'] = NEPTUNE_glomerulus

GPS_table['NEPTUNE_tubulointerstitium'] = NEPTUNE_tubulointerstitium

entry_glo = length(which(NEPTUNE_glomerulus!=0))

enrty_tub = length(which(NEPTUNE_tubulointerstitium!=0))



print("gtex")

####### gtex

gtex = fread("./lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

sqtl = fread("./lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

###duplicate entfernen:
#rsid,variant_id,gene,tissue,locus_id,signal_id
gtex=gtex[-which(duplicated(gtex[,c(1,3,17,18,31,32)])),]
sqtl1=sqtl[-which(duplicated(sqtl[,c(1,3,18,19,32,33)])),]

for(i in 1:nrow(gtex)){

	if(!any(gtex$gene[i]==real_genes$approved_symbol) & any(gtex$gene[i]==real_genes$genes)){
	
		gtex$gene[i]=real_genes$approved_symbol[which(real_genes$genes==gtex$gene[i])[1]]
	}
}
for(i in 1:nrow(sqtl1)){

	if(!any(sqtl1$gene[i]==real_genes$approved_symbol) & any(sqtl1$gene[i]==real_genes$genes)){
	
		sqtl1$gene[i]=real_genes$approved_symbol[which(real_genes$genes==sqtl1$gene[i])[1]]
	}
}

###duplicate entfernen:
#rsid,variant_id,gene,tissue,locus_id,signal_id
gtex=gtex[-which(duplicated(gtex[,c(1,3,17,18,31,32)])),]
sqtl2=sqtl1[-which(duplicated(sqtl[,c(1,3,18,19,32,33)])),]

#write.table(gtex, file= "./lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(sqtl, file= "./lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

gtex_kidney = gtex[which(gtex$tissue == "Kidney_Cortex"),]

gtex_without_kidney = gtex[which(gtex$tissue != "Kidney_Cortex"),]

sqtl_kidney = sqtl2[which(sqtl2$tissue == "Kidney_Cortex"),]

sqtl_without_kidney = sqtl2[-which(sqtl2$tissue == "Kidney_Cortex"),]

write.table(gtex_kidney, file= "./lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup_kidney.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(sqtl_kidney, file= "./lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup_kidney.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(gtex_without_kidney, file= "./lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup_without_kidney.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(sqtl_without_kidney, file= "./lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup_without_kidney.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#gtex_kidney = fread("./lmm/gtex/cred_var_gtex_kidney_eqtl.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
#gtex_without_kidney = fread("./tw/ckd-supptable15_gtex_other_eqtl_210301_nodup.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#sqtl_kidney = fread("./lmm/gtex/cred_var_gtex_kidney_sqtl.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

#sqtl_without_kidney = fread("./tw/ckd-supptable16_gtex_other_sqtl_210301_nodup.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)



GTEX_eQTL_kidney_tissue = c()

GTEX_eQTL_any_other_tissue = c()

GTEX_sQTL_kidney_tissue = c()

GTEX_sQTL_any_other_tissue = c()

for (i in 1:nrow(GPS_table)){

	if(any(GPS_table$Gene[i] == gtex_kidney$gene)) {
	
		a = rep(0,max(cred_var_99$signal_id[which(cred_var_99$region_id==GPS_table$locus_id[i])]))
		
		b = gtex_kidney[which(gtex_kidney$gene==GPS_table$Gene[i] & gtex_kidney$region_id==GPS_table$locus_id[i]),]
		
		for(y in 1:nrow(b)){
		
			a[b$signal_id[y]]=a[b$signal_id[y]]+1
		}
	
		GTEX_eQTL_kidney_tissue[i] =paste(a,collapse="|")
		
	}else{
	
		GTEX_eQTL_kidney_tissue[i] =0
		
	}
	
	if(any(GPS_table$Gene[i] == gtex_without_kidney$gene)) {
	
		a = rep(0,max(cred_var_99$signal_id[which(cred_var_99$region_id==GPS_table$locus_id[i])]))
		
		b = gtex_without_kidney[which(gtex_without_kidney$gene==GPS_table$Gene[i] & gtex_without_kidney$region_id==GPS_table$locus_id[i]),]
		
		b = b[which(!duplicated(b$rs_id_dbSNP151_GRCh38p7)),]
		
		for(y in 1:nrow(b)){
		
			a[b$signal_id[y]]=a[b$signal_id[y]]+1
		}
	
		GTEX_eQTL_any_other_tissue[i] =paste(a,collapse="|")
		
	}else{
	
		GTEX_eQTL_any_other_tissue[i] =0
		
	}
	
	if(any(GPS_table$Gene[i] == sqtl_kidney$Gene)) {
	
		a = rep(0,max(cred_var_99$signal_id[which(cred_var_99$region_id==GPS_table$locus_id[i])]))
		
		b=sqtl_kidney[which(sqtl_kidney$Gene==GPS_table$Gene[i] & sqtl_kidney$region_id==GPS_table$locus_id[i]),]
		
		for(y in 1:nrow(b)){
		
			a[b$signal_id[y]]=a[b$signal_id[y]]+1
		}
	
		GTEX_sQTL_kidney_tissue[i] =paste(a,collapse="|")
		
	}else{
	
		GTEX_sQTL_kidney_tissue[i] =0
		
	}
	
	if(any(GPS_table$Gene[i] == sqtl_without_kidney$gene)) {
	
		a = rep(0,max(cred_var_99$signal_id[which(cred_var_99$region_id==GPS_table$locus_id[i])]))
		
		b=sqtl_without_kidney[which(sqtl_without_kidney$gene==GPS_table$Gene[i] &sqtl_without_kidney$region_id==GPS_table$locus_id[i]),]
		
		b = b[which(!duplicated(b$rs_id_dbSNP151_GRCh38p7)),]
		
		for(y in 1:nrow(b)){
		
			a[b$signal_id[y]]=a[b$signal_id[y]]+1
		}
	
		GTEX_sQTL_any_other_tissue[i] =paste(a,collapse="|")
		
	}else{
	
		GTEX_sQTL_any_other_tissue[i] =0
		
	}
	
}

GPS_table['GTEX_eQTL_kidney_tissue'] <- GTEX_eQTL_kidney_tissue

GPS_table['GTEX_eQTL_any_other_tissue'] <- GTEX_eQTL_any_other_tissue

GPS_table['GTEX_sQTL_kidney_tissue'] <- GTEX_sQTL_kidney_tissue

GPS_table['GTEX_sQTL_any_other_tissue'] <- GTEX_sQTL_any_other_tissue

entries_eQTL_kidney = length(which(GTEX_eQTL_kidney_tissue!=0))

entries_eQTL_other = length(which(GTEX_eQTL_any_other_tissue!=0))

entries_sQTL_kidney = length(which(GTEX_sQTL_kidney_tissue!=0))

entries_sQTL_other = length(which(GTEX_sQTL_any_other_tissue!=0))



print("MGI")
######## MGI

kidney_genes = read.table( "./MGI_lookup/genes_and_pheno/MGI_table_kidney.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill=TRUE, quote="",comment="")

kidney_genes=kidney_genes[which(!is.na(kidney_genes$human_symbol)),]

for(i in 1:nrow(kidney_genes)){

	if(!any(kidney_genes$human_symbol[i]==real_genes$approved_symbol) & any(kidney_genes$human_symbol[i]==real_genes$genes)){
	
		kidney_genes$human_symbol[i]=real_genes$approved_symbol[which(real_genes$genes==kidney_genes$human_symbol[i])[1]]
	}
}

#write.table(kidney_genes, file= "./MGI_lookup/genes_and_pheno/MGI_table_kidney_right_human_gene_name.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


MGI = c()

for (i in 1:nrow(GPS_table)){

	if(any(GPS_table$Gene[i]== kidney_genes$human_symbol)){
	
		MGI[i]=length(which(GPS_table$Gene[i]== kidney_genes$human_symbol))
	}else{
		MGI[i]=0
	}
}

GPS_table['MGI_mouse_kidney_phenotyp'] <- MGI

entries_MGI = length(which(MGI!=0))

print(paste("MGI_entries: ", entries_MGI))

##MGI_results

MGI_known_loci = merge(kidney_genes[,c(5,6)],GPS_table[,c(2,3,5,6,7)], by.x="human_symbol", by.y="Gene", all.x=FALSE, all.y=FALSE)

write.table(MGI_known_loci, file="./lmm/all_loci/stats_and_results/candidate_genes_with_MGI_right_name.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

print("OMIM")

########## OMIM

kidney_pheno = fread("./OMIM/OMIM_Groopman_combined/all_kidney_genes_sorted_nearly_without_duplicates_ks.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

for(i in 1:nrow(kidney_pheno)){

	if(!any(kidney_pheno$gene[i]==real_genes$approved_symbol) & any(kidney_pheno$gene[i]==real_genes$genes)){
	
		kidney_pheno$gene[i]=real_genes$approved_symbol[which(real_genes$genes==kidney_pheno$gene[i])[1]]
	}
}

## remove duplcates:

kidney_pheno=kidney_pheno[-which(duplicated(kidney_pheno)),]

a = which(duplicated(kidney_pheno[,c(1,2)]))

b= which(kidney_pheno$gene%in%kidney_pheno$gene[a])

kidney_pheno$annot[b]="OMIM; Groopman et al"

#write.table(kidney_pheno, file= "./OMIM/OMIM_Groopman_combined/all_kidney_genes_sorted_nearly_without_duplicates_ks_right_gene_name.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

OMIM =c()

for (i in 1:nrow(GPS_table)){

	if(any(GPS_table$Gene[i]==kidney_pheno$gene)){
	
		OMIM[i]=length(which(GPS_table$Gene[i]==kidney_pheno$gene))
	}else{
	
		OMIM[i] = 0
	}
}

GPS_table['OMIM_human_kidney_phenotype'] <- OMIM

entries_omim=length(which(OMIM!=0))

## OMIM_results

kidney_genes = unique(kidney_pheno$gene)

for(i in 1:length(kidney_genes)){

	disease = paste(kidney_pheno$disease[which(kidney_pheno$gene==kidney_genes[i])],collapse=" | ")
	
	if(length(which(kidney_pheno$gene==kidney_genes[i]))==1){

		annot = kidney_pheno$annot[which(kidney_pheno$gene==kidney_genes[i])]
	}else{
	
		if(any(kidney_pheno$annot[which(kidney_pheno$gene==kidney_genes[i])]=="OMIM; Groopman et al") | (any(kidney_pheno$annot[which(kidney_pheno$gene==kidney_genes[i])]=="OMIM") & any(kidney_pheno$annot[which(kidney_pheno$gene==kidney_genes[i])]=="Groopman et al"))){
		
			annot = "OMIM; Groopman et al"
		}else{
		
			annot=kidney_pheno$annot[which(kidney_pheno$gene==kidney_genes[i])[1]]
		}
	}
	
	
	x = data.frame(kidney_genes[i],disease,annot,stringsAsFactors=FALSE)
	
	if(i==1){OMIM_genes = x}
	
	else{OMIM_genes= rbind(OMIM_genes,x, stringsAsFactors=FALSE, make.row.names=FALSE)}
	
}

idx = (OMIM_genes[,1] %in% GPS_table$Gene)

OMIM_results = OMIM_genes[idx,]

names(OMIM_results)=c("Gene","disease","annot")
OMIM_results = merge(OMIM_results,GPS_table[,c(2,3,5,6,7)], by="Gene", all.x=FALSE, all.y=FALSE)

write.table(OMIM_results, file="./lmm/all_loci/stats_and_results/candidate_genes_with_OMIM.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


####### Coloc

coloc_tub=read.table( "./lmm/coloc/coloc_tubulo.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

coloc_glo=read.table( "./lmm/coloc/coloc_glomerulus.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

coloc_tub=coloc_tub[which(coloc_tub$PP_H4 >=0.80),]
coloc_glo=coloc_glo[which(coloc_glo$PP_H4 >=0.80),]
#
#	Gen-Namen anpassen
#

##glom
for(i in 1:nrow(coloc_glo)){

	if(!any(coloc_glo$gene[i]==real_genes$approved_symbol) & any(coloc_glo$gene[i]==real_genes$genes)){
	
		coloc_glo$gene[i]=real_genes$approved_symbol[which(real_genes$genes==coloc_glo$gene[i])]
	}
}
###tub
for(i in 1:nrow(coloc_tub)){

	if(!any(coloc_tub$gene[i]==real_genes$approved_symbol) & any(coloc_tub$gene[i]==real_genes$genes)){
	
		coloc_tub$gene[i]=real_genes$approved_symbol[which(real_genes$genes==coloc_tub$gene[i])]
	}
}


coloc=c()

for (i in 1:nrow(GPS_table)){

	coloc[i]="-"

	if(any(GPS_table$Gene[i]==coloc_tub$gene) & any(GPS_table$Gene[i]==coloc_glo$gene)){
	
		coloc[i]="tubulo-interstitial/glomerulus"
	}else{
	
		if(any(GPS_table$Gene[i]==coloc_tub$gene)) coloc[i]="tubulo-interstitial"
		
		if(any(GPS_table$Gene[i]==coloc_glo$gene)) coloc[i]="glomerulus"
	}
}

GPS_table['Coloc']=coloc


####### Cys/bun

region_table = read.table( "./lmm/all_loci/05_region_table/region_table_MHC_cred_var_cys_bun.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

region_table = region_table[which(region_table$signal_id==1),]

cys_bun = c()

for (i in 1:nrow(GPS_table)){

	x=region_table[which(region_table$locus_id==GPS_table$locus_id[i]),]
	
	if(x$bun_association =="yes" & x$cys_association=="yes") cys_bun[i]="cys/bun"
	
	if(x$bun_association =="yes" & x$cys_association=="no") cys_bun[i]="bun"
	
	if(x$bun_association =="no" & x$cys_association=="yes") cys_bun[i]="cys"
	
	if(x$bun_association =="no" & x$cys_association=="no") cys_bun[i]="-"
	
}

GPS_table['Cys/BUN']=cys_bun

####### credible variants

region_table = read.table( "./lmm/all_loci/05_region_table/region_table_cred_var.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cred_var_99 = read.table("./lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)



##############
#
#  Credible Set Size
#
#
##############
no_cred_var_99 = c()

no_cred_var_99_signals = c()

min_cred_var_99=c()

for (i in 1:nrow(GPS_table)){

	x=cred_var_99[which(cred_var_99$region_id==GPS_table$locus_id[i]),]
	
	a=rep(0,max(x$signal_id))
	
	for(y in 1:length(a)){
	
		a[y]=length(which(x$signal_id==y))
	}
	
	min_cred_var_99[i] = min(a)
	
	no_cred_var_99_signals[i]=paste(a,collapse="|")	
	
	x = x[which(!duplicated(x$rsid)),]
	
	no_cred_var_99[i]=nrow(x)
}

GPS_table['credible_variants_in_locus'] <- no_cred_var_99

GPS_table['credible_variants_per_locus_signal']=no_cred_var_99_signals

GPS_table['num credible variants smallest credible set'] = min_cred_var_99


################
#
#
#	Credible variants in Gene
#
#
#################

cred_var_gene = c()

for (i in 1:nrow(GPS_table)){

	x = cred_var_99[which(cred_var_99$chr==GPS_table[i,5] & cred_var_99$pos<=GPS_table[i,7] & cred_var_99$pos>=GPS_table[i,6]),]
	
	cred_var_gene[i] = nrow(x)
	
}

GPS_table['credible variants in gene'] = cred_var_gene

##################
#
#    format locus name
#
##################

for (i in 1:nrow(GPS_table)){

	GPS_table[i,1]= paste("[",GPS_table[i,1],"]",sep="")
}


#write.table(GPS_table, file="./lmm/known_loci/06_GPS_table/GPS_all_at_once_separate_signals.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(GPS_table, file="./lmm/all_loci/06_GPS_table/GPS_all_at_once_separate_signals_no_dup.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


stats=rbind(entries15,entry_glo,enrty_tub,entries_eQTL_kidney,entries_eQTL_other,entries_sQTL_kidney,entries_sQTL_other,entries_MGI,entries_omim)

#write.table(stats,file="./lmm/known_loci/stats_and_results/stats_GPS_entries.txt", sep="\t" , col.names=FALSE, row.names =TRUE, quote=FALSE)
write.table(stats,file="./lmm/all_loci/stats_and_results/stats_GPS_entries_no_dup.txt", sep="\t" , col.names=FALSE, row.names =TRUE, quote=FALSE)
