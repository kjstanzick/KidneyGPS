library(data.table)
library(shiny)
#install.packages("shinyBS")
library(shinyBS)
#install.packages("DT")
library(DT)
library(shinyjs)
library(shinyWidgets)
library(htmltools)
#library(tippy)
#library(png)

all_sig=fread("Z:/lmm/03_eval/metal_eGFR_meta1.TBL.Indep.500k.5e8.indep_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
indepX=fread("Z:/lmm/03_eval/metal_eGFR_meta1.TBL.Indep.500k.5e8.indepX_edited_MHC-tw.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
cred_var=read.table("Z:/lmm/all_loci/07_cred_var/all_cred_var_99_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
region_table=read.table("Z:/lmm/all_loci/05_region_table/region_table_MHC_cred_var_cys_bun.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
CADD=read.table("Z:/lmm/all_loci/cadd_lmm/results/cadd_in_cred_var_no_dup.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
CADD_wo_filter=read.table("Z:/lmm/all_loci/cadd_lmm/results_no_cutoff/cadd_in_gws_var_ks.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NEPTUNE_glo=fread("Z:/lmm/all_loci/neptune/cred_var_with_nephQTL_glomerular_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
NEPTUNE_tub=fread("Z:/lmm/all_loci/neptune/cred_var_with_nephQTL_tubulointerstitial_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
NEPTUNE_tub=NEPTUNE_tub[which(!(is.na(NEPTUNE_tub$Approved.symbol))),]
susztak_glo=fread("Z:/lmm/all_loci/susztak/cred_var_with_susztak_glomerular_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
susztak_tub=fread("Z:/lmm/all_loci/susztak/cred_var_with_susztak_tubulointerstitial_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)


GTEx_eQTL=fread("Z:/lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup_kidney_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
GTEx_sQTL=fread("Z:/lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup_kidney_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
gtex_without_kidney=fread("Z:/lmm/all_loci/gtex_eqtl/all_cred_var_wakefield_with_gtex_right_gene_name_no_dup_without_kidney_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
sqtl_without_kidney=fread("Z:/lmm/all_loci/gtex_sqtl/all_cred_var_wakefield_with_gtex_sqtl_right_gene_name_no_dup_without_kidney_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

# gtex_without_kidney[,which(lapply(gtex_without_kidney,class)=="numeric")]=lapply(gtex_without_kidney[,which(lapply(gtex_without_kidney,class)=="numeric")],signif,digits=2)
# sqtl_without_kidney[,which(lapply(sqtl_without_kidney,class)=="numeric")]=lapply(sqtl_without_kidney[,which(lapply(sqtl_without_kidney,class)=="numeric")],signif,digits=2)


MGI=read.table( "Z:/lmm/all_loci/stats_and_results/candidate_genes_with_MGI_right_name.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill=TRUE, quote="",comment="")
OMIM=read.table("Z:/lmm/all_loci/stats_and_results/candidate_genes_with_OMIM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GPS=read.table("Z:/lmm/all_loci/06_GPS_table/GPS_signal_based.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE)

DM_decline = read.table("Z:/GPS_App/DM_decline_effects_withStanzickId.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE)
names(CADD_wo_filter)=names(CADD)
CADD=rbind(CADD, CADD_wo_filter)
l=which(duplicated(CADD[,-c(19,20)]))
CADD=CADD[-l,]
names(CADD)[which(names(CADD)=="ppa")]="PPA"


#add susztak eqtls to gps
names(NEPTUNE_glo)[which(names(NEPTUNE_glo)=="ppa.x")]="PPA"
names(NEPTUNE_tub)[which(names(NEPTUNE_tub)=="ppa.x")]="PPA"
names(susztak_glo)[c(1,20,12,13,3,5,6,7,8,9,15)]=names(NEPTUNE_glo)[c(1,16,4,5,14,7,8,9,10,12,20)] 
names(susztak_tub)[c(1,20,12,13,3,5,6,7,8,9,15)]=names(NEPTUNE_tub)[c(1,16,4,5,14,7,8,9,10,12,20)]

glo=rbind(NEPTUNE_glo[,c(1,16,4,5,14,7,8,9,10,12,20,29,30)],susztak_glo[,c(1,20,12,13,3,5,6,7,8,9,15,24,25)])
tub=rbind(NEPTUNE_tub[,c(1,16,4,5,14,7,8,9,10,12,20,29,30)],susztak_tub[,c(1,20,12,13,3,5,6,7,8,9,15,24,25)])

for(i in 1:nrow(GPS)){
  GPS$NEPTUNE_glomerulus[i]=length(unique(glo$rsid[which(glo$Approved.symbol==GPS$Gene[i]&glo$signal_id==GPS$signal_id[i]&glo$region_id==GPS$locus_id[i])]))
  GPS$NEPTUNE_tubulointerstitium[i]=length(unique(tub$rsid[which(tub$Approved.symbol==GPS$Gene[i]&tub$signal_id==GPS$signal_id[i]&tub$region_id==GPS$locus_id[i])]))
  GPS$stop.gained.stop.lost.non.synonymus[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS$Gene[i]&CADD$signal_id==GPS$signal_id[i]&CADD$region_id==GPS$locus_id[i]&CADD$ConsScore>6)]))
  GPS$canonical.splice.noncoding.change.synonymous.splice.site[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS$Gene[i]&CADD$signal_id==GPS$signal_id[i]&CADD$region_id==GPS$locus_id[i]&CADD$ConsScore<=6&CADD$ConsScore>4)]))
  GPS$Cadd15_other[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS$Gene[i]&CADD$signal_id==GPS$signal_id[i]&CADD$region_id==GPS$locus_id[i]&&CADD$ConsScore<=4)]))
  if(i%in%c(100,500,1000,4000,8000,10000,12000,14000))print(i)
  }

#columns DM decline drugability pathways

GPS['DM_effect']=rep("-",nrow(GPS))
GPS['decline_effect']=rep("-",nrow(GPS))
GPS['drugable']=rep("-",nrow(GPS))
GPS['enriched_pathway']=rep("-",nrow(GPS))

DM = DM_decline[grep("DM",DM_decline$Type),]
decline = DM_decline[grep("decline",DM_decline$Type),]

for(i in 1:nrow(DM)){
  a= unlist(strsplit(DM$Stanzick.ID[i],".",fixed=T))
  GPS$DM_effect[which(GPS$locus_id==a[1] & GPS$signal_id==a[2])]=DM$Type[i]
}

for(i in 1:nrow(decline)){
  a= unlist(strsplit(decline$Stanzick.ID[i],".",fixed=T))
  GPS$decline_effect[which(GPS$locus_id==a[1] & GPS$signal_id==a[2])]=decline$Type[i]
}

# create GPS restricted to PPA>0.1
GPS_10=GPS
for(i in 1:nrow(GPS_10)){
  GPS_10$NEPTUNE_glomerulus[i]=length(unique(glo$rsid[which(glo$Approved.symbol==GPS_10$Gene[i]&glo$signal_id==GPS_10$signal_id[i]&glo$region_id==GPS_10$locus_id[i]&glo$PPA>0.1)]))
  GPS_10$NEPTUNE_tubulointerstitium[i]=length(unique(tub$rsid[which(tub$Approved.symbol==GPS_10$Gene[i]&tub$signal_id==GPS_10$signal_id[i]&tub$region_id==GPS_10$locus_id[i]&tub$PPA>0.1)]))
  GPS_10$GTEX_eQTL_kidney_tissue[i]=length(unique(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_eQTL$gene==GPS_10$Gene[i]&GTEx_eQTL$signal_id==GPS_10$signal_id[i]&GTEx_eQTL$region_id==GPS_10$locus_id[i]&GTEx_eQTL$PPA>0.1)]))
  GPS_10$GTEX_eQTL_any_other_tissue[i]=length(unique(gtex_without_kidney_show$RSID[which(gtex_without_kidney_show$`Affected gene`==GPS_10$Gene[i]&gtex_without_kidney_show$signal_id==GPS_10$signal_id[i]&gtex_without_kidney_show$region_id==GPS_10$locus_id[i]&gtex_without_kidney_show$PPA>0.1)]))
  GPS_10$GTEX_sQTL_kidney_tissue[i]=length(unique(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_sQTL$gene==GPS_10$Gene[i]&GTEx_sQTL$signal_id==GPS_10$signal_id[i]&GTEx_sQTL$region_id==GPS_10$locus_id[i]&GTEx_sQTL$PPA>0.1)]))
  GPS_10$GTEX_sQTL_any_other_tissue[i]=length(unique(sqtl_without_kidney_show$RSID[which(sqtl_without_kidney_show$`Affected gene`==GPS_10$Gene[i]&sqtl_without_kidney_show$signal_id==GPS_10$signal_id[i]&sqtl_without_kidney_show$region_id==GPS_10$locus_id[i]&sqtl_without_kidney_show$PPA>0.1)]))
  GPS_10$stop.gained.stop.lost.non.synonymus[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_10$Gene[i]&CADD$signal_id==GPS_10$signal_id[i]&CADD$region_id==GPS_10$locus_id[i]&CADD$PPA>0.1&CADD$ConsScore>6)]))
  GPS_10$canonical.splice.noncoding.change.synonymous.splice.site[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_10$Gene[i]&CADD$signal_id==GPS_10$signal_id[i]&CADD$region_id==GPS_10$locus_id[i]&CADD$PPA>0.1&CADD$ConsScore<=6&CADD$ConsScore>4)]))
  GPS_10$Cadd15_other[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_10$Gene[i]&CADD$signal_id==GPS_10$signal_id[i]&CADD$region_id==GPS_10$locus_id[i]&CADD$PPA>0.1&CADD$ConsScore<=4)]))
  if(i%in%c(100,500,1000,4000,8000,10000,12000,14000))print(i)}

# create GPS restricted to PPA>0.5
GPS_50=GPS
for(i in 1:nrow(GPS_50)){
  GPS_50$NEPTUNE_glomerulus[i]=length(unique(glo$rsid[which(glo$Approved.symbol==GPS_50$Gene[i]&glo$signal_id==GPS_50$signal_id[i]&glo$region_id==GPS_50$locus_id[i]&glo$PPA>0.5)]))
  GPS_50$NEPTUNE_tubulointerstitium[i]=length(unique(tub$rsid[which(tub$Approved.symbol==GPS_50$Gene[i]&tub$signal_id==GPS_50$signal_id[i]&tub$region_id==GPS_50$locus_id[i]&tub$PPA>0.5)]))
  GPS_50$GTEX_eQTL_kidney_tissue[i]=length(unique(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_eQTL$gene==GPS_50$Gene[i]&GTEx_eQTL$signal_id==GPS_50$signal_id[i]&GTEx_eQTL$region_id==GPS_50$locus_id[i]&GTEx_eQTL$PPA>0.5)]))
  GPS_50$GTEX_eQTL_any_other_tissue[i]=length(unique(gtex_without_kidney_show$RSID[which(gtex_without_kidney_show$`Affected gene`==GPS_50$Gene[i]&gtex_without_kidney_show$signal_id==GPS_50$signal_id[i]&gtex_without_kidney_show$region_id==GPS_50$locus_id[i]&gtex_without_kidney_show$PPA>0.5)]))
  GPS_50$GTEX_sQTL_kidney_tissue[i]=length(unique(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_sQTL$gene==GPS_50$Gene[i]&GTEx_sQTL$signal_id==GPS_50$signal_id[i]&GTEx_sQTL$region_id==GPS_50$locus_id[i]&GTEx_sQTL$PPA>0.5)]))
  GPS_50$GTEX_sQTL_any_other_tissue[i]=length(unique(sqtl_without_kidney_show$RSID[which(sqtl_without_kidney_show$`Affected gene`==GPS_50$Gene[i]&sqtl_without_kidney_show$signal_id==GPS_50$signal_id[i]&sqtl_without_kidney_show$region_id==GPS_50$locus_id[i]&sqtl_without_kidney_show$PPA>0.5)]))
  GPS_50$stop.gained.stop.lost.non.synonymus[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_50$Gene[i]&CADD$signal_id==GPS_50$signal_id[i]&CADD$region_id==GPS_50$locus_id[i]&CADD$PPA>0.5&CADD$ConsScore>6)]))
  GPS_50$canonical.splice.noncoding.change.synonymous.splice.site[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_50$Gene[i]&CADD$signal_id==GPS_50$signal_id[i]&CADD$region_id==GPS_50$locus_id[i]&CADD$PPA>0.5&CADD$ConsScore<=6&CADD$ConsScore>4)]))
  GPS_50$Cadd15_other[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_50$Gene[i]&CADD$signal_id==GPS_50$signal_id[i]&CADD$region_id==GPS_50$locus_id[i]&CADD$PPA>0.5&CADD$ConsScore<=4)]))
  if(i%in%c(100,500,1000,4000,8000,10000,12000,14000))print(i)
  }

#format GPS entries:
GPS_show=GPS[,c(1,2,26,3:7,21,23,9:20)]
GPS_show$Gene=paste(c(rep("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=",nrow(GPS_show))),GPS$Gene,c(rep("\" target=\"_blank\">", nrow(GPS_show))),GPS$Gene,c(rep("</a>",nrow(GPS_show))),sep="")
names(GPS_show)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Position gene start","Position gene end", "eGFRcys or BUN validation", "# credible variants in signal","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerulus (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL kidney cortex (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex (GTEx)","sQTL other tissue (GTEx)","# kidney phenotypes in mouse","# kidney phenotypes in human", "Coloc in NEPTUNE tissue")
GPS_show[which(GPS_show[,1]=="[PDILT]"),1]="[UMOD-PDILT]"
for(i in 1:nrow(GPS_show)){
  GPS_show[i,'Score']=11-length(which(GPS_show[i,c(11:21)]==0))
  
}
GPS_show$Score=as.integer(GPS_show$Score)
GPS_show$`eGFRcys or BUN validation`=as.factor(GPS_show$`eGFRcys or BUN validation`)
GPS_show$`Coloc in NEPTUNE tissue`=as.factor(GPS_show$`Coloc in NEPTUNE tissue`)
GPS_show$`Signal ID`=paste("Signal ", GPS_show$`Signal ID`, sep="")
max_ppa=c()
for(i in 1:nrow(GPS_show)){
  u=cred_var[which(cred_var$region_id==GPS$locus_id[i]&cred_var$signal_id==GPS$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 2)
}
GPS_show['max PPA']=max_ppa
GPS_show['Diabetes specific effect']=GPS$DM_effect
GPS_show['decline effect']=GPS$decline_effect
GPS_show['drugable']=GPS$drugable
GPS_show['in enriched pathway']=GPS$enriched_pathway

##same for PPA>0.1
GPS_10_show=GPS_10[,c(1,2,26,3:7,21,23,9:20)]
GPS_10_show$Gene=paste(c(rep("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=",nrow(GPS_10_show))),GPS_10$Gene,c(rep("\" target=\"_blank\">", nrow(GPS_10_show))),GPS_10$Gene,c(rep("</a>",nrow(GPS_10_show))),sep="")
names(GPS_10_show)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Position gene start","Position gene end", "eGFRcys or BUN validation", "# credible variants in signal","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerulus (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL kidney cortex (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex (GTEx)","sQTL other tissue (GTEx)","# kidney phenotypes in mouse","# kidney phenotypes in human", "Coloc in NEPTUNE tissue")
GPS_10_show[which(GPS_10_show[,1]=="[PDILT]"),1]="[UMOD-PDILT]"
for(i in 1:nrow(GPS_10_show)){
  GPS_10_show[i,'Score']=11-length(which(GPS_10_show[i,c(11:21)]==0))
  
}
GPS_10_show$Score=as.integer(GPS_10_show$Score)
GPS_10_show$`eGFRcys or BUN validation`=as.factor(GPS_10_show$`eGFRcys or BUN validation`)
GPS_10_show$`Coloc in NEPTUNE tissue`=as.factor(GPS_10_show$`Coloc in NEPTUNE tissue`)
GPS_10_show$`Signal ID`=paste("Signal ", GPS_10_show$`Signal ID`, sep="")
max_ppa=c()
for(i in 1:nrow(GPS_10_show)){
  u=cred_var[which(cred_var$region_id==GPS_10$locus_id[i]&cred_var$signal_id==GPS_10$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 2)
}
GPS_10_show['max PPA']=max_ppa
GPS_10_show['Diabetes specific effect']=GPS$DM_effect
GPS_10_show['decline effect']=GPS$decline_effect
GPS_10_show['drugable']=GPS$drugable
GPS_10_show['in enriched pathway']=GPS$enriched_pathway

##same for PPA>0.5
GPS_50_show=GPS_50[,c(1,2,26,3:7,21,23,9:20)]
GPS_50_show$Gene=paste(c(rep("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=",nrow(GPS_50_show))),GPS_50$Gene,c(rep("\" target=\"_blank\">", nrow(GPS_50_show))),GPS_50$Gene,c(rep("</a>",nrow(GPS_50_show))),sep="")
names(GPS_50_show)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Position gene start","Position gene end", "eGFRcys or BUN validation", "# credible variants in signal","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerulus (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL kidney cortex (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex (GTEx)","sQTL other tissue (GTEx)","# kidney phenotypes in mouse","# kidney phenotypes in human", "Coloc in NEPTUNE tissue")
GPS_50_show[which(GPS_50_show[,1]=="[PDILT]"),1]="[UMOD-PDILT]"
for(i in 1:nrow(GPS_50_show)){
  GPS_50_show[i,'Score']=11-length(which(GPS_50_show[i,c(11:21)]==0))
  
}
GPS_50_show$Score=as.integer(GPS_50_show$Score)
GPS_50_show$`eGFRcys or BUN validation`=as.factor(GPS_50_show$`eGFRcys or BUN validation`)
GPS_50_show$`Coloc in NEPTUNE tissue`=as.factor(GPS_50_show$`Coloc in NEPTUNE tissue`)
GPS_50_show$`Signal ID`=paste("Signal ", GPS_50_show$`Signal ID`, sep="")
max_ppa=c()
for(i in 1:nrow(GPS_50_show)){
  u=cred_var[which(cred_var$region_id==GPS_50$locus_id[i]&cred_var$signal_id==GPS_50$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 2)
}
GPS_50_show['max PPA']=max_ppa
GPS_50_show['Diabetes specific effect']=GPS$DM_effect
GPS_50_show['decline effect']=GPS$decline_effect
GPS_50_show['drugable']=GPS$drugable
GPS_50_show['in enriched pathway']=GPS$enriched_pathway

#format association tables
all_sig=merge(all_sig,indepX[,c(17,24)], by="Indep.500k.5e8.aLociTag",all.x=T, all.y=F, sort = F )
all_sig=all_sig[,c(2:17,1,18:24)]
all_sig_show=all_sig
all_sig_show$RSID=paste(c(rep("<a href=\"https://www.ncbi.nlm.nih.gov/snp/",nrow(all_sig_show))),all_sig_show$RSID,c(rep("\" target=\"_blank\">", nrow(all_sig_show))),all_sig_show$RSID,c(rep("</a>",nrow(all_sig_show))),sep="")
names(all_sig_show)[which(names(all_sig_show)=="RSID")]="RSID*"
all_sig_show$n=as.integer(all_sig_show$n)
all_sig_show[,which(lapply(all_sig_show,class)=="numeric")]=lapply(all_sig_show[,which(lapply(all_sig_show,class)=="numeric")],signif,digits=2)
cred_var_show=cred_var
cred_var_show$rsid=paste(c(rep("<a href=\"https://www.ncbi.nlm.nih.gov/snp/",nrow(cred_var_show))),cred_var_show$rsid,c(rep("\" target=\"_blank\">", nrow(cred_var_show))),cred_var_show$rsid,c(rep("</a>",nrow(cred_var_show))),sep="")

names(cred_var_show)=c("RSID*","Chr", "Pos","EA**","EAF***", "Effect conditioned","StdErr conditioned", "P-Value conditioned", "N","PPA","cppa","ci95","ci99","Locus ID","Signal ID", "Allele1", "OA","Effect unconditioned", "StdErr unconditioned", "P-Value unconditioned")
cred_var_show=cred_var_show[,c(1,14,15,10,9,4,17,5:8,18:20,2,3)]
cred_var_show[,c(4,8:14)]=lapply(cred_var_show[,c(4,8:14)], signif, digits=2)

#files for locuszoom
lz_files=list.files("Z:/GPS_App/App_test/www")
#lz_files=lz_files[-1]
locus_id_files=c()
for(i in 1:length(lz_files)){
  
  a=unlist(strsplit(lz_files[i],"_",fixed=TRUE))
  locus_id_files[i]=a[2]
}

genes = fread( "Z:/inputs/genes/glist-hg19.tw.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)


real_genes = read.table("Z:/inputs/synonym_table.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

real_genes = real_genes[which(!is.na(real_genes$approved_symbol)),]

real_genes = real_genes[which(!is.na(real_genes$genes)),]

for(i in 1:nrow(genes)){
  
  if(!any(genes$Gene[i]==real_genes$approved_symbol) & any(genes$Gene[i]==real_genes$genes)){
    
    genes$Gene[i]=real_genes$approved_symbol[which(real_genes$genes==genes$Gene[i])[1]]
  }
}


genes_sorted <- genes[order(genes$Chr, genes$Pos1), ]

names(NEPTUNE_glo)[which(names(NEPTUNE_glo)=="ppa.x")]="PPA"
names(NEPTUNE_tub)[which(names(NEPTUNE_tub)=="ppa.x")]="PPA"
names(susztak_glo)[which(names(susztak_glo)=="ppa")]="PPA"
names(susztak_tub)[which(names(susztak_tub)=="ppa")]="PPA"
names(GTEx_eQTL)[which(names(GTEx_eQTL)=="ppa.x")]="PPA"
names(GTEx_sQTL)[which(names(GTEx_sQTL)=="ppa.x")]="PPA"
names(sqtl_without_kidney)[which(names(sqtl_without_kidney)=="ppa.x")]="PPA"
names(gtex_without_kidney)[which(names(gtex_without_kidney)=="ppa.x")]="PPA"
names(cred_var)[which(names(cred_var)=="ppa.x")]="PPA"
names(cred_var)[which(names(cred_var)=="region_id")]="Locus Id"
names(cred_var)[which(names(cred_var)=="signal_id")]="Signal Id"

names(all_sig)=tolower(names(all_sig))



names(all_sig)=c("MarkerName","EA","OA","n","EAF","Effect","StdErr","P.value","Direction","mac","Chr","Pos","P.value.GC","StdErr.GC","RSID","Indep.500k.5e8.aTopHit","Indep.500k.5e8.aLociTag","Indep.500k.5e8.aNumLocusSNPs","Indep.500k.5e8.aLocusCoordinates","Indep.500k.5e8.aLocusSize","Indep.500k.5e8.aLocusAnnot","NearestGene","Nearest Gene Distance","Locus_id")
names(all_sig_show)=c("MarkerName","EA","OA","N","EAF","beta","StdErr","P.value","Direction","mac","Chr","Pos","P.value.GC","StdErr.GC","RSID*","Indep.500k.5e8.aTopHit","Indep.500k.5e8.aLociTag","Indep.500k.5e8.aNumLocusSNPs","Indep.500k.5e8.aLocusCoordinates","Indep.500k.5e8.aLocusSize","Indep.500k.5e8.aLocusAnnot","Nearest gene","Nearest gene distance", "Locus ID")

region_table_show=region_table
names(region_table_show)=c("Locus ID", "Chr", "Start of locus", "End of locus", "Locus size", "# significant SNPs in locus*", "# genes in locus", "Genes in locus", "# credible variants 95", "# credible variants", "Signal ID", "BUN validation of locus", "eGFRcys validation of locus")
region_table_show=region_table_show[,c(1,11,2:6,10,7,8,12,13)]
## replace empty genes in locus FALSE to "-":
region_table_show[which(region_table_show[,9]=="FALSE"),9]="-"


NEPTUNE_glo_show=NEPTUNE_glo[,c(1,16,4,5,14,7,8,9,10,12,20)]
NEPTUNE_tub_show=NEPTUNE_tub[,c(1,16,4,5,14,7,8,9,10,12,20)]
susztak_glo_show=susztak_glo[,c(1,20,12,13,3,5,6,7,8,9,15)]
susztak_tub_show=susztak_tub[,c(1,20,12,13,3,5,6,7,8,9,15)]
GTEx_eQTL_show=GTEx_eQTL[,c(1,19,20,21,17,14,15,9,10,8,23)]
gtex_without_kidney_show=gtex_without_kidney[,c(1,19,20,21,17,18,14,15,9,10,8,23,32,33)]
gtex_without_kidney[,c(2:length(names(gtex_without_kidney)))]<-NULL # save some workspace
GTEx_sQTL_show=GTEx_sQTL[,c(1,20,21,22,18,15,16,10,11,9,24)]
sqtl_without_kidney_show=sqtl_without_kidney[,c(1,20,21,22,18,19,15,16,10,11,9,24,33,34)] 
sqtl_without_kidney[,c(2:length(names(sqtl_without_kidney)))]<-NULL # save some workspace
names(NEPTUNE_glo_show)=c("RSID","PPA","Chr","Pos","Affected gene", "Other allele", "Effect allele*", "Effect on expression", "StdErr", "P-Value","EAF")
names(NEPTUNE_tub_show)=names(NEPTUNE_glo_show)
names(susztak_glo_show)=names(NEPTUNE_glo_show)
names(susztak_tub_show)=names(NEPTUNE_glo_show)
names(GTEx_eQTL_show)=names(NEPTUNE_glo_show)
names(gtex_without_kidney_show)=c("RSID","PPA","Chr","Pos","Affected gene", "Tissue","Other allele", "Effect allele*", "Effect on expression", "StdErr", "P-Value","EAF","region_id","signal_id")
names(GTEx_sQTL_show)=c("RSID","PPA","Chr","Pos","Affected gene", "Other allele", "Effect allele*", "Effect on splicing", "StdErr", "P-Value","EAF")
names(sqtl_without_kidney_show)=c("RSID","PPA","Chr","Pos","Affected gene", "Tissue", "Other allele", "Effect allele*", "Effect on splicing", "StdErr", "P-Value","EAF","region_id","signal_id")
NEPTUNE_glo_show['Effect direction on expression']=ifelse(NEPTUNE_glo_show$`Effect on expression`>0,"upregulating","downregulating")
NEPTUNE_tub_show['Effect direction on expression']=ifelse(NEPTUNE_tub_show$`Effect on expression`>0,"upregulating","downregulating")
susztak_glo_show['Effect direction on expression']=ifelse(susztak_glo_show$`Effect on expression`>0,"upregulating","downregulating")
susztak_tub_show['Effect direction on expression']=ifelse(susztak_tub_show$`Effect on expression`>0,"upregulating","downregulating")
GTEx_eQTL_show['Effect direction on expression']=ifelse(GTEx_eQTL_show$`Effect on expression`>0,"upregulating","downregulating")
# GTEx_sQTL_show['Effect direction']=ifelse(GTEx_sQTL_show$`Effect on splicing`>0,"upregulating","downregulating")
gtex_without_kidney_show['Effect direction on expression']=ifelse(gtex_without_kidney_show$`Effect on expression`>0,"upregulating","downregulating")
# sqtl_without_kidney_show['Effect direction']=ifelse(sqtl_without_kidney_show$`Effect on splicing`>0,"upregulating","downregulating")

NEPTUNE_glo_show=NEPTUNE_glo_show[,c(1,2,5,7,6,11,12,8:10,3,4)]
NEPTUNE_tub_show=NEPTUNE_tub_show[,c(1,2,5,7,6,11,12,8:10,3,4)]
susztak_glo_show=susztak_glo_show[,c(1,2,5,7,6,11,12,8:10,3,4)]
susztak_tub_show=susztak_tub_show[,c(1,2,5,7,6,11,12,8:10,3,4)]
GTEx_eQTL_show=GTEx_eQTL_show[,c(1,2,5,7,6,11,12,8:10,3,4)]
GTEx_sQTL_show=GTEx_sQTL_show[,c(1,2,5,7,6,11,8:10,3,4)]
gtex_without_kidney_show=gtex_without_kidney_show[,c(1,2,5,6,8,7,12,15,9:11,3,4,13,14)]
sqtl_without_kidney_show=sqtl_without_kidney_show[,c(1,2,5,6,8,7,12,9:11,3,4,13,14)]


GTEx_eQTL_show[,which(lapply(GTEx_eQTL_show,class)=="numeric")]=lapply(GTEx_eQTL_show[,which(lapply(GTEx_eQTL_show,class)=="numeric")],signif,digits=2)
GTEx_sQTL_show[,which(lapply(GTEx_sQTL_show,class)=="numeric")]=lapply(GTEx_sQTL_show[,which(lapply(GTEx_sQTL_show,class)=="numeric")],signif,digits=2)
gtex_without_kidney_show[,which(lapply(gtex_without_kidney_show,class)=="numeric")]=lapply(gtex_without_kidney_show[,which(lapply(gtex_without_kidney_show,class)=="numeric")],signif,digits=2)
sqtl_without_kidney_show[,which(lapply(sqtl_without_kidney_show,class)=="numeric")]=lapply(sqtl_without_kidney_show[,which(lapply(sqtl_without_kidney_show,class)=="numeric")],signif,digits=2)
NEPTUNE_glo_show[,which(lapply(NEPTUNE_glo_show,class)=="numeric")]=lapply(NEPTUNE_glo_show[,which(lapply(NEPTUNE_glo_show,class)=="numeric")],signif,digits=2)
NEPTUNE_tub_show[,which(lapply(NEPTUNE_tub_show,class)=="numeric")]=lapply(NEPTUNE_tub_show[,which(lapply(NEPTUNE_tub_show,class)=="numeric")],signif,digits=2)
susztak_glo_show[,which(lapply(susztak_glo_show,class)=="numeric")]=lapply(susztak_glo_show[,which(lapply(susztak_glo_show,class)=="numeric")],signif,digits=2)
susztak_tub_show[,which(lapply(susztak_tub_show,class)=="numeric")]=lapply(susztak_tub_show[,which(lapply(susztak_tub_show,class)=="numeric")],signif,digits=2)

NEPTUNE_glo_show['Source**']="NEPTUNE"
NEPTUNE_tub_show['Source**']="NEPTUNE"
susztak_glo_show['Source**']="Sheng et al, Nat Genet 2021"
susztak_tub_show['Source**']="Sheng et al, Nat Genet 2021"

glo_show=rbind(NEPTUNE_glo_show,susztak_glo_show)
tub_show=rbind(NEPTUNE_tub_show,susztak_tub_show)



MIM_number=read.table("Z:/OMIM/OMIM_gene_pheno/genemap2_edited.txt", header = TRUE,  sep="\t", stringsAsFactors = FALSE, quote="", fill=TRUE)
MIM_number$Gene.Symbols[grep(",", MIM_number$Gene.Symbols)]=toupper(MIM_number$Approved.Symbol[grep(",", MIM_number$Gene.Symbols)])
for(i in 1:nrow(MIM_number)){
  
  if(!any(MIM_number$Gene.Symbols[i]==real_genes$approved_symbol) & any(MIM_number$Gene.Symbols[i]==real_genes$genes)){
    
    MIM_number$Gene.Symbols[i]=real_genes$approved_symbol[which(real_genes$genes==MIM_number$Gene.Symbols[i])[1]]
  }
}

OMIM_MIM=merge(OMIM, MIM_number[,c("MIM.Number", "Gene.Symbols")], by.x="Gene", by.y="Gene.Symbols", all.x=TRUE, all.y=FALSE, sort = FALSE)
OMIM_MIM=OMIM_MIM[-150,]
OMIM_MIM=OMIM_MIM[c(1:222,244,223:243),] #TNXA has no MIM, so it gets sorted to the end. This is the correction
OMIM_MIM$Gene[which(!is.na(OMIM_MIM$MIM.Number))]=paste(c(rep("<a href=\"https://www.omim.org/entry/",nrow(OMIM_MIM)-1)),OMIM_MIM$MIM.Number[which(!is.na(OMIM_MIM$MIM.Number))],c(rep("\" target=\"_blank\">", nrow(OMIM_MIM)-1)),OMIM_MIM$Gene[which(!is.na(OMIM_MIM$MIM.Number))],c(rep("</a>",nrow(OMIM_MIM)-1)),sep="")
names(OMIM_MIM)=c("Gene*","Genetic disorder with kidney phenotype", "Source**", "Locus ID")


MGI_show=MGI[,c(1,3,2)]
names(MGI_show)=c("Gene", "Locus ID", "Phenotype in mouse")

CADD_show=CADD[,c(1,29,15,33,34,10,5,6,16,11,12,17)]
names(CADD_show)=c("RSID", "PPA", "Affected gene", "Locus ID", "Signal ID", "Functional consequence", "REF", "ALT",  "Aminoacid position", "Reference aminoacid", "Alternative aminoacid","CADD Phred-Score")
CADD_show$PPA=signif(CADD_show$PPA, digits=2)

links=scan(file = "Z:/GPS_App/App_test/www/links_plots.txt", what="character")
links_id<-c()
for(i in 1:length(links)){
  z=unlist(strsplit(links[i],"/"))
  links_id[i]=unlist(strsplit(z[6],"_"))[2]
}

save.image(file="Z:/GPS_App/ab_080722/data/workspace.RData")