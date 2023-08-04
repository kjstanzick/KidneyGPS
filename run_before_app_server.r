library(data.table)
#library(shiny)
#install.packages("shinyBS")
#library(shinyBS)
#install.packages("DT")
#library(DT)
#library(shinyjs)
#library(shinyWidgets)
#library(htmltools)
#library(tippy)
#library(png)
rm(list = ls(all = TRUE))

all_sig=fread("./metal_eGFR_meta1.TBL.Indep.500k.5e8.indep_egfr_lowering.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
indepX=fread("./metal_eGFR_meta1.TBL.Indep.500k.5e8.indepX_edited_MHC-tw.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
cred_var=read.table("./05_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
locus_table=read.table("./04_locus_table/locus_table_cys_bun.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
signal_table=read.table("./06_signal_table/signal_table_cred_var.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
CADD=read.table("./10_protein_altering/results_2023-05-22/cadd_in_cred_var.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NEPTUNE_glo=fread("./08_eQTLs/neptune/cred_var_with_nephQTL_glomerular.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
NEPTUNE_tub=fread("./08_eQTLs/neptune/cred_var_with_nephQTL_tubulointerstitial.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
NEPTUNE_tub=NEPTUNE_tub[which(!(is.na(NEPTUNE_tub$Approved.symbol))),]
susztak_glo=fread("./08_eQTLs/susztak/cred_var_with_nephQTL_glomerular.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
susztak_tub=fread("./08_eQTLs/susztak/cred_var_with_nephQTL_tubulointerstitial.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)


GTEx_eQTL=fread("./08_eQTLs/gtex/all_cred_var_with_gtex.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
GTEx_sQTL=fread("./09_sQTLs/gtex/all_cred_var_with_gtex.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)


MGI=read.table( "./MGI_lookup/genes_and_pheno/MGI_table_kidney.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill=TRUE, quote="",comment="") #raw data
OMIM=read.table("./OMIM/OMIM_Groopman_combined/all_kidney_genes_sorted_nearly_without_duplicates_ks.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) #raw data improved
ADTKD=read.table("./all_genes_with_classification.txt", sep = "\t", header = T, stringsAsFactors = F) # raw data from manuscript
coloc_tub=read.table( "./coloc/results_T/coloc_tubulo.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coloc_glo=read.table( "./coloc/results_G/coloc_glomerulus.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GPS=read.table("./07_GPS/table_for_GPS_signal_based.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE)

DM_decline = read.table("./DM_decline_effects_withStanzickId_mg_clean.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE) # manually curated
DM_stats = read.table("./12_DM/DM_and_noDM_statistics_for_cred_vars.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE)
drugability = read.table("./11_drugability/targets_drugs_diseases_in_GPS_marked_kidney_drugs.txt", header=TRUE, sep = "\t", stringsAsFactors=FALSE, comment="", quote="",fill=TRUE)



### Gene Name Reference
genes = fread( "./inputs/glist.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

genes_sorted <- genes[order(genes$Chr, genes$Pos1), ]

real_genes = fread( "./inputs/genes_an_ids.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

print("file load finished")

### prepare CADD:


if(any(duplicated(CADD))) CADD = CADD[-which(duplicated(CADD)),]
CADD = merge(CADD, cred_var[,c('rsid','ppa','region_id','signal_id')], by.x="RSID", by.y="rsid", all.x=T, all.y=F, sort=F)
names(CADD)[which(names(CADD)=="ppa")]="PPA"

####prepare MGI:
# gene names
MGI=MGI[which(!is.na(MGI$human_symbol)),]

for(i in 1:nrow(MGI)){
  
  if(!any(MGI$human_symbol[i]==real_genes$Approved.symbol) & any(MGI$human_symbol[i]==real_genes$genes)){
    
    MGI$human_symbol[i]=real_genes$Approved.symbol[which(real_genes$genes==MGI$human_symbol[i])[1]]
  }
}
#prep table
MGI = merge(MGI[,c('phenotype','human_symbol')],GPS[,c('Gene','locus_id')], by.x="human_symbol", by.y="Gene", all.x=FALSE, all.y=FALSE)
if(any(duplicated(MGI))) MGI=MGI[-which(duplicated(MGI)),]
MGI_show=MGI[,c('human_symbol','phenotype','locus_id')]
names(MGI_show)=c("Gene", "Phenotype in mouse", "Locus ID")


###prepare human pheno:
if(any(duplicated(OMIM))){
  kidney_pheno=OMIM[-which(duplicated(OMIM)),]
  }else{kidney_pheno=OMIM}

a = which(duplicated(kidney_pheno[,c(1,2)]))

b= which(kidney_pheno$gene%in%kidney_pheno$gene[a])

kidney_pheno$annot[b]="OMIM; Groopman et al"

OMIM_genes = unique(kidney_pheno$gene)

for(i in 1:length(OMIM_genes)){
  
  disease = paste(kidney_pheno$disease[which(kidney_pheno$gene==OMIM_genes[i])],collapse=" | ")
  
  if(length(which(kidney_pheno$gene==OMIM_genes[i]))==1){
    
    annot = kidney_pheno$annot[which(kidney_pheno$gene==OMIM_genes[i])]
  }else{
    
    if(any(kidney_pheno$annot[which(kidney_pheno$gene==OMIM_genes[i])]=="OMIM; Groopman et al") | (any(kidney_pheno$annot[which(kidney_pheno$gene==OMIM_genes[i])]=="OMIM") & any(kidney_pheno$annot[which(kidney_pheno$gene==OMIM_genes[i])]=="Groopman et al"))){
      
      annot = "OMIM; Groopman et al"
    }else{
      
      annot=kidney_pheno$annot[which(kidney_pheno$gene==OMIM_genes[i])[1]]
    }
  }
  
  
  x = data.frame(OMIM_genes[i],disease,annot,stringsAsFactors=FALSE)
  
  if(i==1){OMIM_all = x}
  
  else{OMIM_all= rbind(OMIM_all,x, stringsAsFactors=FALSE, make.row.names=FALSE)}
  
}

idx = (OMIM_all[,1] %in% GPS$Gene)

OMIM = OMIM_all[idx,]
names(OMIM)=c("Gene","disease","annot")
OMIM = merge(OMIM,GPS[,c('Gene','locus_id')], by="Gene", all.x=FALSE, all.y=FALSE)
if(any(duplicated(OMIM))) OMIM = OMIM[-which(duplicated(OMIM)),]
#combine OMIM & ADTKD
print("OMIM & ADTKD")
for(i in 1:nrow(ADTKD)){
  if(any(OMIM$Gene==ADTKD$Gene[i])){
    if(length(which(OMIM$Gene==ADTKD$Gene[i]))>1)stop(i)
    OMIM$disease[which(OMIM$Gene==ADTKD$Gene[i])]=paste(OMIM$disease[which(OMIM$Gene==ADTKD$Gene[i])], ADTKD$disease[i], sep= " | ")
    OMIM$annot[which(OMIM$Gene==ADTKD$Gene[i])] = paste(OMIM$annot[which(OMIM$Gene==ADTKD$Gene[i])], "Wopperer et al", sep = "; ")
  }else{
    OMIM[(nrow(OMIM)+1),]="NA"
    OMIM$Gene[nrow(OMIM)]=ADTKD$Gene[i]
    OMIM$disease[nrow(OMIM)]=ADTKD$disease[i]
    OMIM$annot[nrow(OMIM)]="Wopperer et al"
    OMIM$locus_id[nrow(OMIM)]=ADTKD$locus_id[i]
  }
}

# add MIM numberr and OMIM links:
MIM_number=read.table("/stk05236/OMIM/OMIM_gene_pheno/genemap2_edited.txt", header = TRUE,  sep="\t", stringsAsFactors = FALSE, quote="", fill=TRUE)
MIM_number$Gene.Symbols[grep(",", MIM_number$Gene.Symbols)]=toupper(MIM_number$Approved.Symbol[grep(",", MIM_number$Gene.Symbols)])
for(i in 1:nrow(MIM_number)){
  
  if(!any(MIM_number$Gene.Symbols[i]==real_genes$Approved.symbol) & any(MIM_number$Gene.Symbols[i]==real_genes$genes)){
    
    MIM_number$Gene.Symbols[i]=real_genes$Approved.symbol[which(real_genes$genes==MIM_number$Gene.Symbols[i])[1]]
  }
}



OMIM_MIM=merge(OMIM, MIM_number[,c("MIM.Number", "Gene.Symbols")], by.x="Gene", by.y="Gene.Symbols", all.x=TRUE, all.y=FALSE, sort = FALSE, na.last)
if(any(OMIM_MIM$Gene=="NPHP3"& OMIM_MIM$MIM.Number==606995)) OMIM_MIM=OMIM_MIM[-(which(OMIM_MIM$Gene=="NPHP3"& OMIM_MIM$MIM.Number==606995)),] #wrong MIM Number, right one is still included
#genes without MIM get sorted to the end. This is the correction:
pos_OMIM=c()
for(i in 1:nrow(OMIM)){
  if(length(which(OMIM_MIM$Gene==OMIM$Gene[i] & OMIM_MIM$locus_id==OMIM$locus_id[i]))!=1)stop(i)
  pos_OMIM[i]=which(OMIM_MIM$Gene==OMIM$Gene[i] & OMIM_MIM$locus_id==OMIM$locus_id[i])
}

OMIM_MIM=OMIM_MIM[pos_OMIM,]

OMIM_MIM$Gene[which(!is.na(OMIM_MIM$MIM.Number))]=paste(c(rep("<a href=\"https://www.omim.org/entry/",nrow(OMIM_MIM)-1)),OMIM_MIM$MIM.Number[which(!is.na(OMIM_MIM$MIM.Number))],c(rep("\" target=\"_blank\">", nrow(OMIM_MIM)-1)),OMIM_MIM$Gene[which(!is.na(OMIM_MIM$MIM.Number))],c(rep("</a>",nrow(OMIM_MIM)-1)),sep="")
names(OMIM_MIM)=c("Gene*","Genetic disorder with kidney phenotype", "Source**", "Locus ID")


#### prepare eQTL tables:
print("eQTLs")
# # prepare gtex
# if(any(duplicated(GTEx_eQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")]))) GTEx_eQTL=GTEx_eQTL[-which(duplicated(GTEx_eQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")])),]
# if(any(duplicated(GTEx_sQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")]))) GTEx_sQTL=GTEx_sQTL[-which(duplicated(GTEx_sQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")])),]
# 
# for(i in 1:nrow(GTEx_eQTL)){
#   
#   if(!any(GTEx_eQTL$gene[i]==real_genes$Approved.symbol) & any(GTEx_eQTL$gene[i]==real_genes$genes)){
#     
#     GTEx_eQTL$gene[i]=real_genes$Approved.symbol[which(real_genes$genes==GTEx_eQTL$gene[i])[1]]
#   }
# }
# for(i in 1:nrow(GTEx_sQTL)){
#   
#   if(!any(GTEx_sQTL$gene[i]==real_genes$Approved.symbol) & any(GTEx_sQTL$gene[i]==real_genes$genes)){
#     
#     GTEx_sQTL$gene[i]=real_genes$Approved.symbol[which(real_genes$genes==GTEx_sQTL$gene[i])[1]]
#   }
# }
# 
# ###duplicate entfernen:
# #rsid,variant_id,gene,tissue
# if(any(duplicated(GTEx_eQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")]))) GTEx_eQTL=GTEx_eQTL[-which(duplicated(GTEx_eQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")])),]
# if(any(duplicated(GTEx_sQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")]))) GTEx_sQTL=GTEx_sQTL[-which(duplicated(GTEx_sQTL[,c("rs_id_dbSNP151_GRCh38p7","variant_id","gene","tissue")])),]


gtex_without_kidney = GTEx_eQTL[which(GTEx_eQTL$tissue != "Kidney_Cortex"),]
sqtl_without_kidney = GTEx_sQTL[which(GTEx_sQTL$tissue != "Kidney_Cortex"),]
GTEx_eQTL = GTEx_eQTL[which(GTEx_eQTL$tissue == "Kidney_Cortex"),]
GTEx_sQTL = GTEx_sQTL[which(GTEx_sQTL$tissue == "Kidney_Cortex"),]

# add ppa, region_id and signal_id to eqtls and sqtls
gtex_without_kidney = merge(gtex_without_kidney[,c("rs_id_dbSNP151_GRCh38p7","gene","tissue","ref","alt","slope","slope_se","pval_nominal")],cred_var[,c("rsid","ppa","chr","pos","eaf","region_id","signal_id")],by.x="rs_id_dbSNP151_GRCh38p7",by.y="rsid",all.x=T,all.y=F,sort=F)
sqtl_without_kidney = merge(sqtl_without_kidney[,c("rs_id_dbSNP151_GRCh38p7","gene","tissue","ref","alt","slope","slope_se","pval_nominal")],cred_var[,c("rsid","ppa","chr","pos","eaf","region_id","signal_id")],by.x="rs_id_dbSNP151_GRCh38p7",by.y="rsid",all.x=T,all.y=F,sort=F)
GTEx_eQTL = merge(GTEx_eQTL[,c("rs_id_dbSNP151_GRCh38p7","gene","ref","alt","slope","slope_se","pval_nominal")],cred_var[,c("rsid","ppa","chr","pos","eaf","region_id","signal_id")],by.x="rs_id_dbSNP151_GRCh38p7",by.y="rsid",all.x=T,all.y=F,sort=F)
GTEx_sQTL = merge(GTEx_sQTL[,c("rs_id_dbSNP151_GRCh38p7","gene","ref","alt","slope","slope_se","pval_nominal")],cred_var[,c("rsid","ppa","chr","pos","eaf","region_id","signal_id")],by.x="rs_id_dbSNP151_GRCh38p7",by.y="rsid",all.x=T,all.y=F,sort=F)



NEPTUNE_glo = merge(NEPTUNE_glo[,c("rsid","CHROM","POS","Approved.symbol","REF","ALT","BETA","SE","P")],cred_var[,c("rsid","eaf","region_id","signal_id","ppa")],by="rsid", all.x=T, all.y=F, sort=F)
susztak_glo = merge(susztak_glo[,c("rsID","Gene_Symbol","Ref","Alt","Beta","Std","Pvalue")],cred_var[,c("rsid","eaf","region_id","signal_id","ppa","chr","pos")],by.x="rsID",by.y="rsid", all.x=T, all.y=F, sort=F)
NEPTUNE_tub = merge(NEPTUNE_tub[,c("rsid","CHROM","POS","Approved.symbol","REF","ALT","BETA","SE","P")],cred_var[,c("rsid","eaf","region_id","signal_id","ppa")],by="rsid", all.x=T, all.y=F, sort=F)
susztak_tub = merge(susztak_tub[,c("rsID","Gene_Symbol","Ref","Alt","Beta","Std","Pvalue")],cred_var[,c("rsid","eaf","region_id","signal_id","ppa","chr","pos")],by.x="rsID",by.y="rsid", all.x=T, all.y=F, sort=F)

if(any(duplicated(susztak_tub))) susztak_tub = susztak_tub[-which(duplicated(susztak_tub)),]
if(any(duplicated(susztak_glo))) susztak_glo = susztak_glo[-which(duplicated(susztak_glo)),]
#sort columns 

gtex_without_kidney = gtex_without_kidney[,c('rs_id_dbSNP151_GRCh38p7','ppa','chr','pos','gene','tissue','ref','alt','slope','slope_se','pval_nominal','eaf','region_id','signal_id')]
sqtl_without_kidney = sqtl_without_kidney[,c('rs_id_dbSNP151_GRCh38p7','ppa','chr','pos','gene','tissue','ref','alt','slope','slope_se','pval_nominal','eaf','region_id','signal_id')]
GTEx_eQTL = GTEx_eQTL[,c('rs_id_dbSNP151_GRCh38p7','ppa','chr','pos','gene','ref','alt','slope','slope_se','pval_nominal','eaf','region_id','signal_id')]
GTEx_sQTL = GTEx_sQTL[,c('rs_id_dbSNP151_GRCh38p7','ppa','chr','pos','gene','ref','alt','slope','slope_se','pval_nominal','eaf','region_id','signal_id')]
NEPTUNE_glo = NEPTUNE_glo[,c('rsid','ppa','CHROM','POS','Approved.symbol','REF','ALT','BETA','SE','P','eaf','region_id','signal_id')]
NEPTUNE_tub = NEPTUNE_tub[,c('rsid','ppa','CHROM','POS','Approved.symbol','REF','ALT','BETA','SE','P','eaf','region_id','signal_id')]
susztak_glo = susztak_glo[,c('rsID','ppa','chr','pos','Gene_Symbol','Ref','Alt','Beta','Std','Pvalue','eaf','region_id','signal_id')]
susztak_tub = susztak_tub[,c('rsID','ppa','chr','pos','Gene_Symbol','Ref','Alt','Beta','Std','Pvalue','eaf','region_id','signal_id')]

#combining susztak and NEPTUNE
names(susztak_glo)=names(NEPTUNE_glo)
names(susztak_tub)=names(NEPTUNE_tub)

glo=rbind(NEPTUNE_glo,susztak_glo)
tub=rbind(NEPTUNE_tub,susztak_tub)



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

names(cred_var_show)=c("RSID*","Chr", "Pos","EA**","EAF***", "Effect conditioned","StdErr conditioned", "P-Value conditioned", "N","PPA","cppa","ci95","ci99","Locus ID","Signal ID","MarkerName", "Allele1", "OA","Effect unconditioned", "StdErr unconditioned", "P-Value unconditioned")
cred_var_show=cred_var_show[,c('RSID*','Locus ID','Signal ID','PPA','N','EA**','OA','EAF***','Effect conditioned','StdErr conditioned','P-Value conditioned','Effect unconditioned','StdErr unconditioned','P-Value unconditioned','Chr','Pos')]
cred_var_show[,which(lapply(cred_var_show,class)=="numeric")]=lapply(cred_var_show[,which(lapply(cred_var_show,class)=="numeric")], signif, digits=2)
cred_var$chr_pos = paste(cred_var$chr, cred_var$pos, sep=":") #prepare search chr:pos
#files for locuszoom
lz_files=list.files("/stk05236/GPS_App/App_test/www")
#lz_files=lz_files[-1]
locus_id_files=c()
for(i in 1:length(lz_files)){
  
  a=unlist(strsplit(lz_files[i],"_",fixed=TRUE))
  locus_id_files[i]=a[2]
}

names(all_sig)=tolower(names(all_sig))



names(all_sig)=c("MarkerName","EA","OA","n","EAF","Effect","StdErr","P.value","Direction","mac","Chr","Pos","P.value.GC","StdErr.GC","RSID","Indep.500k.5e8.aTopHit","Indep.500k.5e8.aLociTag","Indep.500k.5e8.aNumLocusSNPs","Indep.500k.5e8.aLocusCoordinates","Indep.500k.5e8.aLocusSize","Indep.500k.5e8.aLocusAnnot","NearestGene","Nearest Gene Distance","Locus_id")
names(all_sig_show)=c("MarkerName","EA","OA","N","EAF","beta","StdErr","P.value","Direction","mac","Chr","Pos","P.value.GC","StdErr.GC","RSID*","Indep.500k.5e8.aTopHit","Indep.500k.5e8.aLociTag","Indep.500k.5e8.aNumLocusSNPs","Indep.500k.5e8.aLocusCoordinates","Indep.500k.5e8.aLocusSize","Indep.500k.5e8.aLocusAnnot","Nearest gene","Nearest gene distance", "Locus ID")


### prepare region table and signal table

locus_table$cys_bun = ifelse((locus_table$bun_association =="yes" & locus_table$cys_association=="yes"),"cys/bun",ifelse((locus_table$bun_association =="yes" & locus_table$cys_association=="no"),"bun",ifelse((locus_table$bun_association =="no" & locus_table$cys_association=="yes"), "cys", "-"))) 

region_table = merge(signal_table,locus_table[,c('locus_id','bun_association','cys_association')], by.x="Locus_id", by.y="locus_id", sort=F)
region_table[which(region_table$num_genes_in_region=="FALSE"),'num_genes_in_region']="-"
region_table_show=region_table[,c('Locus_id','chr','region_start_250kb','region_end_250kb','locus_size','num_sig_SNP','num_genes_in_region','genes_in_region','no_cred_var_95','no_cred_var_99','signal_id','bun_association','cys_association')]
names(region_table_show)=c("Locus ID", "Chr", "Start of locus", "End of locus", "Locus size", "# significant SNPs in locus*", "# genes in locus", "Genes in locus", "# credible variants 95", "# credible variants", "Signal ID", "BUN validation of locus", "eGFRcys validation of locus")
region_table_show=region_table_show[,c('Locus ID','Signal ID','Chr','Start of locus','End of locus','Locus size','# significant SNPs in locus*','# credible variants','# genes in locus','Genes in locus','BUN validation of locus','eGFRcys validation of locus')]


### prepare show-tables for eQTLs and sQTLs

NEPTUNE_glo_show=NEPTUNE_glo[,-which(names(NEPTUNE_glo)%in%c("region_id","signal_id"))]
NEPTUNE_tub_show=NEPTUNE_tub[,-which(names(NEPTUNE_tub)%in%c("region_id","signal_id"))]
susztak_glo_show=susztak_glo[,-which(names(susztak_glo)%in%c("region_id","signal_id"))]
susztak_tub_show=susztak_tub[,-which(names(susztak_tub)%in%c("region_id","signal_id"))]
GTEx_eQTL_show=GTEx_eQTL[,-which(names(GTEx_eQTL)%in%c("region_id","signal_id"))]
gtex_without_kidney_show=gtex_without_kidney
#gtex_without_kidney[,c(2:length(names(gtex_without_kidney)))]<-NULL # save some workspace
GTEx_sQTL_show=GTEx_sQTL[,-which(names(GTEx_sQTL)%in%c("region_id","signal_id"))]
sqtl_without_kidney_show=sqtl_without_kidney
#sqtl_without_kidney[,c(2:length(names(sqtl_without_kidney)))]<-NULL # save some workspace
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

NEPTUNE_glo_show=NEPTUNE_glo_show[,c('RSID','PPA','Affected gene','Effect allele*','Other allele','EAF','Effect direction on expression','Effect on expression','StdErr','P-Value','Chr','Pos')]
NEPTUNE_tub_show=NEPTUNE_tub_show[,c('RSID','PPA','Affected gene','Effect allele*','Other allele','EAF','Effect direction on expression','Effect on expression','StdErr','P-Value','Chr','Pos')]
susztak_glo_show=susztak_glo_show[,c('RSID','PPA','Affected gene','Effect allele*','Other allele','EAF','Effect direction on expression','Effect on expression','StdErr','P-Value','Chr','Pos')]
susztak_tub_show=susztak_tub_show[,c('RSID','PPA','Affected gene','Effect allele*','Other allele','EAF','Effect direction on expression','Effect on expression','StdErr','P-Value','Chr','Pos')]
GTEx_eQTL_show=GTEx_eQTL_show[,c('RSID','PPA','Affected gene','Effect allele*','Other allele','EAF','Effect direction on expression','Effect on expression','StdErr','P-Value','Chr','Pos')]
GTEx_sQTL_show=GTEx_sQTL_show[,c('RSID','PPA','Affected gene','Effect allele*','Other allele','EAF','Effect on splicing','StdErr','P-Value','Chr','Pos')]
gtex_without_kidney_show=gtex_without_kidney_show[,c('RSID','PPA','Affected gene','Tissue','Effect allele*','Other allele','EAF','Effect direction on expression','Effect on expression','StdErr','P-Value','Chr','Pos','region_id','signal_id')]
sqtl_without_kidney_show=sqtl_without_kidney_show[,c('RSID','PPA','Affected gene','Tissue','Effect allele*','Other allele','EAF','Effect on splicing','StdErr','P-Value','Chr','Pos','region_id','signal_id')]


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




### prepare CADD show-table
CADD_show=CADD[,c('RSID','PPA','GeneName','region_id','signal_id','ConsDetail','Ref','Alt','protPos','oAA','nAA','PHRED')]
names(CADD_show)=c("RSID", "PPA", "Affected gene", "Locus ID", "Signal ID", "Functional consequence", "REF", "ALT",  "Aminoacid position", "Reference aminoacid", "Alternative aminoacid","CADD Phred-Score")
CADD_show$PPA=signif(CADD_show$PPA, digits=2)
if(any(duplicated(CADD_show))) {
  CADD = CADD[-which(duplicated(CADD_show)),]
  CADD_show = CADD_show[-which(duplicated(CADD_show)),]
}

## prepare DM statistics

DM_stats_show = DM_stats[,c("rsid","Allele1","Allele2","Freq1.dm","Effect.dm","StdErr.dm","P.value.dm","n.dm","Freq1.nodm","Effect.nodm","StdErr.nodm","P.value.nodm","n.nodm","pdiff")]
names(DM_stats_show)=c("RSID", "Effect allele", "Other allele", "EAF in DM", "Effect in DM", "StdErr in DM", "P-Value in DM", "N in DM", "EAF in noDM", "Effect in noDM", "StdErr in noDM", "P-Value in noDM", "N in noDM", "P-Value difference between DM and no DM")


### prepare coloc
#coloc_tub=read.table( "/stk05236/lmm/coloc/coloc_tubulo.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#coloc_glo=read.table( "/stk05236/lmm/coloc/coloc_glomerulus.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

coloc_tub=coloc_tub[which(coloc_tub$PP_H4 >=0.80),]
coloc_glo=coloc_glo[which(coloc_glo$PP_H4 >=0.80),]

#	Gen-Namen anpassen
#glom
for(i in 1:nrow(coloc_glo)){
  
  if(!any(coloc_glo$gene[i]==real_genes$Approved.symbol) & any(coloc_glo$gene[i]==real_genes$genes)){
    
    coloc_glo$gene[i]=real_genes$Approved.symbol[which(real_genes$genes==coloc_glo$gene[i])]
  }
}
#tub
for(i in 1:nrow(coloc_tub)){
  
  if(!any(coloc_tub$gene[i]==real_genes$Approved.symbol) & any(coloc_tub$gene[i]==real_genes$genes)){
    
    coloc_tub$gene[i]=real_genes$Approved.symbol[which(real_genes$genes==coloc_tub$gene[i])]
  }
}

#### prepare drugability table

drugability$INDICATI[which(is.na(drugability$INDICATI))] = "-"
drugability_show=drugability[,c('gps_gene','human_gene_name','TARGNAME','TARGTYPE','DRUGNAME','Highest_status','MOA','INDICATI')]
drugability_show$TARGNAME=paste(c(rep("<a href=\"http://db.idrblab.net/ttd/search/ttd/target?search_api_fulltext=",nrow(drugability_show))),drugability$gps_gene,c(rep("\" target=\"_blank\">", nrow(drugability_show))),drugability_show$TARGNAME,c(rep("</a>",nrow(drugability_show))),sep="")
names(drugability_show)=c("Gene", "All respective target genes", "Target name*", "Targettype", "Drug name", "Highest drug status**", "Mode of action", "Disease/Indication")
print("file preparation finished")

### generate GPS-table

for(i in 1:nrow(GPS)){
  GPS[i,'stop.gained.stop.lost.non.synonymus']=length(unique(CADD$RSID[which(CADD$GeneName==GPS$Gene[i]&CADD$signal_id==GPS$signal_id[i]&CADD$region_id==GPS$locus_id[i]&CADD$ConsScore>6)]))
  GPS[i,'canonical.splice.noncoding.change.synonymous.splice.site']=length(unique(CADD$RSID[which(CADD$GeneName==GPS$Gene[i]&CADD$signal_id==GPS$signal_id[i]&CADD$region_id==GPS$locus_id[i]&CADD$ConsScore<=6&CADD$ConsScore>4)]))
  GPS[i,'Cadd15_other']=length(unique(CADD$RSID[which(CADD$GeneName==GPS$Gene[i]&CADD$signal_id==GPS$signal_id[i]&CADD$region_id==GPS$locus_id[i]&CADD$ConsScore<=4)]))
  GPS[i,'NEPTUNE_glomerulus']=length(unique(glo$rsid[which(glo$Approved.symbol==GPS$Gene[i]&glo$signal_id==GPS$signal_id[i]&glo$region_id==GPS$locus_id[i])]))
  GPS[i,'NEPTUNE_tubulointerstitium']=length(unique(tub$rsid[which(tub$Approved.symbol==GPS$Gene[i]&tub$signal_id==GPS$signal_id[i]&tub$region_id==GPS$locus_id[i])]))
  GPS[i,'GTEX_eQTL_kidney_tissue']=length(unique(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_eQTL$gene==GPS$Gene[i]&GTEx_eQTL$signal_id==GPS$signal_id[i]&GTEx_eQTL$region_id==GPS$locus_id[i])]))
  GPS[i,'GTEX_eQTL_any_other_tissue']=length(unique(gtex_without_kidney_show$RSID[which(gtex_without_kidney_show$`Affected gene`==GPS$Gene[i]&gtex_without_kidney_show$signal_id==GPS$signal_id[i]&gtex_without_kidney_show$region_id==GPS$locus_id[i])]))
  GPS[i,'GTEX_sQTL_kidney_tissue']=length(unique(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_sQTL$gene==GPS$Gene[i]&GTEx_sQTL$signal_id==GPS$signal_id[i]&GTEx_sQTL$region_id==GPS$locus_id[i])]))
  GPS[i,'GTEX_sQTL_any_other_tissue']=length(unique(sqtl_without_kidney_show$RSID[which(sqtl_without_kidney_show$`Affected gene`==GPS$Gene[i]&sqtl_without_kidney_show$signal_id==GPS$signal_id[i]&sqtl_without_kidney_show$region_id==GPS$locus_id[i])]))
  GPS[i,'MGI_mouse_kidney_phenotyp']= length(which(GPS$Gene[i]==MGI$human_symbol))
  GPS[i,'OMIM_human_kidney_phenotype'] = length(unlist(strsplit(OMIM$disease[which(OMIM$Gene==GPS$Gene[i])], split = "|", fixed=T)))
  a=paste(GPS$locus_id[i],GPS$signal_id[i],sep=".")
  GPS[i,'coloc']="-"
  if(any(GPS$Gene[i]==coloc_tub$gene & a==coloc_tub$LocusID) & any(GPS$Gene[i]==coloc_glo$gene& a==coloc_glo$LocusID)){
    GPS[i,'coloc']="tubulo-interstitial/glomerulus"
  }else{
    if(any(GPS$Gene[i]==coloc_tub$gene & a==coloc_tub$LocusID)) GPS[i,'coloc']="tubulo-interstitial"
    if(any(GPS$Gene[i]==coloc_glo$gene& a==coloc_glo$LocusID)) GPS[i,'coloc']="glomerulus"
  }
  GPS[i,'Cys.BUN']=locus_table$cys_bun[which(locus_table$locus_id==GPS$locus_id[i])]
  GPS[i,'credible_variants_per_locus_signal'] = length(unique(cred_var$rsid[which(cred_var$region_id==GPS$locus_id[i]&cred_var$signal_id==GPS$signal_id[i])]))
  if(any(GPS$Gene[i]==drugability$gps_gene[which(drugability$includes_kideny_disease)])){
	GPS[i,'drugable'] = "yes for kidney disease"
  }else {
	if(any(GPS$Gene[i]==drugability$gps_gene[which(!drugability$includes_kideny_disease)])){
		GPS[i,'drugable'] = "yes for other disease"
	}else{GPS[i,'drugable'] = "no"}
  }
  if(!i%%2000)print(i)
  }

#columns DM decline drugability pathways

GPS['DM_effect']=rep("-",nrow(GPS))
GPS['decline_effect']=rep("-",nrow(GPS))
# GPS['drugable']=rep("-",nrow(GPS))
# GPS['enriched_pathway']=rep("-",nrow(GPS))
# 
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
max_ppa=c()
for(i in 1:nrow(GPS)){
  u=cred_var[which(cred_var$region_id==GPS$locus_id[i]&cred_var$signal_id==GPS$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 3)
}
GPS['max_PPA']=max_ppa


# create GPS restricted to PPA>0.1
GPS_10=GPS
for(i in 1:nrow(GPS_10)){
  GPS_10$NEPTUNE_glomerulus[i]=length(unique(glo$rsid[which(glo$Approved.symbol==GPS_10$Gene[i]&glo$signal_id==GPS_10$signal_id[i]&glo$region_id==GPS_10$locus_id[i]&glo$ppa>0.1)]))
  GPS_10$NEPTUNE_tubulointerstitium[i]=length(unique(tub$rsid[which(tub$Approved.symbol==GPS_10$Gene[i]&tub$signal_id==GPS_10$signal_id[i]&tub$region_id==GPS_10$locus_id[i]&tub$ppa>0.1)]))
  GPS_10$GTEX_eQTL_kidney_tissue[i]=length(unique(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_eQTL$gene==GPS_10$Gene[i]&GTEx_eQTL$signal_id==GPS_10$signal_id[i]&GTEx_eQTL$region_id==GPS_10$locus_id[i]&GTEx_eQTL$ppa>0.1)]))
  GPS_10$GTEX_eQTL_any_other_tissue[i]=length(unique(gtex_without_kidney_show$RSID[which(gtex_without_kidney_show$`Affected gene`==GPS_10$Gene[i]&gtex_without_kidney_show$signal_id==GPS_10$signal_id[i]&gtex_without_kidney_show$region_id==GPS_10$locus_id[i]&gtex_without_kidney_show$PPA>0.1)]))
  GPS_10$GTEX_sQTL_kidney_tissue[i]=length(unique(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_sQTL$gene==GPS_10$Gene[i]&GTEx_sQTL$signal_id==GPS_10$signal_id[i]&GTEx_sQTL$region_id==GPS_10$locus_id[i]&GTEx_sQTL$ppa>0.1)]))
  GPS_10$GTEX_sQTL_any_other_tissue[i]=length(unique(sqtl_without_kidney_show$RSID[which(sqtl_without_kidney_show$`Affected gene`==GPS_10$Gene[i]&sqtl_without_kidney_show$signal_id==GPS_10$signal_id[i]&sqtl_without_kidney_show$region_id==GPS_10$locus_id[i]&sqtl_without_kidney_show$PPA>0.1)]))
  GPS_10$stop.gained.stop.lost.non.synonymus[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_10$Gene[i]&CADD$signal_id==GPS_10$signal_id[i]&CADD$region_id==GPS_10$locus_id[i]&CADD$PPA>0.1&CADD$ConsScore>6)]))
  GPS_10$canonical.splice.noncoding.change.synonymous.splice.site[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_10$Gene[i]&CADD$signal_id==GPS_10$signal_id[i]&CADD$region_id==GPS_10$locus_id[i]&CADD$PPA>0.1&CADD$ConsScore<=6&CADD$ConsScore>4)]))
  GPS_10$Cadd15_other[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_10$Gene[i]&CADD$signal_id==GPS_10$signal_id[i]&CADD$region_id==GPS_10$locus_id[i]&CADD$PPA>0.1&CADD$ConsScore<=4)]))
  if(!i%%2000)print(i)
  }

# create GPS restricted to PPA>0.5
GPS_50=GPS
for(i in 1:nrow(GPS_50)){
  GPS_50$NEPTUNE_glomerulus[i]=length(unique(glo$rsid[which(glo$Approved.symbol==GPS_50$Gene[i]&glo$signal_id==GPS_50$signal_id[i]&glo$region_id==GPS_50$locus_id[i]&glo$ppa>0.5)]))
  GPS_50$NEPTUNE_tubulointerstitium[i]=length(unique(tub$rsid[which(tub$Approved.symbol==GPS_50$Gene[i]&tub$signal_id==GPS_50$signal_id[i]&tub$region_id==GPS_50$locus_id[i]&tub$ppa>0.5)]))
  GPS_50$GTEX_eQTL_kidney_tissue[i]=length(unique(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_eQTL$gene==GPS_50$Gene[i]&GTEx_eQTL$signal_id==GPS_50$signal_id[i]&GTEx_eQTL$region_id==GPS_50$locus_id[i]&GTEx_eQTL$ppa>0.5)]))
  GPS_50$GTEX_eQTL_any_other_tissue[i]=length(unique(gtex_without_kidney_show$RSID[which(gtex_without_kidney_show$`Affected gene`==GPS_50$Gene[i]&gtex_without_kidney_show$signal_id==GPS_50$signal_id[i]&gtex_without_kidney_show$region_id==GPS_50$locus_id[i]&gtex_without_kidney_show$PPA>0.5)]))
  GPS_50$GTEX_sQTL_kidney_tissue[i]=length(unique(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7[which(GTEx_sQTL$gene==GPS_50$Gene[i]&GTEx_sQTL$signal_id==GPS_50$signal_id[i]&GTEx_sQTL$region_id==GPS_50$locus_id[i]&GTEx_sQTL$ppa>0.5)]))
  GPS_50$GTEX_sQTL_any_other_tissue[i]=length(unique(sqtl_without_kidney_show$RSID[which(sqtl_without_kidney_show$`Affected gene`==GPS_50$Gene[i]&sqtl_without_kidney_show$signal_id==GPS_50$signal_id[i]&sqtl_without_kidney_show$region_id==GPS_50$locus_id[i]&sqtl_without_kidney_show$PPA>0.5)]))
  GPS_50$stop.gained.stop.lost.non.synonymus[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_50$Gene[i]&CADD$signal_id==GPS_50$signal_id[i]&CADD$region_id==GPS_50$locus_id[i]&CADD$PPA>0.5&CADD$ConsScore>6)]))
  GPS_50$canonical.splice.noncoding.change.synonymous.splice.site[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_50$Gene[i]&CADD$signal_id==GPS_50$signal_id[i]&CADD$region_id==GPS_50$locus_id[i]&CADD$PPA>0.5&CADD$ConsScore<=6&CADD$ConsScore>4)]))
  GPS_50$Cadd15_other[i]=length(unique(CADD$RSID[which(CADD$GeneName==GPS_50$Gene[i]&CADD$signal_id==GPS_50$signal_id[i]&CADD$region_id==GPS_50$locus_id[i]&CADD$PPA>0.5&CADD$ConsScore<=4)]))
  if(!i%%2000)print(i)
  }

#format GPS entries:
GPS_show=GPS[,c('Locus_name','locus_id','signal_id','Gene','distance_to_lead_variant','Chr','Start_of_gene','End_of_gene','Cys.BUN','credible_variants_per_locus_signal','stop.gained.stop.lost.non.synonymus','canonical.splice.noncoding.change.synonymous.splice.site','Cadd15_other','NEPTUNE_glomerulus','NEPTUNE_tubulointerstitium','GTEX_eQTL_kidney_tissue','GTEX_eQTL_any_other_tissue','GTEX_sQTL_kidney_tissue','GTEX_sQTL_any_other_tissue','MGI_mouse_kidney_phenotyp','OMIM_human_kidney_phenotype','coloc','DM_effect','decline_effect','drugable')]
GPS_show$Gene=paste(c(rep("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=",nrow(GPS_show))),GPS$Gene,c(rep("\" target=\"_blank\">", nrow(GPS_show))),GPS$Gene,c(rep("</a>",nrow(GPS_show))),sep="")

names(GPS_show)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Position gene start","Position gene end", "eGFRcys or BUN validation", "# credible variants in signal","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerulus (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL kidney cortex (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex (GTEx)","sQTL other tissue (GTEx)","# kidney phenotypes in mouse","# kidney phenotypes in human", "Coloc in NEPTUNE tissue", "DM_effect", "decline_effect", "known_drug_target")
GPS_show[,1]=paste("[",GPS_show[,1],"]",sep="")
GPS_show[which(GPS_show[,1]=="[PDILT]"),1]="[UMOD-PDILT]"
for(i in 1:nrow(GPS_show)){
  GPS_show[i,'Score']=11-length(which(GPS_show[i,c(11:21)]==0))
  
}
GPS_show$Score=as.integer(GPS_show$Score)
GPS_show$`eGFRcys or BUN validation`=as.factor(GPS_show$`eGFRcys or BUN validation`)
GPS_show$`Coloc in NEPTUNE tissue`=as.factor(GPS_show$`Coloc in NEPTUNE tissue`)
GPS_show$`Signal ID`=as.factor(paste("Signal ", GPS_show$`Signal ID`, sep=""))
GPS_show$DM_effect=as.factor(GPS_show$DM_effect)
GPS_show$decline_effect=as.factor(GPS_show$decline_effect)
GPS_show$`Locus ID`=as.factor(GPS_show$`Locus ID`) 
GPS_show$known_drug_target = as.factor(GPS_show$known_drug_target)


names(cred_var)[which(names(cred_var)=="region_id")]="Locus Id"
names(cred_var)[which(names(cred_var)=="signal_id")]="Signal Id"

max_ppa=c()
for(i in 1:nrow(GPS_show)){
  u=cred_var[which(cred_var$`Locus Id`==GPS$locus_id[i]&cred_var$`Signal Id`==GPS$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 2)
}
GPS_show['max PPA']=max_ppa


##same for PPA>0.1
GPS_10_show=GPS_10[,c('Locus_name','locus_id','signal_id','Gene','distance_to_lead_variant','Chr','Start_of_gene','End_of_gene','Cys.BUN','credible_variants_per_locus_signal','stop.gained.stop.lost.non.synonymus','canonical.splice.noncoding.change.synonymous.splice.site','Cadd15_other','NEPTUNE_glomerulus','NEPTUNE_tubulointerstitium','GTEX_eQTL_kidney_tissue','GTEX_eQTL_any_other_tissue','GTEX_sQTL_kidney_tissue','GTEX_sQTL_any_other_tissue','MGI_mouse_kidney_phenotyp','OMIM_human_kidney_phenotype','coloc','DM_effect','decline_effect','drugable')]
GPS_10_show$Gene=GPS_show$Gene

names(GPS_10_show)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Position gene start","Position gene end", "eGFRcys or BUN validation", "# credible variants in signal","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerulus (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL kidney cortex (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex (GTEx)","sQTL other tissue (GTEx)","# kidney phenotypes in mouse","# kidney phenotypes in human", "Coloc in NEPTUNE tissue", "DM_effect", "decline_effect","known_drug_target")
GPS_10_show[,1]=GPS_show[,1]
GPS_10_show[which(GPS_10_show[,1]=="[PDILT]"),1]="[UMOD-PDILT]"
for(i in 1:nrow(GPS_10_show)){
  GPS_10_show[i,'Score']=11-length(which(GPS_10_show[i,c(11:21)]==0))
  
}
GPS_10_show$Score=as.integer(GPS_10_show$Score)
GPS_10_show$`eGFRcys or BUN validation`=as.factor(GPS_10_show$`eGFRcys or BUN validation`)
GPS_10_show$`Coloc in NEPTUNE tissue`=as.factor(GPS_10_show$`Coloc in NEPTUNE tissue`)
GPS_10_show$`Signal ID`=as.factor(paste("Signal ", GPS_10_show$`Signal ID`, sep=""))
GPS_10_show$DM_effect=as.factor(GPS_10_show$DM_effect)
GPS_10_show$decline_effect=as.factor(GPS_10_show$decline_effect)
GPS_10_show$`Locus ID`=as.factor(GPS_10_show$`Locus ID`)
GPS_10_show$known_drug_target = as.factor(GPS_10_show$known_drug_target)
max_ppa=c()
for(i in 1:nrow(GPS_10_show)){
  u=cred_var[which(cred_var$`Locus Id`==GPS_10$locus_id[i]&cred_var$`Signal Id`==GPS_10$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 2)
}
GPS_10_show['max PPA']=max_ppa


##same for PPA>0.5
GPS_50_show=GPS_50[,c('Locus_name','locus_id','signal_id','Gene','distance_to_lead_variant','Chr','Start_of_gene','End_of_gene','Cys.BUN','credible_variants_per_locus_signal','stop.gained.stop.lost.non.synonymus','canonical.splice.noncoding.change.synonymous.splice.site','Cadd15_other','NEPTUNE_glomerulus','NEPTUNE_tubulointerstitium','GTEX_eQTL_kidney_tissue','GTEX_eQTL_any_other_tissue','GTEX_sQTL_kidney_tissue','GTEX_sQTL_any_other_tissue','MGI_mouse_kidney_phenotyp','OMIM_human_kidney_phenotype','coloc','DM_effect','decline_effect','drugable')]
GPS_50_show$Gene=GPS_show$Gene
names(GPS_50_show)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Position gene start","Position gene end", "eGFRcys or BUN validation", "# credible variants in signal","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerulus (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat Genet, 2021])","eQTL kidney cortex (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex (GTEx)","sQTL other tissue (GTEx)","# kidney phenotypes in mouse","# kidney phenotypes in human", "Coloc in NEPTUNE tissue", "DM_effect", "decline_effect","known_drug_target")
GPS_50_show[,1]=GPS_show[,1]
GPS_50_show[which(GPS_50_show[,1]=="[PDILT]"),1]="[UMOD-PDILT]"
for(i in 1:nrow(GPS_50_show)){
  GPS_50_show[i,'Score']=11-length(which(GPS_50_show[i,c(11:21)]==0))
  
}
GPS_50_show$Score=as.integer(GPS_50_show$Score)
GPS_50_show$`eGFRcys or BUN validation`=as.factor(GPS_50_show$`eGFRcys or BUN validation`)
GPS_50_show$`Coloc in NEPTUNE tissue`=as.factor(GPS_50_show$`Coloc in NEPTUNE tissue`)
GPS_50_show$`Signal ID`=as.factor(paste("Signal ", GPS_50_show$`Signal ID`, sep=""))
GPS_50_show$DM_effect=as.factor(GPS_50_show$DM_effect)
GPS_50_show$decline_effect=as.factor(GPS_50_show$decline_effect)
GPS_50_show$`Locus ID`=as.factor(GPS_50_show$`Locus ID`)
GPS_50_show$known_drug_target = as.factor(GPS_50_show$known_drug_target)
max_ppa=c()
for(i in 1:nrow(GPS_50_show)){
  u=cred_var[which(cred_var$`Locus Id`==GPS_50$locus_id[i]&cred_var$`Signal Id`==GPS_50$signal_id[i]),]
  max_ppa[i]=signif(max(u$ppa) ,digits = 2)
}
GPS_50_show['max PPA']=max_ppa


# prepare filter options
effector_genes_noppa = which(GPS$stop.gained.stop.lost.non.synonymus!=0|GPS$canonical.splice.noncoding.change.synonymous.splice.site!=0|GPS$Cadd15_other!=0|GPS$NEPTUNE_glomerulus!=0|GPS$NEPTUNE_tubulointerstitium!=0|GPS$GTEX_eQTL_kidney_tissue!=0|GPS$GTEX_sQTL_kidney_tissue!=0)
effector_genes_10ppa = which(GPS_10$stop.gained.stop.lost.non.synonymus!=0|GPS_10$canonical.splice.noncoding.change.synonymous.splice.site!=0|GPS_10$Cadd15_other!=0|GPS_10$NEPTUNE_glomerulus!=0|GPS_10$NEPTUNE_tubulointerstitium!=0|GPS_10$GTEX_eQTL_kidney_tissue!=0|GPS_10$GTEX_sQTL_kidney_tissue!=0)
effector_genes_50ppa = which(GPS_50$stop.gained.stop.lost.non.synonymus!=0|GPS_50$canonical.splice.noncoding.change.synonymous.splice.site!=0|GPS_50$Cadd15_other!=0|GPS_50$NEPTUNE_glomerulus!=0|GPS_50$NEPTUNE_tubulointerstitium!=0|GPS_50$GTEX_eQTL_kidney_tissue!=0|GPS_50$GTEX_sQTL_kidney_tissue!=0)
kidney_genes = which(GPS$MGI_mouse_kidney_phenotyp!=0|GPS$OMIM_human_kidney_phenotype!=0)



sourceDir <- getSrcDirectory(function(dummy) {dummy})
links=scan(file = paste(dummy,"www/links_plots.txt",sep="/"), what="character")
links_id<-c()
for(i in 1:length(links)){
  z=unlist(strsplit(links[i],"/"))
  links_id[i]=unlist(strsplit(z[6],"_"))[2]
}



workspace_image=paste(sourceDir,"data/workspace.RData", sep="/")
save.image(file=workspace_image)