library(data.table)

setwd("/stk05236/clean_GPS/")

cadd_new = fread( "./10_protein_altering/results/new_cadd16phred15_short.out",  header = TRUE, sep = " ", stringsAsFactors = FALSE,  data.table = FALSE)


#######
#
#	richtig Allele raussuchen
#
#######


#######
## TW: Vielleicht ist das schneller: 
MarkerName = paste(cadd_new[,1], cadd_new$Pos,sep=":")
isRefFirst = cadd_new$Ref < cadd_new$Alt
MarkerName = ifelse(isRefFirst, paste(MarkerName,":",cadd_new$Ref,"_",cadd_new$Alt,sep=""), paste(MarkerName,":",cadd_new$Alt,"_",cadd_new$Ref,sep=""))
#######

ref = fread( "/stk05236/lmm/03_eval/european/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

SNPs_in_analysis = which (MarkerName %in% ref$MarkerName)

cadd_new['MarkerName'] <- MarkerName

cadd_new_right_SNPs = cadd_new[SNPs_in_analysis,]

#######
#
#	alle ohne annotiertes Gen rausschmeißen
#
#######

cadd_new_right_SNPs = cadd_new_right_SNPs[which(!is.na(cadd_new_right_SNPs$GeneName)),]

#######
#
#	SNP soll in Gen liegen, auf das er sich bezieht
#
########

gen_ref = fread( "./inputs/glist.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

in_gene =c ()

for (i in 1:nrow(cadd_new_right_SNPs)){

	if(!any(cadd_new_right_SNPs$GeneName[i]==real_genes$approved_symbol) & any(cadd_new_right_SNPs$GeneName[i]==real_genes$genes)){
	
		cadd_new_right_SNPs$GeneName[i]=real_genes$approved_symbol[which(real_genes$genes==cadd_new_right_SNPs$GeneName[i])[1]]
	}

	if(any(gen_ref$Chr == cadd_new_right_SNPs[i,1] & gen_ref$Pos1 <= cadd_new_right_SNPs$Pos[i] & gen_ref$Pos2 >= cadd_new_right_SNPs$Pos[i] & gen_ref$Gene == cadd_new_right_SNPs$GeneName[i])){
	
		in_gene[i] = TRUE
		
	}else{ 
		
		in_gene[i] = FALSE
	}
}

#

######## die nächste Zeile schränkt auf SNPs in genes ein. Will man alle mit CADD15 muss man das raus nehmen 

cadd_new_right_SNPs_and_gene =cadd_new_right_SNPs[which(in_gene),]


cred_var_99 = read.table("./05_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cadd_with_rsid = merge(cadd_new_right_SNPs_and_gene, ref[,c(1,15)], by="MarkerName", all.x=TRUE, all.y =FALSE, sort=FALSE)

in_cred_var = which(cadd_with_rsid$RSID %in% cred_var_99$rsid)

cadd_in_cred_var = cadd_with_rsid[in_cred_var,]

write.table(cadd_in_cred_var, file="./10_protein_altering/results/cadd_in_cred_var.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

