library(data.table)

setwd("/stk05236/")

cadd_new = fread( "./lmm/all_loci/cadd_lmm/results/new_cadd16phred15_short.out",  header = TRUE, sep = " ", stringsAsFactors = FALSE,  data.table = FALSE)

rows = c(1:nrow(cadd_new))

cadd_new['rows'] <- rows



#######
#
#	richtig Allele raussuchen
#
#######


# alleles = cadd_new[,c(3,4)] ## REF & ALT
# alleles = t(apply(alleles,1,FUN=function(x) sort(x,decreasing=FALSE))) #REF & ALT alphabetisch sortieren
# alleles = paste(alleles, collapse="_")
# MarkerName = paste(cadd_new[,1], cadd_new$Pos, allels, sep=":") # MarkerName um diesen mit Ref zu vergleichen

#######
## TW: Vielleicht ist das schneller: 
MarkerName = paste(cadd_new[,1], cadd_new$Pos,sep=":")
isRefFirst = cadd_new$Ref < cadd_new$Alt
MarkerName = ifelse(isRefFirst, paste(MarkerName,":",cadd_new$Ref,"_",cadd_new$Alt,sep=""), paste(MarkerName,":",cadd_new$Alt,"_",cadd_new$Ref,sep=""))
#######



#ref = fread( "./lmm/03_eval/metal_eGFR_meta1.TBL.map.annot.gc.gz",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)
ref = fread( "./lmm/03_eval/european/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

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

gen_ref = fread( "./inputs/genes/glist-hg19.tw.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

real_genes = read.table("./inputs/synonym_table.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, fill=TRUE, quote="")

real_genes = real_genes[which(!is.na(real_genes$approved_symbol)),]

real_genes = real_genes[which(!is.na(real_genes$genes)),]

for(i in 1:nrow(gen_ref)){

	if(!any(gen_ref$Gene[i]==real_genes$approved_symbol) & any(gen_ref$Gene[i]==real_genes$genes)){
	
		gen_ref$Gene[i]=real_genes$approved_symbol[which(real_genes$genes==gen_ref$Gene[i])[1]]
	}
}


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

#rows_interesting_cadds = cadd_new_right_SNPs$rows[which(in_gene)]

#write.table(rows_interesting_cadds, file="./cadd_lmm/results/interesting_rows_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

######## die nächste Zeile schränkt auf SNPs in genes ein. Will man alle mit CADD15 muss man das raus nehmen 

cadd_new_right_SNPs_and_gene =cadd_new_right_SNPs[which(in_gene),]

#Consequence=unique(cadd_new_right_SNPs_and_gene$Consequence)

#write.table(Consequence, file="./cadd_lmm/results/all_Consequence_right_SNP_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


Marker_unique = unique(cadd_new_right_SNPs_and_gene$MarkerName)

highest_score= rep(FALSE,nrow(cadd_new_right_SNPs_and_gene))

rows2= c(1:nrow(cadd_new_right_SNPs_and_gene))

cadd_new_right_SNPs_and_gene['rows2'] <- rows2

for(i in 1:length(Marker_unique)){

	x=cadd_new_right_SNPs_and_gene[which(cadd_new_right_SNPs_and_gene$MarkerName==Marker_unique[i]),]
	
	maximum = max(x$ConsScore)
	
	b = which(x$ConsScore == maximum)
	
	
	
	pos = x$rows2[b]
	
	highest_score[pos] = TRUE
}

cadd_right_highest_score = cadd_new_right_SNPs_and_gene[which(highest_score),]

highest_Cons = unique(cadd_right_highest_score$Consequence)

####duplicate entfernen:
#cadd_right_highest_score=cadd_right_highest_score[-which(duplicated(cadd_right_highest_score[,-c(16,17,19)])),]

#write.table(highest_Cons, file="./cadd_lmm/results/highest_Cons_Consequence_right_SNP_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#duplicats = cadd_right_highest_score$rows2[which(duplicated(cadd_right_highest_score$MarkerName))]

#Cons_Table_highest_without_duplicats = table(cadd_right_highest_score$Consequence[-duplicats])

#Cons_Table_all_without_duplicats = table(cadd_new_right_SNPs_and_gene$Consequence[-duplicats])

#write.table(Cons_Table_all_without_duplicats, file="./cadd_lmm/results/Consequence_table_without_duplicats_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(Cons_Table_highest_without_duplicats, file="./cadd_lmm/results/highest_Consequence_table_without_duplicats_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

more_consequences = rep(FALSE,nrow(cadd_new_right_SNPs_and_gene))

for(i in 1:length(Marker_unique)){

	x=cadd_new_right_SNPs_and_gene[which(cadd_new_right_SNPs_and_gene$MarkerName==Marker_unique[i]),]
	
	if (nrow(x)>1 & any(duplicated(x$gene))){
	
		more_consequences[x$rows2]=TRUE
	
	}
}
more_than_one_consequence = cadd_new_right_SNPs_and_gene[which(more_consequences),]

cred_var_99 = read.table("./lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cadd_with_rsid = merge(cadd_new_right_SNPs_and_gene, ref[,c(1,15)], by="MarkerName", all.x=TRUE, all.y =FALSE, sort=FALSE)

in_cred_var = which(cadd_with_rsid$RSID %in% cred_var_99$rsid)

cadd_in_cred_var = cadd_with_rsid[in_cred_var,]

write.table(cadd_in_cred_var, file="./lmm/all_loci/cadd_lmm/results/cadd_in_cred_var.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(cadd_in_cred_var$rows, file="./cadd_lmm/results/rows_cadd_in_cred_var_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#variants_with_more_consequences =cadd_in_cred_var[which(cadd_in_cred_var$RSID%in%cadd_in_cred_var$RSID[which(duplicated(cadd_in_cred_var$RSID))]),c(16,4,5,7,8,9,12,13)]

#write.table(variants_with_more_consequences, file="./cadd_lmm/results/variants_with_more_consequences_in_cred_var_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#cons_table = table(cadd_new_right_SNPs_and_gene$Consequence)
#
#write.table(cons_table, file="./cadd_lmm/results/cons_table_cred_regions_new.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
