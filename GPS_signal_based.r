setwd("/stk05236/")

library(data.table)

GPS=read.table("./lmm/all_loci/06_GPS_table/GPS_all_at_once_separate_signals_no_dup.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

region_table = read.table( "./lmm/all_loci/05_region_table/region_table_cred_var.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

GPS_by_signal=c()



for(i in 1:nrow(GPS)){

	print(i)

	a = max(region_table$signal_id[which(region_table$Locus_id==GPS$locus_id[i])])
	
	signal_id =c(1:a)
	
	x = GPS[c(rep(i,a)),]
	
	x['signal_id']=signal_id
	
	x['num_signals']=rep(a,a)
	
	for(y in 1:nrow(x)){
	
		if(x[y,9]!=0){
		
			x[y,9]=unlist(strsplit(x[y,9], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,10]!=0){
		
			x[y,10]=unlist(strsplit(x[y,10], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,11]!=0){
		
			x[y,11]=unlist(strsplit(x[y,11], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,12]!=0){
		
			x[y,12]=unlist(strsplit(x[y,12], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,13]!=0){
		
			x[y,13]=unlist(strsplit(x[y,13], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,14]!=0){
		
			x[y,14]=unlist(strsplit(x[y,14], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,15]!=0){
		
			x[y,15]=unlist(strsplit(x[y,15], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,16]!=0){
		
			x[y,16]=unlist(strsplit(x[y,16], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,17]!=0){
		
			x[y,17]=unlist(strsplit(x[y,17], split="|", fixed=TRUE))[y]
		}
		
		if(x[y,23]!=0){
		
			x[y,23]=unlist(strsplit(x[y,23], split="|", fixed=TRUE))[y]
		}
	}
	
	if(i==1){
		GPS_by_signal=x
	}else{
		GPS_by_signal=rbind(GPS_by_signal,x)
	}
}


##### correct Coloc

coloc_tub=read.table( "./lmm/coloc/coloc_tubulo.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

coloc_glo=read.table( "./lmm/coloc/coloc_glomerulus.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

coloc_tub=coloc_tub[which(coloc_tub$PP_H4 >=0.80),]
coloc_glo=coloc_glo[which(coloc_glo$PP_H4 >=0.80),]
#
#	Gen-Namen anpassen
#
real_genes = read.table("./inputs/synonym_table.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

real_genes = real_genes[which(!is.na(real_genes$approved_symbol)),]

real_genes = real_genes[which(!is.na(real_genes$genes)),]

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

for (i in 1:nrow(GPS_by_signal)){

	a=paste(GPS_by_signal$locus_id[i],GPS_by_signal$signal_id[i],sep=".")

	coloc[i]="-"

	if(any(GPS_by_signal$Gene[i]==coloc_tub$gene & a==coloc_tub$LocusID) & any(GPS_by_signal$Gene[i]==coloc_glo$gene& a==coloc_glo$LocusID)){
	
		coloc[i]="tubulo-interstitial/glomerulus"
	}else{
	
		if(any(GPS_by_signal$Gene[i]==coloc_tub$gene & a==coloc_tub$LocusID)) coloc[i]="tubulo-interstitial"
		
		if(any(GPS_by_signal$Gene[i]==coloc_glo$gene& a==coloc_glo$LocusID)) coloc[i]="glomerulus"
	}
}

GPS_by_signal$Coloc=coloc

write.table(GPS_by_signal, file="./lmm/all_loci/06_GPS_table/GPS_signal_based.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)