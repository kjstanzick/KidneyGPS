setwd("/stk05236/clean_GPS/")

library(data.table)

region_table = read.table( "./04_locus_table/locus_table.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

leads = read.table( "./02_cys_bun/new_meta_MHC_cys_bun_validation.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

bun_association = c()

cys_association =c()

for (i in 1:nrow(leads)){

	if(leads$effect_direction_bun[i] =="opposite" & leads$nominal_significant_bun[i]=="yes"){
	
		bun_association[i]="yes"
		
	}else{
	
		bun_association[i]="no"
	}
		
	if(leads$effect_direction_cys[i] =="same" & leads$nominal_significant_cys[i]=="yes"){
	
		cys_association[i]="yes"
		
	}else{
	
		cys_association[i]="no"
	}
}

locus_id=leads$Locus_id

association=data.frame(locus_id,bun_association,cys_association,stringsAsFactors=FALSE)

names(region_table)[which(names(region_table)=="Locus_id")]="locus_id"

region_table1=merge(region_table,association,by="locus_id",all.x=TRUE,all.y=FALSE)

write.table(region_table1, file="./04_locus_table/locus_table_cys_bun.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

