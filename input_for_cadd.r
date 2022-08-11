setwd("/stk05236/")

cred_99 = read.table("./lmm/all_loci/07_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ids=unique(cred_99$region_id)

locus_id=c()

chr=c()

pos1=c()

pos2=c()

signal_id=c()

for(i in 1:length(ids)){

	x=cred_99[which(cred_99$region_id==ids[i]),]
	
	for(y in 1:max(x$signal_id)){
	
		a=x[which(x$signal_id==y),]
	
		chr[length(chr)+1]=unique(a$chr)
	
		pos1[length(pos1)+1]=min(a$pos)
	
		pos2[length(pos2)+1]=max(a$pos)
	
		locus_id[length(locus_id)+1]=ids[i]
		
		signal_id[length(signal_id)+1]=y
	}
}

pos=paste(pos1,pos2,sep="-")

chr_pos=paste(chr,pos,sep=":")

input_cadd = data.frame(locus_id,signal_id,chr_pos,stringsAsFactors=FALSE)

write.table(input_cadd, file="./lmm/all_loci/cadd_lmm/inputs/input_cadd.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


