rm(list = ls(all = TRUE))


setwd("/stk05236/clean_GPS/10_protein_altering/")

dir.create("./results")

fileInCredReg = "./inputs/input_cadd.txt"

fileCadd = "/stk05236/cadd/v1.6_GRCh37/whole_genome_SNVs_inclAnno.tsv.gz"
fileCaddIndex = "/stk05236/cadd/v1.6_GRCh37/whole_genome_SNVs_inclAnno.tsv.gz.tbi"
blnWriteFullCadd = FALSE

fileOutBase = "./results/new_cadd16"
	
##########################################################################################################################################################################

fileOutAll = paste(fileOutBase,".out",sep="")
fileOutPhred15 = paste(fileOutBase,"phred15.out",sep="")
fileOutPhred15Short = paste(fileOutBase,"phred15_short.out",sep="")

tcr = read.table(fileInCredReg, sep="\t", header = T, stringsAsFactors=F)
  # locus_id signal_id               chr_pos
#    n100        1 8:134331059-134332960
#    n10         1  14:73388687-73940636
#    n10         2  14:73341942-73372021

for(i in 1:nrow(tcr)) {

	print(paste(i, "of", nrow(tcr)))

	lid = paste(tcr$locus_id[i],tcr$signal_id[i],sep=".")
	cvreg = tcr$chr_pos[i]
	chr = strsplit(cvreg,":",fixed=T)[[1]][1]
	posreg = strsplit(cvreg,":",fixed=T)[[1]][2]
	pos1 = strsplit(posreg,"-",fixed=T)[[1]][1]
	pos2 = strsplit(posreg,"-",fixed=T)[[1]][2]

	fileOutAll = paste(fileOutBase,".out",sep="")
	fileOutPhred15 = paste(fileOutBase,"phred15.out",sep="")

	if(blnWriteFullCadd) {
		if(i==1) {
			if(file.exists(fileOutAll)) file.remove(fileOutAll)
			system(paste("zcat ",fileCadd," | head -n 2 | tail -1 | sed 's/$/\t'lid'/' > ", fileOutAll ,sep=""))
		}
		system(paste("tabix ",fileCadd," ",fileCaddIndex," ", chr,":",pos1,"-",pos2, " | sed 's/$/\t'",lid,"'/' >> ", fileOutAll ,sep=""))
	} 	
	if(i==1) {
		if(file.exists(fileOutPhred15)) file.remove(fileOutPhred15)	
		system(paste("zcat ",fileCadd," | head -n 2 | tail -1 | sed 's/$/\t'lid'/' > ", fileOutPhred15 ,sep=""))
	}
	system(paste("tabix ",fileCadd," ",fileCaddIndex," ", chr,":",pos1,"-",pos2, " | awk '{ if($9 >= 5) { print $0}}' | sed 's/$/\t'",lid,"'/' >> ", fileOutPhred15 ,sep=""))
}

## get relevant columns: 
fileout15ishort = paste(fileOutBase,"phred15_short.out",sep="")
system(paste("awk '{print $1,$2,$3,$4,$5,$8,$9,$10,$17,$18,$19,$20,$21,$29,$116,$117}' ",fileOutPhred15," > ",fileOutPhred15Short,sep=""))
# 1-5, 8-10, 19-21, 116-117, additional coulmns hopefully aminoacid change and position

setwd("/stk05236/")