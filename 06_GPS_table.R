library(data.table)

setwd("/stk05236/")

new = fread( "./lmm/all_loci/04_candidate_genes/lead_genes_annotated.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

genes = fread( "./inputs/genes/glist-hg19.tw.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

real_genes = read.table("./inputs/synonym_table.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

real_genes = real_genes[which(!is.na(real_genes$approved_symbol)),]

real_genes = real_genes[which(!is.na(real_genes$genes)),]

for(i in 1:nrow(genes)){

	if(!any(genes$Gene[i]==real_genes$approved_symbol) & any(genes$Gene[i]==real_genes$genes)){
	
		genes$Gene[i]=real_genes$approved_symbol[which(real_genes$genes==genes$Gene[i])[1]]
	}
}



#GPS = data.frame(Locus/signal_name=character(), locus_signal_no=integer(), Gene=character(), distance_to_lead_variant=numeric(),Chromosome=integer(),Start_of_gene=numeric(),End_of_gene=numeric())

Locus <- c("A")

genes_in_region <-c("A")

locus_id <-0

chr <- 0

Start_of_gene<-0

End_of_gene<-0

GPS <- data.frame(Locus,locus_id,genes_in_region,chr,Start_of_gene,End_of_gene,stringsAsFactors=FALSE)


for (i in 1:nrow(new)){

	genes_in_region <- unlist(strsplit(new$genes_in_region[i],","))
		if(new$genes_in_region[i]==FALSE){
		genes_in_region = new$Indep.500k.5e8.NearestGene[i]
		}
	Locus <- rep(new$Indep.500k.5e8.NearestGene[i],length(genes_in_region))
	locus_id <- rep(new$Locus_id[i],length(genes_in_region))
	Start_of_gene <- as.numeric(unlist((strsplit(new$genes_start[i],","))))
		if(new$genes_in_region[i]==FALSE){
		Start_of_gene <- as.numeric(genes$Pos1[which(genes_in_region==genes$Gene)])
		}
	chr <- unlist(strsplit(new$genes_chr[i],","))
		if(new$genes_in_region[i]==FALSE){
		chr <- genes$Chr[which(genes_in_region==genes$Gene)]
		}
	
	End_of_gene <- as.numeric(unlist(strsplit(new$genes_end[i],",")))
		if(new$genes_in_region[i]==FALSE){
		End_of_gene <- as.numeric(genes$Pos2[which(genes_in_region==genes$Gene)])
		}
		
	x <- data.frame(Locus,locus_id,genes_in_region,chr,Start_of_gene,End_of_gene, stringsAsFactors=FALSE)
	GPS <- rbind(GPS,x)
	genes_in_region=vector()
	print(i)
}

GPS <- GPS[-1,]

names(GPS)<- c("Locus_name","locus_id","Gene","Chr","Start_of_gene","End_of_gene")


distance_to_lead_variant <- vector()

for (i in 1:nrow(GPS)){

			
	if(new$pos[which(new$Locus_id==GPS$locus_id[i])]<GPS$Start_of_gene[i]){
		
		distance_to_lead_variant[i] = GPS$Start_of_gene[i]-new$pos[which(new$Locus_id==GPS$locus_id[i])]
	}else{
		
		if(new$pos[which(new$Locus_id==GPS$locus_id[i])]<=GPS$End_of_gene[i]){
			
			distance_to_lead_variant[i] =0
		}else{
			
			distance_to_lead_variant[i]=GPS$End_of_gene[i]-new$pos[which(new$Locus_id==GPS$locus_id[i])]
		}
	}
	
}

GPS['distance_to_lead_variant'] <- distance_to_lead_variant	

GPS <- GPS[,c(1:3,7,4:6)]

any(duplicated(GPS))

write.table(GPS, file="./lmm/all_loci/06_GPS_table/table_for_GPS.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)