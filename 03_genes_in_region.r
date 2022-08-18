library(data.table)

setwd("/stk05236/clean_GPS/")

## eGFRcrea lead variants
all_meta = fread( "/stk05236/lmm/03_eval/metal_eGFR_meta1.TBL.Indep.500k.5e8.indepX_edited_MHC-tw.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

######################
###Locus start und Locus end

Chr_Region = strsplit (all_meta$Indep.500k.5e8.aLocusCoordinates,":") # aufspalten Chromosom und Region

		#unlist macht einen langen Vektor draus
		#matrix macht aus dem Vektor eine Matrix mit ncol Spalten

Chr_Region <- matrix(data=(unlist(Chr_Region)), ncol=2, byrow=TRUE)
Chr_Region <- as.data.frame(Chr_Region, stringsAsFactors=FALSE)

names(Chr_Region)[2] <- "Region"

Region = strsplit(Chr_Region$Region, "_") #Aufspalten Regionanfang und Regionende
Region <- matrix(data=(unlist(Region)), ncol=2, byrow=TRUE)

all_meta['region.start'] <- as.numeric(Region[,1])-250000
all_meta['region.end'] <- as.numeric(Region[,2])+250000


######################
###Gene per locus +-250kb 
genes_sorted = read.table("./inputs/glist.txt", header=TRUE,sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

genes_in_region <- vector(length=nrow(all_meta))

n_genes_in_region <- rep(0,nrow(all_meta))

genes_chr <- vector(length=nrow(all_meta))

genes_start <- vector(length=nrow(all_meta))

genes_end <- vector(length=nrow(all_meta))

for (i in 1:nrow(all_meta)){

	for (x in 1:nrow(genes_sorted)){
	
		
			if(all_meta$chr[i]==genes_sorted$Chr[x]){
			
				if(all_meta$region.start[i] >= genes_sorted$Pos1[x]){
				
					if(all_meta$region.start[i]<=genes_sorted$Pos2[x]){
					
						if(genes_in_region[i]==FALSE){
							genes_in_region[i]=genes_sorted$Gene[x]
							genes_chr[i]=genes_sorted$Chr[x]
							genes_start[i]=as.numeric(genes_sorted$Pos1[x])
							genes_end[i]=as.numeric(genes_sorted$Pos2[x])
							n_genes_in_region[i] = n_genes_in_region[i]+1
						}else{
							genes_in_region[i]=paste(genes_in_region[i],genes_sorted$Gene[x],sep=",")
							genes_chr[i]=paste(genes_chr[i],genes_sorted$Chr[x],sep=",")
							genes_start[i]=paste(genes_start[i],genes_sorted$Pos1[x],sep=",")
							genes_end[i]=paste(genes_end[i],genes_sorted$Pos2[x],sep=",")
							n_genes_in_region[i] = n_genes_in_region[i]+1
						}
					}
					
				}else{
				
					if(all_meta$region.end[i]>= genes_sorted$Pos1[x]){
					
						if(genes_in_region[i]==FALSE){
							genes_in_region[i]=genes_sorted$Gene[x]
							genes_chr[i]=genes_sorted$Chr[x]
							genes_start[i]=as.numeric(genes_sorted$Pos1[x])
							genes_end[i]=as.numeric(genes_sorted$Pos2[x])
							n_genes_in_region[i] = n_genes_in_region[i]+1
						}else{
							genes_in_region[i]=paste(genes_in_region[i],genes_sorted$Gene[x],sep=",")
							genes_chr[i]=paste(genes_chr[i],genes_sorted$Chr[x],sep=",")
							genes_start[i]=paste(genes_start[i],genes_sorted$Pos1[x],sep=",")
							genes_end[i]=paste(genes_end[i],genes_sorted$Pos2[x],sep=",")
							n_genes_in_region[i] = n_genes_in_region[i]+1
						}
					}
				}
			
			}
		
	}
}



all_meta['genes_in_region'] <- genes_in_region
all_meta['genes_chr'] <- genes_chr
all_meta['genes_start'] <- genes_start
all_meta['genes_end'] <- genes_end
all_meta['n_genes_in_region'] <- n_genes_in_region

all_meta = all_meta[order(all_meta[,8]),] #p.value


write.table(all_meta, file="./03_candidate_genes/lead_genes_annotated.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
			