library(data.table)

setwd("/stk05236/clean_GPS/")

meta = fread( "./03_candidate_genes/lead_genes_annotated.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)



region_start_250kb = meta$region.start

region_end_250kb = meta$region.end

locus_size = region_end_250kb - region_start_250kb

region_table = data.frame(meta$Locus_id,meta$chr,region_start_250kb,region_end_250kb,locus_size, meta$Indep.500k.5e8.aNumLocusSNPs,meta$n_genes_in_region,meta$genes_in_region,stringsAsFactors = FALSE)

names(region_table) = c("Locus_id","chr","region_start_250kb","region_end_250kb","locus_size","num_sig_SNP","num_genes_in_region","genes_in_region")


write.table(region_table, file="./04_locus_table/locus_table.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)