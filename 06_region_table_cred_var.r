setwd("/stk05236/")

region_table = read.table("./04_locus_table/locus_table.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cred_var_95 = read.table("./05_cred_var/all_cred_var_95.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cred_var_99 = read.table("./05_cred_var/all_cred_var_99.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE)


no_cred_var_95 = c()

no_cred_var_99 = c()

region_table_signals = region_table

signal_id_neu = c()

for (i in 1:nrow(region_table)){

	
	if (any(region_table_signals$Locus_id[i]==cred_var_95$region_id & cred_var_95$signal_id >1)){
	
		idx = which(cred_var_95$region_id==region_table_signals$Locus_id[i])
		
		signals_in_region = max(cred_var_95$signal_id[idx])
		
		no_cred_var_95[i] = length(which(region_table_signals$Locus_id[i]==cred_var_95$region_id & cred_var_95$signal_id ==1))

		no_cred_var_99[i] = length(which(region_table_signals$Locus_id[i]==cred_var_99$region_id & cred_var_99$signal_id ==1))
	
		signal_id_neu[i] = 1
		
		for (y in 2:signals_in_region){
		
			region_table_signals[nrow(region_table_signals)+1,]=region_table_signals[i,]
		
			no_cred_var_95[nrow(region_table_signals)] = length(which(region_table_signals$Locus_id[i]==cred_var_95$region_id & cred_var_95$signal_id ==y))
		
			no_cred_var_99[nrow(region_table_signals)] = length(which(region_table_signals$Locus_id[i]==cred_var_99$region_id & cred_var_99$signal_id ==y))
		
			signal_id_neu[nrow(region_table_signals)]=y
		}
		
	}else{
	
	no_cred_var_95[i] = length(which(region_table_signals$Locus_id[i]==cred_var_95$region_id & cred_var_95$signal_id ==1))

	no_cred_var_99[i] = length(which(region_table_signals$Locus_id[i]==cred_var_99$region_id & cred_var_99$signal_id ==1))
	
	signal_id_neu[i] = 1
	}
}

region_table_signals['no_cred_var_95'] <- no_cred_var_95

region_table_signals['no_cred_var_99'] <- no_cred_var_99

region_table_signals['signal_id'] <- signal_id_neu

write.table(region_table_signals, file="./06_signal_table/signal_table_cred_var.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)