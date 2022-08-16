setwd("/stk05236/clean_GPS")

files = list.files("./07_cred_var/cred_var_by_id")

all_cred_95 = c()

all_cred_99 = c()

region_id = c()

signal_id = c()

for (y in 1:length(files)){

	input_file = paste("./05_cred_var/cred_var_by_id/",files[y],sep="")
	
	tin = read.table(as.character(input_file),  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	
	cred_95 = tin[which(tin$ci95==TRUE),]
	
	cred_99 = tin[which(tin$ci99==TRUE),]
	
	name=unlist(strsplit(files[y],".",fixed=TRUE))
	
	
	
	########in eine Tabelle
	
	id=unlist(strsplit(name,"_"))
	
	##### id[4] ist die region_id
	
	region_id[y] = id[4]
	
	names(cred_95) = c("rsid","chr","pos","ea","eaf","beta","se","p","n","ppa","cppa","ci95","ci99")

	names(cred_99) = c("rsid","chr","pos","ea","eaf","beta","se","p","n","ppa","cppa","ci95","ci99")
	
	
	
	if(id[6]=="cond"){
	
		signal_id[y]=id[5]
	}else{signal_id[y] = 1}

	
	cred_95['region_id']= rep(region_id[y],nrow(cred_95))
	
	cred_95['signal_id']= rep(signal_id[y],nrow(cred_95))
	
	cred_99['region_id']= rep(region_id[y],nrow(cred_99))
	
	cred_99['signal_id']= rep(signal_id[y],nrow(cred_99))
	
	if (length(nrow(all_cred_95))==0){ 
		
		all_cred_95 = cred_95
	}else{
	
		all_cred_95=rbind(all_cred_95,cred_95)
		
	}
	
	if (length(nrow(all_cred_99))==0){ 
		
		all_cred_99 = cred_99
	}else{
	
		all_cred_99=rbind(all_cred_99,cred_99)
		
	}
	
	print(y)
	
}

}
write.table(all_cred_95, file="./05_cred_var/all_cred_var_95.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
write.table(all_cred_99, file="./05_cred_var/all_cred_var_99.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
