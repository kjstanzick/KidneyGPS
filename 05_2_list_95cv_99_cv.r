setwd("/stk05236/clean_GPS")

files = list.files("./05_cred_var/cred_var_by_id")

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
	
	##### id[3] ist die region_id
	
	region_id[y] = id[3]
	
	names(cred_95) = c("rsid","chr","pos","ea","eaf","beta","se","p","n","ppa","cppa","ci95","ci99")

	names(cred_99) = c("rsid","chr","pos","ea","eaf","beta","se","p","n","ppa","cppa","ci95","ci99")
	
	
	
	if(id[5]=="cond"){
	
		signal_id[y]=id[4]
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



ref=fread("/stk05236/lmm/03_eval/european/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",  header = TRUE, sep = "\t", stringsAsFactors = FALSE,  data.table = FALSE)

all_cred_95_m=merge(all_cred_95, ref[,c(1,2,3,6,7,8,15)], by.x="rsid", by.y="RSID", all.x=TRUE, all.y=F,  sort=FALSE)
all_cred_95_m$ea=toupper(all_cred_95_m$ea)
all_cred_95_m$Allele1=toupper(all_cred_95_m$Allele1)
all_cred_95_m$Allele2=toupper(all_cred_95_m$Allele2)
#ea und beta_cond auf eGFR lowering drehen
not_lowering=which(all_cred_95_m$beta>0)

for(i in not_lowering){
	all_cred_95_m$beta[i]=(-1)*all_cred_95_m$beta[i]
	all_cred_95_m$eaf[i]=1-all_cred_95_m$eaf[i]
	temp_allele=c(all_cred_95_m$Allele1[i],all_cred_95_m$Allele2[i])
	all_cred_95_m$ea[i]=temp_allele[which(temp_allele!=all_cred_95_m$ea[i])]
}

#Allele1, Allele2 und Effect uncond auf ea anpassen

for(i in 1:nrow(all_cred_95_m)){
	if(all_cred_95_m$Allele1[i]!=all_cred_95_m$ea[i]){
		all_cred_95_m$Effect[i]=(-1)*all_cred_95_m$Effect[i]  
		all_cred_95_m$Allele1[i] <- paste(all_cred_95_m$Allele1[i],all_cred_95_m$Allele2[i], sep = "")
		all_cred_95_m$Allele2[i] <- substr(all_cred_95_m$Allele1[i],0,nchar(all_cred_95_m$Allele1[i]) - nchar(all_cred_95_m$Allele2[i]))
		all_cred_95_m$Allele1[i] <- substr(all_cred_95_m$Allele1[i],nchar(all_cred_95_m$Allele2[i]) + 1, nchar(all_cred_95_m$Allele1[i]))
	}
}

all_cred_99_m=merge(all_cred_99, ref[,c(1,2,3,6,7,8,15)], by.x="rsid", by.y="RSID", all.x=TRUE, all.y=F, sort=FALSE)
all_cred_99_m$ea=toupper(all_cred_99_m$ea)
all_cred_99_m$Allele1=toupper(all_cred_99_m$Allele1)
all_cred_99_m$Allele2=toupper(all_cred_99_m$Allele2)
#ea und beta_cond auf eGFR lowering drehen
not_lowering=which(all_cred_99_m$beta>0)

for(i in not_lowering){
	all_cred_99_m$beta[i]=(-1)*all_cred_99_m$beta[i]
	all_cred_99_m$eaf[i]=1-all_cred_99_m$eaf[i]
	temp_allele=c(all_cred_99_m$Allele1[i],all_cred_99_m$Allele2[i])
	all_cred_99_m$ea[i]=temp_allele[which(temp_allele!=all_cred_99_m$ea[i])]
}

#Allele1, Allele2 und Effect uncond auf ea anpassen

for(i in 1:nrow(all_cred_99_m)){
	if(all_cred_99_m$Allele1[i]!=all_cred_99_m$ea[i]){
		all_cred_99_m$Effect[i]=(-1)*all_cred_99_m$Effect[i]  
		all_cred_99_m$Allele1[i] <- paste(all_cred_99_m$Allele1[i],all_cred_99_m$Allele2[i], sep = "")
		all_cred_99_m$Allele2[i] <- substr(all_cred_99_m$Allele1[i],0,nchar(all_cred_99_m$Allele1[i]) - nchar(all_cred_99_m$Allele2[i]))
		all_cred_99_m$Allele1[i] <- substr(all_cred_99_m$Allele1[i],nchar(all_cred_99_m$Allele2[i]) + 1, nchar(all_cred_99_m$Allele1[i]))
	}
}

names(all_cred_95_m)[c(19:21)]=paste(names(all_cred_95_m[c(19:21)]),".unconditioned",sep="")
names(all_cred_99_m)[c(19:21)]=paste(names(all_cred_99_m[c(19:21)]),".unconditioned",sep="")

write.table(all_cred_95_m, file="./05_cred_var/all_cred_var_95.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
write.table(all_cred_99_m, file="./05_cred_var/all_cred_var_99.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
