########################################################################################################
## Function to calculate Posterior Probablities of Association (PPA)  from betas and SEs              ##
##																									  ##
## Original Publication of method and script:														  ##
## Wakefield J. (2007). A Bayesian measure of the probability of false discovery in genetic 		  ##
## epidemiology studies. American journal of human genetics, 81(2), 208–227. 						  ##
## https://doi.org/10.1086/519024																	  ##
########################################################################################################

# calc the ppa
fnPPAWakefield=function(beta,se,W) { 
	
	### W i
	
	ABF = sqrt((se^2 + W)/(se^2))*exp(-((beta/se)^2)*W/(2*(se^2+W)))
	
	## correct for prior direction (Wakefield has P(b|H0)/P(b|H1))
	## but we need P(b|H1)/P(b|H0)
	ABF = 1/ABF
	
	PPA = ABF/sum(ABF)
	
	return(PPA)
	
}


setwd("/stk05236/clean_GPS")
dir.create("./05_cred_var/cred_var_by_id")

path_to_cond = "/stk05236/credible/inputs_lmm"
files = list.files(path_to_cond)

known = c()

for (i in 1:length(files)){

	a = unlist(strsplit(files[i],"_"))
	
	if(unlist(strsplit(a[5],split=""))[1]=="k") {
		known[i] = TRUE
	}else{
		known[i] = FALSE
	}
}
	
files = files[which(known)]


for (y in 1:length(files)){ 


	input_file = paste(path_to_cond,"/",files[y],sep="")

	tin = read.table(as.character(input_file),  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

	names(tin)=tolower(names(tin))


# remove missings
	if (any(is.na(tin$betax))){
		tin = tin[-which(is.na(tin$betax)),]
	}

	W = 0.005^2

	tin$ppa = fnPPAWakefield(tin$betax,tin$sex,W)
# W is prior for variance of true genetic effect estimates, bg~N(0,W)

# order by decresing PP
	tin = tin[order(tin$ppa, decreasing=T), ]

# sum up PPs

	tin$cppa = rep(NA,nrow(tin))
	tin$ci95 = rep(FALSE,nrow(tin))
	tin$ci99 = rep(FALSE,nrow(tin))

	tin$cppa[1] = tin$ppa[1]
	tin$ci95[1] = TRUE
	tin$ci99[1] = TRUE


	x95=0
	x99=0

	for(i in 2:nrow(tin)) {
			tin$cppa[i] = tin$cppa[i-1] + tin$ppa[i]
		if(tin$cppa[i] < 0.95) tin$ci95[i] = TRUE
		if(tin$cppa[1] > 0.95) x95=1
		if(tin$cppa[i] > 0.95){
			if(x95==0){
				tin$ci95[i] = TRUE
				x95=1
			}
		}

		if(tin$cppa[i] < 0.99) tin$ci99[i] = TRUE
		if(tin$cppa[1] > 0.99) x99=1
		if(tin$cppa[i] > 0.99){
			if(x99==0){		
				tin$ci99[i] = TRUE
				x99=1
			}
		}

	}
	

# write result
	name=unlist(strsplit(files[y],".",fixed=TRUE))
	
	output_file = paste("./05_cred_var/cred_var_by_id/",name[2],"_cred_var.txt",sep="")
	
	
	write.table(tin, file=as.character(output_file), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
