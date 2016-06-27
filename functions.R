### this script contains functions defined for the pipeline ###
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# this function selects typical cells, and output selected typical cells as a matrix
# index is the index of measurements based on which typical cells are selected, should be a vector of integers. 
typical_cell_selection <- function(input,number,index) {
	selected <- matrix(,number,length(index))
	to_use <- input[,index]
	median_data <- apply(scale(to_use),2,mean)
	pos_dis <- apply(scale(to_use),1,function(x)(sum(abs(x-median_data))))
	rank_dist <- rank(pos_dis,ties.method="first")
	to_select <- which(rank_dist<=number)
	selected <- as.matrix(to_use[to_select,])
	colnames(selected) <- colnames(input)[index]
	return(selected)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
appends <- function(lists,to_add) {
	x <- lists
	x[[length(x)+1]] <- to_add
	return(x)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combine_data <- function(directory) {
	library(data.table)
	infiles <- list.files(directory)
	outdata <- list()
	for (i_file in 1:length(infiles)) {
		to_add <- read.table(paste0(directory,infiles[i_file]),header=T,colClasses="numeric",sep="\t")
		outdata <- appends(outdata,to_add)
	}
	outdata <- rbindlist(lapply(outdata,as.data.frame))
	return(outdata)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### this function performs SVM on designated datasets defined by user ######
### this function returns the accuracy of SVM classification, and the direction of classification boundary
### healthy: matrix of typical healthy cells to be classified
### disease: matrix of typical diseased cells to be classified
### feature_index: index of measurements used in classification, should be a vector of integers, length must be larger than 1
SVM <- function(healthy,disease,feature_index=c(1:ncol(healthy))) {
	library(kernlab)
	library(ks)
	AvsB <- matrix(1,nrow(healthy),1)
	healthy <- cbind(healthy,AvsB)
	AvsB <- matrix(-1,nrow(disease),1)
	disease <- cbind(disease,AvsB)
	myTable <- rbind(healthy,disease)
	AvsB <- factor(myTable[,ncol(myTable)])
	top_measures <- feature_index

	n_meas <- length(top_measures)
	x <- as.matrix(myTable[,top_measures])
	center <- apply(x,2,mean)
	std <- apply(x,2,sd)
	x <- scale(x)
	myTable2 <- data.frame(x,AvsB)
	myModel <- ksvm(AvsB ~ ., data=myTable2,type="C-svc", kernel="vanilladot", C = 10, prob.model=TRUE)

	SVM_coeff <- alpha(myModel)[[1]]
	SVM_coeff2 <- coef(myModel)[[1]]
	SVM_index <- alphaindex(myModel)[[1]]
	weight <- rep(0,each=ncol(x))
	for (i in 1:length(SVM_coeff)) {
		weight <- weight + SVM_coeff2[i]*x[SVM_index[i],]
	}
	norm <- sqrt(sum(weight**2))
	weightnorm <- weight/norm
	label <- sign(SVM_coeff/SVM_coeff2)
	SVM_b <- 0
	for (i in 1:length(SVM_coeff)) {	
		SVM_b <- SVM_b + (sum(weight*x[SVM_index[i],])-label[i])
	}
	SVM_b <- SVM_b/length(SVM_coeff)
	SVM_bn <- SVM_b/norm

	dist <- rep(0,each=nrow(x))
	for (i in 1:nrow(x)) {
		dist[i] <- sum(weightnorm*x[i,])-SVM_bn
	}
	accuracy <- 100*sum(dist*as.numeric(as.character(AvsB))>0)/nrow(x)
	# gap size
	n_pat_class_1 <- sum(AvsB == levels(AvsB)[2])
	gap_size <- min(dist[1:n_pat_class_1])-max(dist[(n_pat_class_1+1):nrow(x)])
	results <- list(accuracy, weightnorm, center, std)
	names(results) <- c("accuracy","weightnorm", "center", "std")
	return(results)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# to_be_normed: matrix containing single cells information to be normalized
# center: center of typical cells
# std: standard deviations of typical cells
# weight: direction of boundary plane given by SVM(), weight should be formatted as a n by 1 matrix, n is the number of dimension
normalization <- function(to_be_normed, center, std, weight) {
	scaled <- scale(to_be_normed,center=center,scale=std)
	normalized <- to_be_normed%*%weight
	percent <- 100*sum(normalized>0)/nrow(to_be_normed)
	return(list(normalized,percent))
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
get_normalized <- function(directory,index,center,std,weight) {
	infiles <- list.files(directory)
	output <- c(1:length(infiles))
	cell_count <- c(1:length(infiles))
	name <- c(1:length(infiles))
	for (i in 1:length(infiles)) {
		name[i] <- sub(".txt","",infiles[i])
		indata <- read.table(paste0(directory,infiles[i]),header=T,colClasses="numeric",sep="\t")
		to_use <- indata[,index]
		normed <- normalization(to_use, center, std, weight)
		output[i] <- normed[[2]]
		cell_count[i] <- nrow(indata)
	}
	return(list(output,cell_count,name))
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
FPR <- function(avg,std,threshold,number_of_plates) {
	value <- pt((threshold-avg)*sqrt(number_of_plates)/std,df=number_of_plates-1)
	return(value)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hits_identification <- function(disease_percent,siRNA_percent,std_siRNA_percent,healthy_number,siRNA_number,siRNA_name,number_of_plates) {
	threshold <- mean(disease_percent) + 5*sd(disease_percent)
	number_ratio <- siRNA_number/apply(healthy_number,2,median)
	number_adjusted <- apply(number_ratio,1,mean) + 2*apply(number_ratio,1,sd)
	pre_selected <- which(siRNA_percent>threshold)
	number_selected <- which(siRNA_number>=0.5)
	pre_hit <- intersect(pre_selected,number_selected)
	screened <- rep(0,length(pre_hit))
	index <- 1
	for (i in pre_hit) {
		fdr <- FPR(siRNA_percent[i],std_siRNA_percent[i],threshold,number_of_plates)
		if (fdr<0.2) {screened[index] <- 1}
		index <- index + 1
	}
	selected <- pre_hit[which(screened)]
	hits <- siRNA_name[selected]
	return(hits)
}
#-------------------------------------------------------------- clustering and output siRNAs inside each cluster ----------------------------------------------------------------------------
cluster <- function(avg_siRNA,outdir,siRNA_name,k,master) {
	gene_name <- read.table(master,header=T,sep="\t")
	clusters <- kmeans(avg_siRNA,centers=k,nstart=50)
	output <- paste0(outdir,"siRNA clusters Gene ID.txt")
	output2 <- paste0(outdir,"siRNA clusters symbol.txt")
	for (i_k in 1:k) {
		gene <- c()
		symbol <- c()
		siRNA_name <- rownames(RNA)[which(clusters$cluster==i_k)]
		for (i_RNA in 1:length(siRNA_name)) {
			gene <- c(gene,gene_name$GeneID[grep(paste0(siRNA_name[i_RNA],"$"),as.character(gene_name$GeneSymbol))])
			symbol <- c(symbol,as.character(gene_name$GeneSymbol)[grep(paste0(siRNA_name[i_RNA],"$"),gene_name$GeneSymbol)])
		}
		write(file=output,gene,ncolumns=length(gene),sep="\t",append=T)
		write(file=output2,symbol,ncolumns=length(symbol),sep="\t",append=T)
	}
}
#--------------------------------------------------------------- correlation analysis -------------------------------------------------------------------------------------------------------
pre_correlation <- function(directory, channels, store_dir, centers, stds, weightss) {
	infiles <- list.files(directory)
	multiple <- matrix(,length(infiles),(length(channels)+1))
	num_channel <- length(channels)
	name <- c(1:length(infiles))
	for (i in 1:length(infiles)) {
		indata <- read.table(paste0(directory,infiles[i]),header=T,colClasses="numeric",sep="\t")
		name[i] <- sub(".txt","",infiles[i])
		output <- matrix(,nrow(indata),length(channels))
		for (i_channel in 1:length(channels)) {
			index <- channels[[i_channel]]
			to_use <- indata[,index]
			center <- centers[[i_channel]]
			std <- stds[[i_channel]]
			weight <- weightss[[i_channel]]
			transformed <- normalization(to_use, center, std, weight)
			output[,i_channel] <- transformed[[1]]
		}
		colnames(output) <- names(channels)
		rownames(output) <- name
		write.table(output,paste0(store_dir,infiles[i]),append=F,row.names=F,sep="\t") # store transformed data
		to_calculate <- output>0
		sums <- apply(to_calculate,1,sum)
#--------------------------------------------------------------calculate percentage of cells healthy-like in multiple channels --------------------------------------------------------------
		for (j_channel in 1:length(channels)) {
			to_add <- 100*length(sums==(num_channel-1)&(!to_calculate[,j_channel]))/nrow(indata)
			multiple[i,j_channel] <- to_add
		}
		multiple[i,ncol(multiple)] <- 100*length(sums==num_channel)/nrow(indata)
	}
	colnames(multiple) <- c(paste("healthy in all channels except",names(channels)),"healthy in all channels")
	rownames(multiple) <- name
	return(multiple)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
correlation_analysis <- function(healthy_dir,disease_dir,siRNA_dir,channels,store_dir,center,std,weight) {
	library(GGally)
	healthy <- list()
	disease <- list()
	siRNA <- list()
	for (i in 1:length(healthy_dir)) {
		h_dir <- healthy_dir[i]
		d_dir <- disease_dir[i]
		s_dir <- siRNA_dir[i]
		centers <- center[[i]]
		stds <- std[[i]]
		weightss <- weight[[i]]
		h_add <- pre_correlation(h_dir, channels, store_dir, centers, stds, weightss)
		healthy <- appends(healthy,h_add)
		d_add <- pre_correlation(d_dir, channels, store_dir, centers, stds, weightss)
		disease <- appends(disease,d_add)
		s_add <- pre_correlation(s_dir, channels, store_dir, centers, stds, weightss)
		siRNA <- appends(siRNA,s_add)
	}
	avg_healthy <- apply(simplify2array(healthy), 1:2, mean)
	avg_disease <- apply(simplify2array(disease), 1:2, mean)
	avg_siRNA <- apply(simplify2array(siRNA), 1:2, mean)
	total_data <- rbind(avg_healthy,avg_disease)
	total_data <- rbind(total_data,avg_siRNA)
	correlation <- apply(total_data,2,function(x) cor(x,total_data[,ncol(total_data)]))
	correlation <- correlation[-ncol(correlation)]
	correlation <- as.matrix(correlation)
	rownames(correlation) <- paste(colnames(total_data)[1:length(channels)],colnames(total_data)[ncol(total_data)],sep="_")
	return(list(total_data,correlation))
}