### main script to run the pipeline ###
rm(list=ls())
source("\\code\\functions.R") # define functions used in this script
number_of_plates <- 4
number_of_healthy_controls <- 12
number_of_disease_controls <- 12
number_of_siRNA <- 320
#----------------------------------------------------------------- set directories and global parameters --------------------------------------------------------------------------------------
healthy_control_dir <- paste0("C:\\Data\\plates\\plate ",c(1:number_of_plates),"\\healthy control\\")
disease_control_dir <- paste0("C:\\Data\\plates\\plate ",c(1:number_of_plates),"\\disease control\\")
siRNA_dir <- paste0("C:\\Data\\plates\\plate ",c(1:number_of_plates),"\\siRNA\\")
output_dir <- "C:\\Data\\output\\"
store_dir <- paste0("C:\\Data\\plates\\plate ",c(1:number_of_plates),"\\transformed\\")
hits <- list()
#----------------------------------------------------------------- set parameters for pipeline-------------------------------------------------------------------------------------------------
DAPI <- c(1:12)
lamin <- c(13:15)
progerin <- c(16:18)
gamma <- c(19:21)
channels <- list(DAPI,lamin,progerin,gamma)
names(channels) <- c("DAPI","lamin","progerin","gamma")
number_of_typical_cells <- 300
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------- define output matrices ----------------------------------------------------------------------------------------------------
print("Begin analysis...")
avg_healthy <- matrix(,number_of_healthy_controls,length(channels)) # percentage of healthy-like cells in each healthy control sample, averaged over all replicate plates
avg_disease <- matrix(,number_of_disease_controls,length(channels))
avg_siRNA <- matrix(,number_of_siRNA,length(channels))
std_siRNA <- matrix(,number_of_siRNA,length(channels))
centers <- list()
stds <- list()
weightss <- list()
names(centers) <- names(channels)
names(stds) <- names(channels)
names(weightss) <- names(channels)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for (i_channel in 1:length(channels)) {
	channel <- names(channels)[i_channel]
	percent_healthy <- matrix(,number_of_healthy_controls,number_of_plates) # percentage of healthy-like cells in i_channel th channel for each replicate plate
	percent_disease <- matrix(,number_of_disease_controls,number_of_plates)
	percent_siRNA <- matrix(,number_of_siRNA,number_of_plates)
	number_healthy <- matrix(,number_of_healthy_controls,number_of_plates) # number of cells in each sample for each replicate plate
	number_disease <- matrix(,number_of_disease_controls,number_of_plates)
	number_siRNA <- matrix(,number_of_siRNA,number_of_plates)
	to_center <- list()
	to_std <- list()
	to_weight <- list()
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------ select typical cells ------------------------------------------------------------------------------------------------------
	for (i_plate in 1:number_of_plates) {
		healthy <- combine_data(healthy_control_dir[i_plate])
		disease <- combine_data(disease_control_dir[i_plate])
		typical_healthy <- typical_cell_selection(healthy,number_of_typical_cells,get(channel))
		typical_disease <- typical_cell_selection(disease,number_of_typical_cells,get(channel))
#------------------------------------------------------------------ find SVM classification boundary ----------------------------------------------------------------------------------------
		classificatoin <- SVM(healthy,disease,feature_index)
		boundary_direction <- classification$weightnorm
#-------------------------------------------------------- project individual cell to the direction of classification boundary ---------------------------------------------------------------
		normed_healthy <- get_normalized(healthy_control_dir[i_plate],get(channel),classification$center,classification$std,boundary_direction)
		normed_disease <- get_normalized(disease_control_dir[i_plate],get(channel),classification$center,classification$std,boundary_direction)
		normed_siRNA <- get_normalized(siRNA_dir[i_plate],get(channel),classification$center,classification$std,boundary_direction)
		to_center <- appends(to_center,classification$center)
		to_std <- appends(to_std,classification$std)
		to_weight <- appends(to_weight,boundary_direction)
		percent_healthy[,i_plate] <- normed_healthy[[1]]
		percent_disease[,i_plate] <- normed_disease[[1]]
		percent_siRNA[,i_plate] <- normed_siRNA[[1]]
		number_healthy[,i_plate] <- normed_healthy[[2]]
		number_disease[,i_plate] <- normed_disease[[2]]
		number_siRNA[,i_plate] <- normed_siRNA[[2]]
		siRNA_name <- normed_siRNA[[3]]
	}
	names(to_center) <- paste("plate",c(1:number_of_plates),sep="_")
	names(to_std) <- paste("plate",c(1:number_of_plates),sep="_")
	names(to_weight) <- paste("plate",c(1:number_of_plates),sep="_")
	avg_healthy[,i_channel] <- apply(percent_healthy,1,mean)
	avg_disease[,i_channel] <- apply(percent_disease,1,mean)
	avg_siRNA[,i_channel] <- apply(percent_siRNA,1,mean)
	std_siRNA[,i_channel] <- apply(percent_siRNA,1,sd)
#-------------------------------------------------------------------------------- select hits -----------------------------------------------------------------------------------------------
	to_hits <- hits_identification(avg_disease[,i_channel],avg_siRNA[,i_channel],std_siRNA[,i_channel],number_healthy,number_siRNA,siRNA_name,number_of_plates)
	hits <- appends(hits,to_hits)
	centers <- appends(centers,to_center)
	stds <- appends(stds,to_std)
	weightss <- appends(weightss,to_weight)
}
#-------------------------------------------------------------------------------- output results --------------------------------------------------------------------------------------------
colnames(avg_healthy) <- names(channels)
colnames(avg_disease) <- names(channels)
colnames(avg_siRNA) <- names(channels)
colnames(std_siRNA) <- names(channels)
colnames(number_healthy) <- paste("plate",c(1:number_of_plates),sep="_")
colnames(number_disease) <- paste("plate",c(1:number_of_plates),sep="_")
colnames(number_siRNA) <- paste("plate",c(1:number_of_plates),sep="_")
rownames(avg_siRNA) <- siRNA_name
rownames(std_siRNA) <- siRNA_name
rownames(number_siRNA) <- siRNA_name
names(hits) <- names(channels)
write.table(avg_healthy,paste0(output_dir,"percentage of healthy-like cells in healthy controls.txt"),append=F,row.names=T,sep="\t")
write.table(avg_disease,paste0(output_dir,"percentage of healthy-like cells in disease controls.txt"),append=F,row.names=T,sep="\t")
write.table(avg_siRNA,paste0(output_dir,"percentage of healthy-like cells in siRNA samples.txt"),append=F,row.names=T,sep="\t")
write.table(std_siRNA,paste0(output_dir,"std of percentage of healthy-like cells in siRNA samples.txt"),append=F,row.names=T,sep="\t")
write.table(number_healthy,paste0(output_dir,"number of cells in healthy controls.txt"),append=F,row.names=T,sep="\t")
write.table(number_disease,paste0(output_dir,"number of cells in disease controls.txt"),append=F,row.names=T,sep="\t")
write.table(number_siRNA,paste0(output_dir,"number of cells in siRNA samples.txt"),append=F,row.names=T,sep="\t")
write.table(hits,paste0(output_dir,"selected siRNA hits.txt"),append=F,row.names=F,sep="\t")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------optional sample level analysis----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------clustering analysis--------------------------------------------------------------------------------------------
cluster(avg_siRNA,output_dir,siRNA_name,k,master)
# master is the directory to the list of corresponding Gene ID and symbols. Gene IDs should be stored under colname "GeneID", and Gene symbols should be named "GeneSymbol"
#-----------------------------------------------------------------------------correlation analysis-------------------------------------------------------------------------------------------
correlation <- correlation_analysis(healthy_control_dir,disease_control_dir,siRNA_dir,channels,store_dir,centers,stds,weightss)
write.table(correlation[[1]],paste0(output_dir,"percentage of healthy-like cells in multiple channels.txt"),append=F,row.names=T,sep="\t")
write.table(correlation[[2]],paste0(output_dir,"correlation of percentage of healthy-like cells in multiple channels.txt"),append=F,row.names=T,sep="\t")