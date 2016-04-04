# This script runs normalization and LASSO for training sets of
# .0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 RNASeq and the rest
# microarray data. The testing data is either 100% microarray or
# 100% RNA seq. This script will output graphs of the accuracies
# and kappas

#DIRECTORIES
# Change paths below if you want to read from/write to another directory.
initial_data = "data/"
normalized_data = "normalized_data/"
output = "output/"
write(format(Sys.time(), "%b %d %Y %X"), stderr())



############ functions ##################

# This function takes a file of headers and a file containing a matrix and 
#filters only for columns with headers from the array
filter_matrix = function(old_matrix,header_list) {
	old_train = read.table(old_matrix, sep='\t', row.names=1, header=T,strip.white = TRUE)
	headers <- scan(file = header_list, what = "character", sep="\n")
	return(old_train[, headers])
}

############# Experiment 1 ################## 

## Step 1: Normalize the microArray and RNA Seq Input
target_file = "BRCARNASeq.pcl"
ref_file = "BRCAarray.pcl"
data_folder = initial_data
output_folder = normalized_data
prefix = "BRCA"

# COMMENT THIS PART OUT IF ALREADY NORMALIZED!!
system(paste("Rscript","normalize_data.R",target_file,ref_file,data_folder,output_folder,prefix))

## Step 2: Of the normalized output, cut out the experiments that don't have
## 		   samples in both microArray and RNASeq format

# Get the headers of all the samples that are in both microArray and RNA seq 
# format
array_header=paste0(data_folder,prefix,"_arrayHeaders.txt")
rnaSeq_header=paste0(data_folder,prefix,"_rnaSeqHeaders.txt")
combine_header=paste0(data_folder,prefix,"_headersInCommon.txt")

system(sprintf('head -n1 %s%s | tr "\t" "\n" | awk \'NR>1\'|sort > %s', data_folder,ref_file, array_header))
system(sprintf('head -n1 %s%s | tr "\t" "\n" | awk \'NR>1\'|sort > %s', data_folder,target_file, rnaSeq_header))
system(sprintf('comm -1 -2 %s %s| sed -e "s/-/./g"> %s',array_header,rnaSeq_header,combine_header))

# Get normalized tables for microARRAY and RNA seq containing only samples that
# have both already
input_dir = normalized_data
#1 Zero to one transform reference data. (Reference)
ref_input = 'BRCA_ZEROONE.pcl'
#2 Zero to one scale the TDM normalized data. (targets)
tdm_input = 'BRCA_TDM_ZEROONE.pcl'

#3 Zero to one scale the quantile normalized data (targets)
qn_input = 'BRCA_QN_ZEROONE.pcl'
qn_new_input='BRCA_QN_ZEROONE_filtered'
#4 Zero to one scale the LOG transformed data. (targets)
log_input = 'BRCA_LOG_ZEROONE.pcl'
#5 Zero to one scale the paranormal normalized data.(targets)
npn_input = 'BRCA_NPN_ZEROONE.pcl'
#6 Zero to one scale the paranormal normalized data.(Reference)
ref_npn_input = 'BRCA_REF_NPN_ZEROONE.pcl'
#7 zero-to-one scaling of untransformed data (targets)
un_input = 'BRCA_UN_ZEROONE.pcl'
output_dir = output

ref_new_input = filter_matrix(paste0(input_dir,ref_input),combine_header)
tdm_new_input = filter_matrix(paste0(input_dir,tdm_input),combine_header)
qn_new_input = filter_matrix(paste0(input_dir,qn_input),combine_header)
log_new_input = filter_matrix(paste0(input_dir,log_input),combine_header)
npn_new_input= filter_matrix(paste0(input_dir,npn_input),combine_header)
ref_npn_new_input=filter_matrix(paste0(input_dir,ref_npn_input),combine_header)
un_new_input= filter_matrix(paste0(input_dir,un_input),combine_header)

## Step 3: Seperate the samples into train samples and testing samples

#list the number of samples you plan to use for training and testing
#WARNING: Test_total must be equal or less to the number of total samples
total_samples <- length(readLines(combine_header))

train_total= 300 
test_total=total_samples - train_total

INITIAL_SEED = 2

set.seed(INITIAL_SEED)
seeds = sample(1:total_samples, total_samples, replace = FALSE)
test_seeds = seeds[1:test_total]
train_seeds = seeds[(test_total+1):total_samples]
ref_test_input = cbind(ref_new_input[, test_seeds])
ref_train_input = cbind(ref_new_input[, train_seeds])
tdm_test_input = cbind(tdm_new_input[, test_seeds])
tdm_train_input = cbind(tdm_new_input[, train_seeds])
qn_test_input = cbind(qn_new_input[, test_seeds])
qn_train_input = cbind(qn_new_input[, train_seeds])
log_test_input = cbind(log_new_input[, test_seeds])
log_train_input = cbind(log_new_input[, train_seeds])
npn_test_input = cbind(npn_new_input[, test_seeds])
npn_train_input = cbind(npn_new_input[, train_seeds])
ref_npn_test_input = cbind(ref_npn_new_input[, test_seeds])
ref_npn_train_input = cbind(ref_npn_new_input[, train_seeds])
un_test_input = cbind(un_new_input[, test_seeds])
un_train_input = cbind(un_new_input[, train_seeds])

#write names of testing files (with directory location)
tdm_test_RS_out2 = paste0(normalized_data,prefix,"_tdm_test.txt")
qn_test_RS_out2 = paste0(normalized_data,prefix,"_qn_test.txt")
log_test_RS_out2 = paste0(normalized_data,prefix,"_log_test.txt")
npn_test_RS_out2 = paste0(normalized_data,prefix,"_npn_test.txt")
un_test_RS_out2 = paste0(normalized_data,prefix,"_un_test.txt")
ref_test_MA_out2 = paste0(normalized_data,prefix,"_ref_test.txt")
ref_npn_test_MA_out2 = paste0(normalized_data,prefix,"_ref_npn_test.txt")

#print out testing file tables
write.table(tdm_test_input, file = tdm_test_RS_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
write.table(qn_test_input, file = qn_test_RS_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
write.table(log_test_input, file = log_test_RS_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
write.table(npn_test_input, file = npn_test_RS_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
write.table(un_test_input, file = un_test_RS_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
write.table(ref_test_input, file = ref_test_MA_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
write.table(ref_npn_test_input, file = ref_npn_test_MA_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)

#write names of testing files (without directory location) to use for LASSO script
tdm_test_RS_out1 = paste0(prefix,"_tdm_test.txt")
qn_test_RS_out1 = paste0(prefix,"_qn_test.txt")
log_test_RS_out1 = paste0(prefix,"_log_test.txt")
npn_test_RS_out1 = paste0(prefix,"_npn_test.txt")
un_test_RS_out1 = paste0(prefix,"_un_test.txt")
ref_test_MA_out1 = paste0(prefix,"_ref_test.txt")
ref_npn_test_MA_out1 = paste0(prefix,"_ref_npn_test.txt")


pOfRNASeq = c(.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
write(format(Sys.time(), "%b %d %Y %X"), stderr())
## Step 4: (1) Get pOfRNASeq[x] of the total RNA_Seq train samples and 
## 1 - pOfRNASeq[x] of the total microARRAY train samples and combine
## them into one file (2) Run Lasso on Trial,

for(x in pOfRNASeq){
  write(paste("Lasso",x,":"), stderr())
  write(format(Sys.time(), "%b %d %Y %X"), stderr())
	#initialize matricies
	tdm_train = matrix()
	qn_train = matrix()
	log_train = matrix()
	npn_train = matrix()
	un_train = matrix()

	if (x == 0){
		num_RNASeq_samples = 0
		num_mARRAYx_samples = train_total

		tdm_train = ref_train_input
		qn_train = ref_train_input
		log_train = ref_train_input
		npn_train = ref_npn_train_input
		un_train = ref_train_input
	} 
	if (x == 1){
		num_RNASeq_samples = train_total
		num_mARRAYx_samples = 0

		tdm_train = tdm_train_input
		qn_train = qn_train_input
		log_train = log_train_input
		npn_train = npn_train_input
		un_train = un_train_input
	} 
	if(x!=0 && x!=1){
		num_RNASeq_samples = round(x * train_total)
		num_mARRAYx_samples = train_total - num_RNASeq_samples
		tseeds = sample(1:train_total, train_total, replace = FALSE)
		RNASeq_seeds = tseeds[1:num_RNASeq_samples]
		mARRAY_seeds = tseeds[(num_RNASeq_samples+1):train_total]
		ref_RNA_seq = ref_train_input[,mARRAY_seeds]
		tdm_RNA_seq = tdm_train_input[,RNASeq_seeds]
		qn_RNA_seq = qn_train_input[,RNASeq_seeds]
		log_RNA_seq = log_train_input[,RNASeq_seeds]
		npn_RNA_seq = npn_train_input[,RNASeq_seeds]
		ref_npn_RNA_seq = ref_npn_train_input[,mARRAY_seeds]
		un_RNA_seq = un_train_input[,RNASeq_seeds]

		tdm_train = cbind(ref_RNA_seq,tdm_RNA_seq)
		qn_train = cbind(ref_RNA_seq,qn_RNA_seq)
		log_train = cbind(ref_RNA_seq,log_RNA_seq)
		npn_train = cbind(ref_npn_RNA_seq,npn_RNA_seq)
		un_train = cbind(ref_RNA_seq,un_RNA_seq)
		
	}
	#training outfile names (with directory location)
	tdm_train_out2 = paste0(normalized_data,prefix,"_tdm_train_",(100*x),".txt")
	qn_train_out2 = paste0(normalized_data,prefix,"_qn_train_",(100*x),".txt")
	log_train_out2 = paste0(normalized_data,prefix,"_log_train_",(100*x),".txt")
	npn_train_out2 = paste0(normalized_data,prefix,"_npn_train_",(100*x),".txt")
	un_train_out2 = paste0(normalized_data,prefix,"_un_train_",(100*x),".txt")

	write.table(tdm_train, file = tdm_train_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
	write.table(qn_train, file = qn_train_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
	write.table(log_train, file = log_train_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
	write.table(npn_train, file = npn_train_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
	write.table(un_train, file = un_train_out2, quote= FALSE, sep = "\t",row.names = TRUE, col.names = NA)
  
	#training outfile names (without directory location) to use for LASSO command
	tdm_train_out1 = paste0(prefix,"_tdm_train_",(100*x),".txt")
	qn_train_out1 = paste0(prefix,"_qn_train_",(100*x),".txt")
	log_train_out1 = paste0(prefix,"_log_train_",(100*x),".txt")
	npn_train_out1 = paste0(prefix,"_npn_train_",(100*x),".txt")
	un_train_out1 = paste0(prefix,"_un_train_",(100*x),".txt")
	
	
	prefix1=paste0("test_RS_p",(100*x),"_")
	lasso_command1 = paste("Rscript","lasso_subtype_brca_sk.R")
	lasso_command1 = paste(lasso_command1, un_train_out1) #ref_input
	lasso_command1 = paste(lasso_command1,un_test_RS_out1) #un_input
	lasso_command1 = paste(lasso_command1,tdm_train_out1) #ref_tdm_input
	lasso_command1 = paste(lasso_command1,tdm_test_RS_out1)#tdm_input
	lasso_command1 = paste(lasso_command1, qn_train_out1)#ref_qn_input
	lasso_command1 = paste(lasso_command1, qn_test_RS_out1) #qn_input
	lasso_command1 = paste(lasso_command1, log_train_out1)#ref_log_input
	lasso_command1 = paste(lasso_command1, log_test_RS_out1)#log_input
	lasso_command1 = paste(lasso_command1, npn_train_out1)#ref_npn_input
	lasso_command1 = paste(lasso_command1, npn_test_RS_out1)#npn_input
	lasso_command1 = paste(lasso_command1, prefix1,"&")#pre_fix
	system(lasso_command1)
	# run Lasso script to test on MicroArray data
	prefix2=paste0("test_MA_p",(100*x),"_")
	lasso_command2 = paste("Rscript","lasso_subtype_brca_sk.R")
	lasso_command2 = paste(lasso_command2, un_train_out1) #ref_input
	lasso_command2 = paste(lasso_command2,ref_test_MA_out1) #un_input
	lasso_command2 = paste(lasso_command2,tdm_train_out1) #ref_tdm_input
	lasso_command2 = paste(lasso_command2,ref_test_MA_out1)#tdm_input
	lasso_command2 = paste(lasso_command2, qn_train_out1)#ref_qn_input
	lasso_command2 = paste(lasso_command2, ref_test_MA_out1) #qn_input
	lasso_command2 = paste(lasso_command2, log_train_out1)#ref_log_input
	lasso_command2 = paste(lasso_command2, ref_test_MA_out1)#log_input
	lasso_command2 = paste(lasso_command2, npn_train_out1)#ref_npn_input
	lasso_command2 = paste(lasso_command2, ref_npn_test_MA_out1)#npn_input
	lasso_command2 = paste(lasso_command2, prefix2)#pre_fix
	system(lasso_command2)
	
}

#source("lasso_subtype_brca.R")
#source("normalize_coadread.R")
#source("lasso_subtype_coad.R")
#source("normalize_meta.R")
#source("lasso_subtype_meta.R")
#source("simulated.R")
#source("dist.R")
