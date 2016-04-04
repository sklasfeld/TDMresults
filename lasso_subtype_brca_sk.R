sub# Author: Jeffrey Thompson
# Purpose: Performs multiclass lasso logistic regression to classify cancer subtypes.

# Get command line arguments.

args = commandArgs(trailingOnly = T)
if (length(args) != 11) {
  help_command = "\n\nCommand format:\nRscript lasso_subtype_brca_sk.R "
  help_command = paste(help_command, "$unnormalized_training $unnormalized_testing ")
  help_command = paste(help_command, "$tdm_training $tdm_testing ")
  help_command = paste(help_command, "$quanitile_normalization_training $quanitile_normalization_testing ")
  help_command = paste(help_command, "$log_training $log_testing ")
  help_command = paste(help_command, "$npn_training $npn_testing $prefix \n\n")
  stop(help_command, call.=FALSE)
}else {
  source("package_loader.R")
  source("directories.R")

  NFOLDS = 100
  NSEEDS = 10
  INITIAL_SEED = 2

  input_dir = normalized_data
  #raw datas
  ref_input = args[1]
  ref_clin = 'BRCARNASeqClin.tsv'
  un_input = args[2]
  un_clin = 'BRCARNASeqClin.tsv'
  #tdm
  ref_tdm_input = args[3]
  ref_tdm_clin = 'BRCARNASeqClin.tsv'
  tdm_input = args[4]
  tdm_clin = 'BRCARNASeqClin.tsv'
  #qn
  ref_qn_input = args[5]
  ref_qn_clin = 'BRCARNASeqClin.tsv'
  qn_input = args[6]
  qn_clin = 'BRCARNASeqClin.tsv'
  #log
  ref_log_input = args[7]
  ref_log_clin = 'BRCARNASeqClin.tsv'
  log_input = args[8]
  log_clin = 'BRCARNASeqClin.tsv'
  #npn
  ref_npn_input = args[9]
  ref_npn_clin = 'BRCARNASeqClin.tsv'
  npn_input = args[10]
  npn_clin = 'BRCARNASeqClin.tsv'
  #output
  pre_fix = args[11]
  output_dir = output

  ACCURACY_PLOT_FILENAME = paste0(pre_fix,"brca_accuracies.pdf")
  KAPPA_PLOT_FILENAME = paste0(pre_fix,"brca_kappas.pdf")
  ACCURACY_TABLE_FILENAME = paste0(pre_fix,"brca_accuracies.tsv")
  BALANCED_PLOT_FILENAME = paste0(pre_fix,"brca_balanced_accuracies.pdf")

  load_it(c("glmnet", "caret", "e1071", "stringr", "plyr", "huge"))

  # This function creates dataframes to hold the accuracies for each subtype across runs.
  setup_acc_df = function(df) {
    df = data.frame(matrix(nrow=5, ncol=0))
    rownames(df) = c("Basal", "Her2", "LumA", "LumB", "Normal")
    return (df)
  }

  # Define data.frames to hold classification results.
  refdf = setup_acc_df(refdf)
  tdmdf = setup_acc_df(tdmdf)
  qndf = setup_acc_df(qndf)
  logdf = setup_acc_df(logdf)
  npndf = setup_acc_df(npndf)
  undf = setup_acc_df(undf)
  tdmnotshareddf = setup_acc_df(tdmnotshareddf)
  qnnotshareddf = setup_acc_df(qnnotshareddf)
  lognotshareddf = setup_acc_df(lognotshareddf)
  npnnotshareddf = setup_acc_df(npnnotshareddf)
  unnotshareddf = setup_acc_df(unnotshareddf)

  # Define random seeds.
  set.seed(INITIAL_SEED)
  seeds = sample(1:10000, NSEEDS)
  #seeds = c(7024)

  message(paste("Random seeds:", paste(seeds, collapse=", ")), appendLF=TRUE)

  # Pre-processes test data to ready it for the LASSO.
  # Returns a list with first entry the pre-processed dataset and second entry
  # the classes for those data.
  preprocessTCGA = function(dataFile, clinFile) {
    data = read.table(paste0(input_dir, dataFile), 
      sep='\t', 
      row.names=1, 
      header=T)

    dataT = t(data)
    
    # Remove duplicate samples from the data.
    dataT = dataT[str_detect(substr(rownames(dataT),1,15), "(\\.01)|(\\.11)"),]  

    # Read clinical data.
    clin = read.table(paste0(input_dir, clinFile), sep='\t', row.names=1, header=T)

    # Remove samples without subtype labels.
    clin = subset(clin, clin[,1] != "")
    clin[,1] = droplevels(clin[,1])

    # Filter clinical data to include only those samples with expression data.
    clin.response = clin[,1][match(substr(chartr('-', '.', rownames(dataT)), 1, 15), chartr('-', '.', rownames(clin)))]
    clin.response = na.omit(clin.response)

    # Filter expression data to include only those samples with ids also in clinical data.
    dataT = dataT[substr(chartr('-', '.', rownames(dataT)), 1, 15) %in% chartr('-', '.', rownames(clin)),]
    
    # Filter any data with missing clinical annotation for tumor subtype.
    dataT = dataT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(dataT)),1,15), 
      substr(chartr('-', '.', rownames(clin)),1,15))]),]
    
    return (list(data=dataT, response=clin.response))
  } # end function preprocessTCGA

  # Predict subtypes on TCGA data using previously trained model.
  predictTCGA = function(title, data, model) { 
    # Predict classes based on the trained model.
    lasso.predict = predict(model, data$data, s = "lambda.1se", type="class")
    
    # Make sure all factors are included.
    lasso.predict = factor(lasso.predict, c("Basal", "Her2", "LumA", "LumB", "Normal"))
    
    # Build a contingency table of the results.
    con_table = table(pred = lasso.predict, true = data$response, exclude="none")
    
    # Generate statistics based on the contigency table.
    cm = confusionMatrix(con_table)
    
    return (cm)
  }

  # Load the training data.
  train = read.table(paste0(input_dir, ref_input), sep='\t', row.names=1, header=T)

  # Transpose rows and columns.
  trainT = t(train)

  # Remove duplicate samples from the data.
  trainT = trainT[str_detect(substr(rownames(trainT),1,15), "(\\.01)|(\\.11)"),]  

  # Load the clinical features.
  clin = read.table(paste0(input_dir, ref_clin), sep='\t', row.names=1, header=T)
  #clin[,1] = clin[,3]

  # Remove annotations that have no subtype label.
  clin = subset(clin, clin[,1] != "")
  clin[,1] = droplevels(clin[,1])

  # Filter the expression data to remove those samples without subtype labels.
  trainT = trainT[substr(rownames(trainT),1,15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]

  # Filter clinical data to include only those entries that also have expression data.
  trainclin.subtypes = clin[,1][match(substr(chartr('-', '.', rownames(trainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]
  trainclin.subtypes = na.omit(trainclin.subtypes)
  write.table(trainclin.subtypes, file=paste0(output_dir,pre_fix, "brca_ref_classes.txt"))

  # Filter any expression data with missing clinical annotation for tumor subtype.
  trainT = trainT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(trainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]),]

  # Load the training data.
  tdm_train = read.table(paste0(input_dir, ref_tdm_input), sep='\t', row.names=1, header=T)

  # Transpose rows and columns.
  tdmtrainT = t(tdm_train)

  # Remove duplicate samples from the data.
  tdmtrainT = tdmtrainT[str_detect(substr(rownames(tdmtrainT),1,15), "(\\.01)|(\\.11)"),]

  # Filter the expression data to remove those samples without subtype labels.
  tdmtrainT = tdmtrainT[substr(rownames(tdmtrainT),1,15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]

  # Filter any expression data with missing clinical annotation for tumor subtype.
  tdmtrainT = tdmtrainT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(tdmtrainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]),]

  # Load the training data.
  qn_train = read.table(paste0(input_dir, ref_qn_input), sep='\t', row.names=1, header=T)

  # Transpose rows and columns.
  qntrainT = t(qn_train)

  # Remove duplicate samples from the data.
  qntrainT = qntrainT[str_detect(substr(rownames(qntrainT),1,15), "(\\.01)|(\\.11)"),]

  # Filter the expression data to remove those samples without subtype labels.
  qntrainT = qntrainT[substr(rownames(qntrainT),1,15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]

  # Filter any expression data with missing clinical annotation for tumor subtype.
  qntrainT = qntrainT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(qntrainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]),]

  # Load the training data.
  log_train = read.table(paste0(input_dir, ref_log_input), sep='\t', row.names=1, header=T)

  # Transpose rows and columns.
  logtrainT = t(log_train)

  # Remove duplicate samples from the data.
  logtrainT = logtrainT[str_detect(substr(rownames(logtrainT),1,15), "(\\.01)|(\\.11)"),]

  # Filter the expression data to remove those samples without subtype labels.
  logtrainT = logtrainT[substr(rownames(logtrainT),1,15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]

  # Filter any expression data with missing clinical annotation for tumor subtype.
  logtrainT = qntrainT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(logtrainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]),]

  # Load the training data.
  npn_train = read.table(paste0(input_dir, ref_npn_input), sep='\t', row.names=1, header=T)

  # Transpose rows and columns.
  npntrainT = t(npn_train)

  # Remove duplicate samples from the data.
  npntrainT = npntrainT[str_detect(substr(rownames(npntrainT),1,15), "(\\.01)|(\\.11)"),]  

  # Filter the expression data to remove those samples without subtype labels.
  npntrainT = npntrainT[substr(rownames(npntrainT),1,15) %in% substr(chartr('-', '.', rownames(clin)),1,15),]

  # Filter any expression data with missing clinical annotation for tumor subtype.
  npntrainT = npntrainT[!is.na(clin[,1][match(substr(chartr('-', '.', rownames(npntrainT)),1,15), substr(chartr('-', '.', rownames(clin)),1,15))]),]

  # Pre-process the TCGA data.
  tdm = preprocessTCGA(tdm_input, tdm_clin)
  qn = preprocessTCGA(qn_input, qn_clin)
  log = preprocessTCGA(log_input, log_clin)
  npn = preprocessTCGA(npn_input, npn_clin)
  un = preprocessTCGA(un_input, un_clin)

  # Determine which samples are shared between training and test data.
  tdm_shared = substr(chartr('-', '.', rownames(tdm$data)),1,15) %in% substr(chartr('-', '.', rownames(tdmtrainT)),1,15)  
  qn_shared = substr(chartr('-', '.', rownames(qn$data)),1,15) %in% substr(chartr('-', '.', rownames(qntrainT)),1,15)
  log_shared = substr(chartr('-', '.', rownames(log$data)),1,15) %in% substr(chartr('-', '.', rownames(logtrainT)),1,15)
  npn_shared = substr(chartr('-', '.', rownames(npn$data)),1,15) %in% substr(chartr('-', '.', rownames(npntrainT)),1,15)
  un_shared = substr(chartr('-', '.', rownames(un$data)),1,15) %in% substr(chartr('-', '.', rownames(trainT)),1,15)

  tdm_notshared = tdm
  qn_notshared = qn
  log_notshared = log
  npn_notshared = npn
  un_notshared = un

  # Create datasets for samples that are not shared between training and test.
  tdm_notshared$data = tdm_notshared$data[!tdm_shared,]
  tdm_notshared$response = tdm_notshared$response[!tdm_shared]
  qn_notshared$data = qn_notshared$data[!qn_shared,]
  qn_notshared$response = qn_notshared$response[!qn_shared]
  log_notshared$data = log_notshared$data[!log_shared,]
  log_notshared$response = log_notshared$response[!log_shared]
  npn_notshared$data = npn_notshared$data[!npn_shared,]
  npn_notshared$response = npn_notshared$response[!npn_shared]
  un_notshared$data = un_notshared$data[!un_shared,]
  un_notshared$response = un_notshared$response[!un_shared]

  # Pre-allocate vectors to hold results.
  tdm_accs = vector(mode="list", length=length(seeds))
  qn_accs = vector(mode="list", length=length(seeds))
  log_accs = vector(mode="list", length=length(seeds))
  npn_accs = vector(mode="list", length=length(seeds))
  un_accs = vector(mode="list", length=length(seeds))
  tdmnotshared_accs = vector(mode="list", length=length(seeds))
  qnnotshared_accs = vector(mode="list", length=length(seeds))
  lognotshared_accs = vector(mode="list", length=length(seeds))
  npnnotshared_accs = vector(mode="list", length=length(seeds))
  unnotshared_accs = vector(mode="list", length=length(seeds))

  # Perform a number of iterations equal to the number of random seeds.
  # At each iteration, perform n-fold cross validation to build a model on the training data,
  # then use that model to make predictions on the test data.
  for(seednum in 1:length(seeds)) {
    set.seed(seeds[seednum])
    lasso.model=cv.glmnet(trainT, trainclin.subtypes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
    plot(lasso.model) 

    set.seed(seeds[seednum])
    tdm.lasso.model=cv.glmnet(tdmtrainT, trainclin.subtypes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
    
    set.seed(seeds[seednum])
    qn.lasso.model=cv.glmnet(qntrainT, trainclin.subtypes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
    
    set.seed(seeds[seednum])
    log.lasso.model=cv.glmnet(logtrainT, trainclin.subtypes, family="multinomial", parallel = F, type.measure="class", nfolds=NFOLDS)
    
    set.seed(seeds[seednum])
    npn.lasso.model=cv.glmnet(npntrainT, trainclin.subtypes, family="multinomial", parallel=F, type.measure="class", nfolds=NFOLDS)
    
    ############ TESTING STARTS HERE #############


    acc = predictTCGA("TDM Results", tdm, tdm.lasso.model)
    tdm_accs[[seednum]] = acc
    tdmdf = cbind(tdmdf,as.vector(acc$byClass[,8]))
    
    acc = predictTCGA("QN Results", qn, qn.lasso.model)
    qn_accs[[seednum]] = acc
    qndf = cbind(qndf,as.vector(acc$byClass[,8]))

    acc = predictTCGA("LOG Results", log, log.lasso.model)
    log_accs[[seednum]] = acc
    logdf = cbind(logdf,as.vector(acc$byClass[,8]))

    acc = predictTCGA("NPN Results", npn, npn.lasso.model)
    npn_accs[[seednum]] = acc
    npndf = cbind(npndf, as.vector(acc$byClass[,8]))

    acc = predictTCGA("Untransformed Results", un, lasso.model)
    un_accs[[seednum]] = acc
    undf = cbind(undf, as.vector(acc$byClass[,8]))

    acc = predictTCGA("TDM (not shared) Results", tdm_notshared, lasso.model)
    tdmnotshared_accs[[seednum]] = acc
    tdmnotshareddf = cbind(tdmnotshareddf,as.vector(acc$byClass[,8]))
    
    acc = predictTCGA("QN (not shared) Results", qn_notshared, lasso.model)
    qnnotshared_accs[[seednum]] = acc
    qnnotshareddf = cbind(qnnotshareddf,as.vector(acc$byClass[,8]))

    acc = predictTCGA("LOG(not shared) Results", log_notshared, lasso.model)
    lognotshared_accs[[seednum]] = acc
    lognotshareddf = cbind(lognotshareddf,as.vector(acc$byClass[,8]))

    acc = predictTCGA("NPN (not shared) Results", npn_notshared, npn.lasso.model)
    npnnotshared_accs[[seednum]] = acc
    npnnotshareddf = cbind(npnnotshareddf, as.vector(acc$byClass[,8]))

    acc = predictTCGA("Untransformed (not shared) Results", un_notshared, lasso.model)
    unnotshared_accs[[seednum]] = acc
    unnotshareddf = cbind(unnotshareddf, as.vector(acc$byClass[,8]))
  }

  # Build a table of accuracies across all datasets.
  accuracies = data.frame(Basal=numeric(0), Her2=numeric(0), LumA=numeric(0), LumB=numeric(0), Normal=numeric(0))
  accuracies = rbind(accuracies, apply(tdmdf,1,mean))
  accuracies = rbind(accuracies, apply(qndf,1,mean))
  accuracies = rbind(accuracies, apply(logdf,1,mean))
  accuracies = rbind(accuracies, apply(npndf,1,mean))
  accuracies = rbind(accuracies, apply(undf,1,mean))
  rownames(accuracies) = c("TDM", "QN", "LOG", "NPN", "UNTR")
  colnames(accuracies) = c("Basal", "Her2", "LumA", "LumB", "Normal")

  # Aggregate accuracies:
  tdm_tables = lapply(tdm_accs, function(x) x$table)
  qn_tables = lapply(qn_accs, function(x) x$table)
  log_tables = lapply(log_accs, function(x) x$table)
  npn_tables = lapply(npn_accs, function(x) x$table)
  un_tables = lapply(un_accs, function(x) x$table)

  tdm_reduced = Reduce("+", tdm_tables) / length(tdm_tables)
  qn_reduced = Reduce("+", qn_tables) / length(qn_tables)
  log_reduced = Reduce("+", log_tables) / length(log_tables)
  npn_reduced = Reduce("+", npn_tables) / length(npn_tables)
  un_reduced = Reduce("+", un_tables) / length(un_tables)

  tdm_reduced[1:5,1:5] = t(apply(tdm_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
  qn_reduced[1:5,1:5] = t(apply(qn_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
  log_reduced[1:5,1:5] = t(apply(log_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
  npn_reduced[1:5,1:5] = t(apply(npn_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]
  un_reduced[1:5,1:5] = t(apply(un_reduced, 1, function(x) as.integer(round(x))))[1:5,1:5]

  tdm_cm = confusionMatrix(tdm_reduced)
  qn_cm = confusionMatrix(qn_reduced)
  log_cm = confusionMatrix(log_reduced)
  npn_cm = confusionMatrix(npn_reduced)
  un_cm = confusionMatrix(un_reduced)

  saveRDS(log_cm, paste0(output_dir, pre_fix, "brca_log_cm.RDS"))

  all_accs = data.frame(rbind(tdm_cm$overall, qn_cm$overall, log_cm$overall, npn_cm$overall, un_cm$overall))
  all_accs = cbind(all_accs, method=c("TDM", "QN", "LOG", "NPN", "UNTR"))

  all_accs$method = factor(all_accs$method, levels=all_accs$method)

  nir = all_accs[1,]$AccuracyNull

  # Plot accuracies
  ci = aes(ymin=AccuracyLower, ymax=AccuracyUpper)
  cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(all_accs, aes(x=factor(method), y=Accuracy, color=method)) +
  geom_point(size=3) +
  geom_errorbar(ci, width=.3) +
  ylab("Total Accuracy") +
  xlab("Normalization") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_hline(aes(yintercept = nir), linetype="longdash") +
  scale_color_manual(values=cbPalette) 

  ggsave(paste0(output_dir, pre_fix, ACCURACY_PLOT_FILENAME), plot=last_plot(), width=3, height=2.5)
  cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  #cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(data = all_accs) +
  ylab("Kappa") +
  xlab("Normalization") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_point(aes(x=factor(method), y=Kappa, color=method), shape=18, size=4) +
  scale_color_manual(values=cbPalette) 

  ggsave(paste0(output_dir, pre_fix, KAPPA_PLOT_FILENAME), plot=last_plot(), width=3, height=2.5)

  brca_all = rbind(tdm_cm$byClass[,8],
    qn_cm$byClass[,8],
    log_cm$byClass[,8],
    npn_cm$byClass[,8],
    un_cm$byClass[,8])

  brca_melted = melt(brca_all)
  colnames(brca_melted) = c("Dataset", "Subtype", "Accuracy")
  brca_melted$Dataset = c("TDM", "QN", "LOG", "NPN", "UNTR")
  brca_melted$Dataset = factor(brca_melted$Dataset, levels=c("TDM", "QN", "LOG", "NPN", "UNTR"))

  cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(data=brca_melted, aes(x=Dataset, y=Accuracy, ymax=1)) +
    ylim(0.45,1) +
    coord_flip() + 
    geom_point(size=3, aes(colour=Dataset)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.title=element_text(hjust=0)) +
    ylab("Balanced Accuracy") + theme_bw() + facet_wrap(~Subtype, ncol=5) +
    theme(legend.position="none") +
    scale_color_manual(values=cbPalette)

  ggsave(paste0(output_dir, pre_fix, BALANCED_PLOT_FILENAME), plot=last_plot(), height=2.5, width=6)

  write.table(all_accs, file=paste0(output_dir, pre_fix, ACCURACY_TABLE_FILENAME), sep='\t', row.names=T, col.names=T)
}
