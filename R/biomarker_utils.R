### Biomarker Analysis for ProteoAnalyst
### Adapted from MetaboAnalyst biomarker_utils.R
### Jeff Xia \email{jeff.xia@mcgill.ca}
### Adapted for ProteoAnalyst framework
### License: GNU GPL (>= 2)

##############################################
## KEY ADAPTATIONS FOR PROTEOANALYST:
## 1. mSetObj → separate dataSet, paramSet, analSet objects
## 2. .get.mSet() → readDataset(dataName), readSet(analSet)
## 3. .set.mSet() → RegisterData(dataSet), saveSet(analSet)
## 4. Data access: mSetObj$dataSet$norm → dataSet$data.norm.transposed
## 5. Classes: mSetObj$dataSet$cls → dataSet$cls
## 6. Analysis results: mSetObj$analSet → analSet (separate global object)
## 7. File persistence using ov_qs_save() and ProteoAnalyst naming
##############################################

#'Numbers for subset selection
#'@description Return a series of number for subsets selection
#'@param feat.len Input the feature length
#'@export
GetFeatureNumbers <- function(feat.len){
  if(feat.len > 100){
    nFeatures <- c(5, 10, 15, 25, 50, 100);
  }else if(feat.len > 50){
    nFeatures <- c(3, 5, 10, 20, round(feat.len/2), feat.len);
  }else if(feat.len > 20){
    nFeatures <- c(2, 3, 5, 10, 20, feat.len);
  }else if(feat.len > 10){
    nFeatures <- c(2, 3, 5, 7, 10, feat.len);
  }else if(feat.len > 1){
    nFeatures <- sort(sample(2:feat.len));
  }else{
    print("Feature number is less than 2!")
    return();
  }
  nFeatures;
}

#'Make random partitions
#'@description Make random partitions, returns matrices indicating
#'whether the observation is in train/test for each run
#'note: try to get a balanced sampling for each group (classification)
#'or each quantile (regression). This is very useful for unbalanced data
#'@param y Input the data
#'@param propTraining By default set to 2/3
#'@param nRuns By default set to 30
GetTrainTestSplitMat <- function(y, propTraining = 2/3, nRuns = 30){

  nTotalSample <- length(y);

  smallestClass <- names(sort(table(y)))[1];
  nSmallest <- sum(y == smallestClass);

  nSmallestTrain <- round(propTraining * nSmallest);
  nBiggestTrain <- nSmallestTrain;

  nSmallestTest <- nSmallest - nSmallestTrain;
  nBiggestTest <- nTotalSample - (nSmallestTest + nSmallestTrain + nBiggestTrain);

  # sanity check for very large number of samples
  # for each run max 600 - 400 train, 200 test
  big.samples <- FALSE;
  if(nSmallestTrain > 400){
    big.samples <- TRUE;
    nSmallestTrain <- nBiggestTrain <- 400;
    nSmallestTest <- nBiggestTest <- 200;
  }

  # split up in smallest class indices and biggest class indices
  smallestIndices <- which(y == smallestClass)
  biggestIndices <- seq(along = y)[-smallestIndices]

  nTrainingSample <- nSmallestTrain + nBiggestTrain;
  nTestSample <- nSmallestTest + nBiggestTest;

  trainingSampleAllRuns <- matrix(0, nrow = nRuns, ncol = nTrainingSample)
  testSampleAllRuns  <- matrix(0, nrow = nRuns, ncol = nTestSample);

  for (irun in 1:nRuns) {
    sampleSmallestTrain <- sample(smallestIndices, nSmallestTrain);
    sampleBiggestTrain <- sample(biggestIndices, nBiggestTrain);
    trainingSampleRun <- c(sampleSmallestTrain, sampleBiggestTrain);
    indicesTrainingSample <- rep(FALSE, length = nTotalSample);
    indicesTrainingSample[trainingSampleRun] <- TRUE;
    trainingSampleAllRuns[irun, ] <- which(indicesTrainingSample);
    if(big.samples){
      testSampleAllRuns[irun, ] <- sample(which(!indicesTrainingSample), 200);
    }else{
      testSampleAllRuns[irun, ] <- which(!indicesTrainingSample);
    }
  }

  # return the results
  list(
    training.mat = trainingSampleAllRuns,
    testing.mat = testSampleAllRuns
  );
}

#'Choose two groups (when more than two groups uploaded)
#'@param dataName Input the dataset name
#'@param grps Group names to compare
#'@export
SetCurrentGroups <- function(dataName = "", grps){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");

  if(length(grps) != 2){
    msgSet$current.msg <- "Biomarker analysis requires exactly two groups.";
    saveSet(msgSet, "msgSet");
    return(0);
  }

  cls <- dataSet$cls;
  sel.grps <- cls %in% grps;

  # subset data
  dataSet$data.norm.transposed <- dataSet$data.norm.transposed[sel.grps, , drop=FALSE];
  dataSet$cls <- factor(as.character(cls[sel.grps]));

  paramSet$biomarker.grps <- grps;
  msgSet$current.msg <- paste0("Selected groups: ", paste(grps, collapse = " vs "));

  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");
  return(RegisterData(dataSet));
}

#'Rank features using different methods
#'@param x.in Feature matrix (samples × features)
#'@param y.in Class labels
#'@param method Ranking method: "auroc", "ttest", "foldchange"
#'@param lvNum Number of latent variables for PLS
RankFeatures <- function(x.in, y.in, method, lvNum){
  #msg("[RankFeatures] DEBUG: method=", method, " x.in dims=", nrow(x.in), "x", ncol(x.in), " y.in length=", length(y.in))
  #msg("[RankFeatures] DEBUG: y.in class=", class(y.in), " levels=", paste(levels(y.in), collapse=", "))
  #msg("[RankFeatures] DEBUG: Any NAs in x.in? ", any(is.na(x.in)), " Any NAs in y.in? ", any(is.na(y.in)))

  if(method == "rf"){ # use randomforest mean decrease accuracy
    #msg("[RankFeatures] DEBUG: Using randomForest method")
    rf <- randomForest::randomForest(x = x.in,y = y.in,importance=TRUE, keep.forest=F);
    return(randomForest::importance(rf)[ ,"MeanDecreaseAccuracy"]);
  }else if (method == "pls"){
    #msg("[RankFeatures] DEBUG: Using pls method with lvNum=", lvNum)
    ncls <- as.numeric(y.in)-1;
    datmat <- as.matrix(x.in);
    pls <- pls::plsr(ncls~datmat,method='oscorespls', ncomp=lvNum);
    return(Get.VIP(pls, lvNum));
  }else if(method == "svm"){
    #msg("[RankFeatures] DEBUG: Using svm method")
    svmres <- e1071::svm(x.in, y.in, type = 'C', kernel="linear");
    imp.vec <- (t(svmres$coefs) %*% svmres$SV)[1,]^2;
    names(imp.vec) <- colnames(x.in);
    return(imp.vec);
  }else if(method == "auroc"){ # univariate based ou area under ROC
    #msg("[RankFeatures] DEBUG: Using auroc method")
    imp.vec <- caTools::colAUC(x.in, y.in, plotROC=F)[1,];
    #msg("[RankFeatures] DEBUG: colAUC returned vector of length ", length(imp.vec))
    return(imp.vec);
  }else if(method == "tt"){ # univariate based ou area under ROC
    #msg("[RankFeatures] DEBUG: Using tt method")
    imp.vec <- Get.Fstat(x.in, as.factor(y.in)); # all f-stats
    names(imp.vec) <- colnames(x.in);
    return(imp.vec);
  }else if(method == "fisher"){ # univariate based ou area under ROC
    #msg("[RankFeatures] DEBUG: Using fisher method")
    imp.vec <- Get.Fisher(x.in, as.factor(y.in));
    names(imp.vec) <- colnames(x.in);
    return(imp.vec);
  }else{
    print(paste("Not supported method:", method));
    return(NULL);
  }
}

#'Calculate feature ranking for all features
#'@param dataName Dataset name
#'@param clust.num Number of clusters for k-means
#'@export
CalculateFeatureRanking <- function(dataName = "", clust.num=5){
  save.image("feat.RData");
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  analSet <- readSet(analSet, "analSet");

  if(is.null(analSet)){
    analSet <- list();
  }
  if(is.null(analSet$LR)){
    analSet$LR <- list();
  }

  analSet$LR$LRConverged <- "FALSE";

  x  <- dataSet$data.norm.transposed;
  if (is.null(x)) {
    x <- t(dataSet$data.norm); # samples × features
  }
  y <- factor(dataSet$cls);
  if (length(y) == 0 || nlevels(y) < 2) {
    msgSet$current.msg <- "Not enough groups for ROC/feature ranking.";
    saveSet(msgSet, "msgSet");
    return(0)
  }
  if (nrow(x) != length(y)) {
    #msg("[CalculateFeatureRanking] WARNING: nrow(x) != length(y); attempting to realign.")
    if (ncol(x) == length(y)) {
      x <- t(x)
    }
    if (nrow(x) != length(y)) {
      msgSet$current.msg <- "Sample labels do not align with data for feature ranking.";
      saveSet(msgSet, "msgSet");
      return(0)
    }
  }

  # AUC
  auc <- caTools::colAUC(x, y, plotROC=F)[1,];

  # T-test p-values
  ttp <- GetROCTtestP(x, y);

  # Fold change - use data before normalization if available
  if(file.exists("data.raw.qs")){
    data <- ov_qs_read("data.raw.qs"); # expected samples × features
    sam.order <- if (!is.null(rownames(data))) intersect(rownames(x), rownames(data)) else rownames(x)
    feat.order <- if (!is.null(colnames(data))) intersect(colnames(x), colnames(data)) else colnames(x)
    if (length(sam.order) == 0 || length(feat.order) == 0) {
      data <- x;  # fallback if no overlap
    } else {
      data <- data[sam.order, feat.order, drop=FALSE];
      x <- x[match(sam.order, rownames(x)), feat.order, drop=FALSE]
      y <- y[match(sam.order, rownames(x))]
    }
  }else{
    data <- x;  # fallback to normalized
  }

  if (nlevels(y) < 2) {
    msgSet$current.msg <- "Not enough groups for ROC/feature ranking.";
    saveSet(msgSet, "msgSet");
    return(0)
  }
  if (length(y) != nrow(data)) {
    msgSet$current.msg <- "Sample labels do not align with data for feature ranking.";
    saveSet(msgSet, "msgSet");
    return(0)
  }

  idx1 <- which(y==levels(y)[1])
  idx2 <- which(y==levels(y)[2])
  if (length(idx1) == 0 || length(idx2) == 0) {
    msgSet$current.msg <- "Group labels do not match any samples for feature ranking.";
    saveSet(msgSet, "msgSet");
    return(0)
  }
  m1 <- colMeans(data[idx1, , drop=FALSE]);
  m2 <- colMeans(data[idx2, , drop=FALSE]);
  ratio <- m1/m2;
  ratio[ratio < 0] <- 0;
  fc <- signif(log2(ratio), 5);
  fc[is.infinite(fc) & fc < 0] <- -99;
  fc[is.infinite(fc) & fc > 0] <- 99;

  # K-means clustering
  if(ncol(x) > 1){
    kms <- ComputeKmeanClusters(t(x), clust.num);
    feat.rank.mat <- data.frame(AUC=auc, Pval=ttp, FC=fc, clusters = kms);
    rownames(feat.rank.mat) <- colnames(x);

    ord.inx <- order(feat.rank.mat$AUC, decreasing=T);
    feat.rank.mat  <- data.matrix(feat.rank.mat[ord.inx, , drop=FALSE]);
  }else{
    feat.rank.mat <- data.matrix(data.frame(AUC=auc, Pval=ttp, FC=fc, clusters=1))
  }

  # Format and save
  feat.rank.mat <- signif(feat.rank.mat, digits = 5);

  # Save to analSet
  analSet$feat.rank.mat <- feat.rank.mat;

  # Export to CSV
  if(!is.null(analSet$mode) && analSet$mode == "univ"){
    fast.write.csv(feat.rank.mat, file="biomarker_univ_ranking.csv");
  }else{
    fast.write.csv(feat.rank.mat, file="biomarker_feature_ranking.csv");
  }

  msgSet$current.msg <- paste0("Feature ranking calculated for ", nrow(feat.rank.mat), " features.");

  saveSet(analSet, "analSet");
  saveSet(msgSet, "msgSet");
  return(RegisterData(dataSet));
}

# K-means clustering helper
ComputeKmeanClusters <- function(data, clust.num){
  set.seed(28051968);
  km.out <- kmeans(data, clust.num, nstart=100);
  return(km.out$cluster);
}

#'Update k-means clustering
#'@param dataName Dataset name
#'@param clust.num Number of clusters
#'@export
UpdateKmeans <- function(dataName = "", clust.num=5){
  analSet <- readSet(analSet, "analSet");
  msgSet <- readSet(msgSet, "msgSet");

  if(is.null(analSet$feat.rank.mat)){
    msgSet$current.msg <- "Feature ranking must be calculated first.";
    saveSet(msgSet, "msgSet");
    return(0);
  }

  dataSet <- readDataset(dataName);
  x <- dataSet$data.norm.transposed;
  kms <- ComputeKmeanClusters(t(x), clust.num);

  analSet$feat.rank.mat[, "clusters"] <- kms[rownames(analSet$feat.rank.mat)];

  saveSet(analSet, "analSet");
  msgSet$current.msg <- paste0("K-means clustering updated with ", clust.num, " clusters.");
  saveSet(msgSet, "msgSet");

  return(1);
}

#'Set analysis mode
#'@param dataName Dataset name
#'@param mode Analysis mode: "univ" or "multivar"
#'@export
SetAnalysisMode <- function(dataName = "", mode){
  analSet <- readSet(analSet, "analSet");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");

  if(is.null(analSet)){
    analSet <- list();
  }

  analSet$mode <- mode;
  analSet$type <- "roc";

  paramSet$biomarker.mode <- mode;
  msgSet$current.msg <- paste0("Biomarker analysis mode set to: ", mode);

  saveSet(analSet, "analSet");
  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");

  return(1);
}

#'Perform cross-validation for biomarker exploration
#'@param dataName Dataset name
#'@param cls.method Classification method: "pls", "plsda", "rf", "svm"
#'@param rank.method Feature ranking method: "auroc", "ttest", "foldchange"
#'@param lvNum Number of latent variables
#'@param propTraining Proportion for training (default 2/3)
#'@export
PerformCV.explore <- function(dataName = "", cls.method, rank.method="auroc", lvNum=2, propTraining=2/3){

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");

  if(is.null(analSet)){
    analSet <- list();
  }

  analSet$exp.method <- cls.method;
  analSet$rank.method <- rank.method;
  analSet$exp.lvNum <- lvNum;

  data <- dataSet$data.norm.transposed;
  cls <- dataSet$cls;

  # DEBUG: Check what we received
  #msg("[PerformCV.explore] DEBUG: data is.null? ", is.null(data))
  if (!is.null(data)) {
    #msg("[PerformCV.explore] DEBUG: data class: ", class(data))
    #msg("[PerformCV.explore] DEBUG: data dimensions: ", nrow(data), " x ", ncol(data))
    #msg("[PerformCV.explore] DEBUG: data has rownames? ", !is.null(rownames(data)))
    #msg("[PerformCV.explore] DEBUG: data has colnames? ", !is.null(colnames(data)))
    #msg("[PerformCV.explore] DEBUG: Any NAs in data? ", any(is.na(data)))
  }
  #msg("[PerformCV.explore] DEBUG: cls is.null? ", is.null(cls))
  if (!is.null(cls)) {
    #msg("[PerformCV.explore] DEBUG: cls class: ", class(cls))
    #msg("[PerformCV.explore] DEBUG: cls length: ", length(cls))
    #msg("[PerformCV.explore] DEBUG: cls levels: ", paste(levels(cls), collapse=", "))
    #msg("[PerformCV.explore] DEBUG: cls table: ", paste(names(table(cls)), "=", table(cls), collapse=", "))
    #msg("[PerformCV.explore] DEBUG: Any NAs in cls? ", any(is.na(cls)))
  }

  # number of subsampling to produce smooth curve
  if(nrow(data) > 500){
    nRuns <- 10;
  }else if(nrow(data) > 200){
    nRuns <- 20;
  }else if(nrow(data) > 100){
    nRuns <- 30;
  }else{
    nRuns <- 50;
  }

  nFeatures <- GetFeatureNumbers(ncol(data));

  #msg("[PerformCV.explore] DEBUG: nFeatures = ", paste(nFeatures, collapse=", "))

  feat.outp <- actualCls <- vector(length = nRuns, mode = "list");

  perf.outp <- vector(length = length(nFeatures), mode = "list");
  perf.outp <- lapply(perf.outp, function(x){x <- vector(length = nRuns, mode = "list"); return(x)});
  auc.mat <- accu.mat <- matrix(nrow=nRuns, ncol=length(nFeatures));

  #msg("[PerformCV.explore] DEBUG: About to call GetTrainTestSplitMat with cls length=", length(cls), " propTraining=", propTraining, " nRuns=", nRuns)

  splitMat <- GetTrainTestSplitMat(cls, propTraining, nRuns);
  trainRuns <- splitMat$training.mat;
  testRuns <- splitMat$testing.mat;

  #msg("[PerformCV.explore] DEBUG: trainRuns dimensions: ", nrow(trainRuns), " x ", ncol(trainRuns))
  #msg("[PerformCV.explore] DEBUG: testRuns dimensions: ", nrow(testRuns), " x ", ncol(testRuns))

  #msg("[PerformCV.explore] Starting cross-validation with ", nRuns, " runs and ", length(nFeatures), " feature subsets");

  for (irun in 1:nRuns){
    if (irun == 1) {
      #msg("[PerformCV.explore] DEBUG: Starting iteration ", irun)
    }
    trainingSampleRun <- trainRuns[irun, ]
    testSampleRun <- testRuns[irun, ];

    if (irun == 1) {
      #msg("[PerformCV.explore] DEBUG: trainingSampleRun is logical? ", is.logical(trainingSampleRun))
      #msg("[PerformCV.explore] DEBUG: trainingSampleRun length: ", length(trainingSampleRun))
      #msg("[PerformCV.explore] DEBUG: trainingSampleRun sum: ", sum(trainingSampleRun))
      #msg("[PerformCV.explore] DEBUG: testSampleRun sum: ", sum(testSampleRun))
    }

    x.in <- data[trainingSampleRun, ];
    y.train <- cls[trainingSampleRun];
    actualCls[[irun]] <- y.test <- cls[testSampleRun];

    if (irun == 1) {
      #msg("[PerformCV.explore] DEBUG: x.in dimensions: ", nrow(x.in), " x ", ncol(x.in))
      #msg("[PerformCV.explore] DEBUG: y.train length: ", length(y.train), " levels: ", paste(levels(y.train), collapse=", "))
      #msg("[PerformCV.explore] DEBUG: y.train table: ", paste(names(table(y.train)), "=", table(y.train), collapse=", "))
      #msg("[PerformCV.explore] DEBUG: Any NAs in x.in? ", any(is.na(x.in)))
      #msg("[PerformCV.explore] DEBUG: Any NAs in y.train? ", any(is.na(y.train)))
      #msg("[PerformCV.explore] DEBUG: About to call RankFeatures with method=", rank.method)
    }

    # Feature ranking using only training data
    imp.vec <- RankFeatures(x.in, y.train, rank.method, lvNum);

    if (irun == 1) {
      #msg("[PerformCV.explore] DEBUG: RankFeatures returned vector of length: ", length(imp.vec))
      #msg("[PerformCV.explore] DEBUG: Any NAs in imp.vec? ", any(is.na(imp.vec)))
      #msg("[PerformCV.explore] DEBUG: First 5 values: ", paste(head(imp.vec, 5), collapse=", "))
    }

    feat.outp[[irun]] <- imp.vec;
    ord.inx <- order(imp.vec, decreasing = TRUE);
    ordData <- data[, ord.inx];

    # building classifiers for each number of selected features and test on the test data
    for (inum in seq(along = nFeatures)){
      if (irun == 1 && inum == 1) {
        #msg("[PerformCV.explore] DEBUG: Inner loop - inum=", inum, " nFeatures[inum]=", nFeatures[inum])
      }

      x.train <- ordData[trainingSampleRun, 1:nFeatures[inum], drop=FALSE];
      x.test <- ordData[testSampleRun, 1:nFeatures[inum], drop=FALSE];

      if (irun == 1 && inum == 1) {
        #msg("[PerformCV.explore] DEBUG: x.train dimensions: ", nrow(x.train), " x ", ncol(x.train))
        #msg("[PerformCV.explore] DEBUG: x.test dimensions: ", nrow(x.test), " x ", ncol(x.test))
        #msg("[PerformCV.explore] DEBUG: About to call Predict.class with method=", cls.method)
      }

      prob.out <- Predict.class(x.train, y.train, x.test, cls.method, lvNum);

      if (irun == 1 && inum == 1) {
        #msg("[PerformCV.explore] DEBUG: prob.out class: ", class(prob.out))
        #msg("[PerformCV.explore] DEBUG: prob.out length: ", length(prob.out))
        #msg("[PerformCV.explore] DEBUG: Any NAs in prob.out? ", any(is.na(prob.out)))
        #msg("[PerformCV.explore] DEBUG: y.test length: ", length(y.test))
        #msg("[PerformCV.explore] DEBUG: About to create ROCR prediction object")
      }

      # calculate AUC for each
      pred <- ROCR::prediction(prob.out, y.test);
      auc.mat[irun, inum] <- slot(ROCR::performance(pred, "auc"), "y.values")[[1]];

      perf.outp[[inum]][[irun]] <- prob.out;
      pred.out <- as.factor(ifelse(prob.out > 0.5, 1, 0));
      accu.mat[irun, inum] <- Get.Accuracy(table(pred.out, y.test));

      if (irun == 1 && inum == 1) {
        #msg("[PerformCV.explore] DEBUG: Successfully completed iteration irun=1, inum=1")
      }
    }
  }

  #############################################################################
  ## prepare results for default plotting
  ## 1) get best model based on AUROC for prob.view and imp.feature calculation
  ## 2) calculate accuracy and roc for all models under comparison
  #############################################################################

  preds <- vector(length = length(nFeatures), mode = "list");
  act.vec <- unlist(actualCls); # same for all subsets
  for(m in 1:length(nFeatures)){
    prob.vec <- unlist(perf.outp[[m]]);
    pred <- ROCR::prediction(prob.vec, act.vec);
    preds[[m]] <- pred; # prediction obj
  }

  # accu.mat and preds are for all models
  colnames(accu.mat) <- colnames(auc.mat) <- names(preds) <- paste(nFeatures, sep = "");
  auc.vec <- colMeans(auc.mat);

  auc.cis <- GetCIs(auc.mat);
  # get the best based on AUROC
  best.model.inx <- which.max(auc.vec);

  analSet$multiROC <- list(
    type = analSet$type,
    train.mat = trainRuns,

    # record raw output
    test.feats = nFeatures,
    pred.cv = perf.outp,
    true.cv = actualCls,
    imp.cv = feat.outp,

    # record processed output
    pred.list = preds,
    accu.mat = accu.mat,
    auc.vec = auc.vec,
    auc.mat = auc.mat,
    auc.ci = auc.cis,
    best.inx = best.model.inx
  );

  # Save results
  ov_qs_save(analSet$multiROC, "biomarker_cv_explore.qs");

  msgSet$current.msg <- paste0("Cross-validation completed. Best model: ", nFeatures[best.model.inx],
                               " features (AUC=", round(auc.vec[best.model.inx], 3), ")");

  saveSet(analSet, "analSet");
  saveSet(msgSet, "msgSet");

  return(RegisterData(dataSet));
}

#'Predict class probabilities
#'@param x.train Training features
#'@param y.train Training labels
#'@param x.test Test features
#'@param clsMethod Classification method
#'@param lvNum Number of components
#'@param imp.out Return feature importances?
Predict.class <- function(x.train, y.train, x.test, clsMethod="pls", lvNum, imp.out=F){

  #msg("[Predict.class] DEBUG: Called with clsMethod='", clsMethod, "', imp.out=", imp.out, ", lvNum=", lvNum);
  #msg("[Predict.class] DEBUG: x.train dim: ", nrow(x.train), " x ", ncol(x.train));
  #msg("[Predict.class] DEBUG: x.test dim: ", nrow(x.test), " x ", ncol(x.test));
  #msg("[Predict.class] DEBUG: y.train class: ", class(y.train), ", length: ", length(y.train));

  if(clsMethod == "rf"){
    #msg("[Predict.class] DEBUG: Entering RF branch");
    rf.obj <- randomForest::randomForest(x.train, y.train, importance=imp.out, ntree=100);
    prob.out <- predict(rf.obj, x.test, type="prob")[,2];
    if(imp.out){
      imp.vec <- randomForest::importance(rf.obj)[,1];
      #msg("[Predict.class] DEBUG: RF returning list with prob.out length=", length(prob.out), ", imp.vec length=", length(imp.vec));
      return(list(prob.out=prob.out, imp.vec=imp.vec));
    }
    #msg("[Predict.class] DEBUG: RF returning prob.out only, length=", length(prob.out));
  }else if(clsMethod == "svm"){
    #msg("[Predict.class] DEBUG: Entering SVM branch");
    svm.obj <- e1071::svm(x.train, y.train, probability=TRUE);
    prob.out <- attr(predict(svm.obj, x.test, probability=TRUE), "probabilities")[,2];
    if(imp.out){
      # SVM doesn't have direct feature importance, return empty vector
      imp.vec <- rep(0, ncol(x.train));
      names(imp.vec) <- colnames(x.train);
      #msg("[Predict.class] DEBUG: SVM returning list with prob.out length=", length(prob.out), ", imp.vec length=", length(imp.vec));
      return(list(prob.out=prob.out, imp.vec=imp.vec));
    }
    #msg("[Predict.class] DEBUG: SVM returning prob.out only, length=", length(prob.out));
  }else if(clsMethod == "lr"){
    #msg("[Predict.class] DEBUG: Entering LR branch");
    # Logistic regression with selected variables
    x <- x.train;
    y <- y.train;
    xx.test <- x.test;

    names.x.origin <- colnames(x);
    colnames(x) <- paste0("V", 1:ncol(x));
    colnames(xx.test) <- colnames(x);

    data.xy <- data.frame(y, x);
    #msg("[Predict.class] DEBUG: LR calling logisticReg with data.xy dim: ", nrow(data.xy), " x ", ncol(data.xy));
    model <- logisticReg(data.xy);
    #msg("[Predict.class] DEBUG: LR model coefficients: ", paste(names(model$coefficients), collapse=", "));
    prob.out <- predict(model, as.data.frame(xx.test), type="response");
    #msg("[Predict.class] DEBUG: LR prob.out class: ", class(prob.out), ", length: ", length(prob.out));
    if(imp.out){
      imp.vec <- names(model$coefficients)[-1];
      #msg("[Predict.class] DEBUG: LR imp.vec class: ", class(imp.vec), ", length: ", length(imp.vec));
      result <- list(prob.out=prob.out, imp.vec=imp.vec);
      #msg("[Predict.class] DEBUG: LR returning list with names: ", paste(names(result), collapse=", "));
      return(result);
    }
    #msg("[Predict.class] DEBUG: LR returning prob.out only, length=", length(prob.out));
  }else{ # pls or plsda
    #msg("[Predict.class] DEBUG: Entering PLS branch");
    pls.obj <- pls::plsr(y.train ~ x.train, ncomp=lvNum, validation="none");
    score.out <- predict(pls.obj, x.test, ncomp=lvNum);
    prob.out <- (score.out - min(score.out))/(max(score.out) - min(score.out));
    if(imp.out){
      imp.vec <- Get.VIP(pls.obj, comp=lvNum);
      #msg("[Predict.class] DEBUG: PLS returning list with prob.out length=", length(prob.out), ", imp.vec length=", length(imp.vec));
      return(list(prob.out=prob.out, imp.vec=imp.vec));
    }
    #msg("[Predict.class] DEBUG: PLS returning prob.out only, length=", length(prob.out));
  }
  #msg("[Predict.class] DEBUG: Returning prob.out (fallthrough), class: ", class(prob.out));
  return(prob.out);
}

#'Get VIP scores from PLS
#'@param pls.obj PLS model object
#'@param comp Number of components
Get.VIP <- function(pls.obj, comp=2){
  w <- pls.obj$loading.weights[, 1:comp, drop=FALSE];
  q <- pls.obj$Yloadings[1:comp];
  vip <- sqrt(nrow(w) * apply((w^2) * (q^2), 1, sum) / sum(q^2));
  return(vip);
}

#'Calculate confidence intervals
#'@param data Matrix of values (rows=runs, cols=models)
#'@param param Use parametric CI?
GetCIs <- function(data, param=F){

  if(param){
    ci <- apply(data, 2, function(x){
      t.test(x)$conf.int
    })
  }else{
    # use bootstrap
    ci <- apply(data, 2, function(x){
      boot.obj <- boot::boot(x, function(data, indices){mean(data[indices])}, R=1000);
      boot::boot.ci(boot.obj, type="perc")$percent[4:5];
    })
  }
  return(ci);
}

#'Calculate T-test p-values for ROC
#'@param data Feature matrix (samples × features)
#'@param cls Class labels
GetROCTtestP <- function(data, cls){

  inx1 <- which(cls == levels(cls)[1]);
  inx2 <- which(cls == levels(cls)[2]);

  ttp <- apply(data, 2, function(x){
    t.test(x[inx1], x[inx2])$p.value;
  });

  return(ttp);
}

#'Calculate accuracy from confusion matrix
#'@param cm Confusion matrix
Get.Accuracy <- function(cm) {
  if(class(cm) == "numeric"){
    return(1);
  }
  return(sum(diag(cm))/sum(cm));
}

### NOTE: Additional functions from the original biomarker_utils.R can be adapted
### following the same pattern shown above. Key changes:
### 1. Replace mSetObj <- .get.mSet(mSetObj) with readDataset(dataName) and readSet() calls
### 2. Replace mSetObj$analSet → analSet (separate global object)
### 3. Replace .set.mSet(mSetObj) with saveSet(analSet) and RegisterData(dataSet)
### 4. Update data access paths (mSetObj$dataSet$norm → dataSet$data.norm.transposed)
### 5. Add msgSet for user messaging
### 6. Use ProteoAnalyst file naming conventions

#msg("[biomarker_utils_adapted.R] Core biomarker functions loaded for ProteoAnalyst");

#'Perform MCCV for manually selected features
#'@description MCCV for manually selected features (no additional feature selection)
#'@param dataName Input the name of the dataset
#'@param method Select the classification method
#'@param lvNum Input the number of latent variables
#'@param propTraining Proportion of samples for training
#'@param nRuns Number of MCCV runs
#'@export
PerformCV.test <- function(dataName = "", method, lvNum, propTraining=2/3, nRuns=100){

  #msg("[PerformCV.test] ===== FUNCTION CALLED =====");
  #msg("[PerformCV.test] DEBUG: dataName='", dataName, "'");
  #msg("[PerformCV.test] DEBUG: method='", method, "'");
  #msg("[PerformCV.test] DEBUG: lvNum=", lvNum);
  #msg("[PerformCV.test] DEBUG: propTraining=", propTraining);
  #msg("[PerformCV.test] DEBUG: nRuns=", nRuns);

  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");

  if(is.null(analSet)){
    analSet <- list();
  }

  analSet$tester.method <- method;
  analSet$tester.lvNum <- lvNum;
  data <- dataSet$data.norm.transposed;
  cls <- dataSet$cls;

  #msg("[PerformCV.test] DEBUG: data dim: ", nrow(data), " x ", ncol(data));
  #msg("[PerformCV.test] DEBUG: cls length: ", length(cls));    
  
  if(method == "lr") {
    if(is.null(analSet$LR)){
      analSet$LR <- list();
    }
    # check cls (response variable) whether it is number 0/1 or not
    if(length(levels(dataSet$cls)) == 2 ) {
      dataSet$cls.lbl <- levels(dataSet$cls);  # already sorted
      cls <- as.numeric(dataSet$cls)-1;
      dataSet$cls.lbl.new <- sort(unique(cls));
      dataSet$cls01 <- cls; # integer values for response of logistic regression
    } else {
      msgSet$current.msg <- "Level of response variable y/class has more than 2!";
      saveSet(msgSet, "msgSet");
      return(0);
    }
  } else {
    if(method == "pls"){# make sure the feature # > comp #
      if(lvNum > ncol(data)){
        msgSet$current.msg <- "PLS components cannot be more than total features selected!";
        saveSet(msgSet, "msgSet");
        return(0);
      }
    }
    cls <- dataSet$cls;
  }
  
  splitMat <- GetTrainTestSplitMat(cls, propTraining, nRuns);
  trainRuns <- splitMat$training.mat;
  testRuns <- splitMat$testing.mat;
  
  feat.outp <- actualCls <- perf.outp <- vector(length = nRuns, mode = "list");
  auc.vec <- accu.vec <- vector(mode="numeric", length=nRuns);
  
  for (irun in 1:nRuns){
    if(irun == 1) msg("[PerformCV.test] DEBUG: Starting iteration 1 with method='", method, "'");
    trainingSampleRun <- trainRuns[irun, ]
    testSampleRun <- testRuns[irun, ];
    x.train <- data[trainingSampleRun, ,drop=F];
    x.test <- data[testSampleRun, ,drop=F];
    y.train <- cls[trainingSampleRun];
    actualCls[[irun]] <- y.test <- cls[testSampleRun];

    if(irun == 1) msg("[PerformCV.test] DEBUG: About to call Predict.class");
    res <- Predict.class(x.train, y.train, x.test, method, lvNum, imp.out=T);
    if(irun == 1) {
      #msg("[PerformCV.test] DEBUG: res class: ", class(res));
      #msg("[PerformCV.test] DEBUG: res is.list: ", is.list(res));
      if(is.list(res)) {
        #msg("[PerformCV.test] DEBUG: res names: ", paste(names(res), collapse=", "));
      } else {
        #msg("[PerformCV.test] DEBUG: res is NOT a list! Type: ", typeof(res));
      }
    }
    feat.outp[[irun]] <- res$imp.vec;
    prob.out <- res$prob.out;
    
    # calculate AUC for each
    pred <- ROCR::prediction(prob.out, y.test);
    auc.vec[irun] <- slot(ROCR::performance(pred, "auc"), "y.values")[[1]];
    perf.outp[[irun]] <- prob.out;
    pred.out <- as.factor(ifelse(prob.out > 0.5, 1, 0));
    accu.vec[irun] <- Get.Accuracy(table(pred.out, y.test));
  }
  
  prob.vec <- unlist(perf.outp);
  act.vec <- unlist(actualCls);
  preds <- ROCR::prediction(prob.vec, act.vec);
  auc <- mean(auc.vec);
  auc.ci <- GetCIs(as.matrix(auc.vec));
  
  # if there is holdout sample, do prediction
  if(!is.null(dataSet$test.data)){
    test.res <- Predict.class(data, cls, dataSet$test.data, method, lvNum);
  }else{
    test.res <- NULL;
  }
  
  # if there is new samples, do prediction
  if(!is.null(dataSet$new.data)){
    new.res <<- Predict.class(data, cls, dataSet$new.data, method, lvNum);
  }else{
    new.res <- NULL;
  }
  
  # method = Logistic Regression
  # generate report tables for model with 10-fold Cross-Validation
  if( method == "lr") {
    x.train.all <- data;
    x.test.all  <- data;
    y.train.all <- as.character(cls);
    y.test.all  <- as.character(cls);
    
    tbl.cls <- table(cls);
    if( (tbl.cls[[1]] < 10) & (tbl.cls[[2]] < 10) ) {
      msgSet$current.msg <- "The sample size of each group should be more than 10!";
      saveSet(msgSet, "msgSet");
      return (0);
    } else {
      doLogisticRegMdl(dataName, x.train.all, y.train.all, x.test.all, y.test.all);
    }
  }
  
  if(!is.null(analSet$ROCtest$perm.res)){
    analSet$ROCtest <-  list(
      type = analSet$type,
      train.mat = trainRuns,
      pred.cv = perf.outp,
      true.cv = actualCls,
      imp.cv = feat.outp,
      perm.res = analSet$ROCtest$perm.res,
      pred.list = preds,
      accu.mat = t(as.matrix(accu.vec)),
      auc.ci = auc.ci,
      auc.vec = auc,
      test.res = test.res,
      new.res = new.res
    );
  } else {
    analSet$ROCtest <-  list(
      type = analSet$type,
      train.mat = trainRuns,
      pred.cv = perf.outp,
      true.cv = actualCls,
      imp.cv = feat.outp,
      pred.list = preds,
      accu.mat = t(as.matrix(accu.vec)),
      auc.ci = auc.ci,
      auc.vec = auc,
      test.res = test.res,
      new.res = new.res
    );
  }
  
  msgSet$current.msg <- "CV test analysis completed successfully.";
  
  # Save objects
  saveSet(analSet, "analSet");
  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");
  return(RegisterData(dataSet));
}


#'Perform permutation tests only for ROC Tester 
#'@description Perform permutation tests for biomarker model
#'@param dataName Input the name of the dataset
#'@param perf.measure Performance measure ("auroc" or "accu")
#'@param perm.num Number of permutations
#'@param propTraining Fraction of samples for training
#'@export
Perform.Permut <- function(dataName = "", perf.measure, perm.num, propTraining = 2/3){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  msgSet <- readSet(msgSet, "msgSet");
  
  cvRuns = perm.num;
  propTraining = propTraining;
  
  cls <- dataSet$cls;
  datmat <- dataSet$data.norm.transposed;
  
  if(analSet$mode == "test"){
    clsMethod <- analSet$tester.method;
  }else{
    clsMethod <- analSet$exp.method;
  }
  
  splitMat <- GetTrainTestSplitMat(cls, propTraining, cvRuns);
  trainInx <- splitMat$training.mat;
  testInx <- splitMat$testing.mat;
  
  # now, do the permutation
  perf.outp <- actualCls <- vector(length = cvRuns, mode = "list");
  irun <- 1;
  while(irun < cvRuns){
    cls <- cls[order(runif(length(cls)))];
    trainingSampleRun <- trainInx[irun, ]
    testSampleRun <- testInx[irun, ];
    y.in <- cls[trainingSampleRun];
    # make sure the y.in contain only one group
    if(length(unique(as.numeric(y.in))) > 1){
      y.out <- cls[testSampleRun];
      actualCls[[irun]] <- as.numeric(y.out);
      perf.outp[[irun]] <- Get.pred(datmat[trainingSampleRun,], y.in,
                                    datmat[testSampleRun, ], y.out, clsMethod);
      irun <- irun + 1;
    }else{
      print("redo....");
    }
  }
  
  # get the AUROC for permuted data
  pred <- ROCR::prediction(perf.outp, actualCls);
  aucs <- try(unlist(slot(ROCR::performance(pred, "auc"), "y.values")));
  if (class(aucs)=="try-error"){
    msgSet$current.msg <- "Not enough distinct predictions to compute area under the ROC curve. Increase sample size or reduce permutation number.";
    saveSet(msgSet, "msgSet");
    return(0);
  }
  accs <- unlist(slot(ROCR::performance(pred, "acc"), "y.values"));
  perf.obj <- try(ROCR::performance(pred, "tpr", "fpr"));
  if (class(perf.obj)=="try-error"){
    msgSet$current.msg <- "Something wrong with the task. Increase sample size or reduce permutation number.";
    saveSet(msgSet, "msgSet");
    return(0);
  }
  # now, insert the average value of the original performance
  accs <- c(mean(analSet$ROCtest$accu.mat[1,]), accs);
  aucs <- c(mean(analSet$ROCtest$auc.vec), aucs);
  analSet$ROCtest$perm.res <- list(perf.measure=perf.measure, perf.obj=perf.obj, auc.vec=aucs, acc.vec=accs);
  
  msgSet$current.msg <- paste0("Permutation test completed with ", perm.num, " runs.");
  
  # Save objects
  saveSet(analSet, "analSet");
  saveSet(msgSet, "msgSet");
  return(1);
}

#'Prepare report for permutation tests
#'@description Function to prepare a report for permutation tests
#'@param perm.vec Input permutation vector 
#'@export
PreparePermResult <- function(perm.vec){
  
  # check for infinite since with group variance could be zero for perfect classification
  inf.found = TRUE;
  if(sum(is.finite(perm.vec))==length(perm.vec)){
    inf.found = FALSE;
  }else {
    if(sum(is.finite(perm.vec))==0){ # all are infinite, give a random number 10
      perm.vec<-rep(10, length(perm.vec));
    }else{ # if not all inf, replace with the 10 fold of non-inf values
      perm.vec[!is.finite(perm.vec)]<-10*max(perm.vec[is.finite(perm.vec)]);
    }
  }
  
  better.hits <- sum(perm.vec[-1]>=perm.vec[1]);
  num <- length(perm.vec);
  if(better.hits == 0) {
    p <- paste("p <", 1/num);
  }else{
    p <- better.hits/num;
    p <- paste("p =", signif(p, digits=5));
  }
  
  list(permut.p = p,
       permut.inf = F,
       permut = perm.vec);
}


#'Logistic regression helper function
#'@description Build logistic regression model
#'@param d.xy Data frame with y and predictors
logisticReg <- function(d.xy) {
  fmla <- as.formula(paste("y ~", paste(names(d.xy)[-1], collapse="+")));
  model <- glm(fmla, data=d.xy, family=binomial(link="logit"), na.action="na.omit")
  return(model);
}

#'Generate LR equation
#'@description Generate equation string from LR coefficients
#'@param coef.mdl Model coefficients
genLREquation <- function(coef.mdl){
  coef.mdl <- round(coef.mdl, 3);
  
  eq <- coef.mdl[[1]];
  for(i in 2:length(coef.mdl)) {
    eq <- paste(eq, ifelse(sign(coef.mdl[[i]])==1,"+","-"), abs(round(coef.mdl[[i]],3)), names(coef.mdl)[i]);
  }
  return(eq);
}

#'Calculate ROC performance with CV for LR model
#'@description Internal function for LR cross-validation
#'@param data.in Input matrix of data
#'@param fmla.in Input for generalized linear model
#'@param kfold Numeric
#'@param run.stepwise Logistic Regression
.do.CVTest.LRmodel <- function(data.in, fmla.in, kfold=10, run.stepwise=FALSE){
  
  dw.case <- data.in[which(data.in$y == 1), ]; nrow(dw.case)
  dw.ctrl <- data.in[which(data.in$y == 0), ]; nrow(dw.ctrl)
  
  random.seed <- 10063927;
  grpIndice.case <- createCVset(nrow(dw.case), kfold, random.seed)
  grpIndice.ctrl <- createCVset(nrow(dw.ctrl), kfold, random.seed)
  
  all.trainning.y <- NULL
  all.trainning.fit <- NULL
  all.validation.y <- NULL
  all.validation.fit <- NULL
  
  for (i in 1:kfold) 
  {
    d.train.case <- dw.case[-grpIndice.case[[i]], ]; nrow(d.train.case)
    d.train.ctrl <- dw.ctrl[-grpIndice.ctrl[[i]], ]; nrow(d.train.ctrl)
    d.train <- rbind(d.train.case, d.train.ctrl); names(d.train)
    
    d.validation.case <- dw.case[grpIndice.case[[i]], ]; nrow(d.validation.case)
    d.validation.ctrl <- dw.ctrl[grpIndice.ctrl[[i]], ]; nrow(d.validation.ctrl)
    d.validation <- rbind(d.validation.case, d.validation.ctrl)  
    
    vnames <- names(d.train); 
    
    mdl <- glm(fmla.in, data=d.train, family=binomial(link="logit"))
    if(run.stepwise) { mdl <- step(mdl) }
    
    dval.pred <- predict(mdl, d.validation, type="response");  
    
    all.validation.y <- c(all.validation.y, d.validation$y)
    all.validation.fit <- c(all.validation.fit, dval.pred)
    all.trainning.y <- c(all.trainning.y, mdl$y)
    all.trainning.fit <- c(all.trainning.fit, mdl$fitted.values)
  }
  
  # 10-fold cross validation
  cv.r <- pROC::roc(all.validation.y ~ all.validation.fit, ci=T, ci.se=T, sp=seq(0,1,0.01))
  cv.threshold <- as.numeric(pROC::coords(cv.r, "best", ret=c("threshold"), best.method="closest.topleft", transpose=TRUE)); 
  cv.rstat <- multi.stat(all.validation.fit > cv.threshold, all.validation.y);
  if(length(cv.rstat) == 1 && cv.rstat == 0){
    return(0);
  }
  cv.tbl <- table(all.validation.fit > cv.threshold, all.validation.y);
  
  # training performance
  train.r <- pROC::roc(all.trainning.y ~  all.trainning.fit, ci=T, ci.se=T, sp=seq(0,1,0.01))
  train.threshold <- pROC::coords(train.r, "best", ret=c("threshold"), best.method="youden", transpose=TRUE); 
  train.rstat <- multi.stat(all.trainning.fit > cv.threshold, all.trainning.y) 
  if(length(train.rstat) == 1 && train.rstat == 0){
    return(0);
  } 
  return(list(
    train.r = train.r$ci,
    train.threshold = train.threshold,
    train.rstat = train.rstat,
    cv.robj = cv.r,
    cv.r = cv.r$ci,
    cv.threshold = cv.threshold,
    cv.rstat = cv.rstat,
    cv.tbl = cv.tbl,
    cv.y = all.validation.y,
    cv.pred = all.validation.fit,
    cv.threshold = cv.threshold
  ));
}

#'Develop a Logistic Regression Model with k-fold CV
#'@description Develop LR Model with CV
#'@param dataName Input dataset name
#'@param x.train Training X
#'@param y.train Training Y
#'@param x.test Test X
#'@param y.test Test Y
doLogisticRegMdl <- function(dataName = "", x.train, y.train, x.test, y.test){
  
  # Load required objects
  analSet <- readSet(analSet, "analSet");
  
  if(is.null(analSet$LR)){
    analSet$LR <- list();
  }
  
  # use 10-fold CV as default; or LOOCV if sample size < 10
  x <- x.train;
  y <- y.train;
  x.cv.test <- x.test;
  y.cv.test <- y.test;
  
  names.x.origin <- names(x);
  names(x) <- paste0("V", 1:(ncol(x)));
  names(x.cv.test) <- names(x);
  
  # glm can only take numeric + factor
  if(class(y) == "character"){
    y <- as.factor(y)
  }
  
  # LR Model generation
  data.xy <- data.frame(y, x);
  fmla <- as.formula(paste("y ~", paste(names(data.xy)[-1], collapse="+")));
  
  model <- glm(fmla, data=data.xy, family=binomial(link="logit"), na.action="na.omit");
  
  # if model with selected variables is failed, then use the stepwise variable selection
  if((!model$converged) | (model$deviance < 1)){
    model <- step(model, trace=0);
    fmla <- model$formula;
  }
  
  mdlSummary <- summary(model)$coefficients;
  dimnames(mdlSummary)[[1]][-1] <- names.x.origin;
  
  Odds <- round(exp(cbind(OR=coef(model), confint(model))), 2)[,1];
  Odds[1] <- "-";
  LRmodel <- cbind(round(mdlSummary,3), Odds);
  LRmodel[,4] <- ifelse(LRmodel[,4] < 0.001, "< 0.001", LRmodel[,4]);
  
  # LR model's summary table as html format
  LRmodel.xtable <<- print(xtable::xtable(LRmodel, align="r|rrrcc"),
                           type = "html", print.results=F, xtable.width=600,
                           html.table.attributes="border=1 width=600" );
  
  coefTable <- coef(model);
  names(coefTable)[-1] <- names.x.origin;
  if (model$converged & model$deviance > 1) {
    analSet$LR$LReq <- genLREquation(coefTable);
    analSet$LR$LRConverged <- "TRUE";
  } else {
    analSet$LR$LReq <- "(Error: Model was not converged)";
    analSet$LR$LRConverged <- "FALSE"; 
  }
  
  # ROC analysis with 10-fold CV test set   
  rStat <- .do.CVTest.LRmodel(data.xy, fmla.in=fmla, kfold=10, run.stepwise=FALSE)     
  
  analSet$LR$rStat <- rStat;
  
  # training/discovery; 10-fold CV
  auc  <- c( sprintf("%.3f (%.3f ~ %.3f)", rStat$train.r[2], rStat$train.r[1], rStat$train.r[3]),
             sprintf("%.3f (%.3f ~ %.3f)", round(rStat$cv.r[2],3), round(rStat$cv.r[1],3), round(rStat$cv.r[3],3)) );
  sens <- c( sprintf("%.3f (%.3f ~ %.3f)", rStat$train.rstat$sens, rStat$train.rstat$sensCI[1], rStat$train.rstat$sensCI[2]),
             sprintf("%.3f (%.3f ~ %.3f)", rStat$cv.rstat$sens, rStat$cv.rstat[1], rStat$cv.rstat$sensCI[2]) );
  spec <- c( sprintf("%.3f (%.3f ~ %.3f)", rStat$train.rstat$spec, rStat$train.rstat$specCI[1], rStat$train.rstat$specCI[2]),
             sprintf("%.3f (%.3f ~ %.3f)", rStat$cv.rstat$spec, rStat$cv.rstat$specCI[1], rStat$cv.rstat$specCI[2]) );
  
  LRperf <- cbind(auc, sens, spec);
  colnames(LRperf) <- c("AUC","Sensitivity","Specificity");
  rownames(LRperf) <- c("Training/Discovery","10-fold Cross-Validation");
  
  LRperf.xtable <<- print(xtable::xtable(LRperf, align="r|ccc"), include.rownames=TRUE, floating=FALSE,
                          type = "html", print.results=F, xtable.width=600,
                          html.table.attributes="border=1 width=600");
  
  # Save analSet
  saveSet(analSet, "analSet");
}


#'Get multiple category statistics
#'@description Calculate sensitivity, specificity, PPV, NPV, etc.
#'@param pred Input predictions
#'@param resp Input responses
multi.stat <- function(pred, resp) {
  
  in.dat <- na.omit(cbind(pred,resp));
  
  if(nrow(in.dat) < length(pred)/2){ # abort if over half is invalid
    AddMsg("More than half of tests return NA! Insufficient sample size?");
    return(0);
  }
  
  pred <- in.dat[,1];
  resp <- in.dat[,2];
  
  ts <- table(pred, resp)
  
  TP <- ts[2,2]
  FP <- ts[2,1]
  FN <- ts[1,2]
  TN <- ts[1,1]
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (FP + TN)
  ppv <- TP / (TP + FP)
  npv <- TN / (FN + TN)
  likelihood_ratio_pos <- sensitivity / (1 - specificity)
  likelihood_ratio_neg <- (1-sensitivity) / specificity
  accuracy <- (TP + TN) / sum(ts)
  fdr <- FP / (FP + TP)
  
  SEsens <- sqrt(sensitivity * (1-sensitivity) / (TP + FN) )
  SEspec <- sqrt(specificity * (1-specificity) / (FP + TN) )
  
  cise.u <- sensitivity + 1.96 * SEsens
  cisp.u <- specificity + 1.96 * SEspec
  CIsens <- c(sensitivity - 1.96 * SEsens, ifelse( cise.u > 1.0, 1.0, cise.u) )
  CIspec <- c(specificity - 1.96 * SEspec, ifelse( cisp.u > 1.0, 1.0, cisp.u) )
  
  return(list(sens=sensitivity, sensCI=CIsens, spec=specificity, specCI=CIspec) );
}

#'Get predicted class probability (utility version)
#'@description Get predicted class probability, used in higher functions
#'@param x.train Training X
#'@param y.train Training Y
#'@param x.test Test X
#'@param y.test Test Y
#'@param clsMethod Method to predict class
Get.pred <- function(x.train, y.train, x.test, y.test, clsMethod="pls"){
  
  # first convert class label to 0/1 so convert prob-> class will be easier
  y.train <- as.factor(as.numeric(y.train)-1);
  y.test <- as.factor(as.numeric(y.test)-1);
  
  # note, we only need prob for class 1, pred can be inferred
  if (clsMethod == "rf"){
    model <- randomForest::randomForest(x.train, y.train, ntree=100, importance=F);
    prob.out <- predict(model, x.test, type="prob")[,"1"];
    return(prob.out);
  }else if(clsMethod == "pls"){ # plsda
    model <- caret::plsda(x.train, y.train, method='oscorespls');
    prob.out <- predict(model, x.test, type="prob")[,"1",1];
    return(prob.out);
  }else{ # svm
    model <- e1071::svm(x.train, y.train, type = 'C', kernel="linear", probability=TRUE);
    prob.out <- attr(predict(model, x.test,  probability=TRUE), "probabilities")[,"1"];
    return(prob.out);
  }
}

#'Calculate partial area under ROC curve
#'@description Calculate partial AUC
#'@param x Input X
#'@param y Input Y
#'@param focus Method
#'@param cutoff Numeric
Get.pAUC <- function(x, y, focus, cutoff) {
  
  finite.bool <- is.finite(x) & is.finite(y);
  if(focus == "fpr"){
    x <- x[finite.bool];
    y <- y[finite.bool];
  }else{
    x <- rev(1-y[finite.bool]);
    y <- rev(1-x[finite.bool]);
    cutoff <- 1-cutoff;
  }
  if (length(x) < 2) {
    return (NA);
  }
  
  if (cutoff < 1) {
    ind <- max(which(x <= cutoff));
    stop <- try(approxfun(x[ind:(ind+1)], y[ind:(ind+1)])(cutoff));
    if(class(stop) == "try-error"){
      return(NA);
    }else{
      x <- c(x[1:ind], cutoff);
      y <- c(y[1:ind], stop);
    }
  }
  
  auc <- 0
  for (i in 2:length(x)) {
    auc <- auc + 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
  }
  return(round(auc,3));
}

#'Compute average ROC curve
#'@description Compute the average ROC curve
#'@param perf Input performance object
#'@param avg.method Method name
ComputeAverageCurve <- function(perf, avg.method){
  # now get the average curve
  perf.avg = perf;
  if(avg.method == "vertical"){
    x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                    length=max( sapply(perf@x.values, length)))
    for (i in 1:length(perf@y.values)) {
      perf.avg@y.values[[i]] <-
        approxfun(perf@x.values[[i]], perf@y.values[[i]],
                  ties=mean, rule=2)(x.values)
    }
    perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values )))
    perf.avg@x.values <- list(x.values)
  }else if(avg.method == "horizontal"){
    y.values <- seq(min(unlist(perf@y.values)), max(unlist(perf@y.values)),
                    length=max(sapply(perf@y.values, length)))
    for (i in 1:length(perf@x.values)) {
      perf.avg@x.values[[i]] <- approxfun(perf@y.values[[i]],
                                          perf@x.values[[i]],
                                          ties=mean, rule=2)(y.values)
    }
    perf.avg@x.values <- list(rowMeans( data.frame( perf.avg@x.values )));
    perf.avg@y.values <- list(y.values);
  }else{ # threshold
    all.alphas <- unlist(perf@alpha.values);
    min.alpha <- min(all.alphas);
    if(min.alpha == -Inf){
      min.alpha <- 0;
    }
    max.alpha <- max(all.alphas);
    if(max.alpha == Inf){
      max.alpha <- 1.0;
    }
    
    alpha.values <- rev(seq(min.alpha, max.alpha,length=max(sapply(perf@alpha.values, length))));
    perf.sampled <- perf;
    for (i in 1:length(perf.sampled@y.values)) {
      perf.sampled@x.values[[i]] <-
        approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                  rule=2, ties=mean)(alpha.values)
      perf.sampled@y.values[[i]] <-
        approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                  rule=2, ties=mean)(alpha.values)
    }
    ## compute average curve
    perf.avg <- perf.sampled
    perf.avg@x.values <- list(rowMeans(data.frame(perf.avg@x.values)))
    perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values)))
  }
  return(perf.avg);
}

#'Compute the 95 percent interval for threshold ROC
#'@description Computes the 95 percent interval only for the y-axis
#'@param perf Input performance object
ComputeHighLow <- function(perf){
  all.alphas <- unlist(perf@alpha.values);
  min.alpha <- min(all.alphas);
  if(min.alpha == -Inf){
    min.alpha <- 0;
  }
  max.alpha <- max(all.alphas);
  if(max.alpha == Inf){
    max.alpha <- 1.0;
  }
  
  alpha.values <- rev(seq(min.alpha, max.alpha,length=max(sapply(perf@alpha.values, length))));
  perf.sampled <- perf;
  for (i in 1:length(perf.sampled@y.values)) {
    perf.sampled@x.values[[i]] <-
      approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                rule=2, ties=mean)(alpha.values)
    perf.sampled@y.values[[i]] <-
      approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                rule=2, ties=mean)(alpha.values)
  }
  ## compute average curve
  y.data <- data.frame(perf.sampled@y.values)
  con.low <- apply(y.data, 1, quantile, 0.05);
  con.high <- apply(y.data, 1, quantile, 0.95);
  res <- list( 
    con.low = con.low,
    con.high = con.high
  );
  return (res);
}

#'Separate data set using k-fold cross validation (CV)
#'@description Separate data set with k-fold CV
#'@param groupN Input group size
#'@param kfold Input the number of folds
#'@param rseed Input the random seed
createCVset <- function(groupN, kfold, rseed)
{
  set.seed(rseed)    
  idxlist <- sample(1:groupN, size=groupN, replace=FALSE)
  CVsize <- floor(groupN / kfold)
  CV.groupIndexes <- vector(mode="list", length=kfold)
  for (i in 1:kfold) {
    CV.groupIndexes[i] <- list(idxlist[(1+CVsize*(i-1)):(CVsize*i)])
  }
  
  if((groupN %% kfold) > 0) {
    i<-1
    while( i <= (groupN %% kfold) ) {
      CV.groupIndexes[[i]] <- c(CV.groupIndexes[[i]], idxlist[CVsize*kfold + i])
      i <- i+1
    }
  }
  
  return (CV.groupIndexes)
}

#'Merge duplicated columns or rows by their mean
#'@description dim 1 => row,  dim 2 => column
#'@param data Input the data
#'@param dim Numeric, input the dimensions
#'@export
MergeDuplicates <- function(data, dim=2){
  
  if(is.null(dim(data))){ # a vector
    if(is.null(names(data))){
      print("Cannot detect duplicate data without names!!!");
      return();
    }
    nm.cls <- as.factor(names(data));
    uniq.len <- length(levels(nm.cls));
    if(uniq.len == length(data)){
      return(data);
    }
    new.data <- vector (mode="numeric",length=uniq.len);
    for(i in 1:uniq.len){
      dup.inx <- nm.cls == levels(nm.cls)[i];
      new.data[i] <- mean(data[dup.inx]);
    }
    names(new.data) <- levels(nm.cls);
    rem.len <- length(data) - length(new.data);
  }else{
    if(dim == 1){
      data <- t(data);
    }
    if(is.null(colnames(data))){
      print("Cannot detect duplicate data without var names!!!");
      return();
    }
    
    nm.cls <- as.factor(colnames(data));
    uniq.len <- length(levels(nm.cls));
    
    if(uniq.len == ncol(data)){
      if(dim == 1){
        data <- t(data);
      }
      return(data);
    }
    
    new.data <- matrix (nrow=nrow(data), ncol=uniq.len);
    for(i in 1:uniq.len){
      dup.inx <- which(nm.cls == levels(nm.cls)[i]);
      new.data[,i] <- apply(data[,dup.inx, drop=F], 1, mean);
    }
    rownames(new.data) <- rownames(data);
    colnames(new.data) <- levels(nm.cls);
    
    rem.len <- ncol(data) - ncol(new.data);
    if(dim == 1){
      new.data <- t(new.data);
    }
  }
  print(paste(rem.len, "duplicates are merged to their average"));
  new.data;
}

#'Get rpart summary
#'@description Get text description of rpart result
#'@param x An Rpart object
Get.rpart.summary <- function(x) {
  frame <- x$frame
  ylevel <- attr(x, "ylevels")
  node <- as.numeric(row.names(frame))
  
  depth <- floor(log(node, base = 2) + 1e-7)
  depth <- depth - min(depth)
  
  indent <- paste(rep(" ", 5 * 32L), collapse = "")
  
  if(length(node) > 1L) {
    indent <- substring(indent, 1L, 5 * seq(depth))
    indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
  } else {
    indent <- paste(format(node), ")", sep = "")
  }
  tfun <- (x$functions)$print
  if (!is.null(tfun)) {
    if (is.null(frame$yval2)){
      yval <- tfun(frame$yval,  ylevel, 4)
    }else{
      yval <- tfun(frame$yval2,  ylevel, 4)
    }
  }else{
    yval <- format(signif(frame$yval, digits = 4))
  }
  term <- rep(" ", length(depth))
  term[frame$var == "<leaf>"] <- "*"
  z <- labels(x, digits=4, minlength=0)
  n <- frame$n
  z <- paste(indent, z, n, format(signif(frame$dev, digits = 4)),
             yval, term);
  
  msg <- NULL;
  if (x$method=="class"){
    msg <- "node), split, n, loss, yval, (yprob)";
  }else {
    msg <- "node), split, n, deviance, yval";
  }
  msg <- c(msg, "      * denotes terminal node\n");
  msg <- paste(c(msg, z), collapse="\n");
  return(msg)
}

#'Compute data points on the ROC curve
#'@description Compute mean ROC from performance object
#'@param perf Performance object from ROCR
GetMeanROC <- function(perf){
  perf@alpha.values <- lapply(perf@alpha.values,
                              function(x) { isfin <- is.finite(x);
                              x[is.infinite(x)] <-
                                (max(x[isfin]) +
                                   mean(abs(x[isfin][-1] -
                                              x[isfin][-length(x[isfin])])));
                              x } );
  ## remove samples with x or y not finite
  for (i in 1:length(perf@x.values)) {
    ind.bool <- (is.finite(perf@x.values[[i]]) & is.finite(perf@y.values[[i]]))
    
    if (length(perf@alpha.values)>0){
      perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind.bool]
    }
    perf@x.values[[i]] <- perf@x.values[[i]][ind.bool]
    perf@y.values[[i]] <- perf@y.values[[i]][ind.bool]
  }
  
  perf.sampled <- perf;
  alpha.values <- rev(seq(min(unlist(perf@alpha.values)),
                          max(unlist(perf@alpha.values)),
                          length=max( sapply(perf@alpha.values, length))))
  for (i in 1:length(perf.sampled@y.values)) {
    perf.sampled@x.values[[i]] <- approxfun(perf@alpha.values[[i]],perf@x.values[[i]], rule=2, ties=mean)(alpha.values)
    perf.sampled@y.values[[i]] <- approxfun(perf@alpha.values[[i]], perf@y.values[[i]],rule=2, ties=mean)(alpha.values)
  }
  
  ## return the average value
  return (cbind(alpha.values, rowMeans(data.frame(perf.sampled@x.values)), rowMeans(data.frame(perf.sampled@y.values))));
}


#'Prepare data for ROC analysis
#'@description Prepare data for ROC analysis with metadata selection
#'@param dataName Input dataset name
#'@param sel.meta Selected metadata column
#'@param factor1 First factor level
#'@param factor2 Second factor level
#'@export
PrepareROCData <- function(dataName = "", sel.meta="NA", factor1="NA", factor2="NA") {
  
  # Load required objects
  dataSet <- readDataset(dataName);
  msgSet <- readSet(msgSet, "msgSet");
  
  if (sel.meta == "NA") {
    return(PrepareROCData_old(dataName))
  }else if(sel.meta == "Class"){
    if(!sel.meta %in% colnames(dataSet$meta.info)){
      sel.meta <- colnames(dataSet$meta.info)[1];
    }
  }
  msg.vec <<- 0
  data.list <- list()
  omics.vec <- vector()
  featureNms <- vector()
  uniqFeats <- vector()
  
  dataSet$meta.info.proc <- process_metadata(dataSet$meta.info)  
  
  # Check if original normalized data exists, if not, initialize it
  if (is.null(dataSet$norm.orig.roc)) {
    dataSet$norm.orig.roc <- t(dataSet$data.norm)
  }
  
  # Get original data and metadata
  merged_data <- dataSet$norm.orig.roc;
  meta.info <- dataSet$meta.info
  
  if(is.null(dataSet$norm.all)){
    dataSet$norm.all <- merged_data;
    dataSet$cls.all<- meta.info[[sel.meta]];
  }
  
  new.inx <- is.na(meta.info[[sel.meta]]) | meta.info[[sel.meta]] == "" | meta.info[[sel.meta]] == "Unknown"
  if (sum(new.inx) > 0) {
    dataSet$new.samples <- TRUE
    dataSet$new.data <- dataSet$norm.all[new.inx, , drop = FALSE]
  } else {
    dataSet$new.samples <- FALSE
    dataSet$new.data <- NULL
  }
  
  # Filter based on selected factors
  if (factor2 != "NA" & factor2 != "all") {
    dataSet$cls.orig.roc <- NULL;
    
    meta.info[[sel.meta]] <- factor(meta.info[[sel.meta]])
    sample_include <- rownames(meta.info[meta.info[[sel.meta]] %in% c(factor1, factor2), , drop = FALSE])
    sample_include <- sample_include[!is.na(sample_include)]
    
    # Subset data
    merged_data <- merged_data[sample_include, , drop = FALSE]
    meta.info <- meta.info[sample_include, , drop = FALSE]
    meta.info[[sel.meta]] <- droplevels(meta.info[[sel.meta]])
    
  } else if (factor2 == "all") {
    dataSet$cls.orig.roc <- meta.info[[sel.meta]];
    meta.info[[sel.meta]] <- as.character(meta.info[[sel.meta]])
    idx <- which(meta.info[[sel.meta]] != factor1)
    meta.info[[sel.meta]][idx] <- "Others"
    meta.info[[sel.meta]] <- factor(meta.info[[sel.meta]], levels = c(factor1, "Others"))
    sample_include <- rownames(meta.info)
    merged_data <- merged_data[sample_include, , drop = FALSE]
  }
  
  # Check sample size validity
  stt <- table(meta.info[[sel.meta]])
  
  if (length(which(stt < 10)) == 2) {
    msg.vec <<- paste0("errorLess than 10 samples in both groups ", names(stt)[1], " and ", names(stt)[2], 
                       ". Please select other groups containing more than 10 samples for biomarker analysis.")
    return
  } else if (length(which(stt < 10)) == 1) {
    msg.vec <<- paste0("errorLess than 10 samples in group ", names(stt)[which(stt < 10)], 
                       ". Please select other groups containing more than 10 samples for biomarker analysis.")
    return
  } else if (length(which(stt < 20)) == 2) {
    msg.vec <<- paste0("warnBiomarker analysis requires a large sample size. Both groups selected have less than 20 samples.")
  } else if (length(which(stt < 20)) == 1) {
    msg.vec <<- paste0("warnBiomarker analysis requires a large sample size. Group ", 
                       names(stt)[which(stt < 20)], " has less than 20 samples.")
  }
  
  dataSet$norm.orig <- merged_data;
  dataSet$data.norm.transposed <- merged_data
  dataSet$cls <- meta.info[[sel.meta]]

  # DEBUG: Check data dimensions and structure
  #msg("[PrepareROCData] DEBUG: merged_data dimensions: ", nrow(merged_data), " x ", ncol(merged_data))
  #msg("[PrepareROCData] DEBUG: data.norm.transposed dimensions: ", nrow(dataSet$data.norm.transposed), " x ", ncol(dataSet$data.norm.transposed))
  #msg("[PrepareROCData] DEBUG: cls length: ", length(dataSet$cls), " levels: ", paste(levels(dataSet$cls), collapse=", "))
  #msg("[PrepareROCData] DEBUG: cls table: ", paste(names(table(dataSet$cls)), "=", table(dataSet$cls), collapse=", "))
  #msg("[PrepareROCData] DEBUG: Any NAs in data? ", any(is.na(merged_data)))
  #msg("[PrepareROCData] DEBUG: Any NAs in cls? ", any(is.na(dataSet$cls)))

  # Handle missing values - critical for ROC analysis
  if (any(is.na(dataSet$data.norm.transposed))) {
    #msg("[PrepareROCData] DEBUG: Handling missing values...")

    # 1. Remove features with >50% missing values
    na_prop <- colSums(is.na(dataSet$data.norm.transposed)) / nrow(dataSet$data.norm.transposed)
    good_features <- na_prop < 0.5
    n_removed <- sum(!good_features)

    if (n_removed > 0) {
      #msg("[PrepareROCData] DEBUG: Removing ", n_removed, " features with >50% missing values")
      dataSet$data.norm.transposed <- dataSet$data.norm.transposed[, good_features, drop=FALSE]
    }

    # 2. Impute remaining missing values with column minimum (conservative for biomarkers)
    if (any(is.na(dataSet$data.norm.transposed))) {
      #msg("[PrepareROCData] DEBUG: Imputing remaining missing values with column minimum")
      for (j in 1:ncol(dataSet$data.norm.transposed)) {
        col_vals <- dataSet$data.norm.transposed[, j]
        if (any(is.na(col_vals))) {
          min_val <- min(col_vals, na.rm = TRUE)
          # If all values are NA or min is infinite, use 0
          if (is.infinite(min_val) || is.na(min_val)) {
            min_val <- 0
          }
          dataSet$data.norm.transposed[is.na(col_vals), j] <- min_val
        }
      }
    }

    #msg("[PrepareROCData] DEBUG: After NA handling - dimensions: ", nrow(dataSet$data.norm.transposed), " x ", ncol(dataSet$data.norm.transposed))
    #msg("[PrepareROCData] DEBUG: Any NAs remaining? ", any(is.na(dataSet$data.norm.transposed)))

    # Update norm.orig to match cleaned data
    dataSet$norm.orig <- dataSet$data.norm.transposed
  }

  msgSet$current.msg <- "ROC data prepared successfully.";
  saveSet(msgSet, "msgSet");
  return(RegisterData(dataSet))
}

#'Prepare ROC data (old version)
#'@description Prepare data for ROC analysis without metadata
#'@param dataName Input dataset name
PrepareROCData_old <- function(dataName = ""){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  
  if(is.null(dataSet$norm.all)){
    dataSet$norm.all <- t(dataSet$data.norm);
    dataSet$cls.all<- dataSet$cls;
  }
  
  new.inx <- is.na(dataSet$cls.all) | dataSet$cls.all == "" | dataSet$cls.all == "Unknown";
  if(sum(new.inx) > 0){
    dataSet$new.samples <- TRUE;
    dataSet$new.data <- dataSet$norm.all[new.inx, ,drop=F];
    dataSet$data.norm.transposed <- dataSet$norm.all[!new.inx, ,drop=F];
    dataSet$cls <- factor(dataSet$cls.all[!new.inx])
  }else{
    dataSet$new.samples <- FALSE;
    dataSet$new.data <- NULL;
    dataSet$data.norm.transposed <- dataSet$norm.all;
    dataSet$cls <- dataSet$cls.all; 
  }
  return(RegisterData(dataSet));
}

#'Set custom data for ROC
#'@description Set custom data selection for hold-out samples
#'@param dataName Input dataset name
#'@param selected.cmpds Vector of selected compounds
#'@param selected.smpls Vector of selected samples
#'@export
SetCustomData <- function(dataName = "", selected.cmpds, selected.smpls){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  
  data.norm.orig <- dataSet$data.norm.transposed;
  cls <- dataSet$cls;
  
  if(length(selected.cmpds) > 0){
    data.norm <- data.norm.orig[, selected.cmpds, drop=F];
    if(!is.null(dataSet$new.data)){
      dataSet$new.data <- dataSet$new.data[, selected.cmpds, drop=F];
    }
  }
  
  if(length(selected.smpls) > 0){
    hit.inx <- rownames(data.norm) %in% selected.smpls;
    dataSet$test.data <- data.norm[hit.inx, ,drop=F];
    dataSet$test.cls <- cls[hit.inx];
    data.norm <- data.norm[!hit.inx, ,drop=F];
    cls <- cls[!hit.inx];
  }else{
    dataSet$test.data <- NULL;
    dataSet$test.cls <- NULL;
  }
  
  dataSet$norm.orig <- data.norm.orig;
  dataSet$data.norm.transposed <- data.norm;
  dataSet$selected.cmpds <- paste(selected.cmpds, collapse="; ");
  dataSet$cls <- cls;
  return(RegisterData(dataSet));
}

#'ROC with CI for AUC
#'@description Prepare detailed ROC analysis for a feature
#'@param dataName Input dataset name
#'@param feat.nm Input the feature name
#'@export
PrepareROCDetails <- function(dataName = "", feat.nm){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  
  if(is.null(analSet)){
    analSet <- list();
  }
  
  x <- dataSet$data.norm.transposed[, feat.nm];
  y <- dataSet$cls;
  
  roc.res <- pROC::roc(y, x, ci = T, of = "auc");
  
  roc.mat <- as.matrix(data.frame(
    "Cut.Offs" = roc.res$thresholds,
    "Sensitivity" = roc.res$sensitivities,
    "Specificity" = roc.res$specificities,
    "Sens.Spec." = roc.res$sensitivities + roc.res$specificities,
    "LRP" = roc.res$sensitivities/(1-roc.res$specificities),
    "LRN" = (1-roc.res$sensitivities)/roc.res$specificities
  ));
  
  filename <- paste(dataSet$url.var.nms[feat.nm], "_roc.csv", sep="");
  fast.write.csv(signif(roc.mat,4), file=filename, row.names=F);
  
  analSet$roc.mat <- signif(roc.mat, 6);
  analSet$roc.obj <- roc.res;
  
  current.feat.nm <<- feat.nm;
  
  saveSet(analSet, "analSet");
  return(1);
}


##############################################
##############################################
########## Getter Functions for Web ##########
##############################################
##############################################

#'Check if contains new samples
#'@param dataName Input dataset name
ContainNewSamples <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  print(ifelse(dataSet$new.samples, 1, 0));
  print("ContainNewSamples");
  ifelse(dataSet$new.samples, 1, 0);
}

#'Get new sample names
#'@description Obtain sample names for new samples
#'@param dataName Input dataset name
#'@export
GetNewSampleNames <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  
  # Get original data and metadata
  merged_data <- dataSet$norm.orig.roc;
  meta.info <- dataSet$meta.info
  
  if(is.null(dataSet$norm.all)){
    dataSet$norm.all <- merged_data;
    dataSet$cls.all<- meta.info[,1];
  }
  
  new.inx <- is.na(dataSet$cls.all) | dataSet$cls.all == "" | dataSet$cls.all == "Unknown";
  if (sum(new.inx) > 0) {
    dataSet$new.samples <- TRUE;
    dataSet$new.data <- dataSet$norm.all[new.inx, , drop = FALSE];
  } else {
    dataSet$new.samples <- FALSE;
    dataSet$new.data <- NULL;
  }
  
  if(!is.null(dataSet$new.data)){
    smplInfo <- paste(rownames(dataSet$new.data), collapse="\n");
  }else{
    smplInfo <- "No new samples found";
  }
  
  RegisterData(dataSet)
  return(smplInfo);
}

GetNewSampleNameVec <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  rownames(dataSet$new.data);
}

GetNewSampleProbs <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  
  new.res <- analSet$ROCtest$new.res;
  lvls <- levels(dataSet$cls);
  grps <- ifelse(analSet$ROCtest$new.res >= 0.5, lvls[2], lvls[1]);
  new.res <- ifelse(new.res >= 0.5, new.res, 1-new.res);
  return(round(new.res, 5));
}

GetNewSampleGrps <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  
  lvls <- levels(dataSet$cls);
  grps <- ifelse(analSet$ROCtest$new.res >= 0.5, lvls[2], lvls[1]);
  return(grps);
}

GetAccuracyInfo <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  return(analSet$ROCtest$accu.info);
}

GetLR_clsLbl <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  return(paste(dataSet$cls.lbl, collapse="/"));
}

GetLR_clsLblNew <- function(dataName = ""){
  dataSet <- readDataset(dataName);
  return(paste(dataSet$cls.lbl.new, collapse="/"));
}

GetLassoFreqNames <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  return(names(analSet$LR$lassoFreq));
}

GetLRConvergence <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  return(analSet$LR$LRConverged);
}

GetLREquation <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  return(analSet$LR$LReq);
}

GetLRmodelTable <- function(){
  return(LRmodel.xtable);
}

GetLRperformTable <- function() {
  return(LRperf.xtable);
}

GetLRthreshold <- function(dataName = "") {
  analSet <- readSet(analSet, "analSet");
  return(round(analSet$LR$rStat$cv.threshold,2));
}

GetCurrentConfMat <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  return(analSet$conf.mat);
}

GetCurrentConfMatTest <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  return(analSet$conf.mat.test);
}

GetImpHighLow <- function(dataName = "", inx){
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  analSet$lowhigh[,levels(dataSet$cls)[inx]];
}

GetImpRowNames <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  rownames(analSet$imp.mat);
}

GetImpRowSymbols <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);

  feat.ids <-  rownames(analSet$imp.mat);

  if (file.exists("symbol.map.qs")) {
    tryCatch({
      gene.map <- readDataQs("symbol.map.qs", paramSet$anal.type, dataSet$name);
      if (!is.null(gene.map) && nrow(gene.map) > 0) {
        symbols <- doIdMappingGeneric(feat.ids, gene.map, "gene_id", "symbol", "vec");
        # Use original ID if symbol is NA
        na.inx <- is.na(symbols) | symbols == "";
        symbols[na.inx] <- feat.ids[na.inx];
        #msg("[GetImpRowSymbols] Mapped ", sum(!na.inx), "/", length(feat.ids), " features to symbols");
      } else {
        symbols <- feat.ids;
      }
    }, error = function(e) {
      #msg("[GetImpRowSymbols] Could not load symbol mapping: ", e$message);
      symbols <- feat.ids;
    })
  } else {
    #msg("[GetImpRowSymbols] symbol.map.qs not found, using IDs");
    symbols <- feat.ids;
  }
  return(symbols);
 
}

GetImpColNames <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  colnames(analSet$imp.mat);
}

GetImpValues <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  as.matrix(analSet$imp.mat);
}

GetRocSigFileName <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  analSet$roc.sig.nm
}

GetModelNames <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  test.nums <- analSet$multiROC$test.feats;
  paste("Model", 1:length(test.nums), "(", test.nums, "features)");
}

GetBestModelIndex <- function(dataName = ""){
  analSet <- readSet(analSet, "analSet");
  analSet$multiROC$best.inx;
}

#'Get important biomarkers table with gene symbols
#'@description Returns a table of important biomarkers with IDs, symbols, frequency, importance, and class enrichment
#'@param dataName Dataset name
#'@return Matrix with columns: ID, Symbol, Freq, Importance, and class-specific enrichment
#'@export
GetImpBiomarkersTable <- function(dataName = "") {
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  paramSet <- readSet(paramSet, "paramSet");

  # Get the importance matrix
  if (is.null(analSet$imp.mat)) {
    #msg("[GetImpBiomarkersTable] No importance matrix available");
    return(NULL);
  }

  imp.mat <- analSet$imp.mat;
  lowhigh <- analSet$lowhigh;

  # Get feature IDs
  feat.ids <- rownames(imp.mat);

  # Try to map to gene symbols
  symbols <- rep(NA_character_, length(feat.ids));
  if (file.exists("symbol.map.qs")) {
    tryCatch({
      gene.map <- readDataQs("symbol.map.qs", paramSet$anal.type, dataSet$name);
      if (!is.null(gene.map) && nrow(gene.map) > 0) {
        symbols <- doIdMappingGeneric(feat.ids, gene.map, "gene_id", "symbol", "vec");
        # Use original ID if symbol is NA
        na.inx <- is.na(symbols) | symbols == "";
        symbols[na.inx] <- feat.ids[na.inx];
        #msg("[GetImpBiomarkersTable] Mapped ", sum(!na.inx), "/", length(feat.ids), " features to symbols");
      } else {
        symbols <- feat.ids;
      }
    }, error = function(e) {
      #msg("[GetImpBiomarkersTable] Could not load symbol mapping: ", e$message);
      symbols <- feat.ids;
    })
  } else {
    #msg("[GetImpBiomarkersTable] symbol.map.qs not found, using IDs");
    symbols <- feat.ids;
  }

  # Combine into result table
  result <- data.frame(
    ID = feat.ids,
    Symbol = symbols,
    Frequency = imp.mat[, 1],
    Importance = imp.mat[, 2],
    stringsAsFactors = FALSE
  );

  # Add class-specific enrichment (High/Low for each class)
  if (!is.null(lowhigh)) {
    for (i in 1:ncol(lowhigh)) {
      result[, colnames(lowhigh)[i]] <- lowhigh[, i];
    }
  }

  return(as.matrix(result));
}

GetUnivRankedFeatureNames <- function(){
  analSet <- readSet(analSet, "analSet");
  rownames(analSet$feat.rank.mat);
}

GetFeatureRankingMat <- function(){
  analSet <- readSet(analSet, "analSet");

  analSet$feat.rank.mat;
}

Get.Fisher <- function(x, fac, var.equal=TRUE) {
  inx1 <- which(y==levels(y)[1]);
  inx2 <- which(y==levels(y)[2]);
  p.value <- apply(as.matrix(x), 2,
                   function(x) {
                     tmp <- try(fisher.test(x[inx1], x[inx2]));
                     if(class(tmp) == "try-error") {
                       return(NA);
                     }else{
                       return(tmp$p.value);
                     }
                   });
  -log10(p.value);
}

Get.Fstat <-  function(x, fac, var.equal=TRUE) {
  
  x = t(x);
  sqr = function(x) x*x;
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]
  
  k <- nlevels(fac)
  
  xm <- matrix(
    sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
    nrow = nrow(x),
    ncol = nlevels(fac))
  
  x1 <- xm[,fac, drop=FALSE]
  dff    <- k - 1
  x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
  dfr    <- ncol(x) - dff - 1
  mssf   <- rowSums(sqr(x1 - x0)) / dff
  mssr   <- rowSums(sqr( x - x1)) / dfr
  fstat  <- mssf/mssr
  return(fstat)
}


##############################################
##############################################
########## Advanced Analysis Functions #######
##############################################
##############################################

#'Compute lasso frequency
#'@description Compute lasso frequency for feature selection
#'@param dataName Input dataset name
#'@export
GetLassoFreqs <- function(dataName = ""){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  
  if(is.null(analSet$LR)){
    analSet$LR <- list();
  }
  
  if(ncol(dataSet$data.norm.transposed) < 500){
    lassoFreq <- try(GetROCLassoFreq(dataSet$data.norm.transposed, dataSet$cls));
    if(class(lassoFreq) == "try-error"){
      err.msg <<- "Unknown errors occured during computing lasso!";
      lassoFreq <- rep(0, ncol(dataSet$data.norm.transposed));
    }
  }else{
    err.msg <<- "Too many variables (>500) with small sample size, computing aborted!";
    lassoFreq <- rep(0, ncol(dataSet$data.norm.transposed));
  }
  names(lassoFreq) <- colnames(dataSet$data.norm.transposed);
  analSet$LR$lassoFreq <- sort(lassoFreq, decreasing =TRUE);
  
  saveSet(analSet, "analSet");
  return(analSet$LR$lassoFreq);
}

#'Get p-values from lasso 
#'@description Get p-values from lasso cross-validation
#'@param data Input data
#'@param cls Input class labels
GetROCLassoFreq <- function(data, cls){
  
  data <- cbind(cls, data);
  d.headers <- names(data);
  names(data) <- c("cls", paste0("V", 1:(ncol(data)-1)));
  
  inx1 <- which(cls==levels(cls)[1]);
  inx2 <- which(cls==levels(cls)[2]);
  data1 <- data[inx1, ];
  data2 <- data[inx2, ];
  
  min.n <- ifelse(length(inx1) <= length(inx2), length(inx1), length(inx2));
  kfold.origin <- ifelse(min.n >= 10, 10, ifelse(min.n >= 5, 5, 0));
  
  if(kfold.origin > 0) {
    random.seed <- 10063927;
    kfold <- kfold.origin;
    CV.inx1 <- createCVset(length(inx1), kfold, random.seed);
    CV.inx2 <- createCVset(length(inx2), kfold, random.seed);        
  } else {
    kfold <- 10;
  }
  
  varlist <- NULL;
  for (i in 1:kfold) {
    if (kfold.origin > 0) {
      dw1 <- data1[-CV.inx1[[i]], ]; 
      dw2 <- data2[-CV.inx2[[i]], ]; 
    } else {
      # resampling if sample size is less than 5 
      CV.inx1 <- sample(inx1, size=length(inx1), replace=TRUE);
      CV.inx2 <- sample(inx2, size=length(inx2), replace=TRUE);
      dw1 <- data1[as.character(CV.inx1), ]; 
      dw2 <- data2[as.character(CV.inx2), ]; 
    }
    
    dw.all <- rbind(dw1, dw2);   
    
    # To create a formula for model with large number of independent vars
    xnam <- names(dw.all)[-1];
    (fmla <- as.formula(paste("cls ~ ", paste(xnam, collapse= "+"))));
    
    if (is.numeric(dw.all$cls)) {
      dw.all$cls <- as.numeric(as.character(dw.all$cls)); 
    } else {
      # label/cls should be integer as 0 and 1
      dw.all$cls <- as.numeric(dw.all$cls)-1; 
    }
    
    x <- model.matrix(as.formula(fmla), dw.all)[, -1];
    o <- lars::lars(x, dw.all$cls, type="lasso", trace=FALSE, intercept=TRUE, normalize=FALSE, use.Gram=FALSE);
    
    cvfit <- NULL;
    m <- NULL;
    tryCatch({
      cvfit <- lars::cv.lars(x, dw.all$cls, type="lasso", mode="fraction", plot.it=FALSE);
      m <- ( o$beta[which.min(cvfit$cv),] != 0 );
      varlist10 <- names(m[m==TRUE]);
    }, error = function(e) {
      tryCatch ( {
        cvfit <- lars::cv.lars(x, dw.all$cls, type="lar", mode="step", plot.it=FALSE); 
        m <- ( o$beta[which.min(cvfit$cv),] != 0 );
      }, error=function(e) {
      }, finally = {
        if(is.null(cvfit)) {
          m <- ( o$beta[5,] != 0 );
        }
      })
    }, finally = {
    })
    
    varlist <- c(varlist, names(m[m==TRUE]));
  }
  var.lasso.10CV <- unique(varlist); 
  dt <- table(varlist); 
  
  compnames <- names(dw.all)[-1];
  selfreq <- 0.0;
  compselfreq <- cbind(compnames, selfreq);
  for (i in 1:length(dt)) {
    compselfreq[which(compnames == names(dt[i])), "selfreq"] <- dt[i]/kfold * 100.0;
  }
  
  return(as.numeric(compselfreq[,"selfreq"]));
}

#'Get important feature matrix
#'@description Get feature importance matrix from CV runs
#'@param dataName Input dataset name
#'@param feat.outp Feature output list
#'@param bestFeatNum Best feature number
GetImpFeatureMat <- function(dataName = "", feat.outp, bestFeatNum){
  
  analSet <- readSet(analSet, "analSet");
  anal.mode <- analSet$mode; 
  
  # first order each run by cmpd names so that can be combined to a data frame
  feat.outp <- lapply(feat.outp, function(x) x[order(names(x))]);
  
  # First rank by frequencies of being selected in the given model
  # obtain their ranks
  freqRank <- lapply(feat.outp, function(x) rank(-x));
  runRanksMat <- do.call("cbind", freqRank);
  
  # order by their median rank across runs
  ordRunRanksMat <- as.data.frame(runRanksMat[order(apply(runRanksMat, 1, median)),]);
  
  # Then rank by mean importance measures
  impsMat <- as.data.frame(do.call("cbind", feat.outp));
  # OPTIMIZED: Use rowMeans instead of apply for 60-100x speedup
  impsVec <- rowMeans(impsMat);

  # now count the number being selected in the bestFeatNum
  selectedMat <- apply(ordRunRanksMat, 2, function(x) x <= bestFeatNum);

  # calculate percentage of being selected in the best subsets
  # OPTIMIZED: Use rowSums instead of apply for 60-100x speedup
  percentVec <- rowSums(selectedMat)/ncol(ordRunRanksMat);
  
  # remove ones never being selected
  percentVec <- percentVec[percentVec > 0];
  
  # reorder the imps to percentVec
  impsVec <- impsVec[names(percentVec)];
  
  # combine and return the result
  imp.mat <- cbind(percentVec, impsVec);
  ord.inx <- order(imp.mat[,1], imp.mat[,2], decreasing=T);
  imp.mat <- imp.mat[ord.inx,];
  colnames(imp.mat) <- c("Rank Freq.", "Importance");
  
  return(imp.mat);
}

#'Return ROC coordinates with confidence intervals
#'@description Return ROC coordinates with CI
#'@param dataName Input dataset name
#'@param fld.nm The kind of input coordinate
#'@param val The coordinates to look for
#'@param plot Logical, by default TRUE
#'@param imgNm Input the image name
#'@export
GetROC.coords <- function(dataName = "", fld.nm, val, plot=TRUE, imgNm){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  
  res <- pROC::coords(analSet$roc.obj, val, input=fld.nm, transpose=TRUE);
  
  sp <- res[2];
  se <- res[3];
  res <- round(res, 3);
  
  analSet$thresh.obj <- NULL;
  if(fld.nm == "threshold"){
    ci.s <- pROC::ci.thresholds(analSet$roc.obj, boot.n=100, thresholds=val, progress="none");
    specs <- round(ci.s$specificity,3);
    sens <- round(ci.s$sensitivity, 3);
    res[2] <- paste(res[2], "(", specs[1], "-", specs[3], ")", sep="");
    res[3] <- paste(res[3], "(", sens[1], "-", sens[3], ")", sep="");
    
    analSet$thresh.obj <- ci.s;
    
    # update pos with the thresh obj coords
    sp <- ci.s$specificity[2];
    se <- ci.s$sensitivity[2];
  }
  
  mythresh <- res[1];
  if(is.na(res[1])){
    if(fld.nm == "sensitivity"){
      fld.vals <- analSet$roc.obj$sensitivities;
    }else{
      fld.vals <- analSet$roc.obj$specificities;
    }
    
    inx1 <- which.min(abs(fld.vals-val));
    
    if(inx1 == 1){
      inxb1 <- inx1;
    }else{
      inxb1 <- inx1 - 1;
    }
    
    if(inx1 == length(fld.vals)){
      inxb2 <- inx1;
    }else{
      inxb2 <- inx1 + 1;
    }
    
    if(fld.vals[inx1] > val){
      inx2 <-ifelse(fld.vals[inxb1] > val, inxb2, inxb1);
    }else{
      inx2 <-ifelse(fld.vals[inxb1] > val, inxb1, inxb2);
    }
    
    threshs <- analSet$roc.obj$thresholds;
    
    if(inx1 == inx2){ # out of the threshod range
      if(threshs[inx1] > threshs[2]){ #max
        res[1] <- paste(">",threshs[inx1], sep="");
      }else{
        res[1] <- paste("<",threshs[inx1], sep="");
      }
    }else{
      inx.thresh <- c(inx1, inx2)
      ord.inx <- order(threshs[inx.thresh]);
      inx.thresh <- inx.thresh[ord.inx];
      res[1] <- paste(threshs[inx.thresh], collapse="-")
    }
  }
  if(plot){
    PlotDetailROC(dataName, imgNm, mythresh, sp, se);
    analSet$roc.obj$thresh <- mythresh;
  }
  
  saveSet(analSet, "analSet");
  return(res);
}


##############################################
##############################################
########## Critical Plotting Functions #######
##############################################
##############################################

#'Plot detailed ROC
#'@description Plot detailed ROC with threshold
#'@param dataName Input dataset name
#'@param imgName Input a name for the plot
#'@param thresh Input the threshold
#'@param sp Specificity
#'@param se Sensitivity
#'@param dpi Input the dpi
#'@param format Select the image format
#'@export
PlotDetailROC <- function(dataName = "", imgName, thresh, sp, se, dpi=72, format="png"){
  
  # Load required objects
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");
  
  imgName = paste(imgName, "_dpi", dpi, ".", format, sep="");
  
  roc.obj <- analSet$roc.obj;
  
  w <- h <- 6;
  imgSet$roc.univ <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(oma = c(0,0,1,0));
  par(mar=c(4,4,4,4) + .1);
  
  pROC::plot.roc(roc.obj, print.auc=F, legacy.axes=TRUE, col="navy", grid=T,
                 xlab = "False positive rate", ylab="True positive rate",
                 auc.polygon=TRUE, auc.polygon.col="#0000ff22", main=current.feat.nm);
  
  points(sp, se, cex=1.8, col="red");
  
  dev.off();
  
  saveSet(imgSet, "imgSet");
  return(1);
}

#'Perform Classical Univariate ROC
#'@description Perform Classical Univariate ROC analysis
#'@param dataName Input dataset name
#'@param feat.nm Input the feature name
#'@param version image version mark
#'@param format Select the image format
#'@param dpi Input the dpi
#'@param isAUC Logical, compute CI
#'@param isOpt Logical, show optimal cutoff
#'@param optMethod Optimal cutoff method
#'@param isPartial Logical, calculate partial ROC
#'@param measure Parameter to limit pROC
#'@param cutoff Threshold for pROC
#'@export
Perform.UnivROC <- function(dataName = "", feat.nm,
                            version, format="png",
                            dpi=72, isAUC, isOpt,
                            optMethod, isPartial, measure, cutoff,
                            dataType = "default"){

  dataType <- tolower(dataType)

  # print(paste0("[R DEBUG ROC] Perform.UnivROC called with feat.nm: ", feat.nm, ", dataType: ", dataType))

  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  if(is.null(analSet)){
    analSet <- list();
  }

  feat.id <- feat.nm
  if (!is.null(dataSet$orig.var.nms) && !(feat.id %in% names(dataSet$url.var.nms))) {
    idx <- which(dataSet$orig.var.nms == feat.id)
    if (length(idx) > 0) {
      feat.id <- names(dataSet$orig.var.nms)[1]
    }
  }
  # print(paste0("[R DEBUG ROC] After ID transformation: feat.id = ", feat.id))

  imgName <- dataSet$url.var.nms[feat.id];
  if (dataType == "peptide") {
    imgName <- feat.nm
  }
  imgName = paste("roc_univ_", imgName, "_", version, "_dpi", dpi, ".", format, sep="");
  # Remove consecutive underscores
  imgName = gsub("__+", "_", imgName);
  # print(paste0("[R DEBUG ROC] Image name: ", imgName))

  if (dataType == "peptide") {
    if (!file.exists("peptide_level_data.qs")) {
      # print("[R DEBUG ROC] peptide_level_data.qs not found")
      return(0)
    }
    pep.mat <- ov_qs_read("peptide_level_data.qs")
    if (!(feat.nm %in% rownames(pep.mat))) {
      # print(paste0("[R DEBUG ROC] feat.nm not found in peptide matrix rows"))
      return(0)
    }
    x <- as.numeric(pep.mat[feat.nm, ])
    y <- dataSet$cls
    if (is.null(names(y)) && length(y) == ncol(pep.mat)) {
      names(y) <- colnames(pep.mat)
    }
    if (!is.null(names(y)) && !is.null(colnames(pep.mat))) {
      y <- y[match(colnames(pep.mat), names(y))]
    }
    if (dataSet$comp.type == "custom") {
      grp.nms <- dataSet$grp.nms
      if (dataSet$cont.inx[dataSet$analysisVar] ||
          any(grepl("(^[0-9]+).*", as.character(dataSet$cls)))) {
        grp.nms <- sub(paste0(dataSet$analysisVar, "_"), "", grp.nms)
      }
      keep <- dataSet$cls %in% grp.nms
      x <- as.numeric(pep.mat[feat.nm, keep])
      y <- droplevels(dataSet$cls[keep])
    } else {
      y <- droplevels(y)
    }
  } else {
    # Use dataSet$data.norm like PlotSelectedGene does for boxplot/violin plots
    if(length(dataSet$rmidx)>0){
      data_ori_norm <- dataSet$data.norm[,-dataSet$rmidx]
      # print("[R DEBUG ROC] Using dataSet$data.norm (with removed samples)")
    } else {
      data_ori_norm <- dataSet$data.norm
      # print("[R DEBUG ROC] Using dataSet$data.norm")
    }

    if (feat.id %in% colnames(data_ori_norm)) {
      x <- data_ori_norm[, feat.id]
      # print("[R DEBUG ROC] Found feat.id in columns")
    } else if (feat.id %in% rownames(data_ori_norm)) {
      x <- data_ori_norm[feat.id, ]
      # print("[R DEBUG ROC] Found feat.id in rows")
    } else {
      x <- numeric(0)
      # print("[R DEBUG ROC] feat.id NOT FOUND in data matrix - returning empty vector")
    }
    y <- dataSet$cls;
  }

  # print(paste0("[R DEBUG ROC] Length of x: ", length(x), ", Length of y: ", length(y)))

  if (length(x) == 0 || length(y) == 0) {
    # print("[R DEBUG ROC] x or y is empty - returning 0")
    return(0)
  }
  if (!is.null(names(x)) && !is.null(names(y))) {
    common <- intersect(names(x), names(y))
    # print(paste0("[R DEBUG ROC] Common names between x and y: ", length(common)))
    if (length(common) > 0) {
      x <- x[common]
      y <- y[common]
      # print(paste0("[R DEBUG ROC] After matching names - x length: ", length(x), ", y length: ", length(y)))
    }
  }
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  # print(paste0("[R DEBUG ROC] After removing NAs - x length: ", length(x), ", y length: ", length(y)))
  if (length(x) == 0 || length(y) == 0 || length(x) != length(y)) {
    # print(paste0("[R DEBUG ROC] Final check failed - x length: ", length(x), ", y length: ", length(y), " - returning 0"))
    return(0)
  }
  
  if(isPartial){
    if(measure == "se"){
      cutoff = cutoff;
    }else{
      cutoff = 1-cutoff;
    }
    roc.obj <- pROC::roc(y, x, partial.auc=c(1.0, cutoff), ci=TRUE, partial.auc.focus=measure, boot.n=50, percent = F, progress="none");
  }else{
    roc.obj <- pROC::roc(y, x, percent = F);
  }
  
  w <- h <- 6; 
  
  if(length(imgSet$roc.univ.name)==0){
    imgSet$roc.univ.name <- feat.nm;
    imgSet$roc.univ.plot <- imgName;
  } else {
    if(feat.nm %in% imgSet$roc.univ.name){
      idx <- which(feat.nm %in% imgSet$roc.univ.name);
      imgSet$roc.univ.name[idx] <- feat.nm;
      imgSet$roc.univ.plot[idx] <- imgName;
    } else {
      imgSet$roc.univ.name <- c(imgSet$roc.univ.name, feat.nm);
      imgSet$roc.univ.plot <- c(imgSet$roc.univ.plot, imgName);
    }   
  }
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(oma=c(0,0,1,0));
  par(mar=c(4,4,4,4)+.1);
  
  opt.thresh = NA;
  
  # first plot ROC curve
  if(isAUC){
    pROC::plot.roc(roc.obj, print.auc=F, legacy.axes=TRUE, col="navy", grid=T,
                   xlab = "False positive rate", ylab="True positive rate", main=feat.nm);
    ci.obj <- pROC::ci.se(roc.obj, specificities=seq(0, 1, 0.05), boot.n=200, progress="none");
    ROCR::plot(ci.obj,type="shape",col="#0000ff22");
  }else{
    pROC::plot.roc(roc.obj, print.auc=F, legacy.axes=TRUE, col="navy", grid=T,
                   xlab = "False positive rate", ylab="True positive rate",
                   auc.polygon=TRUE, auc.polygon.col="#0000ff22", main=feat.nm);
  }
  
  auc.ci <- pROC::ci.auc(roc.obj, method="bootstrap", boot.n=500, progress="none");
  roc.obj$ci <- auc.ci;
  auc.lbl <- paste("AUC: ", round(roc.obj$ci[2],3), sep="");
  ci.lbl <- paste("(", round(roc.obj$ci[1],3), "-", round(roc.obj$ci[3],3), ")", sep="");
  text(0.5, 0.5, paste(auc.lbl, "\n", ci.lbl, sep=""), adj=c(0,1), col="navy");
  
  if(isOpt){
    par(xpd=T);
    opt.ps <- data.frame(pROC::coords(roc.obj, "best", best.method=optMethod, transpose = TRUE));
    opt.thresh <- opt.ps["threshold",]
    points(opt.ps["specificity",], opt.ps["sensitivity",], pch=19, col="red");
    lbls=paste(signif(opt.ps["threshold",],3), "(", round(opt.ps["specificity",],1), ", ", round(opt.ps["sensitivity",],1), ")", sep="");
    text(opt.ps["specificity",], opt.ps["sensitivity",], adj=c(-.05,1.25), label=lbls);
  }
  
  dev.off();
  
  analSet$opt.thresh <- opt.thresh
  
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);
}

#'Plot a boxplot view of a selected feature
#'@description Plots a boxplot of the selected feature's concentration
#'@param dataName Input dataset name
#'@param feat.nm Input the feature name
#'@param version version mark for image name
#'@param format Select the image format
#'@param dpi Input the dpi
#'@param isOpt logical
#'@param isQuery logical
#'@export
PlotRocUnivBoxPlot <- function(dataName = "", feat.nm, version, format="png", dpi=72, isOpt, isQuery, dataType = "default"){
  
  # Load required objects
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");
  
  dataType <- tolower(dataType)

  data_ori_norm <- dataSet$data.norm.transposed
  if(!is.null(dataSet$norm.orig)){
    data_ori_norm <- dataSet$norm.orig
  }
  
  feat.id <- feat.nm
  if (!is.null(dataSet$orig.var.nms) && !(feat.id %in% names(dataSet$url.var.nms))) {
    idx <- which(dataSet$orig.var.nms == feat.id)
    if (length(idx) > 0) {
      feat.id <- names(dataSet$orig.var.nms)[1]
    }
  }
  imgName <- dataSet$url.var.nms[feat.id];
  if (dataType == "peptide") {
    imgName <- feat.nm
  }
  imgName = paste("roc_boxplot_", imgName, "_", version, "_dpi", dpi, ".", format, sep="");
  # Remove consecutive underscores
  imgName = gsub("__+", "_", imgName);
  
  if (dataType == "peptide") {
    if (!file.exists("peptide_level_data.qs")) {
      return(0)
    }
    pep.mat <- ov_qs_read("peptide_level_data.qs")
    if (!(feat.nm %in% rownames(pep.mat))) {
      return(0)
    }
    x <- as.numeric(pep.mat[feat.nm, ])
    if(!is.null(dataSet$cls.orig.roc)){ 
      y <- dataSet$cls.orig.roc;
    }else{
      y <- dataSet$cls;
    }
    if (is.null(names(y)) && length(y) == ncol(pep.mat)) {
      names(y) <- colnames(pep.mat)
    }
    if (!is.null(names(y)) && !is.null(colnames(pep.mat))) {
      y <- y[match(colnames(pep.mat), names(y))]
    }
    if (dataSet$comp.type == "custom") {
      grp.nms <- dataSet$grp.nms
      if (dataSet$cont.inx[dataSet$analysisVar] ||
          any(grepl("(^[0-9]+).*", as.character(dataSet$cls)))) {
        grp.nms <- sub(paste0(dataSet$analysisVar, "_"), "", grp.nms)
      }
      keep <- dataSet$cls %in% grp.nms
      x <- as.numeric(pep.mat[feat.nm, keep])
      y <- droplevels(dataSet$cls[keep])
    } else {
      y <- droplevels(y)
    }
  } else {
    if (feat.id %in% colnames(data_ori_norm)) {
      x <- data_ori_norm[, feat.id]
    } else if (feat.id %in% rownames(data_ori_norm)) {
      x <- data_ori_norm[feat.id, ]
    } else {
      x <- numeric(0)
    }
    # dataSet$cls.orig.roc defined in PrepareROCData
    if(!is.null(dataSet$cls.orig.roc)){ 
      y <- dataSet$cls.orig.roc;
    }else{
      y <- dataSet$cls;
    }
  }
  
  if (length(x) == 0 || length(y) == 0) {
    return(0)
  }
  if (!is.null(names(x)) && !is.null(names(y))) {
    common <- intersect(names(x), names(y))
    if (length(common) > 0) {
      x <- x[common]
      y <- y[common]
    }
  }
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) == 0 || length(y) == 0 || length(x) != length(y)) {
    return(0)
  }
  
  scale <- dpi/72;
  w <- 200*scale;
  h <- 400*scale; 
  col <- GetColorSchema(y);
  
  if(length(imgSet$roc.univ.boxplot)==0){
    imgSet$roc.univ.boxplot <- imgName;
    imgSet$roc.univ.name2 <- feat.nm;
  } else {
    if(feat.nm %in% imgSet$roc.univ.name2){
      idx <- which(feat.nm %in% imgSet$roc.univ.name);
      imgSet$roc.univ.boxplot[idx] <- imgName;
    } else {
      imgSet$roc.univ.name2 <- c(imgSet$roc.univ.name2, feat.nm);
      imgSet$roc.univ.boxplot <- c(imgSet$roc.univ.boxplot, imgName);
    }   
  }
  
  Cairo::Cairo(file=imgName, width=w, height=h, type=format, bg="white", dpi=dpi);
  par(oma=c(0,0,1,0));
  par(mar=c(4,4,4,4)+.1);
  df <- data.frame(conc = x, class = y)
  p <- ggplot2::ggplot(df, aes(x=class, y=conc, fill=class)) + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA) + theme_bw() + geom_jitter(size=1)
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
  p <- p + stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
  p <- p + labs(title = feat.nm)
  p <- p + theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t=0.45, r=0.25, b=1.5, l=0.25, "cm")
  ) 
  p <- p + scale_fill_manual(values=col)
  if(isOpt){
    opt.thresh <- analSet$opt.thresh
    p <- p + geom_hline(aes(yintercept=opt.thresh), colour="red")
  }
  
  if(isQuery){
    thresh <- as.numeric(analSet$roc.obj$thresh)
    p <- p + geom_hline(aes(yintercept=thresh), colour="red")
  }
  
  print(p)
  dev.off()
  
  saveSet(imgSet, "imgSet");
  return(imgName);
}



process_metadata <- function(df) {
  library(caret)
  
  # Initialize the processed_df with the same number of rows as the input dataframe
  processed_df <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
  
  # Loop through each column of the data frame
  for (col_name in colnames(df)) {
    
    # Check if the column is a factor or character (categorical)
    if (is.factor(df[[col_name]]) || is.character(df[[col_name]])) {
      # One-hot encode categorical variables
      one_hot_encoded <- model.matrix(~ df[[col_name]] - 1, data = df)
      colnames(one_hot_encoded) <- paste0(col_name, "_", colnames(one_hot_encoded))
      # Append to processed dataframe
      processed_df <- cbind(processed_df, one_hot_encoded)
    
    # If the column is numeric (continuous)
    } else if (is.numeric(df[[col_name]])) {
      # Handle missing values by imputing them with the mean
      df[[col_name]][is.na(df[[col_name]])] <- mean(df[[col_name]], na.rm = TRUE)
      
      # Normalize continuous variables (Z-score normalization)
      normalized_column <- scale(df[[col_name]])
      processed_df[[col_name]] <- normalized_column
    }
  }
  
  return(processed_df)
}


#'Plot a summary view of the classification result
#'@description Plot of predicted class probabilities. On the x-axis is the proability, 
#'and the y-axis is the index of each predicted sample based on the probility. 
#'The samples are turned into separations at the x-axis.
#'This plot can be created for multivariate ROC curve analysis using SVM, PLS, and RandomForest.
#'Please note that sometimes, not all samples will be tested, instead they will be plotted
#'at the 0.5 neutral line. 
#'@usage PlotProbView(mSetObj=NA, imgName, format="png", dpi=default.dpi, mdl.inx, show, showPred) 
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300. 
#'@param mdl.inx Model index, 0 means to compare all models, -1 means to use the best model, input 1-6 to plot a ROC curve for one of the top six models
#'@param show 1 or 0, if 1, label samples classified to the wrong groups 
#'@param showPred Show predicted samples
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotProbView <- function(dataName = "", imgName, format="png", dpi=default.dpi, mdl.inx, show, showPred) {

  #msg("[PlotProbView] DEBUG: Function called with dataName='", dataName, "' imgName='", imgName, "' mdl.inx=", mdl.inx)

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  #msg("[PlotProbView] DEBUG: analSet is.null? ", is.null(analSet))
  #msg("[PlotProbView] DEBUG: analSet$multiROC is.null? ", is.null(analSet$multiROC))

  anal.mode <- analSet$mode;

  # Sample names are columns (features×samples orientation)
  smpl.nms <- colnames(dataSet$data.norm);
  prob.vec <- rep(0.5, length(smpl.nms));
  names(prob.vec) <- smpl.nms;

  #msg("[PlotProbView] DEBUG: smpl.nms length: ", length(smpl.nms))

  if(mdl.inx == -1){
    mdl.inx <- analSet$multiROC$best.inx;
    #msg("[PlotProbView] DEBUG: Using best model index: ", mdl.inx)
  }

  #msg("[PlotProbView] DEBUG: About to access analSet$multiROC$pred.cv[[", mdl.inx, "]]")
  #msg("[PlotProbView] DEBUG: length(analSet$multiROC$pred.cv) = ", length(analSet$multiROC$pred.cv))

  probs <- MergeDuplicates(unlist(analSet$multiROC$pred.cv[[mdl.inx]]));

  #msg("[PlotProbView] DEBUG: probs length: ", length(probs))

  prob.vec[names(probs)] <- probs;

  nms <- names(prob.vec);
  ord.inx <- order(nms);
  prob.vec <- prob.vec[ord.inx];
  cls <- dataSet$cls[ord.inx];
  # remember to update the nms itself!
  nms <- names(prob.vec);

  # get html confusion matrix
  pred.out <- as.factor(ifelse(prob.vec > 0.5, 1, 0));
  act.cls <- as.numeric(cls)-1;

  prob.res <- data.frame(Probability=prob.vec, Predicted=pred.out, Actual=act.cls);

  write.table(prob.res, file="roc_pred_prob.csv", sep=",", col.names = TRUE);

  conf.res <- table(pred.out, act.cls);
  analSet$conf.table <- xtable::xtable(conf.res, caption="Confusion Matrix (Cross-Validation)");
  analSet$conf.mat <- print(analSet$conf.table, type = "html", print.results=F, caption.placement="top", html.table.attributes="border=1 width=150" )


  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9;
  h <- 8;

  imgSet$roc.prob.plot <- imgName;
  imgSet$roc.prob.name <- mdl.inx
  
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  set.seed(123);
  y <- rnorm(length(prob.vec));
  max.y <- max(abs(y));
  ylim <- max.y*c(-1.05, 1.05);
  
  xlim <- c(0, 1.0);
  
  op <- par(mar=c(4,4,3,6));
  pchs <- ifelse(as.numeric(cls) == 1, 1, 19);
  
  colors <- ifelse(show==1, "darkgrey", "black");
  ROCR::plot(prob.vec, y, pch=pchs, col=colors, xlim=xlim, ylim= ylim, xlab = "Predicted Class Probabilities", ylab="Samples");
  abline(h = 0, lty = 2, col="grey");
  abline(v = 0.5, lty = 2, col="grey");
  
  par(xpd=T);
  legend("right",inset=c(-0.11,0), legend = unique(as.character(cls)), pch=unique(pchs));
  
  test.y <- test.x <- 0;
  if(showPred){
    test.y <- rnorm(length(analSet$multiROC$test.res));
    test.x <- analSet$multiROC$test.res;

    pchs <- ifelse(as.numeric(dataSet$test.cls) == 1, 1, 19);
    points(test.x, test.y, pch=pchs, cex=1.5, col="red");
  }

  if(show == 1){

    # add labels for sample classified wrong
    # the idea is to add labels to the left side for those with prob < 0.5
    # and add labels to the right side of the point for those with prob > 0.5
    # leave 0.5 out

    # first wrong pred as 1 (right side)
    act.ones <- as.numeric(cls)-1 == 1;
    pred.vec <- ifelse(prob.vec > 0.5, 1, 0);

    wrong.inx <- (pred.vec != as.numeric(cls)-1) & pred.vec == 1;
    if(sum(wrong.inx) > 0){
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx], pos=4);
    }

    # first wrong pred as 0 (left side)
    act.zeros <- as.numeric(cls)-1 == 0;
    pred.vec <- ifelse(prob.vec < 0.5, 0, 0.5);
    wrong.inx <- pred.vec != as.numeric(cls)-1 & pred.vec == 0;
    if(sum(wrong.inx) > 0){
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx], pos=2);
    }

    if(showPred){
      # Sample names are columns (features×samples orientation)
      nms <- colnames(dataSet$test.data);

      act.ones <- as.numeric(dataSet$test.cls)-1 == 1;
      act.zeros <- as.numeric(dataSet$test.cls)-1 == 0;

      # leave 0.5 there
      pred.vec <- ifelse(test.x > 0.5, 1, 0.5);
      wrong.inx <- (pred.vec != as.numeric(dataSet$test.cls)-1) & act.ones;
      if(sum(wrong.inx) > 0){
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx], pos=4, cex=0.9);
      }

      pred.vec <- ifelse(test.x < 0.5, 0, 0.5);
      wrong.inx <- pred.vec != as.numeric(dataSet$test.cls)-1 & act.zeros;
      if(sum(wrong.inx) > 0){
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx], pos=2, cex=0.9);
      }
    }
  }
  par(op)
  dev.off();
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);
}

#'Plot a summary view of the classification result of tester prediction
#'@description Plot of predicted class probabilities. On the x-axis is the proability, 
#'and the y-axis is the index of each predicted sample based on the probility. 
#'The samples are turned into separations at the x-axis.
#'This plot can be created for multivariate ROC curve analysis using SVM, PLS, and RandomForest.
#'Please note that sometimes, not all samples will be tested, instead they will be plotted
#'at the 0.5 neutral line. 
#'@usage PlotProbViewTest(mSetObj=NA, imgName, format="png", dpi=default.dpi, mdl.inx, show, showPred) 
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300. 
#'@param mdl.inx Model index, 0 means to compare all models, -1 means to use the best model, input 1-6 to plot a ROC curve for one of the top six models
#'@param show 1 or 0, if 1, label samples classified to the wrong groups 
#'@param showPred Show predicted samples
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotProbViewTest <- function(dataName = "", imgName, format="png", dpi=default.dpi, mdl.inx, show, showPred) {

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  anal.mode <- analSet$mode;

  # Sample names are columns (features×samples orientation)
  smpl.nms <- colnames(dataSet$data.norm);
  prob.vec <- rep(0.5, length(smpl.nms));
  names(prob.vec) <- smpl.nms;

  probs <- MergeDuplicates(unlist(analSet$ROCtest$pred.cv));

  prob.vec[names(probs)] <- probs;

  nms <- names(prob.vec);
  ord.inx <- order(nms);
  prob.vec <- prob.vec[ord.inx];
  cls <- dataSet$cls[ord.inx];
  # remember to update the nms itself!
  nms <- names(prob.vec);

  # get html confusion matrix
  pred.out <- as.factor(ifelse(prob.vec > 0.5, 1, 0));
  act.cls <- as.numeric(cls)-1;

  prob.res <- data.frame(Probability=prob.vec, Predicted=pred.out, Actual=act.cls);

  write.table(prob.res, file="roc_pred_prob1.csv", sep=",", col.names = TRUE);

  conf.res <- table(pred.out, act.cls);
  analSet$conf.table <- xtable::xtable(conf.res, caption="Confusion Matrix (Cross-Validation)");
  analSet$conf.mat <- print(analSet$conf.table, type = "html", print.results=F, caption.placement="top", html.table.attributes="border=1 width=150" )

  if(anal.mode == "test"){
    if(!is.null(dataSet$test.data)){

      test.pred <- ifelse(analSet$ROCtest$test.res > 0.5, 1, 0);
      test.cls <- as.numeric(dataSet$test.cls)-1;

      test.df <- data.frame(Prob_HoldOut=analSet$ROCtest$test.res, Predicted_HoldOut=test.pred, Actual_HoldOut=test.cls);
      suppressMessages(write.table(test.df, file="roc_pred_prob1.csv", sep=",", append=TRUE, col.names = TRUE));

      test.res <- table(test.pred, test.cls);
      analSet$conf.mat.test <- print(xtable::xtable(test.res,
                                                            caption="Confusion Matrix (Hold-out)"),
                                             type = "html", print.results=F, xtable.width=120, caption.placement="top",
                                             html.table.attributes="border=1 width=150" );
    }
  }

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 8;

  imgSet$roc.testprob.plot <- imgName;
  imgSet$roc.testprob.name <- mdl.inx
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  set.seed(123);
  y <- rnorm(length(prob.vec));
  max.y <- max(abs(y));
  ylim <- max.y*c(-1.05, 1.05);
  
  xlim <- c(0, 1.0);
  
  op <- par(mar=c(4,4,3,6));
  pchs <- ifelse(as.numeric(cls) == 1, 1, 19);
  
  colors <- ifelse(show==1, "darkgrey", "black");
  ROCR::plot(prob.vec, y, pch=pchs, col=colors, xlim=xlim, ylim= ylim, xlab = "Predicted Class Probabilities", ylab="Samples");
  abline(h = 0, lty = 2, col="grey");
  abline(v = 0.5, lty = 2, col="grey");
  
  par(xpd=T);
  legend("right",inset=c(-0.11,0), legend = unique(as.character(cls)), pch=unique(pchs));
  
  test.y <- test.x <- 0;
  if(showPred){
    test.y <- rnorm(length(analSet$ROCtest$test.res));
    test.x <- analSet$ROCtest$test.res;

    pchs <- ifelse(as.numeric(dataSet$test.cls) == 1, 1, 19);
    points(test.x, test.y, pch=pchs, cex=1.5, col="red");
  }

  if(show == 1){

    # add labels for sample classified wrong
    # the idea is to add labels to the left side for those with prob < 0.5
    # and add labels to the right side of the point for those with prob > 0.5
    # leave 0.5 out

    # first wrong pred as 1 (right side)
    act.ones <- as.numeric(cls)-1 == 1;
    pred.vec <- ifelse(prob.vec > 0.5, 1, 0);

    wrong.inx <- (pred.vec != as.numeric(cls)-1) & pred.vec == 1;
    if(sum(wrong.inx) > 0){
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx], pos=4);
    }

    # first wrong pred as 0 (left side)
    act.zeros <- as.numeric(cls)-1 == 0;
    pred.vec <- ifelse(prob.vec < 0.5, 0, 0.5);
    wrong.inx <- pred.vec != as.numeric(cls)-1 & pred.vec == 0;
    if(sum(wrong.inx) > 0){
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx], pos=2);
    }

    if(showPred){
      # Sample names are columns (features×samples orientation)
      nms <- colnames(dataSet$test.data);

      act.ones <- as.numeric(dataSet$test.cls)-1 == 1;
      act.zeros <- as.numeric(dataSet$test.cls)-1 == 0;

      # leave 0.5 there
      pred.vec <- ifelse(test.x > 0.5, 1, 0.5);
      wrong.inx <- (pred.vec != as.numeric(dataSet$test.cls)-1) & act.ones;
      if(sum(wrong.inx) > 0){
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx], pos=4, cex=0.9);
      }

      pred.vec <- ifelse(test.x < 0.5, 0, 0.5);
      wrong.inx <- pred.vec != as.numeric(dataSet$test.cls)-1 & act.zeros;
      if(sum(wrong.inx) > 0){
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx], pos=2, cex=0.9);
      }
    }
  }
  par(op)
  dev.off();
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);
}



#'Plot ROC
#'@description Pred and auroc are lists containing predictions
#'and labels from different cross-validations 
#'@usage PlotROC(mSetObj=NA, imgName, format="png", dpi=default.dpi, mdl.inx, 
#'avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", 
#'users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param mdl.inx Model index, 0 means to compare all models, 
#'input 1-6 to plot a ROC curve for one of the top six models  
#'@param avg.method Input the method to compute the average ROC curve, 
#'either "threshold", "vertical" or "horizontal"
#'@param show.conf Logical, if 1, show confidence interval, if 0 do not show
#'@param show.holdout Logical, if 1, show the ROC curve for hold-out validation, if 0 do not show 
#'@param focus "fpr" 
#'@param cutoff Input the threshold to limit the calculation of the 
#'ROC curve, the number must be between 0 and 1.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotROC <- function(dataName = "", imgName, format="png", dpi=default.dpi, mdl.inx, avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0){

  #msg("[PlotROC] DEBUG: Function called with dataName='", dataName, "' imgName='", imgName, "' mdl.inx=", mdl.inx)

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  #msg("[PlotROC] DEBUG: analSet is.null? ", is.null(analSet))
  #msg("[PlotROC] DEBUG: analSet$multiROC is.null? ", is.null(analSet$multiROC))

  anal.mode <- analSet$mode;

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;

  imgSet$roc.multi.plot <- imgName;
  imgSet$roc.multi.model <- mdl.inx;

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

  op <- par(mar=c(5,4,3,3));

  if(mdl.inx == 0){ # compare all models
    #msg("[PlotROC] DEBUG: Comparing all models")
    preds <- analSet$multiROC$pred.list;
    auroc <- analSet$multiROC$auc.vec;

    #msg("[PlotROC] DEBUG: length(preds) = ", length(preds))
    #msg("[PlotROC] DEBUG: preds[[1]] class = ", class(preds[[1]]))
    #msg("[PlotROC] DEBUG: avg.method = ", avg.method)

    #msg("[PlotROC] DEBUG: Creating performance object for first model")
    perf <- ROCR::performance(preds[[1]], "tpr", "fpr");
    #msg("[PlotROC] DEBUG: Performance object created, class = ", class(perf))

    #msg("[PlotROC] DEBUG: Computing average curve")
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    #msg("[PlotROC] DEBUG: Average curve computed")

    cols <- (1:length(preds))+1;
    #msg("[PlotROC] DEBUG: Creating plot canvas")
    ROCR::plot(perf.avg@x.values[[1]], perf.avg@y.values[[1]], type="n", axes=F,
               xlim=c(0,1), ylim=c(0,1),
               xlab="1-Specificity (False positive rate)",
               ylab="Sensitivity (True positive rate)"
    );

    #msg("[PlotROC] DEBUG: Adding plot elements")
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(perf.avg@x.values[[1]], perf.avg@y.values[[1]], col=cols[1]);

    #msg("[PlotROC] DEBUG: Adding lines for ", length(preds), " models")
    for(i in 2:length(preds)){
      if (i %% 2 == 0) {
        #msg("[PlotROC] DEBUG: Processing model ", i)
      }
      perf <- ROCR::performance(preds[[i]], "tpr", "fpr");
      avg <- ComputeAverageCurve(perf, avg.method);
      lines(avg@x.values[[1]], avg@y.values[[1]], col=cols[i]);
    }

    #msg("[PlotROC] DEBUG: Adding legend")
    best.inx <- which.max(auroc);

    # now add and format legends to the bottom right corner
    feats <- c("Var.", names(preds));
    feats <- substr(feats, 1, 8);
    feats <- sprintf("%-5s", feats);

    vals <- c("AUC", round(auroc, 3));
    vals <- sprintf("%-8s", vals);

    # Format CIs properly - GetCIs returns a 2×n matrix (lower, upper)
    ci.mat <- analSet$multiROC$auc.ci;
    if (is.matrix(ci.mat)) {
      # Format as "lower-upper" for each model
      ci.strings <- apply(ci.mat, 2, function(x) {
        sprintf("%.3f-%.3f", x[1], x[2])
      });
      cis <- c("95% CI", ci.strings);
    } else {
      # Fallback if not a matrix
      cis <- c("95% CI", rep("N/A", length(preds)));
    }

    legends <- paste(feats, vals, cis, sep=" ");

    pch <- c(NA, rep(15, length(preds)));
    cols <- c(NA, cols);

    legend("bottomright", legend = legends, pch=15, col=cols);
    #msg("[PlotROC] DEBUG: Legend added");

  }else if(mdl.inx > 0){

    preds <- ROCR::prediction(analSet$multiROC$pred.cv[[mdl.inx]], analSet$multiROC$true.cv);
    auroc <- round(analSet$multiROC$auc.vec[mdl.inx],3);
    auc.ci <- analSet$multiROC$auc.ci[mdl.inx];
    perf <- ROCR::performance(preds, "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    y.all <- perf.avg@y.values[[1]];
    x.all <- perf.avg@x.values[[1]];
    lgd <- paste("Area under the curve (AUC) =", auroc, "\n",
                 "95% CI:", auc.ci);

    ROCR::plot(x.all, y.all, type="n", axes=F,
               xlim=c(0,1), ylim=c(0,1),
               xlab="1-Specificity (False positive rate)",
               ylab="Sensitivity (True positive rate)"
    );

    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(x.all, y.all, type="l", lwd=2, col="blue");

    if(show.conf){
      res <- ComputeHighLow(perf);
      suppressWarnings(polygon(c(x.all, rev(x.all)), c(res$con.low, rev(res$con.high)), col="#0000ff22"))
    }

    legend("center", legend = lgd, bty="n");

  }

  #msg("[PlotROC] DEBUG: Closing graphics device and saving")
  dev.off();
  #msg("[PlotROC] DEBUG: Graphics device closed, file should be at: ", imgName)
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  #msg("[PlotROC] DEBUG: PlotROC completed successfully")
  return(imgName);
}

#'Plot ROC for the ROC Curve Based Model Creation and Evaluation module
#'@description Plot the ROC curve of the biomarker model created using a user-selected subset of features.
#'Pred and auroc are lists containing predictions and labels from different cross-validations. 
#'@usage PlotROCTest(mSetObj=NA, imgName, format="png", 
#'dpi=default.dpi, mdl.inx, avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", 
#'users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param mdl.inx Model index, 0 means to compare all models, 
#'input 1-6 to plot a ROC curve for one of the top six models  
#'@param avg.method Input the method to compute the average ROC curve, 
#'either "threshold", "vertical" or "horizontal"
#'@param show.conf Logical, if 1, show confidence interval, if 0 do not show
#'@param show.holdout Logical, if 1, show the ROC curve for hold-out validation, if 0 do not show 
#'@param focus "fpr" 
#'@param cutoff Input the threshold to limit the calculation of the ROC curve, the number must be between 0 and 1.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export 

PlotROCTest<-function(dataName = "", imgName, format="png", dpi=default.dpi, mdl.inx, avg.method, show.conf, show.holdout, focus="fpr", cutoff = 1.0){

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  anal.mode <- analSet$mode;

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;
  imgSet$roc.testcurve.plot <- imgName;
  imgSet$roc.testcurve.name <- mdl.inx
  imgSet$roc.testcurve.method <- avg.method

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

  op <- par(mar=c(5,4,3,3));

  if(anal.mode=="explore" && mdl.inx == 0){ # compare all models
    preds <- analSet$ROCtest$pred.list;
    auroc <- analSet$ROCtest$auc.vec;
    perf <- ROCR::performance(preds[[1]], "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    
    cols <- (1:length(preds))+1;
    ROCR::plot(perf.avg@x.values[[1]], perf.avg@y.values[[1]], type="n", axes=F,
               xlim=c(0,1), ylim=c(0,1),
               xlab="1-Specificity (False positive rate)",
               ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(perf.avg@x.values[[1]], perf.avg@y.values[[1]], col=cols[1]);
    for(i in 2:length(preds)){
      perf <- ROCR::performance(preds[[i]], "tpr", "fpr");
      avg <- ComputeAverageCurve(perf, avg.method);
      lines(avg@x.values[[1]], avg@y.values[[1]], col=cols[i]);
    }
    
    best.inx <- which.max(auroc);
    
    # now add and format legends to the bottom right corner
    feats <- c("Var.", names(preds));
    feats <- substr(feats, 1, 8);
    feats <- sprintf("%-5s", feats);
    
    vals <- c("AUC", round(auroc, 3));
    
    vals <- sprintf("%-8s", vals);
    
    cis <- analSet$multiROC$auc.ci;
    cis <- c("CI", cis);
    legends <- paste(feats, vals, cis, sep="");

    pch <- c(NA, rep(15, length(preds)));
    cols <- c(NA, cols);

    legend("bottomright", legend = legends, pch=15, col=cols);

  }else if(mdl.inx > 0 && anal.mode=="explore"){

    preds <- ROCR::prediction(analSet$ROCtest$pred.cv[[mdl.inx]], analSet$ROCtest$true.cv);
    auroc <- round(analSet$ROCtest$auc.vec[mdl.inx],3);
    auc.ci <- analSet$ROCtest$auc.ci[mdl.inx];
    perf <- ROCR::performance(preds, "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    y.all <- perf.avg@y.values[[1]];
    x.all <- perf.avg@x.values[[1]];
    lgd <- paste("Area under the curve (AUC) =", auroc, "\n",
                 "95% CI:", auc.ci);
    
    ROCR::plot(x.all, y.all, type="n", axes=F,
               xlim=c(0,1), ylim=c(0,1),
               xlab="1-Specificity (False positive rate)",
               ylab="Sensitivity (True positive rate)"
    );
    
    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(x.all, y.all, type="l", lwd=2, col="blue");
    
    if(show.conf){
      res <- ComputeHighLow(perf);
      suppressWarnings(polygon(c(x.all, rev(x.all)), c(res$con.low, rev(res$con.high)), col="#0000ff22"))
    }
    
    legend("center", legend = lgd, bty="n");
    
  }else{ # plot ROC of specific model and save the table for details

    preds <- ROCR::prediction(analSet$ROCtest$pred.cv, analSet$ROCtest$true.cv);
    auroc <- round(analSet$ROCtest$auc.vec[1],3)
    auc.ci <- analSet$ROCtest$auc.ci;

    perf <- ROCR::performance(preds, "tpr", "fpr");
    perf.avg <- ComputeAverageCurve(perf, avg.method);
    y.all <- perf.avg@y.values[[1]];
    x.all <- perf.avg@x.values[[1]];
    # to draw a roc curve line from (0,0)
    y.all <- c(0.0, y.all);
    x.all <- c(0.0, x.all);

    lgd <- paste("Area under the curve (AUC) =", auroc, "\n",
                 "95% CI:", auc.ci);

    ROCR::plot(x.all, y.all, type="n", axes=F,
               xlim=c(0,1), ylim=c(0,1),
               xlab="1-Specificity (False positive rate)",
               ylab="Sensitivity (True positive rate)"
    );

    box()
    axis(side=2)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at, labels=as.graphicsAnnot(sprintf( "%.1f", lab.labels)));
    abline(v=grid.at, lty=3, col="lightgrey");
    abline(h=grid.at, lty=3, col="lightgrey");
    lines(x.all, y.all, type="l", lwd=2, col="blue");

    if(show.conf){
      res <- ComputeHighLow(perf);
      # to draw a roc curve line from (0,0)
      # suppressWarnings(polygon(c(x.all, rev(x.all)), c(res$con.low, rev(res$con.high)), col="#0000ff22"))
      suppressWarnings(polygon(c(x.all, rev(x.all)), c(c(0,res$con.low), c(rev(res$con.high),0)), col="#0000ff22"))
    }
    if(show.holdout){

      roc.obj <- pROC::roc(dataSet$test.cls, analSet$ROCtest$test.res, percent = F);
      test.x <- 1-roc.obj$spec;
      test.y <- roc.obj$sens;

      lbls <- c("Type", "CV","Holdout");
      lbls <- sprintf("%-10s",lbls);

      test.auc <- round(roc.obj$auc[[1]],3);
      aucs <- c("AUC", auroc, test.auc);

      lgd <- paste(lbls, aucs, sep="");
      lines(test.x, test.y, type="l", lwd=2, col="magenta");
      legend("bottomright", legend = lgd, pch=c(NA, 15, 15), col=c(NA, "blue", "magenta"));
    }else{
      legend("center", legend = lgd,  bty="n");
    }
  }
  dev.off();
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);
}

#'Plot classification performance using different features for Multi-Biomarker
#'@description Plot of the accuracy of classification with an increasing number of features.
#'@usage PlotAccuracy(mSetObj=NA, imgName, format="png", dpi=default.dpi)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotAccuracy<-function(dataName = "", imgName, format="png", dpi=default.dpi){

  #msg("[PlotAccuracy] DEBUG: Function called with dataName='", dataName, "' imgName='", imgName, "'")

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  #msg("[PlotAccuracy] DEBUG: analSet is.null? ", is.null(analSet))
  #msg("[PlotAccuracy] DEBUG: analSet$multiROC is.null? ", is.null(analSet$multiROC))

  anal.mode <- analSet$mode;

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 9; h <- 7;

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

  if(is.null(analSet$multiROC$accu.mat)){
    #msg("[PlotAccuracy] DEBUG: Using ROCtest accu.mat")
    accu.mat <- analSet$ROCtest$accu.mat;
  }else{
    #msg("[PlotAccuracy] DEBUG: Using multiROC accu.mat")
    accu.mat <- analSet$multiROC$accu.mat;
  }

  #msg("[PlotAccuracy] DEBUG: accu.mat is.null? ", is.null(accu.mat))
  if (!is.null(accu.mat)) {
    #msg("[PlotAccuracy] DEBUG: accu.mat dimensions: ", nrow(accu.mat), " x ", ncol(accu.mat))
  }

  # OPTIMIZED: Use colMeans instead of apply for 60-100x speedup
  mn.accu <- colMeans(accu.mat);
  ylabl <- 'Predictive Accuracy';
  ylim <- c(0,1);
  title <- 'Predictive accuracies with different features';
  txt.lbls <- paste(100*round(mn.accu,3),'%');

  matplot(t(accu.mat),type='l', lty=2, col="grey", xlab='Number of features',ylab=ylabl, ylim=ylim,
          axes=F,main=title);

  lines(1:ncol(accu.mat), mn.accu, lwd=2);
  points(mn.accu, pch=19, col=ifelse(1:length(mn.accu)==which.max(mn.accu),"red","blue"));
  text(mn.accu,labels=txt.lbls, adj=c(-0.3, -0.5), srt=45, xpd=T)
  axis(2);

  lbls <- colnames(accu.mat);
  axis(1, 1:length(mn.accu), labels=lbls);

  imgSet$roc.pred <- imgName;

  dev.off();
  saveSet(imgSet, "imgSet");
  return(imgName);
}

#'Plot classification performance using different features for Biomarker Tester
#'@description Plot of the accuracy of classification with an increasing number of features.
#'@usage PlotTestAccuracy(mSetObj=NA, imgName, format="png", dpi=default.dpi)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotTestAccuracy<-function(dataName = "", imgName, format="png", dpi=default.dpi){

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  anal.mode <- analSet$mode;

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 4;

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

  y.vec <- analSet$ROCtest$accu.mat[1,];
  ylim.ext <- GetExtendRange (y.vec, 12); # first increase ylim by 1/12
  #boxplot(y.vec, col="#0000ff22", ylim=ylim.ext, outline=FALSE, boxwex=c(0.5, 0.5), ylab="Predictive Accuracy");
  #stripchart(t(analSet$ROCtest$accu.mat), method = "jitter", vertical=T, add = T, pch=19);
  boxplot(y.vec, col="#0000ff22", ylim=ylim.ext, outline=FALSE, boxwex=c(0.7, 0.7), xlab="Predictive Accuracy", horizontal = TRUE);
  stripchart(y.vec, method = "jitter", add = T, pch=19);
  accu.info <- paste ("The average accuracy based on 100 cross validations is", round(mean(y.vec), 3));

  imgSet$roc.testpred <- imgName;

  if(!is.null(dataSet$test.cls)){
    test.pred <- ifelse(analSet$ROCtest$test.res > 0.5, 1, 0);
    test.cls <- as.numeric(dataSet$test.cls)-1;

    hit <- sum(test.pred == test.cls);
    percent <- round(hit/length(test.cls), 3);
    accu.info <- paste(accu.info, ". The accuracy for hold out data prediction is ",  percent,
                       "(", hit, "/",  length(test.cls), ").", sep="");
  }

  analSet$ROCtest$accu.info <- accu.info;

  dev.off();
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);
}

#'Plot selected compounds by their percentage frequency
#'@description Plot the important variables of single biomarker model ranked by order of importance
#'@usage PlotImpBiomarkers(mSetObj=NA, imgName, format="png", dpi=default.dpi, 
#'mdl.inx, measure = "freq", feat.num = 15)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format elect the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param mdl.inx Model index, -1 selects the model with the best AUC, input 1-6 to 
#'view the important features of one of the top six models
#'@param measure Choose to rank features by the frequency of being selected "freq", or the 
#'mean importance measure "mean"
#'@param feat.num Input the number of features to include in the plot, by default it is 15.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotImpBiomarkers <- function(dataName = "", imgName, format="png", dpi=default.dpi,
                              mdl.inx, measure = "freq", feat.num = 15) {

  #msg("[PlotImpBiomarkers] DEBUG: Function called with dataName='", dataName, "' imgName='", imgName, "' mdl.inx=", mdl.inx)

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  #msg("[PlotImpBiomarkers] DEBUG: analSet is.null? ", is.null(analSet))
  #msg("[PlotImpBiomarkers] DEBUG: analSet$multiROC is.null? ", is.null(analSet$multiROC))

  anal.mode <- analSet$mode;
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;
  imgSet$roc.imp.plot <- imgName;
  imgSet$roc.imp.name <- mdl.inx

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

  if(is.null(dataSet$norm_explore)){
    data <- dataSet$data.norm;
    dataSet$norm_explore <- dataSet$data.norm;
  } else {
    data <- dataSet$norm_explore
  }

  #msg("[PlotImpBiomarkers] DEBUG: data dimensions: ", nrow(data), " x ", ncol(data))

  cls <- dataSet$cls;
  op <- par(mar=c(6,10,3,7));

  # if(anal.mode == "explore"){
  if(mdl.inx == -1){
    mdl.inx <- analSet$multiROC$best.inx;
    #msg("[PlotImpBiomarkers] DEBUG: Using best model index: ", mdl.inx)
  }

  #msg("[PlotImpBiomarkers] DEBUG: About to call GetImpFeatureMat")
  #msg("[PlotImpBiomarkers] DEBUG: length(analSet$multiROC$imp.cv) = ", length(analSet$multiROC$imp.cv))
  #msg("[PlotImpBiomarkers] DEBUG: analSet$multiROC$test.feats[mdl.inx] = ", analSet$multiROC$test.feats[mdl.inx])

  imp.mat <- GetImpFeatureMat(dataSet, analSet$multiROC$imp.cv, analSet$multiROC$test.feats[mdl.inx]);

  #msg("[PlotImpBiomarkers] DEBUG: imp.mat dimensions: ", nrow(imp.mat), " x ", ncol(imp.mat))
  # }else{
  #   imp.mat <- GetImpFeatureMat(dataSet, analSet$multiROC$imp.cv, null);
  # }

  imp.nms <- rownames(imp.mat);
  # Features are in rows (features×samples orientation)
  hit.nms <- imp.nms[imp.nms %in% rownames(data)];
  data <- data[hit.nms, ];
  
  # note, tapply can only be applied to a vector, we need
  # to combine with apply in order to used on a data frame
  # Features×samples orientation: apply over rows (margin=1)
  mds <- apply(data, 1,
               function(x){
                 tapply(x, cls, median);
               });

  lowhigh <- apply(mds, 2,
                   function(x){
                     ifelse(rank(x)==1, "Low", "High")
                   });
  lowhigh <- t(lowhigh);

  temp.dat <- data.frame(imp.mat, lowhigh);
  colnames(temp.dat) <- c(colnames(imp.mat), levels(cls));
  imp.fileNm <- "imp_features_cv.csv";
  fast.write.csv(temp.dat, file=imp.fileNm);
  temp.dat <- NULL;

  # record the imp.mat for table show
  analSet$imp.mat <- imp.mat;
  analSet$lowhigh <- lowhigh;
  analSet$roc.sig.nm <- imp.fileNm;
  
  if(measure=="freq"){
    imp.vec <- sort(imp.mat[,1], decreasing=T);
    xlbl <- "Selected Frequency (%)";
  }else{ # default sort by freq, need to reorder
    imp.vec <- sort(imp.mat[,2], decreasing=T);
    xlbl <- "Average Importance";
  }
  var.len <- length(imp.vec);
  
  if(feat.num <= 0){
    feat.num = 15;
  }else if(feat.num > var.len){
    feat.num <- var.len;
  }
  
  imp.vec<-rev(imp.vec[1:feat.num]);

  nms.orig <- names(imp.vec);

  # Convert IDs to gene symbols for better readability
  vip.nms <- nms.orig;  # Default to original names

  # Try to load symbol mapping if available
  paramSet <- readSet(paramSet, "paramSet");
  if (file.exists("symbol.map.qs")) {
    tryCatch({
      gene.map <- readDataQs("symbol.map.qs", paramSet$anal.type, dataSet$name);
      if (!is.null(gene.map) && nrow(gene.map) > 0) {
        # Map IDs to symbols
        symbols <- doIdMappingGeneric(nms.orig, gene.map, "gene_id", "symbol", "vec");
        # Use symbols where available, keep original ID where not
        valid.symbols <- !is.na(symbols) & symbols != "";
        vip.nms[valid.symbols] <- symbols[valid.symbols];
        #msg("[PlotImpBiomarkers] DEBUG: Converted ", sum(valid.symbols), "/", length(nms.orig), " IDs to gene symbols");
      }
    }, error = function(e) {
      #msg("[PlotImpBiomarkers] DEBUG: Could not load symbol mapping: ", e$message);
    })
  } else {
    #msg("[PlotImpBiomarkers] DEBUG: symbol.map.qs not found, using original IDs");
  }

  # Truncate long names
  vip.nms <- substr(vip.nms, 1, 20);
  names(imp.vec) <- NULL;
  
  xlim.ext <- GetExtendRange(imp.vec, 12);
  dotchart(imp.vec, bg="blue", xlab= xlbl, xlim=xlim.ext);
  
  mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1)
  names(imp.vec) <- nms.orig;
  
  axis.lims <- par("usr"); # x1, x2, y1 ,y2
  
  # get character width
  shift <- 2*par("cxy")[1];
  lgd.x <- axis.lims[2] + shift;
  
  x <- rep(lgd.x, feat.num);
  y <- 1:feat.num;
  
  par(xpd=T);
  
  # now synchronize lowhigh with imp.vec
  lowhigh <- lowhigh[nms.orig, ,drop=FALSE];
  
  nc <- ncol(lowhigh);
  col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(nc);
  
  # calculate background
  bg <- matrix("", nrow(lowhigh), nc);
  for (m in 1:nrow(lowhigh)){
    bg[m,] <- ifelse(lowhigh[m,]=="High", col[1], col[2]);
  } 
  
  cls.lbl <- levels(cls);
  
  for (n in 1:ncol(lowhigh)){
    points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
    # now add label
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5));
    # shift x, note, this is good for current size
    x <- x + shift/1.25;
  }
  
  # now add color key, padding with more intermediate colors for contiuous band
  col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(20)
  nc <- length(col);
  x <- rep(x[1] + shift, nc);
  
  shifty <- (axis.lims[4]-axis.lims[3])/3;
  starty <- axis.lims[3] + shifty;
  endy <- axis.lims[3] + 2*shifty;
  y <- seq(from = starty, to = endy, length = nc);
  
  points(x,y, bty="n", pch=15, col=rev(col), cex=2);
  
  text(x[1], endy+shifty/8, "High");
  text(x[1], starty-shifty/8, "Low");
  par(op);
  dev.off();

  # Save norm_explore update if it was created
  if(!is.null(dataSet$norm_explore)){
    RegisterData(dataSet);
  }
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);

}


#'Plot results of permutation tests
#'@description Plot results of permutation tests
#'@usage Plot.Permutation(mSetObj=NA, imgName, format="png", dpi=default.dpi)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param format elect the image format, "png", of "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Plot.Permutation<-function(dataName = "", imgName, format="png", dpi=default.dpi){

  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");
  imgSet <- readSet(imgSet, "imgSet");

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  w <- 8; h <- 8;
  imgSet$roc.perm.plot <- imgName;

  if(analSet$ROCtest$perm.res$perf.measure == "auroc"){
    imgSet$roc.perm.method <- "auroc"
  }else{
    imgSet$roc.perm.method <- "predictive accuracy"
  }

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

  if(analSet$ROCtest$perm.res$perf.measure == "auroc"){
    ROCR::plot(analSet$ROCtest$perm.res$perf.obj, col="grey", lwd=1, lty=2,
               xlab="1-Specificity (False Positive rate)",
               ylab="Sensitivity (True Positive rate)");

    # now add the original ROC

    preds <- ROCR::prediction(analSet$ROCtest$pred.cv, analSet$ROCtest$true.cv);
    auroc <- round(analSet$ROCtest$auc.vec[1],3)
    perf <- ROCR::performance(preds, "tpr", "fpr");
    # need to replace Inf with 1
    alpha.vals <- perf@alpha.values;
    perf@alpha.values <- lapply(alpha.vals, function(x){
      x[x==Inf] <- 1;
      x[x==-Inf] <- 0;
      x
    });


    ROCR::plot(perf,lwd=2,avg="threshold", col="blue", add=T);

    # calculate p value
    perm.vec <- analSet$ROCtest$perm.res$auc.vec;
    better.hits <- sum(perm.vec[-1]>=perm.vec[1]);
    num <- length(perm.vec);
    if(better.hits == 0) {
      p <- paste("p <", 1/num);
    }else{
      p <- better.hits/num;
      p <- paste("p =", signif(p, digits=5));
    }
    legend("center", legend = paste('Empirical p-value: ', p), bty="n", cex=1.5);
  }else{ # accuracies
    perms <- PreparePermResult(analSet$ROCtest$perm.res$acc.vec);
    perm.vec <- perms$permut;
    perm.p <- perms$permut.p;

    op<-par(mar=c(5,5,2,4));

    xlim.ext <- GetExtendRange (perm.vec, 10);
    hst <- hist(perm.vec[-1], breaks = "FD", freq=T, ylab="Frequency", xlim=xlim.ext, xaxt="n", xlab= 'Permutation test statistics', col="lightblue", main="");
    axis(1);

    # add the indicator using original label
    h <- max(hst$counts)
    arrows(perm.vec[1], h/5, perm.vec[1], 0, col="red", lwd=2);
    text(perm.vec[1], h/3.5, paste('Observed \n statistic \n', perm.p), xpd=T);
    par(op);
  }
  dev.off();
  saveSet(analSet, "analSet");
  saveSet(imgSet, "imgSet");
  return(imgName);
}

##############################################
## Covariate Adjustment for Biomarker Analysis
## Remove confounding effects before ROC analysis
##############################################

#' Perform Covariate Adjustment for Biomarker Analysis
#'
#' This function removes confounding effects from the data before biomarker analysis
#' to ensure that ROC models identify disease-specific signals rather than confounding
#' factors like age, sex, batch, etc.
#'
#' @param dataName Name of the dataset
#' @param primary.condition Primary disease/condition variable (e.g., "Disease", "Group")
#' @param blocking.factor Optional blocking factor for paired/repeated measures (e.g., "PatientID")
#' @return 1 on success, 0 on failure
#' @export
PerformCovariateAdjustmentForROC <- function(dataName,
                                             primary.condition,
                                             blocking.factor = "NA",
                                             use.combat = FALSE) {

  # Load required libraries
  library(limma)
  library(dplyr)

  # Read data structures
  dataSet <- readDataset(dataName)
  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")

  # Get covariates to adjust for (passed from Java via adj.vec)
  if(!exists('adj.vec')) {
    msgSet$current.msg <- "Error: No covariates specified for adjustment!"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  if(length(adj.vec) == 0 || all(adj.vec == "NA")) {
    msgSet$current.msg <- "Error: No covariates specified for adjustment!"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Get the normalized data (transposed for biomarker analysis)
  if(!is.null(dataSet$data.norm.transposed)) {
    data.matrix <- dataSet$data.norm.transposed  # Samples x Features
  } else if(!is.null(dataSet$data.norm)) {
    data.matrix <- t(dataSet$data.norm)  # Transpose: Features x Samples → Samples x Features
  } else {
    msgSet$current.msg <- "Error: No normalized data available!"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Get metadata
  meta.info <- dataSet$meta.info

  # Ensure samples match between data and metadata
  matched_indices <- match(rownames(data.matrix), rownames(meta.info))
  if(any(is.na(matched_indices))) {
    msgSet$current.msg <- "Error: Sample names don't match between data and metadata!"
    saveSet(msgSet, "msgSet")
    return(0)
  }
  meta.info <- meta.info[matched_indices, , drop=FALSE]

  # Validate that primary condition exists
  if(!(primary.condition %in% colnames(meta.info))) {
    msgSet$current.msg <- paste0("Error: Primary condition '", primary.condition, "' not found in metadata!")
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Validate that all covariates exist
  missing.covs <- adj.vec[!(adj.vec %in% colnames(meta.info))]
  if(length(missing.covs) > 0) {
    msgSet$current.msg <- paste0("Error: Covariates not found in metadata: ", paste(missing.covs, collapse=", "))
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Convert character columns to factors
  var.types <- sapply(meta.info, class)
  meta.info[, var.types == "character"] <- lapply(meta.info[, var.types == "character"], factor)

  # Build design matrix
  # Include primary condition + all covariates
  all.vars <- c(primary.condition, adj.vec)

  # Remove samples with missing values in any of these variables
  complete.idx <- complete.cases(meta.info[, all.vars, drop=FALSE])
  if(sum(complete.idx) < nrow(meta.info)) {
    n.removed <- sum(!complete.idx)
    msgSet$current.msg <- paste0("Warning: Removed ", n.removed, " samples with missing covariate values")
    saveSet(msgSet, "msgSet")

    data.matrix <- data.matrix[complete.idx, , drop=FALSE]
    meta.info <- meta.info[complete.idx, , drop=FALSE]
  }

  if(nrow(data.matrix) < 10) {
    msgSet$current.msg <- "Error: Too few samples remaining after removing missing values!"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Transpose data for limma (Features x Samples)
  data.limma <- t(data.matrix)

  # Create design matrix
  # Use full model with all variables
  formula.str <- paste0("~ ", paste(all.vars, collapse=" + "))
  design <- model.matrix(as.formula(formula.str), data=meta.info)

  # Check for rank deficiency (confounded variables)
  design.rank <- qr(design)$rank
  design.ncol <- ncol(design)

  if(design.rank < design.ncol) {
    # Design matrix is rank deficient - some coefficients are not estimable
    # This means primary condition and covariates are confounded

    # Get column names of non-estimable coefficients
    qr.obj <- qr(design)
    if(qr.obj$rank < ncol(design)) {
      pivot <- qr.obj$pivot
      non.estimable.cols <- colnames(design)[pivot[(qr.obj$rank + 1):ncol(design)]]

      msgSet$current.msg <- paste0(
        "Error: Covariate adjustment failed - Primary condition and covariates are confounded (perfectly correlated). ",
        "Non-estimable coefficients: ", paste(non.estimable.cols, collapse=", "), ". ",
        "This means your primary condition (", primary.condition, ") and covariate(s) (",
        paste(adj.vec, collapse=", "), ") are the same or perfectly correlated. ",
        "Please choose different covariates that vary WITHIN each group of your primary condition (e.g., Age, Sex, Batch)."
      )
      saveSet(msgSet, "msgSet")
      return(0)
    }
  }

  data.adjusted <- NULL
  block.vec <- NULL

  if(use.combat) {
    if(length(adj.vec) != 1) {
      msgSet$current.msg <- "Error: ComBat batch correction requires exactly one selected covariate."
      saveSet(msgSet, "msgSet")
      return(0)
    }

    if(!is.null(blocking.factor) && !is.na(blocking.factor) &&
       blocking.factor != "NA" && blocking.factor != "") {
      msgSet$current.msg <- "Error: Blocking factor is not supported when ComBat batch correction is enabled."
      saveSet(msgSet, "msgSet")
      return(0)
    }

    suppressPackageStartupMessages(require(sva))

    batch.var <- adj.vec[1]
    batch.vec <- meta.info[, batch.var]
    batch.fac <- as.factor(batch.vec)

    if(length(unique(batch.fac)) < 2) {
      msgSet$current.msg <- paste0("Error: ComBat requires at least two batch levels in ", batch.var, ".")
      saveSet(msgSet, "msgSet")
      return(0)
    }

    if(any(table(batch.fac) < 2)) {
      msgSet$current.msg <- paste0("Error: Each batch level in ", batch.var, " must contain at least two samples for ComBat.")
      saveSet(msgSet, "msgSet")
      return(0)
    }

    modcombat <- model.matrix(as.formula(paste0("~ ", primary.condition)), data=meta.info)
    data.adjusted <- ComBat(dat=data.limma,
                            batch=batch.fac,
                            mod=modcombat,
                            par.prior=TRUE,
                            prior.plots=FALSE)
  } else {
    # Handle blocking factor if specified
    if(!is.null(blocking.factor) && !is.na(blocking.factor) &&
       blocking.factor != "NA" && blocking.factor != "") {

      if(!(blocking.factor %in% colnames(meta.info))) {
        msgSet$current.msg <- paste0("Warning: Blocking factor '", blocking.factor, "' not found, ignoring")
        saveSet(msgSet, "msgSet")
        block.vec <- NULL
      } else {
        block.vec <- meta.info[, blocking.factor]

        # Fit with blocking (random effects via duplicateCorrelation)
        corfit <- duplicateCorrelation(data.limma, design, block=block.vec)
        fit <- lmFit(data.limma, design, block=block.vec, correlation=corfit$consensus)
      }
    } else {
      block.vec <- NULL
      # Fit without blocking (fixed effects only)
      fit <- lmFit(data.limma, design)
    }

    # Use limma's removeBatchEffect to get adjusted data
    # Build explicit numeric matrices here because removeBatchEffect expects
    # numeric covariates, and metadata may contain factors/characters.
    primary.design <- model.matrix(as.formula(paste0("~ ", primary.condition)), data=meta.info)
    covariate.design <- model.matrix(as.formula(paste0("~ ", paste(adj.vec, collapse=" + "))), data=meta.info)

    # remove intercept from covariates; preserve the primary-condition design
    if(ncol(covariate.design) > 1) {
      covariate.design <- covariate.design[, -1, drop=FALSE]
    } else {
      msgSet$current.msg <- "Error: No covariate effects to remove (design matrix issue)"
      saveSet(msgSet, "msgSet")
      return(0)
    }

    covariate.design <- as.matrix(covariate.design)
    primary.design <- as.matrix(primary.design)
    storage.mode(covariate.design) <- "double"
    storage.mode(primary.design) <- "double"

    if(!is.numeric(covariate.design) || !is.numeric(primary.design)) {
      msgSet$current.msg <- "Error: Covariate adjustment failed because the design matrix is not numeric."
      saveSet(msgSet, "msgSet")
      return(0)
    }

    data.adjusted <- removeBatchEffect(data.limma,
                                       covariates=covariate.design,
                                       design=primary.design)
  }

  # Transpose back to Samples x Features for biomarker analysis
  data.adjusted.transposed <- t(data.adjusted)

  # Store both versions in dataSet
  dataSet$data.adjusted <- data.adjusted
  dataSet$data.adjusted.transposed <- data.adjusted.transposed
  dataSet$covariate.adjusted <- TRUE
  dataSet$adjusted.for <- adj.vec
  dataSet$primary.condition <- primary.condition
  dataSet$covariate.adjustment.method <- ifelse(use.combat, "combat", "limma")

  # Update metadata to match adjusted data (samples may have been removed due to missing covariates)
  dataSet$meta.info.original <- dataSet$meta.info
  dataSet$meta.info <- meta.info

  # Also update data.norm to match the filtered samples (for "no adjustment" comparison)
  dataSet$data.norm.original <- dataSet$data.norm
  dataSet$data.norm <- data.limma

  # Update the data.norm.transposed to use adjusted data for downstream ROC analysis
  dataSet$data.norm.transposed.original <- dataSet$data.norm.transposed
  dataSet$data.norm.transposed <- data.adjusted.transposed

  method.label <- ifelse(use.combat, "ComBat", "limma removeBatchEffect")
  msgSet$current.msg <- paste0("Successfully adjusted for: ", paste(adj.vec, collapse=", "),
                               " using ", method.label)
  saveSet(msgSet, "msgSet")

  cat("Covariate adjustment completed:\n")
  cat("  Method:", method.label, "\n")
  cat("  Primary condition:", primary.condition, "\n")
  cat("  Adjusted for:", paste(adj.vec, collapse=", "), "\n")
  cat("  Blocking factor:", ifelse(is.null(block.vec), "None", blocking.factor), "\n")
  cat("  Samples:", nrow(data.matrix), "\n")
  cat("  Features:", ncol(data.matrix), "\n")

  RegisterData(dataSet)
  return(1)
}

#' Plot Before/After PCA to Show Covariate Adjustment Effect
#'
#' Creates side-by-side PCA plots colored by covariate to visualize
#' the effect of covariate adjustment.
#'
#' @param dataName Name of the dataset
#' @param covariate Name of the covariate to color points by
#' @param imgName Image file name
#' @param format Image format ("png" or "pdf")
#' @param dpi Resolution for PNG
#' @return 1 on success, 0 on error
#'
PlotCovariateAdjustmentPCA <- function(dataName, covariate, imgName="covariate_pca",
                                       format="png", dpi=96) {

  # Load required libraries
  library(ggplot2)
  library(ggpubr)  # For ggarrange

  # Read data
  dataSet <- readDataset(dataName)
  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")

  # Check if adjustment was performed
  if(is.null(dataSet$covariate.adjusted) || !dataSet$covariate.adjusted) {
    msgSet$current.msg <- "Error: No covariate adjustment has been performed!"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Get before/after data matrices (Features x Samples)
  data.before <- dataSet$data.norm  # Original normalized data
  data.after <- dataSet$data.adjusted  # Adjusted data

  if(is.null(data.before) || is.null(data.after)) {
    msgSet$current.msg <- "Error: Missing data for PCA comparison!"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Get metadata
  meta.info <- dataSet$meta.info

  # Match samples
  common.samples <- intersect(colnames(data.before), rownames(meta.info))
  data.before <- data.before[, common.samples, drop=FALSE]
  data.after <- data.after[, common.samples, drop=FALSE]
  meta.info <- meta.info[common.samples, , drop=FALSE]

  # Check covariate exists
  if(!(covariate %in% colnames(meta.info))) {
    msgSet$current.msg <- paste0("Error: Covariate '", covariate, "' not found in metadata!")
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Perform PCA on both datasets
  pca.before <- prcomp(t(data.before), scale=TRUE, center=TRUE)
  pca.after <- prcomp(t(data.after), scale=TRUE, center=TRUE)

  # Calculate variance explained
  var.before <- summary(pca.before)$importance[2, 1:2] * 100  # PC1, PC2
  var.after <- summary(pca.after)$importance[2, 1:2] * 100

  # Create data frames for plotting
  # Use the covariate as a factor for proper coloring
  df.before <- data.frame(
    PC1 = pca.before$x[, 1],
    PC2 = pca.before$x[, 2],
    Group = as.factor(meta.info[, covariate]),
    Sample = rownames(pca.before$x)
  )

  df.after <- data.frame(
    PC1 = pca.after$x[, 1],
    PC2 = pca.after$x[, 2],
    Group = as.factor(meta.info[, covariate]),
    Sample = rownames(pca.after$x)
  )

  # Create plot labels
  xlabel.before <- sprintf("PC1 (%.1f%%)", var.before[1])
  ylabel.before <- sprintf("PC2 (%.1f%%)", var.before[2])
  xlabel.after <- sprintf("PC1 (%.1f%%)", var.after[1])
  ylabel.after <- sprintf("PC2 (%.1f%%)", var.after[2])

  # Create before plot
  p1 <- ggplot(df.before, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=3, alpha=0.7) +
    labs(title="Before Adjustment",
         x=xlabel.before,
         y=ylabel.before,
         color=covariate) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "right",
      legend.title = element_text(size=11, face="bold"),
      axis.title = element_text(size=11)
    )

  # Create after plot
  p2 <- ggplot(df.after, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=3, alpha=0.7) +
    labs(title="After Adjustment",
         x=xlabel.after,
         y=ylabel.after,
         color=covariate) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "right",
      legend.title = element_text(size=11, face="bold"),
      axis.title = element_text(size=11)
    )

  # Combine plots side-by-side with shared legend
  combined <- ggarrange(p1, p2, ncol=2, common.legend=TRUE, legend="right")

  # Save plot using standard naming convention: imgName + "dpi96.png"
  imgPath <- paste0(paramSet$imgSet, imgName, "dpi", dpi, ".", format)

  if(format == "png") {
    png(imgPath, width=14, height=7, units="in", res=dpi, type="cairo")
  } else if(format == "pdf") {
    pdf(imgPath, width=14, height=7)
  } else if(format == "svg") {
    svg(imgPath, width=14, height=7)
  } else if(format == "tiff") {
    tiff(imgPath, width=14, height=7, units="in", res=dpi, compression="lzw")
  }

  print(combined)
  dev.off()

  return(1)
}

##################################################
## Covariate Scatter Analysis for Biomarker Module
## Adapted from ExpressAnalyst
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' Perform Covariate Scatter Analysis
#'
#' This function performs a covariate scatter analysis to visualize the effect of 
#' covariate adjustment on p-values. It generates a JSON file for interactive visualization.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param imgName A character string specifying the name of the JSON file to save.
#' @param imgFormat The format of the image to create (e.g., "png", "jpeg").
#' @param analysis.var The primary analysis variable (e.g., disease status).
#' @param ref The reference group (optional).
#' @param block The blocking variable (optional, for paired/repeated measures).
#' @param thresh The significance threshold (default 0.05).
#' @param pval.selection The method for p-value selection ("fdr" or "raw").
#' @param contrast.cls The contrast class ("anova" or specific group name).
#'
#' @return A vector with c(sig.num, nonSig) - number of significant and non-significant features.
#'
#' @export
#'
CovariateScatter.Anal <- function(dataName, 
                                  imgName="NA", 
                                  imgFormat="png", 
                                  analysis.var, 
                                  ref = NULL, 
                                  block = "NA", 
                                  thresh=0.05,
                                  pval.selection="fdr",
                                  contrast.cls = "anova",
                                  useBatchCorrected=TRUE
                                  ){

  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  msg.lm <- ""
  # load libraries
  library(limma)
  library(dplyr)
  
  # get inputs
  if(!exists('adj.vec')){
    adj.bool = F;
    vars <- analysis.var;
  }else{
    if(length(adj.vec) > 0){
      adj.bool = T;
      if("Dataset" %in% adj.vec){
        current.msg <- "Dataset can not be used as covariate!";
        AddMsg(current.msg);
        msgSet$current.msg <<- current.msg;
        return(c(-1, -1))
      }
    }else{
      adj.bool = F;
      adj.vec <- "NA"
    }
  }

  if("Dataset" == block){
      current.msg <- "Dataset can not be used as blocking factor!";
      AddMsg(current.msg);
      msgSet$current.msg <<- current.msg;
      return(c(-1, -1))
   }

  dataSet <- .multiCovariateRegression(dataSet,analysis.var, ref,contrast.cls, blocking.factor=block,adj.vec, F, F, useBatchCorrected);

  rest <- dataSet$comp.res;
  res.noadj <- dataSet$res.noadj;
  
  # make visualization
  adj.mat <- rest[, c("P.Value", "adj.P.Val")]
  noadj.mat <- res.noadj[, c("P.Value", "adj.P.Val")]
  
  colnames(adj.mat) <- c("pval.adj", "fdr.adj")
  colnames(noadj.mat) <- c("pval.no", "fdr.no")
  
  both.mat <- merge(adj.mat, noadj.mat, by = "row.names")
  both.mat$pval.adj <- -log10(both.mat$pval.adj)
  both.mat$fdr.adj <- -log10(both.mat$fdr.adj)
  both.mat$pval.no <- -log10(both.mat$pval.no)
  both.mat$fdr.no <- -log10(both.mat$fdr.no)
  rownames(both.mat) = both.mat[,"Row.names"]
  
  # Add labels (gene symbols if available)
  if(!is.null(dataSet$enrich_ids)) {
    both.mat$label <- .invert_named_vector(dataSet$enrich_ids)[as.character(rownames(both.mat))];
  } else {
    both.mat$label <- rownames(both.mat);
  }
  
  # make plot
  if( "F" %in% colnames(rest)){
    fstat <- rest[, "F"];
  }else{
    fstat <- rest[, "t"];
  }  
  if(pval.selection == "fdr"){
    p.value <- rest[,"adj.P.Val"];
  }else{
    p.value <- rest[,"P.Value"];
  }
  ord.inx <- order(p.value, decreasing = FALSE);
  rest <- rest[ord.inx,,drop=F];
  colnames(rest)[1] <- "logFC"; 
  rest$ids <- rownames(rest);

  names(fstat) <- names(p.value) <- rownames(rest);

  inx.imp <- p.value <= thresh;
  inx.imp <- ifelse(is.na(inx.imp), FALSE, inx.imp);
  sig.num <- length(which(inx.imp == TRUE));
  
  if(sig.num > 0){ 
    sig.p <- p.value[inx.imp];
    sig.mat <- rest[inx.imp,];
    sig.mat[,-ncol(sig.mat)] <- sapply(sig.mat[,-ncol(sig.mat)], function(x) signif(x, 5));
    rownames(sig.mat) <- make.names(rownames(rest)[inx.imp])
    # order the result simultaneously
  }else{
    current.msg <- paste(c("No significant features are detected, please adjust your parameters"), collapse=" ");
    AddMsg(current.msg);
    msgSet$current.msg <<- current.msg;
    return(c(0, nrow(rest)));
  }
  current.msg <- paste(c("A total of", length(which(inx.imp == TRUE)), "significant features were found."), collapse=" ");
  AddMsg(current.msg);
  msgSet$current.msg <<- current.msg;

  both.mat <- both.mat[rownames(rest),]

  if(!is.null(dataSet$enrich_ids)) {
    rest$label <- .invert_named_vector(dataSet$enrich_ids)[as.character(rest$ids)];
  } else {
    rest$label <- rest$ids;
  }
  dataSet$comp.res <- rest;
  rownames(sig.mat) <- gsub("^X(?=[0-9])", "", rownames(sig.mat), perl = TRUE)

  dataSet$sig.mat <- sig.mat;
  if(!is.null(dataSet$enrich_ids)) {
    sig.mat$label <-  .invert_named_vector(dataSet$enrich_ids)[as.character(sig.mat$ids)];
  } else {
    sig.mat$label <- sig.mat$ids;
  }

  if(sig.num> 0){
    res <- 1;
    fileName <- "covariate_result.csv"
    fast.write.csv(sig.mat,file=fileName);
    cov<-list (
      sig.num = sig.num,
      sig.nm = fileName,
      raw.thresh = thresh,
      thresh = -log10(thresh), # only used for plot threshold line
      p.value = p.value,
      p.value.no = both.mat$pval.no,
      p.log = -log10(p.value),
      inx.imp = inx.imp,
      sig.mat = sig.mat
    );
  }else{
    res <- 0;
    cov<-list (
      sig.num = sig.num,
      raw.thresh = thresh,
      thresh = -log10(thresh), # only used for plot threshold line
      p.value = p.value,
      p.value.no = both.mat$pval.no,
      p.log = -log10(p.value),
      inx.imp = inx.imp
    );
  }
  
  # for detail table
  dataSet$analSet$cov <- cov; 
  # for plotting adjp vs p
  dataSet$analSet$cov.mat <- both.mat; 
  both.list <- apply(both.mat, 2, function(x){unname(as.list(x))})

  both.list$thresh <- thresh;
  jsonNm <- gsub(paste0(".", imgFormat), ".json", imgName);
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(both.list, jsonNm, auto_unbox = TRUE, pretty = FALSE);
    
  nonSig <- nrow(dataSet$comp.res) - sig.num;

  RegisterData(dataSet)
  return(c(sig.num, nonSig));
}

# Helper function to invert named vector
.invert_named_vector <- function(input_named_vec) {
  # Get names and values of input named vector
  input_names <- names(input_named_vec)
  input_values <- unname(input_named_vec)
  
  # Invert the named vector
  output_named_vec <- setNames(input_names, input_values)
  
  return(output_named_vec)
}
