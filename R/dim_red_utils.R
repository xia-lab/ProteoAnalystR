##################################################
## R script for ProteoAnalyst
## Description: Computing PCA coordinates
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

SaveClusterJSONLoading <- function(dataName="", fileNm, clustOpt, nb){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    .saveExpressClusterLoadingJSON(dataSet, fileNm, clustOpt,paramSet, nb);
  }else{
    .saveMetaClusterLoadingJSON(dataSet, fileNm, clustOpt,paramSet, nb);
  }
}

SaveClusterJSON <- function(dataName="", fileNm, clustOpt, opt){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    .saveExpressClusterJSON(dataSet, fileNm, clustOpt, paramSet ,opt);
  }else{
    .saveMetaClusterJSON(dataSet, fileNm, clustOpt,paramSet , opt);
  }
}

.preparePcaFeatureMatrix <- function(dat, scaleData=FALSE){
  dat <- as.matrix(dat);
  keep <- apply(dat, 1, function(x) all(is.finite(x)));
  if(scaleData){
    variances <- apply(dat, 1, stats::var);
    keep <- keep & is.finite(variances) & variances > .Machine$double.eps;
  }
  dat <- dat[keep, , drop=FALSE];
  if(nrow(dat) == 0 || ncol(dat) == 0){
    AddErrMsg("PCA requires at least one feature with finite, non-constant values.");
    stop("No valid features for PCA");
  }
  return(dat);
}

.pcaExplainedVariance3D <- function(pca){
  variance <- rep(0, 3);
  imp.pca <- summary(pca)$importance;
  if(!is.null(dim(imp.pca)) && nrow(imp.pca) >= 2 && ncol(imp.pca) > 0){
    dim.count <- min(3, ncol(imp.pca));
    variance[seq_len(dim.count)] <- as.numeric(imp.pca[2, seq_len(dim.count)]);
  }
  variance[!is.finite(variance)] <- 0;
  return(variance);
}

.pcaAxisLabels3D <- function(pca, includeVariance=TRUE){
  if(!includeVariance){
    return(paste("PC", 1:3, sep=""));
  }
  variance <- .pcaExplainedVariance3D(pca);
  return(paste("PC", 1:3, " (", 100*round(variance, 3), "%)", sep=""));
}

.pcaComponentMatrix3D <- function(component.mat){
  component.mat <- as.matrix(component.mat);
  coords <- matrix(0, nrow=nrow(component.mat), ncol=3);
  if(nrow(component.mat) > 0 && ncol(component.mat) > 0){
    dim.count <- min(3, ncol(component.mat));
    coords[, seq_len(dim.count)] <- signif(component.mat[, seq_len(dim.count), drop=FALSE], 5);
  }
  rownames(coords) <- rownames(component.mat);
  colnames(coords) <- paste("Dim", 1:3, sep="");
  return(coords);
}

.pcaJsonCoords3D <- function(component.mat){
  coords <- data.frame(t(.pcaComponentMatrix3D(component.mat)), check.names=FALSE);
  colnames(coords) <- NULL;
  return(coords);
}

.saveMetaClusterJSON <- function(dataSet, fileName, clustOpt,paramSet, opt){
    
    msgSet <- readSet(msgSet, "msgSet");
    paramSet <- readSet(paramSet, "paramSet");
    analSet <- readSet(analSet, "analSet");

    mdata.all <- paramSet$mdata.all;

    inmex.meta <- ov_qs_read("inmex_meta.qs");
    datanm.vec <- names(mdata.all)[mdata.all==1];

    dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
    dat <- inmex.meta$data[, dat.inx, drop=F]; 

    # need to deal with missing values 
    dat <- na.omit(dat);

    pca3d <- list();
    if(clustOpt == "pca"){
        dat <- .preparePcaFeatureMatrix(dat, scaleData=TRUE);
        if(opt == "all"){
            pca <- prcomp(t(dat), center=T, scale=T);
            }else{
            dat <- dat[which(rownames(dat) %in% analSet$loadEntrez),]
            dat <- .preparePcaFeatureMatrix(dat, scaleData=TRUE);
            pca <- prcomp(t(dat), center=T, scale=T);
            }

        pca3d$score$axis <- .pcaAxisLabels3D(pca);
        coords <- .pcaJsonCoords3D(pca$x);
    }else if(clustOpt == "umap"){
        require('uwot');
        # OPTIMIZED: Calculate ncol once
        n_cols <- ncol(dat)
        if(n_cols < 100){
            neighbor_num <- n_cols
        }else{
            neighbor_num <- 100;
        }

        ndat <- as.matrix(t(dat));
        res <- umap(ndat, n_components=3, n_neighbors=neighbor_num);
        pca3d$score$axis <- paste("UMAP dim ", 1:3, sep="");
        coords <- data.frame(t(signif(res, 5)));
    }else{
        require('Rtsne');
        ndat <- as.matrix(t(dat));
        max.perx <- floor((nrow(ndat)-1)/3);
        if(max.perx > 30){
            max.perx <- 30;
        }
        res <- Rtsne(ndat, dims = 3, perplexity=max.perx);
        pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
        coords <- data.frame(t(signif(res$Y, 5)));
    }

    colnames(coords) <- NULL;
    pca3d$score$xyz <- coords;
    pca3d$score$name <- colnames(dat);

    facA <- as.character(inmex.meta$cls.lbl[dat.inx]);
    if(all.numeric(facA)){
        facA <- paste("Group", facA);
    }
    pca3d$score$facA <- facA;

    facB <-  as.character(inmex.meta$data.lbl[dat.inx]);
    if(all.numeric(facB)){
        facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;

    # now set color for each group
    cols <- unique(GetColorSchema(facB));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;

    # add shape sphere, triangles, square, pentagon (first two)
    pca3d$score$shapes <- c("sphere", "triangle");

    if(clustOpt == "pca"){
      mypos <- .pcaComponentMatrix3D(pca$x);
    }else{
      mypos <- t(coords);
      colnames(mypos) <- paste("Dim", 1:3, sep="");
    }
    coords <- data.frame(Class=facA, Data=facB, mypos);

    pos.xyz <- mypos;
    pos.xyz <- unitAutoScale(pos.xyz);
    rownames(pos.xyz) = pca3d$score$name;
    ov_qs_save(pos.xyz, "score_pos_xyz.qs");

    fast.write(coords, file="proteoanalyst_3d_pos.csv");

    pca3d$org <- paramSet$data.org
    pca3d$analType <- paramSet$anal.type
    pca3d$naviString <- "Scatter 3D"
    ov_qs_save(pca3d, "pca3d.qs");

    paramSet$jsonNms$pcascore <- fileName
    paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
    # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
    jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
    msgSet$current.msg <- "Annotated data is now ready for 3D visualization!";
    saveSet(msgSet, "msgSet");
    saveSet(paramSet, "paramSet");

    return(1);
}


.saveMetaClusterLoadingJSON <- function(dataSet, fileName, clustOpt, paramSet, nb){
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  mdata.all <- paramSet$mdata.all;

  inmex.meta <- ov_qs_read("inmex_meta.qs");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  nb <- as.numeric(5000) # set to max 5000 datapoints
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$data[, dat.inx, drop=F]; 
  
  # need to deal with missing values 
  dat <- na.omit(dat);
  dat <- .preparePcaFeatureMatrix(dat, scaleData=TRUE);
  variances <- apply(dat,1, function(x){var(x)})
  df <- data.frame(var = variances, inx = seq.int(1,length(variances)))
  df <- df[order(-df$var),];

  #do not take subset of loading data points now
  if(nb < length(df$inx)){
    inx <- df$inx[c(1:nb)];
  }else{
    inx <- df$inx;
  }
  dat <- dat[inx,];
  
  pca3d <- list();
  
  pca <- prcomp(t(dat), center=T, scale=T);    
  pca3d$score$axis <- .pcaAxisLabels3D(pca, includeVariance=FALSE);
  coords <- .pcaJsonCoords3D(pca$rotation);
  
  pca3d$score$xyz <- coords;
  pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation), paramSet$data.org, paramSet$data.idType);
  pca3d$score$entrez <- rownames(pca$rotation);
  
  analSet$loadEntrez <- pca3d$score$entrez
  mypos <- .pcaComponentMatrix3D(pca$rotation);
  rownames(mypos) <- analSet$loadEntrez;
  mypos <- unitAutoScale(mypos);
  ov_qs_save(mypos, "loading_pos_xyz.qs");
  
  coords <- data.frame(mypos);
  fast.write(coords, file="proteoanalyst_loadings_3d_pos.csv");
  
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
  paramSet$jsonNms$pcaload <- fileName;
  ov_qs_save(pca3d, "pca3d.qs");
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
  msgSet$current.msg <- "Annotated data is now ready for 3D visualization!";
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  saveSet(analSet, "analSet");

  return(1);
}


.saveExpressClusterLoadingJSON <- function(dataSet, fileName, clustOpt, paramSet, nb){  
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  nb <- as.numeric(25000) # set to max 5000 datapoints
  if(clustOpt == "pca"){
    dat <- .preparePcaFeatureMatrix(dat, scaleData=TRUE);
    pca <- prcomp(t(dat), center=T, scale=T);
    pca3d$score$axis <- .pcaAxisLabels3D(pca);
    coords <- .pcaJsonCoords3D(pca$rotation);
    
    pca3d$score$xyz <- coords;
    pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation), paramSet$data.org, paramSet$data.idType);
    pca3d$score$entrez <-rownames(pca$rotation);
    weights <- .pcaExplainedVariance3D(pca);
    if(!any(weights > 0)){
      weights <- rep(1, 3);
    }
    mypos <- .pcaComponentMatrix3D(pca$rotation);
    meanpos <- apply(abs(mypos),1, function(x){weighted.mean(x, weights)})
    df <- data.frame(pos = meanpos, inx = seq.int(1,length(meanpos)))
    df <- df[order(-df$pos),]
    
    if(nrow(df) > nb){
      inx <- df$inx[c(1:nb)]
      mypos <- mypos[inx,, drop=FALSE];
      pca3d$score$xyz <- coords[, inx, drop=FALSE]
      pca3d$score$name <- pca3d$score$name[inx]
      pca3d$score$entrez <- pca3d$score$entrez[inx]
    }
  }

  pca3d$cls <- dataSet$meta.info;
  # see if there is secondary
  analSet$loadEntrez <- pca3d$score$entrez
  rownames(mypos) <- pca3d$score$name;
  rownames(mypos) <- analSet$loadEntrez;
  mypos <- unitAutoScale(mypos);
  ov_qs_save(mypos, "loading_pos_xyz.qs");
  
  fast.write(mypos, file="proteoanalyst_3d_load_pos.csv");
  ov_qs_save(pca3d, "pca3d.qs");
  paramSet$jsonNms$pcaload <- fileName
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
  msgSet$current.msg <- "Annotated data is now ready for PCA 3D visualization!";
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  saveSet(analSet, "analSet");

  return(1);
}

# single expression data
.saveExpressClusterJSON <- function(dataSet, fileName, clustOpt,paramSet, opt){
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  
  if(clustOpt == "pca"){
    dat <- .preparePcaFeatureMatrix(dat, scaleData=FALSE);
    #if(opt == "all"){
      pca <- prcomp(t(dat));
   #}else{
    #  dat <- dat[which(rownames(dat) %in% analSet$loadEntrez),]
    #  pca <- prcomp(t(dat), center=T, scale=T);
    #}
    pca3d$score$axis <- .pcaAxisLabels3D(pca);
    coords <- .pcaJsonCoords3D(pca$x);
  }else if(clustOpt == "umap"){
    require('uwot');
    # OPTIMIZED: Calculate ncol once
    n_cols <- ncol(dat)
    if(n_cols < 100){
      neighbor_num <- n_cols
    }else{
      neighbor_num <- 100;
    }
    dat <- as.matrix(t(dat));
    res <- umap(dat, n_components=3, n_neighbors=neighbor_num);
    pca3d$score$axis <- paste("UMAP dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res, 5)));
    
  }else{ # tsne
    require('Rtsne');
    dat <- as.matrix(t(dat));
    max.perx <- floor((nrow(dat)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    res <- Rtsne(dat, dims = 3, perplexity=max.perx);
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res$Y, 5)));
  }
  
  colnames(coords) <- NULL;
  pca3d$score$xyz <- coords;
  pca3d$score$name <- colnames(dataSet$data.norm);
  
  if(clustOpt == "pca"){
    pos.xyz <- .pcaComponentMatrix3D(pca$x);
  }else{
    pos.xyz <- t(coords);
    colnames(pos.xyz) <- paste("Dim", 1:3, sep="");
  }
  pos.xyz <- unitAutoScale(pos.xyz);
  rownames(pos.xyz) = colnames(dataSet$data.norm);
  ov_qs_save(pos.xyz, "score_pos_xyz.qs");
  
  facA <- as.character(dataSet$fst.cls);
  if(all.numeric(facA)){
    facA <- paste("Group", facA);
  }
  pca3d$score$facA <- facA;
  
  if(clustOpt == "pca"){
    mypos <- .pcaComponentMatrix3D(pca$x);
  }else{
    mypos <- t(coords);
    colnames(mypos) <- paste("Dim", 1:3, sep="");
  }
  # see if there is secondary
  if(length(dataSet$sec.cls) > 1){
    facB <- as.character(dataSet$sec.cls);
    if(all.numeric(facB)){
      facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;
    
    # set shape based on the first group
    pca3d$score$shapes <- c("sphere", "triangle");
    
    # now set color based on 2nd group
    cols <- unique(GetColorSchema(facB));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
    pca3d$score$colors <- cols;
    
    mypos <- data.frame(factorA=facA, factorB=facB, mypos);
  }else{
    # now set color based on first group
    cols <- unique(GetColorSchema(facA));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;
    mypos <- data.frame(factorA=facA, mypos);
  }
  
  pca3d$cls <- dataSet$meta.info;
  pca3d$org <- paramSet$data.org
  pca3d$analType <- paramSet$anal.type
  pca3d$naviString <- "Scatter 3D"
  ov_qs_save(pca3d, "pca3d.qs");
  paramSet$jsonNms$pcascore <- fileName
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(fileName))
  rownames(mypos) <- colnames(dataSet$data.norm);
  
  fast.write(mypos, file="proteoanalyst_3d_pos.csv");
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(pca3d, fileName, auto_unbox = TRUE, pretty = FALSE);
  msgSet$current.msg <- "Annotated data is now ready for PCA 3D visualization!";
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  return(1);
}


ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA"){
  Sys.setenv(RGL_USE_NULL = TRUE)
  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;
  level <- as.numeric(level)
  names = strsplit(names.vec, "; ")[[1]]
  pos.xyz <-ov_qs_read("score_pos_xyz.qs");

  tryCatch({
    if (!file.exists("score_pos_xyz.qs")) {
      sink(filenm); cat("{}"); sink()
      return(filenm)
    }

    pos.xyz <- ov_qs_read("score_pos_xyz.qs")

    inx <- rownames(pos.xyz) %in% names
    coords <- as.matrix(pos.xyz[inx, c(1:3)])

    if (nrow(coords) < 4) {
      sink(filenm); cat(RJSONIO::toJSON(list())); sink()
      return(filenm)
    }

    # Only rgl::ellipse3d in subprocess
    mesh <- rsclient_isolated_exec(
      func_body = function(input_data) {
        Sys.setenv(RGL_USE_NULL = TRUE)
        pos <- cov(input_data$coords, y = NULL, use = "everything")
        center <- colMeans(input_data$coords)
        t_val <- sqrt(qchisq(input_data$level, 3))
        mesh <- list()
        mesh[[1]] <- rgl::ellipse3d(x = as.matrix(pos), centre = center, t = t_val)
        mesh
      },
      input_data = list(coords = coords, level = level),
      packages = c("rgl", "qs"),
      timeout = 120,
      output_type = "qs"
    )

    # Write JSON in master (subprocess may have different wd)
    if (!is.list(mesh) || !isFALSE(mesh$success)) {
      sink(filenm); cat(RJSONIO::toJSON(mesh)); sink()
    }
  }, error = function(e) {
    message("[ComputeEncasing] ", e$message)
    sink(filenm); cat("{}"); sink()
  })
  return(filenm)
}

ComputeEncasingBatch <- function(filenm, type, groups_json, level = 0.95, omics = "NA") {
  tryCatch({
    level <- as.numeric(level)

    if (!file.exists("score_pos_xyz.qs")) {
      sink(filenm); cat("{}"); sink()
      return(filenm)
    }
    pos.xyz <- qs::qread("score_pos_xyz.qs")

    groups_list <- RJSONIO::fromJSON(groups_json)
    if (is.data.frame(groups_list)) {
      groups_list <- split(groups_list, seq_len(nrow(groups_list)))
    }

    group_data <- lapply(seq_along(groups_list), function(i) {
      g <- groups_list[[i]]
      if (is.character(g)) { gn <- unname(g["grpName"]); ids <- unname(g["names"]) }
      else if (is.data.frame(g)) { gn <- g$grpName[1]; ids <- g$names[1] }
      else { gn <- g$grpName; ids <- g$names }
      nms <- strsplit(ids, "; ")[[1]]
      inx <- rownames(pos.xyz) %in% nms
      list(group = gn, coords = as.matrix(pos.xyz[inx, c(1:3)]))
    })

    result_list <- rsclient_isolated_exec(
      func_body = function(input_data) {
        Sys.setenv(RGL_USE_NULL = TRUE)
        lapply(input_data$groups, function(g) {
          coords <- g$coords
          if (nrow(coords) < 4) return(list(group = g$group, mesh = list(), error = "Insufficient points"))
          tryCatch({
            pos <- cov(coords, y = NULL, use = "everything")
            center <- colMeans(coords)
            t_val <- sqrt(qchisq(input_data$level, 3))
            mesh <- list()
            mesh[[1]] <- rgl::ellipse3d(x = as.matrix(pos), centre = center, t = t_val)
            list(group = g$group, mesh = mesh, error = NULL)
          }, error = function(e) {
            list(group = g$group, mesh = list(), error = e$message)
          })
        })
      },
      input_data = list(groups = group_data, level = level),
      packages = c("rgl", "qs"),
      timeout = 120,
      output_type = "qs"
    )

    if (is.list(result_list) && !isFALSE(result_list$success)) {
      sink(filenm); cat(RJSONIO::toJSON(result_list)); sink()
    } else {
      sink(filenm); cat("[]"); sink()
    }
  }, error = function(e) {
    message("[ComputeEncasingBatch] ", e$message)
    sink(filenm); cat("[]"); sink()
  })
  return(filenm)
}

unitAutoScale <- function(df){
    df <- as.data.frame(df)
    row.nms <- rownames(df);
    col.nms <- colnames(df);
    df <- as.data.frame(lapply(df, function(x){
      sd.x <- sd(x, na.rm=T);
      if(!is.finite(sd.x) || sd.x == 0){
        return(rep(0, length(x)));
      }
      res <- (x - mean(x, na.rm=T))/sd.x;
      res[!is.finite(res)] <- 0;
      return(res);
    }), check.names=FALSE);
    df <- as.matrix(df);
    rownames(df) <- row.nms;
    colnames(df) <- col.nms;
    maxVal <- max(abs(df), na.rm=T)
    if(!is.finite(maxVal) || maxVal == 0){
      df[!is.finite(df)] <- 0;
      return(df);
    }
    df<- df/maxVal
    df[!is.finite(df)] <- 0;
    return(df)
}
