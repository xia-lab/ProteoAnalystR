##################################################
## R script for ProteoAnalyst
## Description: functions for DE analysis
##
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' Select Metadata Factors for DE Analysis
#'
#' This function is used to select metadata factors for a Differential Expression (DE) analysis.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param meta0 The primary metadata factor
#' @param meta1 The secondary metadata factor
#' @param block1 The blocking factor for
#'
#' @author Guangyan Zhou, \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @export
#'
SetSelectedMetaInfo <- function(dataName="", meta0, meta1, block1){
  dataSet <- readDataset(dataName);
  if(meta0 == "NA"){
    RegisterData(dataSet, 0);
  }else{
    rmidx <- which(dataSet$meta.info[, meta0]=="NA")
    if(meta1 != "NA"){
        rmidx <- c(rmidx,which(dataSet$meta.info[, meta1]=="NA"))
    }
    if(length(rmidx)>0){
       meta<- dataSet$meta.info[-rmidx,]
     #print(meta);
       for(col in 1:ncol(meta)){
        meta[,col]<- droplevels(meta[,col])
       }
       dataSet$rmidx <- rmidx
    }else{
        meta<- dataSet$meta.info
    }
    cls <- meta[, meta0];
    dataSet$fst.cls <- cls; # for PCA plotting
    block <- NULL;
    dataSet$sec.cls <- "NA";
    if(meta1 != "NA"){
      if(block1){
        block <- meta[, meta1];
      }else{ # two factor
        cls <- interaction(meta[, c(meta0, meta1)], sep = "_", lex.order = TRUE);
      }
      dataSet$sec.cls <- meta[, meta1]; # for pca coloring
    }
    dataSet$analysisVar <- meta0 
    dataSet$secondVar <- meta1
    dataSet$cls <- cls; # record main cls;
    dataSet$block <- block;
    RegisterData(dataSet, levels(cls)[levels(cls)!="NA"]);
  }
}

#' Perform Differential Analysis
#'
#' This function performs differential analysis based on different types of comparisons.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param anal.type The type of analysis to perform. Options: "default", "custom", "time", "reference", "nested".
#' @param par1 Parameter 1, depending on the analysis type.
#' @param par2 Parameter 2, depending on the analysis type.
#' @param nested.opt Option for nested analysis. Options: "intonly" (default), "all".
#' @param robustTrend Logical. If TRUE, use robust trend test.
#'
#' @return Results of the differential analysis.
#' @details default: all pair-wise comparison (A-B) + (B-C) + (A-C), custom: only compare two groups (A-C), time: compare consecutive groups (B-A) + (C-B), reference: all others against common reference (A-C) + (B-C), nested: (A-B)+(C-D) 
#' @author Guangyan Zhou, \email{guangyan.zhou@mail.mcgill.ca}
#' @export
#' @license MIT
#'
PerformDEAnal<-function (dataName="", anal.type = "default", par1 = NULL, par2 = NULL, nested.opt = "intonly", robustTrend=F){

  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");

  #msg("[DE] PerformDEAnal entry: data=", dataName, " anal.type=", anal.type, " par1=", par1, " par2=", par2)
  #msg("[DE] dataSet dims: data.norm=", nrow(dataSet$data.norm), "x", ncol(dataSet$data.norm), " method=", dataSet$de.method)
  #msg("[DE] data.norm rownames (first 10): ", paste(utils::head(rownames(dataSet$data.norm), 10), collapse=", "))
  #print(head(dataSet$meta.info));

  # Check for duplicate rownames before annotation
  if (any(duplicated(rownames(dataSet$data.norm)))) {
    dup_ids <- rownames(dataSet$data.norm)[duplicated(rownames(dataSet$data.norm))]
    #msg("[DE] WARNING: Found ", length(dup_ids), " duplicate rownames BEFORE annotation remapping")
    #msg("[DE] Duplicate IDs: ", paste(utils::head(dup_ids, 20), collapse=", "))
  }

  # Ensure DE uses annotated IDs if available
  if (file.exists("annotation.qs")) {
    anot.id <- ov_qs_read("annotation.qs")
    #msg("[DE] Found annotation.qs with ", length(anot.id), " entries")
    if (!is.null(anot.id)) {
      if (!is.null(names(anot.id))) {
        mapped <- anot.id[rownames(dataSet$data.norm)]
      } else {
        mapped <- anot.id
      }
      hit <- !is.na(mapped)
      #msg("[DE] Annotation mapping: ", sum(hit), " hits out of ", length(hit), " rows")
      if (any(hit)) {
        old_names <- rownames(dataSet$data.norm)[hit]
        new_names <- mapped[hit]
        #msg("[DE] Remapping examples (old -> new):")
        for (i in 1:min(5, sum(hit))) {
          idx <- which(hit)[i]
          #msg("  ", old_names[i], " -> ", new_names[i])
        }

        # Check for duplicates in mapped IDs BEFORE assignment
        if (any(duplicated(new_names))) {
          dup_ids <- unique(new_names[duplicated(new_names)])
          #msg("[DE] WARNING: Mapped IDs contain ", length(dup_ids), " duplicates")
          #msg("[DE] Duplicate IDs: ", paste(utils::head(dup_ids, 20), collapse=", "))
          #msg("[DE] This means multiple protein IDs map to the same gene ID")
          #msg("[DE] Aggregating duplicates by mean BEFORE remapping...")

          # Build a dataframe with protein IDs as a column (to avoid duplicate rowname error)
          df <- as.data.frame(dataSet$data.norm)
          # Create new ID column with mapped and unmapped IDs
          new_id_vector <- rownames(dataSet$data.norm)
          new_id_vector[hit] <- new_names  # Replace only the mapped ones
          df$Protein <- new_id_vector

          # Aggregate duplicates by mean
          df <- stats::aggregate(. ~ Protein, data = df, FUN = function(x) mean(x, na.rm = TRUE))
          rownames(df) <- df$Protein
          df$Protein <- NULL
          dataSet$data.norm <- as.matrix(df)
          #msg("[DE] After aggregation: ", nrow(dataSet$data.norm), "x", ncol(dataSet$data.norm))
          #msg("[DE] Aggregated data head IDs: ", paste(utils::head(rownames(dataSet$data.norm), 10), collapse=", "))
        } else {
          # No duplicates, safe to assign
          rownames(dataSet$data.norm)[hit] <- new_names
          #msg("[DE] data.norm remapped to annotated IDs (no duplicates); head now: ", paste(utils::head(rownames(dataSet$data.norm), 5), collapse=", "))
        }
      }
    }
  }

  #msg("[DE] PerformDEAnal start data=", dataName, " method=", dataSet$de.method, " annotated=", dataSet$annotated, " id.current=", dataSet$id.current,
  #        " norm.opt=", paramSet$norm.opt)
  #msg("[DE] data.norm head IDs=", paste(utils::head(rownames(dataSet$data.norm), 5), collapse=", "),
  #        " cols=", paste(utils::head(colnames(dataSet$data.norm), 5), collapse=", "))
  if (file.exists("data.anot.qs")) {
    da <- try(ov_qs_read("data.anot.qs"), silent = TRUE)
    if (!inherits(da, "try-error")) {
      #msg("[DE] data.anot.qs head IDs=", paste(utils::head(rownames(da), 5), collapse=", "))
      if (any(duplicated(rownames(da)))) {
        #msg("[DE] WARNING: data.anot.qs has duplicate rownames!")
      }
    }
  }
  if (dataSet$de.method == "deseq2") {
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .run.deseq(dataSet, anal.type, par1, par2, nested.opt);
  }else if (dataSet$de.method == "limma"){
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_limma_edger(dataSet, robustTrend);
  }else if (dataSet$de.method == "deqms"){
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_deqms(dataSet, robustTrend);
  }else if (dataSet$de.method == "edger"){
    dataSet <- prepareEdgeRContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_limma_edger(dataSet, robustTrend);
  }else{ #dataSet$de.method == "wtt"
    dataSet <- prepareContrast(dataSet, anal.type, par1, par2, nested.opt);
    dataSet <- .perform_williams_trend(dataSet, robustTrend);
  }

  # Perform peptide-level DE analysis if peptide data exists
  if (file.exists("peptide_level_data.qs")) {
    PerformPeptideLevelDEAnal(dataName)
  }

  return(RegisterData(dataSet));
}

.run.deseq <- function(dataSet, anal.type, par1, par2, nested.opt) {

  # Read annotated data in parent (not available in child)
  data.anot <- ov_qs_read("data.anot.qs")

  bridge_in <- paste0(tempdir(), "/bridge_", paste0(sample(letters,6,replace=TRUE), collapse=""), "_in.qs")
  bridge_out <- sub("_in.qs", "_out.qs", bridge_in)
  ov_qs_save(list(
    data.anot = data.anot,
    rmidx = dataSet$rmidx,
    fst.cls = dataSet$fst.cls,
    analysisVar = dataSet$analysisVar,
    block = dataSet$block,
    anal.type = anal.type,
    par1 = par1
  ), bridge_in, preset = "fast")
  on.exit(unlink(c(bridge_in, bridge_out)), add = TRUE)

  run_func_via_rsclient(
    func = function(wd, bridge_in, bridge_out) {
      setwd(wd)
      require(DESeq2)
      input <- ov_qs_read(bridge_in)

      rmidx <- input$rmidx
      fst.cls <- input$fst.cls
      analysisVar <- input$analysisVar
      block <- input$block
      anal.type <- input$anal.type
      par1 <- input$par1
      data.anot <- input$data.anot

      formatLevel <- function(x) {
        if (grepl("^[0-9]", x)) paste0(analysisVar, "_", x) else x
      }
      parse_contrast_groups <- function(cstr) {
        comps <- strsplit(cstr, " vs\\.?\\s*")[[1]]
        if (length(comps) != 2) stop(paste("Invalid contrast format:", cstr))
        return(comps)
      }

      if (length(rmidx) > 0) {
        data.anot <- data.anot[, -rmidx]
      }

      if (any(grepl("(^[0-9]+).*", fst.cls))) {
        fst.cls <- paste0(analysisVar, "_", fst.cls)
      }
      fst.cls <- as.character(fst.cls)

      all_conditions <- unique(fst.cls)
      contrast_list <- list()

      colData <- data.frame(condition = factor(fst.cls, levels = all_conditions))

      if (!is.null(block)) {
        colData$block <- factor(block)
        design <- ~ block + condition
      } else {
        design <- ~ condition
      }

      if (anal.type == "default") {
        n_conditions <- length(all_conditions)
        for (i in 1:(n_conditions - 1))
          for (j in (i + 1):n_conditions) {
            contrast_name <- paste0(all_conditions[i], " vs ", all_conditions[j])
            contrast_list[[contrast_name]] <-
              c("condition", all_conditions[j], all_conditions[i])
          }
      } else if (anal.type == "reference") {
        ref <- formatLevel(par1)
        if (!(ref %in% all_conditions)) stop("Reference level not found: ", ref)
        for (cond in setdiff(all_conditions, ref)) {
          contrast_name <- paste0(ref, " vs ", cond)
          contrast_list[[contrast_name]] <- c("condition", cond, ref)
        }
      } else if (anal.type == "custom") {
        comps <- parse_contrast_groups(par1)
        comps <- vapply(comps, formatLevel, "")
        if (!all(comps %in% all_conditions)) stop("Invalid custom contrast: ", par1)
        contrast_name <- paste0(comps[1], " vs ", comps[2])
        contrast_list[[contrast_name]] <- c("condition", comps[2], comps[1])
      }

      count_mat <- data.anot
      if (any(!is.finite(count_mat))) count_mat[!is.finite(count_mat)] <- 0
      dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                                    colData = colData, design = design)
      set.seed(123)
      dds <- DESeq(dds, betaPrior = FALSE)
      ov_qs_save(dds, "deseq.res.obj.rds")

      results_list <- list()
      if (length(contrast_list) > 0) {
        for (contrast_name in names(contrast_list)) {
          res <- results(dds, contrast = contrast_list[[contrast_name]],
                         independentFiltering = FALSE, cooksCutoff = Inf)
          topFeatures <- data.frame(res@listData)
          rownames(topFeatures) <- rownames(res)
          colnames(topFeatures) <- sub("padj", "adj.P.Val", colnames(topFeatures))
          colnames(topFeatures) <- sub("pvalue", "P.Value", colnames(topFeatures))
          colnames(topFeatures) <- sub("log2FoldChange", "logFC", colnames(topFeatures))
          topFeatures <- topFeatures[c("logFC","baseMean","lfcSE","stat","P.Value","adj.P.Val")]
          topFeatures <- topFeatures[order(topFeatures$P.Value), ]
          results_list[[contrast_name]] <- topFeatures
        }
      } else {
        # Inline .get.interaction.results() — not available in child
        interaction_name <- grep("factorA.*factorB.*", resultsNames(dds), value = TRUE)
        if (length(interaction_name) == 0) {
          stop("No interaction term found in model.")
        }
        res <- results(dds, name = interaction_name[1], independentFiltering = FALSE, cooksCutoff = Inf)
        topFeatures <- data.frame(res@listData)
        rownames(topFeatures) <- rownames(res)
        colnames(topFeatures) <- sub("padj", "adj.P.Val", colnames(topFeatures))
        colnames(topFeatures) <- sub("pvalue", "P.Value", colnames(topFeatures))
        colnames(topFeatures) <- sub("log2FoldChange", "logFC", colnames(topFeatures))
        topFeatures <- topFeatures[c("logFC", "baseMean", "lfcSE", "stat", "P.Value", "adj.P.Val")]
        topFeatures <- topFeatures[order(topFeatures$P.Value), ]
        results_list[[1]] <- topFeatures
      }

      ov_qs_save(results_list, bridge_out, preset = "fast")
    },
    args = list(wd = getwd(), bridge_in = bridge_in, bridge_out = bridge_out),
    timeout_sec = 300
  )

  results_list <- if (file.exists(bridge_out)) ov_qs_read(bridge_out) else NULL

  dataSet$comp.res.list <- results_list
  dataSet$comp.res <- results_list[[1]]
  ov_qs_save(results_list, file = "dat.comp.res.qs")
  return(dataSet)
}


prepareContrast <-function(dataSet, anal.type = "reference", par1 = NULL, par2 = NULL, nested.opt = "intonly"){
  
  msgSet <- readSet(msgSet, "msgSet");
  cat(anal.type, par1, par2, nested.opt, "\n")
  set.seed(1337);
  myargs <- list();
  cls <- dataSet$cls;
  dataSet$comp.type <- anal.type;
  grp.nms <- levels(cls);
  analysisVar <- dataSet$analysisVar
  if(dataSet$cont.inx[analysisVar] |  any(grepl("(^[0-9]+).*", grp.nms))){
    if(grepl( "vs",par1)){
      par1 <- strsplit(par1, " vs. ")[[1]]
      par1 <- paste0(analysisVar,"_",par1[1]," vs. ",analysisVar,"_",par1[2])
    }else{
      par1<- paste0(analysisVar,"_",par1)
    }
    if(par2 != "NA"){
      if(grepl( "vs",par2)){
        par2 <- strsplit(par2, " vs. ")[[1]]
        par2 <- paste0(analysisVar,"_",par2[1]," vs. ",analysisVar,"_",par2[2])
      }else{
        par2<- paste0(analysisVar,"_",par2)
      }
    }

    if(any(grepl("(^[0-9]+).*",  colnames(dataSet$design)))){
      colnames(dataSet$design) = as.character(sapply( colnames(dataSet$design),function(x) paste0(analysisVar,"_",x)))
    }
    grp.nms <- paste0(analysisVar,"_",grp.nms)
    
  }

  dataSet$par1 <- par1;
  
  if (anal.type == "default") {
    inx <- 0
    for (m in 1:(length(grp.nms) - 1)) {
      for (n in (m + 1):length(grp.nms)) {
        inx <- inx + 1
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep = "");
      }
    }
    filename <- "sig_features_pairwise";
  } else if (anal.type == "time") {
    for (i in 2:length(grp.nms)) {
      myargs[[i - 1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep = "")
    }
    filename <- "sig_features_time_series";
  } else if (anal.type == "custom") {
    grp.nms <- strsplit(par1, " vs. ")[[1]]
    myargs[[1]] <- paste(grp.nms, collapse = "-")
    dataSet$grp.nms <- grp.nms;
    filename <- paste("sig_features_", paste(grp.nms, collapse = "_vs_"), sep = "")
    dataSet$contrast <- paste(grp.nms, collapse = "_vs_");
  } else if (anal.type == "reference") {
    ref <- par1;
    cntr.cls <- grp.nms[grp.nms != ref]
    myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""));
    dataSet$ref <- ref; 
    filename <- paste("sig_features_reference_", ref, sep = "");
  } else if (anal.type == "nested") {
    grp.nms1 <- strsplit(par1, " vs. ")[[1]]
    grp.nms2 <- strsplit(par2, " vs. ")[[1]]
    if (all(grp.nms1 == grp.nms2)) {
      msgSet$current.msg <-"The two nested groups are the same. Please choose two different groups."
      saveSet(msgSet, "msgSet");      
      return(0)
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2))
    if (nested.opt == "intonly") {
      dataSet$nested.int.opt <- "True";
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    } else {
      dataSet$nested.int.opt <- "False";
      myargs[[1]] <- paste(grp.nms1, collapse = "-")
      myargs[[2]] <- paste(grp.nms2, collapse = "-")
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    }
    dataSet$contrast <- paste(paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
    filename <- paste("sig_features_nested_", paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
  } else {
    print(paste("Not supported: ", anal.type))
  }
    
  dataSet$filename <- paste0(filename, "_", dataSet$de.method);
  require(limma);
  design <- dataSet$design;
  myargs[["levels"]] <- design;
  dataSet$contrast.type <- anal.type;
  contrast.matrix <- do.call(makeContrasts, myargs);
  dataSet$contrast.matrix <- contrast.matrix;
  RegisterData(dataSet);
  return(dataSet);
}

.perform_limma_edger <- function(dataSet, robustTrend = FALSE) {
  ## ------------------------------------------------------------------ ##
  ## 1 · Input checks & dependencies                                    ##
  ## ------------------------------------------------------------------ ##
  require(limma)
  require(edgeR)

  if (is.null(dataSet$design) || is.null(dataSet$contrast.matrix)) {
    stop("design and/or contrast.matrix missing in dataSet. Run prepareEdgeRContrast() first.")
  }

  design           <- dataSet$design
  contrast.matrix  <- dataSet$contrast.matrix
  contrast.names   <- colnames(contrast.matrix)

  paramSet <- readSet(paramSet, "paramSet")
  msgSet   <- readSet(msgSet,   "msgSet")

  data.norm <- if (length(dataSet$rmidx) > 0) {
    dataSet$data.norm[, -dataSet$rmidx, drop = FALSE]
  } else {
    dataSet$data.norm
  }

  if (dataSet$de.method == "limma") {
    if (is.null(dataSet$block)) {
      fit <- lmFit(data.norm, design)
    } else {
      corfit <- duplicateCorrelation(data.norm, design, block = dataSet$block)
      fit <- lmFit(data.norm, design, block = dataSet$block, correlation = corfit$consensus)
    }
    
    if (!is.fullrank(design)) {
      msgSet$current.msg <- "This metadata combination is not full rank! Please use other combination.";
      saveSet(msgSet, "msgSet");  
      return(0)
    }
    
    df.residual <- fit$df.residual
    if (all(df.residual == 0)) {
      msgSet$current.msg <- "All residuals equal 0. There is not enough replicates in each group (no residual degrees of freedom)!";
      saveSet(msgSet, "msgSet");  
      return(0);
    }
    fit2 <- contrasts.fit(fit, contrast.matrix);
    fit2 <- eBayes(fit2, trend=robustTrend, robust=robustTrend);

result.list <- list()
for (nm in colnames(contrast.matrix)) {
  tbl <- topTable(fit2, coef = nm, number = Inf, adjust.method = "fdr")
  colnames(tbl)[colnames(tbl) == "FDR"] <- "adj.P.Val"
  result.list[[nm]] <- tbl
}

dataSet$comp.res.list      <- result.list     
  dataSet$comp.res <- result.list[[1]]

  } else if (dataSet$de.method == "edger") {

    set.seed(1)

    # Retrieve the raw (un‑normalised) count matrix
    cnt.mat <- .get.annotated.data()
    if (length(dataSet$rmidx) > 0)
      cnt.mat <- cnt.mat[, -dataSet$rmidx, drop = FALSE]

    grp.fac <- factor(dataSet$cls)

    if (!is.null(dataSet$block)) {
      blk.fac <- factor(dataSet$block)
      design  <- model.matrix(~ grp.fac + blk.fac)
    } else {
      # Use the stored design if created with ~0+grp.fac; otherwise rebuild
      if (is.null(attr(design, "assign")))
        design <- model.matrix(~ 0 + grp.fac)
    }

    y <- DGEList(counts = cnt.mat, group = grp.fac)
    y <- calcNormFactors(y)

    ## Dispersions
    y <- estimateGLMCommonDisp(y, design, verbose = FALSE)
    y <- tryCatch(
      estimateGLMTrendedDisp(y, design),
      error   = function(e) { msgSet$current.msg <- e$message ; saveSet(msgSet, "msgSet"); return(0) },
      warning = function(w) { msgSet$current.msg <- c(msgSet$current.msg, w$message); saveSet(msgSet, "msgSet"); }
    )
    y <- estimateGLMTagwiseDisp(y, design)

    fit <- glmFit(y, design)

    result.list <- vector("list", length(contrast.names))
    names(result.list) <- contrast.names
    for (nm in contrast.names) {
      lrt <- glmLRT(fit, contrast = contrast.matrix[, nm])
      tbl <- topTags(lrt, n = Inf)$table
      colnames(tbl)[colnames(tbl) == "FDR"]    <- "adj.P.Val"
      colnames(tbl)[colnames(tbl) == "PValue"] <- "P.Value"
        # Order tbl by P.Value
        if ("P.Value" %in% colnames(tbl)) {
          tbl <- tbl[order(tbl$P.Value), ]
        }
      result.list[[nm]] <- tbl
    }
  dataSet$comp.res.list <- result.list
  dataSet$comp.res <- result.list[[1]]

}

  return(dataSet);

}

.perform_deqms <- function(dataSet, robustTrend = FALSE) {
  ## ------------------------------------------------------------------ ##
  ## DEqMS: proteomics-aware differential expression using PSM counts   ##
  ## ------------------------------------------------------------------ ##
  require(limma)

  if (!requireNamespace("DEqMS", quietly = TRUE)) {
    stop("DEqMS requires the 'DEqMS' package. Install with: BiocManager::install('DEqMS')")
  }
  require(DEqMS)

  if (is.null(dataSet$design) || is.null(dataSet$contrast.matrix)) {
    stop("design and/or contrast.matrix missing in dataSet. Run prepareContrast() first.")
  }

  design          <- dataSet$design
  contrast.matrix <- dataSet$contrast.matrix
  contrast.names  <- colnames(contrast.matrix)

  paramSet <- readSet(paramSet, "paramSet")
  msgSet   <- readSet(msgSet,   "msgSet")

  data.norm <- if (length(dataSet$rmidx) > 0) {
    dataSet$data.norm[, -dataSet$rmidx, drop = FALSE]
  } else {
    dataSet$data.norm
  }

  ## ------------------------------------------------------------------ ##
  ## 1 · Standard limma workflow                                        ##
  ## ------------------------------------------------------------------ ##
  if (is.null(dataSet$block)) {
    fit <- lmFit(data.norm, design)
  } else {
    corfit <- duplicateCorrelation(data.norm, design, block = dataSet$block)
    fit <- lmFit(data.norm, design, block = dataSet$block, correlation = corfit$consensus)
  }

  if (!is.fullrank(design)) {
    msgSet$current.msg <- "This metadata combination is not full rank! Please use other combination.";
    saveSet(msgSet, "msgSet");
    return(0)
  }

  df.residual <- fit$df.residual
  if (all(df.residual == 0)) {
    msgSet$current.msg <- "All residuals equal 0. There is not enough replicates in each group (no residual degrees of freedom)!";
    saveSet(msgSet, "msgSet");
    return(0);
  }

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend = robustTrend, robust = robustTrend)

  ## ------------------------------------------------------------------ ##
  ## 2 · Extract PSM/peptide/spectra count information                  ##
  ## ------------------------------------------------------------------ ##
  # Try to find PSM count column from original data or annotations
  psm.count <- NULL
  count.col.name <- NULL

  # Check if pepcount was saved during data loading (from MaxQuant, DIA-NN, etc.)
  # First try to load from file if not in memory
  if (is.null(dataSet$pepcount) && file.exists("pepcount.qs")) {
    #msg("[DEqMS] Loading pepcount from pepcount.qs file")
    dataSet$pepcount <- ov_qs_read("pepcount.qs")
    #msg("[DEqMS] DEBUG: Loaded pepcount length=", length(dataSet$pepcount))
    #msg("[DEqMS] DEBUG: Loaded pepcount range: ", min(dataSet$pepcount, na.rm=TRUE), " - ", max(dataSet$pepcount, na.rm=TRUE))
  }

  if (!is.null(dataSet$pepcount)) {
    #msg("[DEqMS] Aligning pepcount with fit2 rownames")
    #msg("[DEqMS] DEBUG: fit2 rownames (first 10): ", paste(head(rownames(fit2), 10), collapse=", "))
    #msg("[DEqMS] DEBUG: pepcount names (first 10): ", paste(head(names(dataSet$pepcount), 10), collapse=", "))
    psm.count <- dataSet$pepcount[rownames(fit2)]
    na_before <- sum(is.na(psm.count))
    #msg("[DEqMS] DEBUG: Aligned psm.count length=", length(psm.count))
    #msg("[DEqMS] DEBUG: NAs in aligned psm.count: ", na_before)

    # If the data were re-annotated (e.g., UniProt -> Entrez), try to remap
    # pepcount to the annotated IDs so DEqMS sees the correct PSM counts.
    if (na_before > 0 && file.exists("annotation.qs")) {
      anot.id <- ov_qs_read("annotation.qs")
      if (!is.null(names(anot.id))) {
        lvl.opt <- if (!is.null(dataSet$lvl.opt)) dataSet$lvl.opt else if (!is.null(paramSet$lvl.opt)) paramSet$lvl.opt else "sum"
        lvl.opt <- tolower(lvl.opt)
        agg.fun <- switch(lvl.opt,
                          mean = function(x) mean(x, na.rm = TRUE),
                          median = function(x) stats::median(x, na.rm = TRUE),
                          max = function(x) max(x, na.rm = TRUE),
                          sum = function(x) sum(x, na.rm = TRUE),
                          function(x) sum(x, na.rm = TRUE))

        pep.df <- data.frame(
          orig = names(dataSet$pepcount),
          mapped = anot.id[names(dataSet$pepcount)],
          count = as.numeric(dataSet$pepcount),
          stringsAsFactors = FALSE
        )
        pep.df <- pep.df[!is.na(pep.df$mapped), , drop = FALSE]

        if (nrow(pep.df) > 0) {
          remapped <- tapply(pep.df$count, pep.df$mapped, agg.fun)
          remapped <- remapped[rownames(fit2)]
          na_after <- sum(is.na(remapped))
          #msg("[DEqMS] DEBUG: Remapped pepcount via annotation (lvl.opt=", lvl.opt,
          #        "); NAs after remap: ", na_after)
          if (na_after < na_before) {
            psm.count <- remapped
            count.col.name <- paste0("pepcount (", lvl.opt, " via annotation)")
            #msg("[DEqMS] Using pepcount remapped to annotated IDs")
          }
        }
      }
    }

    if (is.null(count.col.name)) {
      count.col.name <- "pepcount (from reader)"
      #msg("[DEqMS] Using pepcount from data reader")
    }
  }

  # Check various possible column names for PSM/peptide counts
  count.col.candidates <- c("PSM.count", "PSMs", "Peptides", "Unique.Peptides",
                            "PeptideCount", "SpectrumCount", "MS.MS.Count",
                            "count", "psm_count", "peptide_count")

  # Try to get from data annotations first
  if (is.null(psm.count)) {
    data.annotated <- tryCatch({
      readDataQs("data.annotated.qs", paramSet$anal.type, dataSet$name)
    }, error = function(e) NULL)

    if (!is.null(data.annotated)) {
      # Check if there's a count attribute
      if (!is.null(attr(data.annotated, "psm.count"))) {
        psm.count <- attr(data.annotated, "psm.count")[rownames(fit2)]
        count.col.name <- "PSM.count"
      } else {
        # Try to find count column in annotations
        for (col.name in count.col.candidates) {
          if (col.name %in% colnames(data.annotated)) {
            psm.count <- data.annotated[[col.name]][rownames(fit2)]
            count.col.name <- col.name
            break
          }
        }
      }
    }
  }

  # If not found in annotations, try original data structure
  if (is.null(psm.count) && !is.null(dataSet$data_orig_full)) {
    for (col.name in count.col.candidates) {
      if (col.name %in% colnames(dataSet$data_orig_full)) {
        psm.count <- dataSet$data_orig_full[[col.name]][rownames(fit2)]
        count.col.name <- col.name
        break
      }
    }
  }

  # Fallback: use minimum detections (number of non-NA values) as proxy
  if (is.null(psm.count)) {
    #msg("[DEqMS] PSM/peptide count not found. Using minimum detections as proxy.")
    psm.count <- apply(data.norm, 1, function(x) sum(!is.na(x) & x > 0))
    count.col.name <- "min.count"
  }

  # Ensure PSM count matches fit2 dimensions and is numeric
  psm.names <- names(psm.count)
  psm.count <- as.numeric(psm.count)
  # Force names to match the model matrix so downstream tables stay aligned
  names(psm.count) <- rownames(fit2)
  psm.count[is.na(psm.count)] <- 1  # Replace NA with 1 to avoid issues
  psm.count[psm.count < 1] <- 1     # Minimum count of 1

  #msg("[DEqMS] Using ", count.col.name, " for variance modeling (range: ",
  #        min(psm.count), " - ", max(psm.count), ")")

  # Check if PSM counts have sufficient variation for DEqMS
  count.range <- max(psm.count) - min(psm.count)
  count.unique <- length(unique(psm.count))

  if (count.range == 0 || count.unique < 3) {
    msg <- paste0("[DEqMS] WARNING: Insufficient variation in PSM counts (range=", count.range,
                  ", unique=", count.unique, "). DEqMS requires variable PSM counts to model variance properly. ",
                  "This typically happens when: (1) using DIA data without PSM counts, (2) all proteins detected in all samples, ",
                  "or (3) using min.count as proxy. Falling back to standard limma.")
    #msg(msg)
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")

    # Fall back to standard limma results (already computed in fit2)
    result.list <- list()
    contrast.names <- colnames(contrast.matrix)

    for (i in seq_along(contrast.names)) {
      nm <- contrast.names[i]
      tbl <- topTable(fit2, coef = i, number = Inf, adjust.method = "BH", sort.by = "none")

      # Ensure all required columns are present and properly formatted
      if ("P.Value" %in% colnames(tbl)) {
        tbl$P.Value <- as.numeric(as.character(tbl$P.Value))
      }
      if ("adj.P.Val" %in% colnames(tbl)) {
        tbl$adj.P.Val <- as.numeric(as.character(tbl$adj.P.Val))
      }
      if ("logFC" %in% colnames(tbl)) {
        tbl$logFC <- as.numeric(as.character(tbl$logFC))
      }

      # Keep row order/names synchronized with the model matrix
      rownames(tbl) <- rownames(fit2)
      result.list[[nm]] <- tbl
    }

    dataSet$comp.res.list <- result.list
    dataSet$comp.res <- result.list[[1]]
    return(dataSet)
  }

  ## ------------------------------------------------------------------ ##
  ## 3 · Apply DEqMS variance adjustment                                ##
  ## ------------------------------------------------------------------ ##
  fit2$count <- psm.count

  fit2 <- tryCatch({
    spectraCounteBayes(fit2)
  }, error = function(e) {
    msg <- paste0("[DEqMS] spectraCounteBayes failed: ", e$message,
                  ". This often happens when PSM counts lack sufficient variation. Falling back to standard limma.")
    #msg(msg)
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(NULL)
  })

  # If DEqMS failed, fall back to limma
  if (is.null(fit2)) {
    # Re-run standard eBayes (fit2 was overwritten above)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, trend = robustTrend, robust = robustTrend)

    result.list <- list()
    contrast.names <- colnames(contrast.matrix)

    for (i in seq_along(contrast.names)) {
      nm <- contrast.names[i]
      tbl <- topTable(fit2, coef = i, number = Inf, adjust.method = "BH", sort.by = "none")

      if ("P.Value" %in% colnames(tbl)) {
        tbl$P.Value <- as.numeric(as.character(tbl$P.Value))
      }
      if ("adj.P.Val" %in% colnames(tbl)) {
        tbl$adj.P.Val <- as.numeric(as.character(tbl$adj.P.Val))
      }
      if ("logFC" %in% colnames(tbl)) {
        tbl$logFC <- as.numeric(as.character(tbl$logFC))
      }

      result.list[[nm]] <- tbl
    }

    dataSet$comp.res.list <- result.list
    dataSet$comp.res <- result.list[[1]]
    return(dataSet)
  }

  ## ------------------------------------------------------------------ ##
  ## 4 · Extract adjusted results                                       ##
  ## ------------------------------------------------------------------ ##
  result.list <- list()
  contrast.names <- colnames(contrast.matrix)

  for (i in seq_along(contrast.names)) {
    nm <- contrast.names[i]

    # 1. Get the standard Limma table first (for structure, logFC, AveExpr, etc.)
    tbl <- topTable(fit2, coef = i, number = Inf, adjust.method = "BH", sort.by = "none")

    # 2. CRITICAL FIX: Overwrite t, P.Value and adj.P.Val with DEqMS adjusted values
    # DEqMS stores results in $sca.t and $sca.p matrices
    if (!is.null(fit2$sca.t) && !is.null(fit2$sca.p)) {
        # Extract specific column for this contrast
        # Note: sca.t/p are matrices aligned with the contrasts
        tbl$t <- fit2$sca.t[, i] 
        tbl$P.Value <- fit2$sca.p[, i]
        
        # Recalculate FDR based on the new P-values
        tbl$adj.P.Val <- p.adjust(tbl$P.Value, method = "BH")
    }

    # Ensure all required columns are present and properly formatted
    required_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

    # Check for missing required columns
    missing_cols <- setdiff(required_cols, colnames(tbl))
    if (length(missing_cols) > 0) {
      #msg("[DEqMS] Warning: Missing columns: ", paste(missing_cols, collapse = ", "))
    }

    # Ensure numeric types
    if ("P.Value" %in% colnames(tbl)) {
      tbl$P.Value <- as.numeric(as.character(tbl$P.Value))
    }
    if ("adj.P.Val" %in% colnames(tbl)) {
      tbl$adj.P.Val <- as.numeric(as.character(tbl$adj.P.Val))
    }
    if ("logFC" %in% colnames(tbl)) {
      tbl$logFC <- as.numeric(as.character(tbl$logFC))
    }

    # Add count information to results
    tbl[[count.col.name]] <- psm.count[rownames(tbl)]

    # Keep row order/names synchronized with the model matrix
    rownames(tbl) <- rownames(fit2)

    # Order tbl by P.Value (This will now sort by the DEqMS P-value)
    if ("P.Value" %in% colnames(tbl)) {
      tbl <- tbl[order(tbl$P.Value), ]
    }
    
    result.list[[nm]] <- tbl
  }

  dataSet$comp.res.list <- result.list
  dataSet$comp.res <- result.list[[1]]
  print(head(dataSet$comp.res));

  msgSet$current.msg <- c(msgSet$current.msg,
                          paste0("[DEqMS] Variance adjusted using ", count.col.name,
                                 ". Mean-variance trend fitted successfully."))
  saveSet(msgSet, "msgSet")

  return(dataSet)
}

.perform_williams_trend <- function(dataSet,
                                    robustTrend = FALSE,
                                    verbose     = TRUE)
{

  expr <- if (length(dataSet$rmidx) > 0)
    dataSet$data.norm[, -dataSet$rmidx, drop = FALSE] else
      dataSet$data.norm

  cls_vals <- if (length(dataSet$rmidx) > 0)
    dataSet$cls[-dataSet$rmidx] else
      dataSet$cls

  unique_cls <- unique(cls_vals)
  dose_order <- suppressWarnings(as.numeric(as.character(unique_cls)))

  if (all(!is.na(dose_order))) {
    ord_levels <- unique_cls[order(dose_order)]
  } else {
    ord_levels <- sort(unique_cls)
  }

  grp <- factor(cls_vals, levels = ord_levels, ordered = TRUE)
  grp <- droplevels(grp)

  if (nlevels(grp) < 3L)
    stop("Williams trend test requires ≥ 3 ordered doses.")

  grp_levels <- levels(grp)
  grp_index  <- lapply(grp_levels, function(lv) which(grp == lv))
  n_per_grp  <- vapply(grp_index, length, integer(1))

  total_n  <- sum(n_per_grp)
  df_resid <- total_n - length(grp_levels)
  if (df_resid <= 0)
    stop("Williams trend test requires replication (df <= 0).")

  gene_ids   <- rownames(expr)
  gene_count <- nrow(expr)

  group_means <- matrix(NA_real_, nrow = gene_count,
                        ncol = length(grp_levels),
                        dimnames = list(gene_ids, grp_levels))
  SSE <- numeric(gene_count)

  for (i in seq_along(grp_levels)) {
    idx <- grp_index[[i]]
    mat <- expr[, idx, drop = FALSE]
    gm  <- rowMeans(mat)
    group_means[, i] <- gm

    if (length(idx) > 1L) {
      centered <- sweep(mat, 1, gm, "-", check.margin = FALSE)
      SSE <- SSE + rowSums(centered^2)
    }
  }

  MSE <- SSE / df_resid

  control_mean <- group_means[, 1, drop = TRUE]
  dose_means   <- group_means[, -1, drop = FALSE]
  pair_logFC   <- sweep(dose_means, 1, control_mean, "-")
  colnames(pair_logFC) <- paste0("Dose_", grp_levels[1],
                                 ".Dose_", grp_levels[-1])

  ave_expr <- rowMeans(expr)

  combo_idx <- lapply(seq_len(length(grp_levels) - 1L),
                      function(k) seq(from = length(grp_levels) - k + 1L,
                                      to   = length(grp_levels)))
  combo_n  <- vapply(combo_idx, function(idx) sum(n_per_grp[idx]), numeric(1))

  diff_mat <- matrix(NA_real_, nrow = gene_count,
                     ncol = length(combo_idx))
  for (j in seq_along(combo_idx)) {
    idx      <- combo_idx[[j]]
    weights  <- n_per_grp[idx] / sum(n_per_grp[idx])
    submeans <- group_means[, idx, drop = FALSE]
    weighted <- sweep(submeans, 2, weights, `*`)
    diff_mat[, j] <- rowSums(weighted) - control_mean
  }

  sqrt_terms <- sqrt(1 / n_per_grp[1] + 1 / combo_n)
  denom_mat  <- outer(sqrt(MSE), sqrt_terms, `*`)
  t_mat      <- diff_mat / denom_mat

  zero_mask <- denom_mat == 0
  if (any(zero_mask, na.rm = TRUE)) {
    t_mat[zero_mask & diff_mat > 0]  <- Inf
    t_mat[zero_mask & diff_mat < 0]  <- -Inf
    t_mat[zero_mask & diff_mat == 0] <- 0
  }

  t_pos <- apply(t_mat, 1, max, na.rm = TRUE)
  t_neg <- apply(t_mat, 1, min, na.rm = TRUE)
  t_pos[is.infinite(t_pos) & t_pos < 0] <- NA_real_
  t_neg[is.infinite(t_neg) & t_neg > 0] <- NA_real_

  both_na <- is.na(t_pos) & is.na(t_neg)
  abs_pos <- abs(t_pos); abs_pos[is.na(abs_pos)] <- -Inf
  abs_neg <- abs(t_neg); abs_neg[is.na(abs_neg)] <- -Inf
  use_pos <- abs_pos >= abs_neg
  t_stat  <- ifelse(use_pos, t_pos, t_neg)
  t_stat[both_na] <- NA_real_

  P.Value   <- 2 * pt(-abs(t_stat), df = df_resid)
  adj.P.Val <- p.adjust(P.Value, "fdr")

  topFeatures <- data.frame(
    pair_logFC,
    AveExpr   = ave_expr,
    t         = t_stat,
    P.Value   = P.Value,
    adj.P.Val = adj.P.Val,
    check.names = FALSE)

  ord <- order(topFeatures$adj.P.Val, na.last = TRUE)
  topFeatures <- topFeatures[ord, , drop = FALSE]

  dataSet$comp.res <- topFeatures
  dataSet$comp.res.list <- make_comp_res_list(
    topFeatures,
    stat.cols = c("AveExpr", "t", "P.Value", "adj.P.Val"))
  dataSet$de.method <- "wtt"

  if (verbose) {
    cat("Williams trend test (BMDExpress-style):",
        nrow(topFeatures), "features processed;",
        "df =", df_resid, "\n")
  }

  return(dataSet)
}



SetupDesignMatrix<-function(dataName="", deMethod){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  cls <- dataSet$cls; 
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls);
  dataSet$design <- design;
  dataSet$de.method <- deMethod;
  dataSet$pval <- 0.05;
  dataSet$fc.val <- 1;

  saveSet(paramSet, "paramSet");
  return(RegisterData(dataSet));
}


# perform limma on given two groups selected 
# used by integarative analysis

#'Perform differential analysis using Limma method (meta-analysis)
#'@description Detect number of DE features in individual data set for quality checking purposes 
#'@param dataName File name of data set.
#'@param grps Selected two groups for comparison
#'@param p.lvl P-value threshold
#'@param fc.lvl Fold-change threshold
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformLimmaDE<-function(dataName="", grps, p.lvl, fc.lvl=NULL){
  #save.image("limma.RDAta");
  dataSet <- readDataset(dataName);
  dataSet$pval <- p.lvl;
    dataSet$fc.lvl <- 0;

  if(length(levels(dataSet$cls))>2){ 
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
  }else{
    sel.inx <- rep(T, ncol(dataSet$data.norm));
  }
  
  group <- factor(dataSet$cls[sel.inx]); # note regenerate factor to drop levels 
  data <- dataSet$data.norm[, sel.inx];
  
  res.limma <- PerformLimma(data, group);
  res.all <- GetLimmaResTable(res.limma$fit.obj);
  
  if(!is.null(fc.lvl)){
    hit.inx <- abs(res.all$logFC)>= fc.lvl & res.all$adj.P.Val <= p.lvl
    dataSet$fc.lvl <- fc.lvl;

  }else{
    hit.inx <- res.all$adj.P.Val <= p.lvl
    dataSet$fc.lvl <- 0;

  }
  if(sum(hit.inx) == 0){
    return (c(1, 0, nrow(res.all)));
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx);
  res <- res.all[hit.inx, , drop=F];
  
  # rm .txt suffix for new names
  shortNm <- substring(dataName, 0, nchar(dataName)-4);
  fast.write(signif(res[,-1],5), file=paste("Sigfeatures_", shortNm, ".csv",sep=""));
  
  sig.count <- nrow(res);
  de.features <- rownames(res);
  non.sig.count <- nrow(data)-sig.count;
  rm(res.all);
  
  gc();
  
  # record the sig gene vec
  output <- c(1, sig.count, non.sig.count);
  return(RegisterData(dataSet, output));
}

# perfor differential analysis for array/RNA seq data
# for two groups only (used for meta-analysis)
PerformLimma<-function(data, group){
  #print(identical(colnames(data), group))
  #print(colnames(data))
  #print(group)

  #print("limma");
  require(limma);
  design <- model.matrix(~-1 + group);
  fit <- lmFit(data, design);
  
  grps.cmp <- paste("group", levels(group)[2], " - ", "group", levels(group)[1], sep="");
  myargs <- list(grps.cmp, levels = design);
  contrast.matrix <- do.call(makeContrasts, myargs);
  fit <- contrasts.fit(fit, contrast.matrix);
  fit <- eBayes(fit);
  gc();
  return (list(fit.obj=fit));
}

# get result table from eBayes fit object
GetLimmaResTable<-function(fit.obj){
  require(limma);
  resTable <- topTable(fit.obj, number=Inf, adjust.method="BH");
  if(!is.null(resTable$ID)){ # for older version
    rownames(resTable) <- resTable$ID;
    resTable$ID <- NULL;
  }
  return (resTable);
}

MultiCovariateRegression <- function(fileName,
                                     analysis.var, # metadata variable name
                                     ref = NULL, # reference class from analysis.var metadata (only if categorical)
                                     contrast = "anova",  # comparison class from analysis.var (only if categorical)
                                     blocking.factor = NULL, 
                                     robustTrend = F, 
                                     internal=F,
                                     useBatchCorrected=TRUE){ # whether returns 0/1 or dataset object

  dataSet <- readDataset(fileName);
  if(!exists('adj.vec')){ ## covariate to adjust for
    adj.vec <- "NA";
  }else{
    if(length(adj.vec) > 0){

    }else{
      adj.vec <- "NA"
    }
  }
  interim <- .multiCovariateRegression(dataSet, analysis.var, ref, contrast, blocking.factor, adj.vec, robustTrend, F, useBatchCorrected)
  if(is.list(interim)){
    res <- 1;
  }else{
    res <- interim;
  }
  return(res);
}

.multiCovariateRegression <- function(dataSet,
                                     analysis.var, # metadata variable name
                                     ref = NULL, # reference class from analysis.var metadata (only if categorical)
                                     contrast = "anova",  # comparison class from analysis.var (only if categorical)
                                     # fixed.effects = NULL,
                                     blocking.factor = NULL,
                                     adj.factors="NA",# metadata variables to adjust for
                                     robustTrend = F,
                                     internal=F, # whether returns 0/1 or dataset object, T for metaanal covariate
                                     useBatchCorrected=TRUE){
  cat("DEBUG [.multiCovariateRegression]: ENTRY - class(dataSet) =", class(dataSet), "\n");
  cat("DEBUG [.multiCovariateRegression]: analysis.var =", analysis.var, "\n");
  cat("DEBUG [.multiCovariateRegression]: ref =", ref, "\n");
  cat("DEBUG [.multiCovariateRegression]: is.null(ref) =", is.null(ref), "\n");
  cat("DEBUG [.multiCovariateRegression]: is.na(ref) =", is.na(ref), "\n");
  cat("DEBUG [.multiCovariateRegression]: contrast =", contrast, "\n");
  cat("DEBUG [.multiCovariateRegression]: useBatchCorrected =", useBatchCorrected, "\n");
  flush.console();

  # load libraries
  library(limma);
  library(dplyr);
  
  # need a line for read dataSet
  msgSet <- readSet(msgSet, "msgSet");
  dataSet$rmidx <- NULL;

  # for embedded inside tools (ProteoAnalyst etc)
  paramSet <- readSet(paramSet, "paramSet");
  useMeta <- !is.null(paramSet$performedBatch) && isTRUE(paramSet$performedBatch);
  if(internal){
    if(useMeta){
    inmex.meta <- ov_qs_read("inmex_meta.qs");
    }else{

    inmex.meta <- ov_qs_read("inmex.meta.orig.qs");

    }
    feature_table <- inmex.meta$data[, colnames(inmex.meta$data) %in% colnames(dataSet$data.norm)];
  }else{
    # Use covariate-adjusted data if available and useBatchCorrected=TRUE
    if(useBatchCorrected && !is.null(dataSet$data.adjusted)) {
      feature_table <- dataSet$data.adjusted;
    } else {
      feature_table <- dataSet$data.norm;
    }
  }
  covariates <- dataSet$meta.info;

  matched_indices <- match(colnames(feature_table), rownames(covariates));
  covariates <- covariates[matched_indices, ,drop=F ];
  dataSet$meta.info <- covariates;
  #fixed.effects <- adj.vec
  # process covariates
  var.types <- lapply(covariates, class) %>% unlist();
  covariates[,c(var.types == "character")] <- lapply(covariates[,c(var.types == "character")], factor);

  # aggregate vars
  all.vars <- c(analysis.var);
  vars <- c(analysis.var);

  # If using covariate-adjusted data, don't include adjustment factors in the model again
  # The effects have already been removed via removeBatchEffect
  if(useBatchCorrected && !is.null(dataSet$covariate.adjusted) && dataSet$covariate.adjusted) {
    # Using adjusted data - only include analysis variable, NOT the adjustment covariates
    if(all(adj.factors != "NA")){
      # Log that we're skipping the adjustment factors since data is already adjusted
      cat("Note: Using covariate-adjusted data. Adjustment factors already removed.\n");
    }
  } else {
    # Not using adjusted data - include adjustment factors in the model
    if(all(adj.factors != "NA")){
      vars = c(vars, adj.factors);
    }
  }
  
  if(!is.null(blocking.factor) && !is.na(blocking.factor) && blocking.factor!="NA" && blocking.factor!="" ){
    all.vars = c(all.vars, blocking.factor);
  }

  all.vars<- unique(all.vars);
    
  covariates <- covariates[,unique(c(vars, all.vars)),drop=F];
  rmidx <-which(apply(covariates, 1, function(x) "NA" %in% x))
  
  if(length(rmidx)>0){
    covariates <- covariates[-rmidx,,drop=F];
    dataSet$rmidx <- rmidx;
  }
  feature_table <- feature_table[,colnames(feature_table) %in% rownames(covariates)];
  
  if(!identical(colnames(feature_table), rownames(covariates))){
    msgSet$current.msg <- "Error - order of samples got mixed up between tables";
    saveSet(msgSet, "msgSet");
    cat("DEBUG [.multiCovariateRegression]: RETURNING 0 - sample order mismatch\n");
    return(0)
  }
  
  
  # get analysis type
  analysis.type = ifelse(dataSet$disc.inx[analysis.var],"disc","cont")
  if(is.na(analysis.type)){
    msgSet$current.msg <- "Analysis var not found in our database!";
    saveSet(msgSet, "msgSet");
    cat("DEBUG [.multiCovariateRegression]: RETURNING 0 - analysis var not found\n");
    return(0)
  }
  
  if(analysis.type == "disc"){
    cat("DEBUG: Discrete analysis type\n");
    cat("DEBUG: analysis.var =", analysis.var, "\n");
    cat("DEBUG: ref =", ref, "\n");
    cat("DEBUG: contrast =", contrast, "\n");
    cat("DEBUG: vars =", paste(vars, collapse=", "), "\n");
    flush.console();

    # build design and contrast matrix
    #  covariates[, analysis.var] <- covariates[, analysis.var] %>% make.names() %>% factor();
    #str(covariates)
    grp.nms <- levels(covariates[, analysis.var]);
    cat("DEBUG: Original grp.nms from levels:", paste(grp.nms, collapse=", "), "\n");
    cat("DEBUG: Number of samples per group:\n");
    print(table(covariates[, analysis.var]));
    flush.console();

    if(any(grepl("(^[0-9]+).*", grp.nms))){
      grp.nms <- paste0(analysis.var,"_",grp.nms);
      if(!(is.null(ref))& ref!="NA"){
        ref <- paste0(analysis.var,"_", ref)
      }
      if(contrast!="anova"){
        contrast <- paste0(analysis.var,"_", contrast)
      }
      cat("DEBUG: Modified grp.nms (starts with number):", paste(grp.nms, collapse=", "), "\n");
      cat("DEBUG: Modified ref:", ref, "\n");
      cat("DEBUG: Modified contrast:", contrast, "\n");
    }
    flush.console();

    for(col in 1:ncol(covariates)){
      if(dataSet$cont.inx[colnames(covariates)[col]]){
        covariates[,col] <- as.numeric(as.character(covariates[,col]))
      }
    }

    cat("DEBUG: Building design matrix with formula: ~ 0", paste0(" + ", vars, collapse = ""), "\n");
    flush.console();
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    cat("DEBUG: Design matrix created, dimensions:", nrow(design), "x", ncol(design), "\n");
    cat("DEBUG: Original design colnames:", paste(colnames(design), collapse=", "), "\n");
    flush.console();

    colnames(design)[1:length(grp.nms)] <- grp.nms;
    cat("DEBUG: After renaming first", length(grp.nms), "columns to grp.nms\n");
    cat("DEBUG: New design colnames:", paste(colnames(design), collapse=", "), "\n");
    flush.console();

    # Check for NA values in design matrix
    if(any(is.na(design))) {
      na.cols <- which(colSums(is.na(design)) > 0);
      msgSet$current.msg <- paste0("Error: Design matrix contains NA values in columns: ",
                                   paste(colnames(design)[na.cols], collapse=", "),
                                   ". Please check your covariate data.");
      saveSet(msgSet, "msgSet");
      cat("DEBUG [.multiCovariateRegression]: RETURNING 0 - design matrix has NA values\n");
      return(0)
    }

    cat("DEBUG: Building contrast matrix\n");
    cat("DEBUG: contrast =", contrast, "\n");
    cat("DEBUG: ref =", ref, "\n");
    cat("DEBUG: is.null(ref):", is.null(ref), "\n");
    cat("DEBUG: is.na(ref):", is.na(ref), "\n");
    cat("DEBUG: grp.nms =", paste(grp.nms, collapse=", "), "\n");
    cat("DEBUG: design matrix dimensions:", nrow(design), "x", ncol(design), "\n");
    cat("DEBUG: design colnames:", paste(colnames(design), collapse=", "), "\n");
    flush.console();

    # Handle missing or invalid ref value
    if(is.null(ref) || length(ref) == 0 || is.na(ref) || ref == "NA" || ref == ""){
      cat("DEBUG: ref is NULL/NA/empty, using first group as reference\n");
      ref <- grp.nms[1];
      cat("DEBUG: Set ref to:", ref, "\n");
      flush.console();
    }

    myargs <- list();
    if(contrast == "anova"){
      contrasts <- grp.nms[grp.nms != ref];
      myargs <- as.list(paste(contrasts, "-", ref, sep = ""));
      cat("DEBUG: ANOVA mode - contrasts:", paste(unlist(myargs), collapse=", "), "\n");
    } else {
      myargs <- as.list(paste(contrast, "-", ref, sep = ""));
      cat("DEBUG: T-test mode - contrast:", paste(unlist(myargs), collapse=", "), "\n");
    }
    flush.console();

    myargs[["levels"]] <- design;
    cat("DEBUG: Calling makeContrasts with myargs:\n");
    print(str(myargs, max.level = 1));
    flush.console();

    contrast.matrix <- do.call(makeContrasts, myargs);

    cat("DEBUG: contrast.matrix created, dimensions:", nrow(contrast.matrix), "x", ncol(contrast.matrix), "\n");
    cat("DEBUG: contrast.matrix:\n");
    print(contrast.matrix);
    cat("DEBUG: any NA in contrast.matrix:", any(is.na(contrast.matrix)), "\n");
    flush.console();

    # Check for NA values in contrast matrix
    if(any(is.na(contrast.matrix))) {
      msgSet$current.msg <- "Error: Contrast matrix contains NA values. This may occur when groups have insufficient samples after covariate adjustment.";
      saveSet(msgSet, "msgSet");
      cat("DEBUG [.multiCovariateRegression]: RETURNING 0 - contrast matrix has NA values\n");
      flush.console();
      return(0)
    }

    # handle blocking factor
    if (is.null(blocking.factor) | is.na(blocking.factor) | blocking.factor=="NA" | blocking.factor == "") {
      fit <- lmFit(feature_table, design);
    } else {
      block.vec <- covariates[,blocking.factor];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec);
      fit <- lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus);
    }

    # Check for NA values in fit coefficients before contrasts.fit
    if(any(is.na(fit$coefficients))) {
      na.features <- rowSums(is.na(fit$coefficients));
      msgSet$current.msg <- paste0("Error: Linear model fit contains ", sum(na.features > 0),
                                   " features with NA coefficients. This indicates numerical issues with the model.");
      saveSet(msgSet, "msgSet");
      cat("DEBUG [.multiCovariateRegression]: RETURNING 0 - fit coefficients have NA values\n");
      return(0)
    }

    # get results
    fit <- contrasts.fit(fit, contrast.matrix);
    fit <- eBayes(fit, trend=robustTrend, robust=robustTrend);
    rest <- topTable(fit, number = Inf);
    
    if(contrast != "anova"){    
      colnames(rest)[1] <- myargs[[1]];
      grp.nms<-c(ref,contrast)
      
    }
    #for meta-anal or when using covariate-adjusted data
    if(internal || (useBatchCorrected && !is.null(dataSet$covariate.adjusted) && dataSet$covariate.adjusted)){
        ### get results with no adjustment using original unadjusted data
        design.noadj <- model.matrix(formula(paste0("~ 0", paste0(" + ", analysis.var, collapse = ""))), data = covariates);
        colnames(design.noadj)[1:length(grp.nms)] <- grp.nms;
        myargs.noadj <- myargs;
        myargs.noadj[["levels"]] <- design.noadj;
        contrast.matrix.noadj <- do.call(makeContrasts, myargs.noadj);

        # Use original unadjusted data for comparison
        feature_table.noadj <- dataSet$data.norm;

        fit.noadj <- lmFit(feature_table.noadj, design.noadj)
        fit.noadj <- contrasts.fit(fit.noadj, contrast.matrix.noadj);
        fit.noadj <- eBayes(fit.noadj, trend=robustTrend, robust=robustTrend);
        res.noadj <- topTable(fit.noadj, number = Inf);
        dataSet$res.noadj <- res.noadj;
    }

    dataSet$contrast.matrix <- contrast.matrix;
    dataSet$par1 <-  myargs[[1]];
    dataSet$grp.nms <- ifelse(any(grepl("(^[0-9]+).*", grp.nms)), paste0(analysis.var,"_",grp.nms),grp.nms);
  } else { 
    
    # build design matrix
    types <- dataSet$cont.inx[vars];
    contIdx <- as.numeric(which(types))
    covariates[,contIdx] <- unlist(lapply(covariates[,contIdx], function(x) as.numeric(as.character((x)))));
    
    if (all(types)) {
      design <- model.matrix(formula(paste0("~", paste0(" + ", vars, collapse = ""))), data = covariates);
    } else {
      design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    }
    
    # handle blocking factor
    if (is.null(blocking.factor) | is.na(blocking.factor) | blocking.factor=="NA") {
      fit <- lmFit(feature_table, design);
    } else {
      block.vec <- covariates[, blocking.factor];
      corfit <- duplicateCorrelation(feature_table, design, block = block.vec);
      fit <- lmFit(feature_table, design, block = block.vec, correlation = corfit$consensus);
    }
    
    # get results
    fit <- eBayes(fit, trend=robustTrend, robust=robustTrend);
    rest <- topTable(fit, number = Inf, coef = analysis.var);
    colnames(rest)[1] <- dataSet$par1 <- analysis.var;
    
    ### get results with no adjustment
    if(internal || (useBatchCorrected && !is.null(dataSet$covariate.adjusted) && dataSet$covariate.adjusted)){
    # Use original unadjusted data for comparison
    feature_table.noadj <- dataSet$data.norm;

    design.noadj <- model.matrix(formula(paste0("~", analysis.var)), data = covariates);
    fit.noadj <- eBayes(lmFit(feature_table.noadj, design.noadj), trend=robustTrend, robust=robustTrend);
    res.noadj <- topTable(fit.noadj, number = Inf);
    dataSet$res.noadj <- res.noadj;
    }
  }
  
  cat("DEBUG [.multiCovariateRegression]: Setting up dataSet before return\n");
  cat("DEBUG [.multiCovariateRegression]: nrow(rest) =", nrow(rest), "\n");
  cat("DEBUG [.multiCovariateRegression]: class(rest) =", class(rest), "\n");

  dataSet$design <- design;
  dataSet$contrast.type <- analysis.type;
  dataSet$comp.res <- rest;
  cat("DEBUG [.multiCovariateRegression]: dataSet$comp.res assigned, nrow =", nrow(dataSet$comp.res), "\n");

  dataSet$comp.res.list <- make_comp_res_list(rest)
  dataSet$de.method <- "limma"
  dataSet$comp.type <- "default"
  dataSet$fit.obj <- fit;

  dataSet$pval <- 0.05;
  dataSet$fc.val <- 1;
  dataSet$analysis.var <- analysis.var;
  dataSet$de.adj <- adj.factors;

  if(all(adj.factors != "NA")){
    dataSet$de.adj <- "NA"
  }else{
    dataSet$de.adj <- adj.factors;
  }

  cat("DEBUG [.multiCovariateRegression]: About to RegisterData\n");
  RegisterData(dataSet);
  cat("DEBUG [.multiCovariateRegression]: RegisterData complete, returning dataSet\n");
  cat("DEBUG [.multiCovariateRegression]: class(dataSet) =", class(dataSet), "\n");
  return(dataSet);
}
make_comp_res_list <- function(resTab,
                               stat.cols = c("AveExpr", "F", "t",
                                             "P.Value", "adj.P.Val",
                                             "B", "FDR"))
{
  ## detect logFC columns automatically -------------------------------
  lfc.cand <- setdiff(colnames(resTab), stat.cols)

  ## strip common suffix/prefix for the test
  strip_logfc <- function(x)
      sub("\\.logFC$","",
      sub("^logFC\\.","", x, ignore.case = TRUE), ignore.case = TRUE)

  uniq.core <- unique(strip_logfc(lfc.cand))

  if (length(uniq.core) == 0L)
      stop("Could not detect any log-fold-change columns automatically. ",
           "Pass lfc.cols explicitly.")

  ## build list --------------------------------------------------------
  out <- lapply(uniq.core, function(core) {

            ## match any of the naming patterns for this contrast
            pat <- paste0("^(", core,
                          "|logFC\\.", core,
                          "|", core, "\\.logFC)$")
            this.lfc <- grep(pat, colnames(resTab), value = TRUE)

            if (length(this.lfc) == 0L)
                stop("Unexpected: no column matched for contrast ", core)

            df <- resTab[ , c(this.lfc[1], stat.cols[stat.cols %in% colnames(resTab)]),
                          drop = FALSE]
            colnames(df)[1] <- "logFC"
            df
         })
  names(out) <- uniq.core
  out
}

parse_contrast_groups <- function(contrast_str) {
  comps <- strsplit(contrast_str, " vs\\.?\\s*")[[1]]
  if (length(comps) != 2) stop(paste("Invalid contrast format:", contrast_str))
  return(comps)
}

.get.interaction.results <- function(dds.path = "deseq.res.obj.rds") {
  dds <- ov_qs_read(dds.path)
  cat("Available result names:\n")
  
  # Automatically detect the interaction term
  interaction_name <- grep("factorA.*factorB.*", resultsNames(dds), value = TRUE)
  if (length(interaction_name) == 0) {
    stop("No interaction term found in model.")
  }

  cat("Extracting interaction term:", interaction_name, "\n")
  res <- results(dds, name = interaction_name[1], independentFiltering = FALSE, cooksCutoff = Inf)

  # Format results
  topFeatures <- data.frame(res@listData)
  rownames(topFeatures) <- rownames(res)
  colnames(topFeatures) <- sub("padj", "adj.P.Val", colnames(topFeatures))
  colnames(topFeatures) <- sub("pvalue", "P.Value", colnames(topFeatures))
  colnames(topFeatures) <- sub("log2FoldChange", "logFC", colnames(topFeatures))
  topFeatures <- topFeatures[c("logFC", "baseMean", "lfcSE", "stat", "P.Value", "adj.P.Val")]
  topFeatures <- topFeatures[order(topFeatures$P.Value), ]

  return(topFeatures)
}

prepareEdgeRContrast <- function(dataSet,
                                 anal.type  = "reference",
                                 par1       = NULL,
                                 par2       = NULL,
                                 nested.opt = "intonly") {
  msgSet <- readSet(msgSet, "msgSet")
  set.seed(1337)
  require(limma)

  cls_raw      <- factor(dataSet$cls)
  raw_levels   <- levels(cls_raw)
  syn_levels   <- make.names(raw_levels)

  cls_syn      <- factor(cls_raw, levels = raw_levels, labels = syn_levels)

  dataSet$cls_raw   <- cls_raw
  dataSet$cls       <- cls_syn               # <- modeling uses sanitized
  dataSet$grp.nms   <- syn_levels
  dataSet$grp.nms_raw <- raw_levels
  dataSet$comp.type <- anal.type
  dataSet$par1      <- par1

  design <- model.matrix(~ 0 + cls_syn)
  colnames(design) <- syn_levels

  # helper: map a UI label (raw or syn) into sanitized
  to_syn <- function(x) {
    if (is.null(x)) return(NULL)
    x <- trimws(as.character(x))
    if (!nzchar(x)) return(NULL)
    if (x %in% raw_levels) return(make.names(x))
    x
  }

  # parse "A vs. B" (accept raw or syn on either side)
  parse_vs <- function(x) {
    if (is.null(x)) return(c(NA_character_, NA_character_))
    parts <- strsplit(x, "\\s*vs\\.?\\s*", perl = TRUE)[[1]]
    if (length(parts) != 2) return(c(NA_character_, NA_character_))
    c(to_syn(parts[1]), to_syn(parts[2]))
  }

  ## build contrasts using sanitized tokens
  if (anal.type == "reference") {
    ref_syn <- to_syn(par1)
    if (is.null(ref_syn) || !(ref_syn %in% syn_levels)) {
      stop("`par1` must specify a valid reference level. You gave '",
           if (is.null(par1)) "NULL" else par1,
           "'. Valid (raw): ", paste(raw_levels, collapse = ", "))
    }
    others <- setdiff(syn_levels, ref_syn)
    conts  <- setNames(lapply(others, \(g) paste0(g, " - ", ref_syn)),
                       paste0(others, "_vs_", ref_syn))

  } else if (anal.type == "default") {
    combs <- combn(syn_levels, 2, simplify = FALSE)
    conts <- setNames(lapply(combs, \(x) paste0(x[1], " - ", x[2])),
                      sapply(combs, \(x) paste0(x[1], "_vs_", x[2])))

  } else if (anal.type == "time") {
    tm <- syn_levels
    conts <- setNames(lapply(seq_len(length(tm) - 1L),
                             \(i) paste0(tm[i + 1L], " - ", tm[i])),
                      paste0(tm[-1], "_vs_", tm[-length(tm)]))

  } else if (anal.type == "custom") {
    grp <- parse_vs(par1)
    if (anyNA(grp) || !all(grp %in% syn_levels)) {
      stop("`par1` must be 'A vs. B'. Valid (raw): ",
           paste(raw_levels, collapse = ", "))
    }
    conts <- setNames(list(paste0(grp[2], " - ", grp[1])),
                      paste0(grp[2], "_vs_", grp[1]))

  } else if (anal.type == "nested") {
    g1 <- parse_vs(par1); g2 <- parse_vs(par2)
    if (anyNA(g1) || anyNA(g2) || !all(c(g1, g2) %in% syn_levels)) {
      stop("`par1` and `par2` must be 'A vs. B'. Valid (raw): ",
           paste(raw_levels, collapse = ", "))
    }
    if (identical(nested.opt, "intonly")) {
      expr  <- paste0("(", g1[1], " - ", g1[2], ") - (", g2[1], " - ", g2[2], ")")
      nm    <- paste0(g1[1], g1[2], "_vs_", g2[1], g2[2], "_interaction")
      conts <- setNames(list(expr), nm)
    } else {
      expr1 <- paste0(g1[2], " - ", g1[1])
      expr2 <- paste0(g2[2], " - ", g2[1])
      expr3 <- paste0("(", g1[2], " - ", g1[1], ") - (", g2[2], " - ", g2[1], ")")
      conts <- c(setNames(list(expr1), paste0(g1[2], "_vs_", g1[1])),
                 setNames(list(expr2), paste0(g2[2], "_vs_", g2[1])),
                 setNames(list(expr3), paste0("int_", g1[2], g1[1], "_vs_", g2[2], g2[1])))
    }

  } else {
    stop("Unsupported `anal.type`: ", anal.type)
  }

  contrast.matrix <- do.call(makeContrasts, c(conts, list(levels = design)))

  dataSet$design          <- design
  dataSet$contrast.matrix <- contrast.matrix
  dataSet$contrast.names  <- colnames(contrast.matrix)
  dataSet$contrast.type   <- anal.type
  dataSet$filename        <- paste0("edgeR_", anal.type, "_", dataSet$de.method)

  RegisterData(dataSet)
  dataSet
}

###################################################
## Peptide-level differential analysis functions
###################################################

#' Check if peptide-level data exists
#' @description Checks if peptide-level data and results are available
#' @export
HasPeptideLevelData <- function() {
  pep_data_exists <- file.exists("peptide_level_data.qs")
  pep_de_exists <- file.exists("peptide_de_results.qs")
  result <- pep_data_exists && pep_de_exists
  return(result)
}

#' Perform peptide-level differential analysis
#' @description Performs DE analysis on peptide-level data after protein-level analysis is complete
#' @param dataName Dataset name
#' @export
PerformPeptideLevelDEAnal <- function(dataName = "") {
  dataSet <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  msgSet <- readSet(msgSet, "msgSet")

  # Check if peptide-level data exists
  if (!file.exists("peptide_level_data.qs")) {
    msgSet$current.msg <- "Peptide-level DE analysis skipped: peptide_level_data.qs not found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Check if peptide-to-protein map exists
  if (!file.exists("peptide_to_protein_map.qs")) {
    msgSet$current.msg <- "Peptide-level DE analysis skipped: peptide_to_protein_map.qs not found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Load peptide-level data
  pep.mat <- ov_qs_read("peptide_level_data.qs")
  pep.map <- ov_qs_read("peptide_to_protein_map.qs")

  # Create a temporary dataset for peptide-level analysis
  peptide.dataSet <- dataSet
  peptide.dataSet$data.norm <- pep.mat

  # Perform DE analysis based on the method used for protein-level
  if (dataSet$de.method == "limma") {
    peptide.dataSet <- .perform_limma_edger(peptide.dataSet, robustTrend = FALSE)
  } else if (dataSet$de.method == "deqms") {
    peptide.dataSet <- .perform_deqms(peptide.dataSet, robustTrend = FALSE)
  } else if (dataSet$de.method == "edger") {
    peptide.dataSet <- .perform_limma_edger(peptide.dataSet, robustTrend = FALSE)
  } else {
    msgSet$current.msg <- paste0("Peptide-level DE analysis: method '", dataSet$de.method, "' not yet supported for peptides.")
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Save peptide-level DE results (shadow save for Arrow zero-copy access)
  shadow_save(peptide.dataSet$comp.res, "peptide_de_results.qs")
  shadow_save(peptide.dataSet$sig.mat, "peptide_sig_mat.qs")

  msgSet$current.msg <- paste0("Peptide-level DE analysis completed: ", nrow(peptide.dataSet$comp.res), " peptides analyzed.")
  saveSet(msgSet, "msgSet")

  return(1)
}

#' Get protein-peptide mapping with DE results
#' @description Returns a data frame with protein and associated peptides with DE stats
#' @param dataName Dataset name
#' @param proteinID Protein ID to get peptides for
#' @export
GetProteinPeptideMapping <- function(dataName = "", proteinID = "") {
  msg("[R DEBUG] GetProteinPeptideMapping called with proteinID: ", proteinID)

  normalize_protein_id <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_character_)
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(NA_character_)
    first <- strsplit(x[1], ";", fixed = TRUE)[[1]][1]
    first <- trimws(first)
    if (grepl("\\|", first)) {
      parts <- strsplit(first, "\\|")[[1]]
      if (length(parts) >= 2 && nzchar(parts[2])) {
        first <- parts[2]
      }
    }
    first <- sub("_.*$", "", first)
    first <- sub("-\\d+$", "", first)
    first
  }

  # Load protein-level DE results from dataset
  dataSet <- readDataset(dataName)
  if (is.null(dataSet)) {
    msg("[R DEBUG] Failed to load dataset: ", dataName)
    return(NULL)
  }

  if (is.null(dataSet$comp.res)) {
    msg("[R DEBUG] dataSet$comp.res is NULL - DE analysis may not have been run")
    return(NULL)
  }

  prot.res <- dataSet$comp.res
  msg("[R DEBUG] Loaded protein DE results from dataSet with ", nrow(prot.res), " rows")

  # Load peptide-level DE results
  if (!file.exists("peptide_de_results.qs")) {
    msg("[R DEBUG] peptide_de_results.qs does not exist")
    return(NULL)
  }

  pep.res <- try(ov_qs_read("peptide_de_results.qs"), silent = TRUE)
  if (inherits(pep.res, "try-error")) {
    msg("[R DEBUG] Error reading peptide_de_results.qs")
    return(NULL)
  }
  msg("[R DEBUG] Loaded peptide_de_results.qs with ", nrow(pep.res), " rows")

  # Load peptide-to-protein map
  if (!file.exists("peptide_to_protein_map.qs")) {
    msg("[R DEBUG] peptide_to_protein_map.qs does not exist")
    return(NULL)
  }

  pep.map <- try(ov_qs_read("peptide_to_protein_map.qs"), silent = TRUE)
  if (inherits(pep.map, "try-error")) {
    msg("[R DEBUG] Error reading peptide_to_protein_map.qs")
    return(NULL)
  }
  msg("[R DEBUG] Loaded peptide_to_protein_map.qs with ", nrow(pep.map), " rows")
  msg("[R DEBUG] Peptide map column names: ", paste(colnames(pep.map), collapse=", "))
  msg("[R DEBUG] First few unique protein IDs in map: ", paste(head(unique(pep.map$Protein), 10), collapse=", "))

  # Get peptides for this protein
  original.protein.id <- proteinID
  protein.id.norm <- normalize_protein_id(proteinID)
  pep.map$Protein <- as.character(pep.map$Protein)
  pep.map$Protein.norm <- vapply(pep.map$Protein, normalize_protein_id, character(1))

  matched.protein.id <- proteinID
  peptides <- pep.map$Peptide[pep.map$Protein == matched.protein.id]
  if (length(peptides) == 0 && !is.na(protein.id.norm)) {
    peptides <- pep.map$Peptide[pep.map$Protein.norm == protein.id.norm]
    if (length(peptides) > 0) {
      match.idx <- which(pep.map$Protein.norm == protein.id.norm)[1]
      matched.protein.id <- pep.map$Protein[match.idx]
      msg("[R DEBUG] Matched via normalized protein ID: ", protein.id.norm, " -> ", matched.protein.id)
    }
  }
  msg("[R DEBUG] Found ", length(peptides), " peptides for protein ", proteinID)

  # If no exact match, try to find similar protein IDs
  if (length(peptides) == 0) {
    msg("[R DEBUG] No exact match for protein ", proteinID)
    # Check if protein ID exists with different formatting
    matching_proteins <- unique(pep.map$Protein[
      grepl(proteinID, pep.map$Protein, fixed = TRUE) |
      (!is.na(protein.id.norm) & grepl(protein.id.norm, pep.map$Protein, fixed = TRUE))
    ])
    if (length(matching_proteins) > 0) {
      msg("[R DEBUG] Found similar protein IDs: ", paste(matching_proteins, collapse=", "))
      msg("[R DEBUG] Using first match: ", matching_proteins[1])
      matched.protein.id <- matching_proteins[1]
      peptides <- pep.map$Peptide[pep.map$Protein == matched.protein.id]
      if (matched.protein.id != original.protein.id) {
        msg("[R DEBUG] Using matched ID for peptides only: ", matched.protein.id, " (input: ", original.protein.id, ")")
      }
    } else {
      msg("[R DEBUG] No similar protein IDs found containing: ", proteinID)
      return(NULL)
    }
  }

  # Extract protein DE stats
  prot.lookup.id <- original.protein.id
  if (!(prot.lookup.id %in% rownames(prot.res))) {
    rn <- rownames(prot.res)
    rn.norm <- vapply(rn, normalize_protein_id, character(1))
    if (!is.na(protein.id.norm) && protein.id.norm %in% rn.norm) {
      prot.lookup.id <- rn[match(protein.id.norm, rn.norm)]
      msg("[R DEBUG] Using normalized protein ID match in DE results: ", prot.lookup.id)
    }
    similar.prot <- rownames(prot.res)[grepl(prot.lookup.id, rownames(prot.res), fixed = TRUE)]
    if (length(similar.prot) == 0 && grepl("\\|", prot.lookup.id)) {
      core.id <- strsplit(prot.lookup.id, "\\|", fixed = FALSE)[[1]][2]
      similar.prot <- rownames(prot.res)[grepl(core.id, rownames(prot.res), fixed = TRUE)]
    }
    if (length(similar.prot) > 0) {
      prot.lookup.id <- similar.prot[1]
      msg("[R DEBUG] Using protein ID match in DE results: ", prot.lookup.id)
    }
  }

  if (prot.lookup.id %in% rownames(prot.res)) {
    prot.fc <- prot.res[prot.lookup.id, "logFC"]
    prot.pval <- prot.res[prot.lookup.id, "P.Value"]
    msg("[R DEBUG] Protein DE stats: FC=", prot.fc, " P=", prot.pval)
  } else {
    prot.fc <- NA
    prot.pval <- NA
    msg("[R DEBUG] Protein ", matched.protein.id, " not found in DE results")
  }

  # Extract peptide DE stats
  # OPTIMIZED: Use list to collect peptide stats, then combine once
  pep.stats.list <- vector("list", length(peptides))
  pep.count <- 0

  for (pep in peptides) {
    if (pep %in% rownames(pep.res)) {
      pep.count <- pep.count + 1
      pep.stats.list[[pep.count]] <- data.frame(
        type = "Peptide",
        id = pep,
        name = pep,
        logFC = pep.res[pep, "logFC"],
        pValue = pep.res[pep, "P.Value"],
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all peptide stats at once
  if (pep.count > 0) {
    pep.stats <- do.call(rbind, pep.stats.list[1:pep.count])
  } else {
    pep.stats <- data.frame(
      type = character(0),
      id = character(0),
      name = character(0),
      logFC = numeric(0),
      pValue = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  msg("[R DEBUG] Found DE stats for ", nrow(pep.stats), " peptides")

  # Create combined data frame with protein first, then peptides
  result <- data.frame(
    type = c("Protein", rep("Peptide", nrow(pep.stats))),
    id = c(original.protein.id, pep.stats$id),
    name = c(original.protein.id, pep.stats$name),
    logFC = c(prot.fc, pep.stats$logFC),
    pValue = c(prot.pval, pep.stats$pValue),
    stringsAsFactors = FALSE
  )
  msg("[R DEBUG] Returning result with ", nrow(result), " rows")

  return(result)
}

#' Get peptide mappings for multiple proteins in a single pass
#' @description Returns a named list keyed by protein ID, each containing peptide rows
#' @param dataName Dataset name
#' @param proteinIDs Character vector of protein IDs
#' @export
GetProteinPeptideMappingBatch <- function(dataName = "", proteinIDs = character(0)) {
  if (length(proteinIDs) == 0) {
    return(list())
  }

  normalize_protein_id <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_character_)
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(NA_character_)
    first <- strsplit(x[1], ";", fixed = TRUE)[[1]][1]
    first <- trimws(first)
    if (grepl("\\|", first)) {
      parts <- strsplit(first, "\\|")[[1]]
      if (length(parts) >= 2 && nzchar(parts[2])) {
        first <- parts[2]
      }
    }
    first <- sub("_.*$", "", first)
    first <- sub("-\\d+$", "", first)
    first
  }

  dataSet <- readDataset(dataName)
  if (is.null(dataSet) || is.null(dataSet$comp.res)) {
    return(setNames(vector("list", length(proteinIDs)), as.character(proteinIDs)))
  }

  if (!file.exists("peptide_de_results.qs") || !file.exists("peptide_to_protein_map.qs")) {
    return(setNames(vector("list", length(proteinIDs)), as.character(proteinIDs)))
  }

  pep.res <- try(ov_qs_read("peptide_de_results.qs"), silent = TRUE)
  pep.map <- try(ov_qs_read("peptide_to_protein_map.qs"), silent = TRUE)
  if (inherits(pep.res, "try-error") || inherits(pep.map, "try-error")) {
    return(setNames(vector("list", length(proteinIDs)), as.character(proteinIDs)))
  }

  pep.map$Protein <- as.character(pep.map$Protein)
  pep.map$Peptide <- as.character(pep.map$Peptide)
  pep.map$Protein.norm <- vapply(pep.map$Protein, normalize_protein_id, character(1))
  pep.rownms <- rownames(pep.res)

  out <- setNames(vector("list", length(proteinIDs)), as.character(proteinIDs))

  for (pid in as.character(proteinIDs)) {
    pid.norm <- normalize_protein_id(pid)
    peptides <- pep.map$Peptide[pep.map$Protein == pid]
    if (length(peptides) == 0 && !is.na(pid.norm)) {
      peptides <- pep.map$Peptide[pep.map$Protein.norm == pid.norm]
    }
    peptides <- unique(peptides)

    if (length(peptides) == 0) {
      out[[pid]] <- data.frame(
        type = character(0),
        id = character(0),
        name = character(0),
        logFC = numeric(0),
        pValue = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    pep.fc <- rep(NA_real_, length(peptides))
    pep.p <- rep(NA_real_, length(peptides))

    idx <- match(peptides, pep.rownms)
    valid <- !is.na(idx)
    if (any(valid)) {
      pep.fc[valid] <- pep.res[idx[valid], "logFC"]
      pep.p[valid] <- pep.res[idx[valid], "P.Value"]
    }

    out[[pid]] <- data.frame(
      type = rep("Peptide", length(peptides)),
      id = peptides,
      name = peptides,
      logFC = pep.fc,
      pValue = pep.p,
      stringsAsFactors = FALSE
    )
  }

  out
}
