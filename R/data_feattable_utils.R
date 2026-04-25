#' Get Significant features from Analysis
#'
#' This function retrieves significant features from the DE analysis based on the specified parameters.
#'
#' @param dataName A character string specifying the name of the dataset.
#' @param res.nm A character string specifying the name of the output result table
#' @param p.lvl The significance threshold for p-values.
#' @param fc.lvl The fold change threshold.
#' @param inx The index for comparison (e.g., in case of multiple comparisons).
#'
#' @return A list of information including the filename, significant gene details, and counts.
#'
#' @author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @examples
#' \dontrun{
#' GetSigfeatures(dataName = "MyData", res.nm = "result_A", p.lvl = 0.05,
#'             fc.lvl = 1, inx = 1)
#' }
#'
#' @export
#' @license MIT License
#'
GetSigfeatures <-function(dataName="", res.nm="nm", p.lvl=0.05, fc.lvl=1, inx=1, FDR = "true"){
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  analSet <- readSet(analSet, "analSet");
  dataSet <- readDataset(dataName);

  paramSet$use.fdr <- as.logical(FDR);
  total <- nrow(dataSet$comp.res);
  resTable <- dataSet$comp.res;
  filename <- dataSet$filename;
  if(is.null(dataSet$fc.lvl)){
      dataSet$fc.lvl <- 0;
  }
  filename <- paste(filename, "_", res.nm, "_fc_" , dataSet$fc.lvl, ".csv", sep="");
  current.msg <- "";
  
  if (is.null(resTable) || nrow(resTable) == 0){
    msgSet$current.msg <- paste(msgSet$current.msg, "No significant features were identified using the given design and cutoff.");
  }

  # Add diagnostic logging
  #msg("[GetSigfeatures] ===========================================")
  #msg("[GetSigfeatures] DIAGNOSTIC: Initial resTable structure")
  #msg("[GetSigfeatures]   class: ", class(resTable))
  #msg("[GetSigfeatures]   is.data.frame: ", is.data.frame(resTable))
  #msg("[GetSigfeatures]   is.matrix: ", is.matrix(resTable))
  #msg("[GetSigfeatures]   typeof: ", typeof(resTable))
  if (!is.null(dim(resTable))) {
    #msg("[GetSigfeatures]   dim: ", paste(dim(resTable), collapse=" x "))
  } else {
    #msg("[GetSigfeatures]   dim: NULL")
  }
  if (!is.null(colnames(resTable))) {
    #msg("[GetSigfeatures]   colnames: ", paste(colnames(resTable), collapse=", "))
  } else {
    #msg("[GetSigfeatures]   colnames: NULL")
  }
  #msg("[GetSigfeatures]   DE method: ", dataSet$de.method)
  #msg("[GetSigfeatures] ===========================================")

  # Ensure resTable is a proper data frame/matrix
  if (!is.data.frame(resTable) && !is.matrix(resTable)) {
    #msg("[GetSigfeatures] ERROR: resTable is not a data frame or matrix. Converting...")
    resTable <- as.data.frame(resTable)
  }
  # Normalize column names to a simple character vector (strip dimensions/attributes)
  cn0 <- colnames(resTable)
  if (!is.null(cn0)) {
    dim(cn0) <- NULL
    colnames(resTable) <- as.character(cn0)
  }

  # Check dimensions
  if (is.null(dim(resTable))) {
    #msg("[GetSigfeatures] ERROR: resTable has no dimensions!")
    msgSet$current.msg <- paste(msgSet$current.msg, "Error: DE results table has incorrect structure.")
    saveSet(msgSet, "msgSet")
    return(c("error", "0", "0", "0", "0", "0", "0"))
  }

  # now rank by logFC, note, the logFC for each comparisons
  # are returned in resTable before the AveExpr columns
  # for two-class, only one column, multiple columns can be involved
  # for > comparisons - in this case, use the largest logFC among all comparisons
  # further filter by logFC
  if (dataSet$de.method=="deseq2"){
    hit.inx <- which(colnames(resTable) == "baseMean");
    dataSet$comp.res <- dataSet$comp.res.list[[inx]];
    resTable <- dataSet$comp.res;
   } else if (dataSet$de.method=="limma" || dataSet$de.method=="deqms" || dataSet$de.method=="wtt" ){
    #msg("[GetSigfeatures] Processing limma/deqms/wtt method...")

    # Safer column name matching with detailed diagnostics
    #msg("[GetSigfeatures] About to check colnames...")
    #msg("[GetSigfeatures]   colnames(resTable) class: ", class(colnames(resTable)))

    cn <- colnames(resTable)
    # Flatten any dimensional attributes to avoid match() errors (e.g., array-like colnames)
    if (!is.null(cn)) {
      dim(cn) <- NULL
      cn <- as.character(cn)
    }

    if (!is.null(cn) && length(cn) > 0) {
      #msg("[GetSigfeatures] Attempting match for AveExpr...")
      tryCatch({
        hit.inx <- match("AveExpr", cn);
        #msg("[GetSigfeatures] Match successful: hit.inx = ", hit.inx)
      }, error = function(e) {
        #msg("[GetSigfeatures] ERROR in match(): ", e$message)
        hit.inx <<- NA
      })
    } else {
      #msg("[GetSigfeatures] WARNING: resTable has no column names")
      hit.inx <- NA
    }

    dataSet$comp.res <- dataSet$comp.res.list[[inx]];
    resTable <- dataSet$comp.res;

    if (is.na(hit.inx) && dataSet$de.method == "wtt") {
      ave.expr <- rowMeans(
        dataSet$data.norm[rownames(resTable), , drop = FALSE],
        na.rm = TRUE)
      resTable$AveExpr <- ave.expr
      # Flatten colnames to avoid match() error
      cn <- colnames(resTable)
      dim(cn) <- NULL
      cn <- as.character(cn)
      hit.inx <- match("AveExpr", cn)
      dataSet$comp.res <- resTable
      dataSet$comp.res.list[[inx]] <- resTable
    }
  } else {
    hit.inx <- which(colnames(resTable) == "logCPM");
    dataSet$comp.res <- dataSet$comp.res.list[[inx]];
    resTable <- dataSet$comp.res;
  }

  if(length(hit.inx) == 0){
    hit.inx <- 1;
  }

  #msg("[GetSigfeatures] DEBUG: About to filter NA rows...")
  resTable <- resTable[!is.na(resTable[,1]),]
  #msg("[GetSigfeatures] DEBUG: After NA filter, nrow=", nrow(resTable))
  orig.resTable <- resTable;
  # select based on p-value
  #msg("[GetSigfeatures] DEBUG: About to filter by p-value (FDR=", FDR, ", p.lvl=", p.lvl, ")...")
  if(FDR == "true"){
      hit.inx.p <- resTable$adj.P.Val <= p.lvl;
  } else {
      hit.inx.p <- resTable$P.Value <= p.lvl;
  }
  #msg("[GetSigfeatures] DEBUG: hit.inx.p: class=", class(hit.inx.p), " length=", length(hit.inx.p), " sum=", sum(hit.inx.p, na.rm=TRUE))

  resTable<-resTable[hit.inx.p,];
  #msg("[GetSigfeatures] DEBUG: After p-value filter, nrow=", nrow(resTable))
  
  if (is.na(hit.inx) || hit.inx < 2) {
    maxFC.inx <- 1
  } else {
    maxFC.inx <- hit.inx - 1
  }

  if (ncol(resTable) == 0) {
    logfc.mat <- matrix(0, nrow = nrow(resTable), ncol = 1,
                        dimnames = list(rownames(resTable), "logFC"))
  } else {
    cols_to_take <- seq_len(min(maxFC.inx, ncol(resTable)))
    logfc.mat <- resTable[, cols_to_take, drop = FALSE];
  }
  #msg("[GetSigfeatures] DEBUG: About to filter by fold-change (fc.lvl=", fc.lvl, ")...")
  if(paramSet$oneDataAnalType == "dose"){
    pos.mat <- abs(logfc.mat);
    fc.vec <- apply(pos.mat, 1, max);   # for > comparisons - in this case, use the largest logFC among all comparisons
    hit.inx.fc <- fc.vec >= fc.lvl;
    resTable <- resTable[hit.inx.fc,];
  } else if (dataSet$de.method=="deseq2" || dataSet$de.method=="edger" || dataSet$de.method=="limma" || dataSet$de.method=="deqms" || dataSet$de.method=="wtt"){
    #msg("[GetSigfeatures] DEBUG: Using limma/deqms/wtt FC filter...")
    pos.mat <- abs(logfc.mat);
    fc.vec <- pos.mat[,1];
    #msg("[GetSigfeatures] DEBUG: fc.vec: class=", class(fc.vec), " length=", length(fc.vec))
    hit.inx.fc <- fc.vec >= fc.lvl;
    #msg("[GetSigfeatures] DEBUG: hit.inx.fc: class=", class(hit.inx.fc), " length=", length(hit.inx.fc), " sum=", sum(hit.inx.fc, na.rm=TRUE))
    resTable <- resTable[hit.inx.fc,];
  }else {
    pos.mat <- abs(logfc.mat[,inx]);
    fc.vec <- pos.mat;
    hit.inx.fc <- fc.vec >= fc.lvl;
    resTable <- resTable[hit.inx.fc,];
  }
  #msg("[GetSigfeatures] DEBUG: After FC filter, nrow=", nrow(resTable))
  
  if (nrow(resTable) == 0){
    msgSet$current.msg <- paste(msgSet$current.msg, "No significant features were identified using the given design and cutoff.");
  }
  
  ### Note, rowname of resTable must be entrez ID
  # calculate differently if dose-response
  de.Num <- nrow(resTable);

  non.de.Num <- nrow(dataSet$data.norm) - de.Num;
  
  # may need to update data, class and meta.info
  #msg("[GetSigfeatures] DEBUG: Loading data.norm, cls, meta.info...")
  data <- dataSet$data.norm;
  cls <- dataSet$cls;
  meta.info <- dataSet$meta.info;
  #msg("[GetSigfeatures] DEBUG: cls: class=", class(cls), " length=", length(cls))
  #msg("[GetSigfeatures] DEBUG: cls dimensions: ", paste(dim(cls), collapse="x"))
  grp.nms <- levels(cls);
  #msg("[GetSigfeatures] DEBUG: grp.nms: ", paste(grp.nms, collapse=", "))
  #msg("[GetSigfeatures] DEBUG: grp.nms: class=", class(grp.nms), " length=", length(grp.nms))
  #msg("[GetSigfeatures] DEBUG: grp.nms dimensions: ", paste(dim(grp.nms), collapse="x"))

  #msg("[GetSigfeatures] DEBUG: About to perform cls %in% grp.nms...")
  # Flatten cls to avoid match() dimension errors in %in% operator
  cls_vec <- as.vector(cls)
  dim(cls_vec) <- NULL
  #msg("[GetSigfeatures] DEBUG: cls_vec: class=", class(cls_vec), " length=", length(cls_vec))
  hit.inx <- cls_vec %in% grp.nms;
  #msg("[GetSigfeatures] DEBUG: hit.inx from %in%: class=", class(hit.inx), " length=", length(hit.inx), " sum=", sum(hit.inx))

  if (sum(hit.inx) < length(hit.inx)){
    msgSet$current.msg <- paste(msgSet$current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
    cls <- factor(cls[hit.inx]);
    cls.lvls <- levels(cls);
    data <- data[,hit.inx];
    meta.info <- dataSet$meta.info[hit.inx,];
  }
  #msg("[GetSigfeatures] DEBUG: About to save data.stat.qs...")
  ov_qs_save(data, file="data.stat.qs");
  #msg("[GetSigfeatures] DEBUG: About to order comp.res by P.Value...")
  o <- with(dataSet$comp.res, order(P.Value, -abs(logFC), na.last = TRUE))
  dataSet$comp.res <- dataSet$comp.res[o, , drop = FALSE]

  #msg("[GetSigfeatures] DEBUG: About to filter comp.res using %in% on rownames...")
  rn_comp <- rownames(dataSet$comp.res)
  rn_res <- rownames(resTable)
  #msg("[GetSigfeatures] DEBUG: rn_comp: class=", class(rn_comp), " length=", length(rn_comp))
  #msg("[GetSigfeatures] DEBUG: rn_res: class=", class(rn_res), " length=", length(rn_res))
  # Flatten rownames to avoid dimension errors in %in%
  dim(rn_comp) <- NULL
  dim(rn_res) <- NULL
  rn_comp <- as.character(rn_comp)
  rn_res <- as.character(rn_res)
  #msg("[GetSigfeatures] DEBUG: After flattening - rn_comp: class=", class(rn_comp), " length=", length(rn_comp))
  dataSet$comp.res <- dataSet$comp.res[!(rn_comp %in% rn_res), , drop = FALSE]
  #msg("[GetSigfeatures] DEBUG: After filtering, comp.res nrow=", nrow(dataSet$comp.res))

  #msg("[GetSigfeatures] DEBUG: About to rbind resTable with comp.res...")
  dataSet$comp.res <- rbind(resTable, dataSet$comp.res)
  #msg("[GetSigfeatures] DEBUG: After rbind, comp.res nrow=", nrow(dataSet$comp.res))

  # Keep every contrast result aligned with the updated comp.res row order
  if (!is.null(dataSet$comp.res.list) && length(dataSet$comp.res.list) > 0) {
    for (i in seq_along(dataSet$comp.res.list)) {
      res_i <- dataSet$comp.res.list[[i]]
      # Reorder/extend rows to match the master comp.res rownames
      res_i <- res_i[rownames(dataSet$comp.res), colnames(res_i), drop = FALSE]
      dataSet$comp.res.list[[i]] <- res_i
    }
  }
  
  
  #msg("[GetSigfeatures] DEBUG: Setting sig.mat...")
  dataSet$sig.mat <- resTable;

  # PHOSPHO DATA CHECK: Phosphosite IDs should NOT be annotated as Entrez IDs
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  #msg("[GetSigfeatures] DEBUG: About to write annotation file (is_phospho=", is_phospho, ", annotated=", dataSet$annotated, ")...")
  if (is_phospho) {
    # For phospho data, write results without Entrez annotation
    #msg("[GetSigfeatures] PHOSPHO DATA: Exporting phosphosite results without Entrez annotation")
    gene.anot <- NULL;
    fast.write(signif(resTable,5), file=filename);
  } else if (dataSet$annotated){ # annotated to entrez
    #msg("[GetSigfeatures] DEBUG: Annotated path - extracting rownames from comp.res...")
    anot.id <- rownames(dataSet$comp.res);
    # Flatten anot.id to avoid dimension errors
    dim(anot.id) <- NULL
    anot.id <- as.character(anot.id)
    #msg("[GetSigfeatures] DEBUG: anot.id: class=", class(anot.id), " length=", length(anot.id))
    gene.anot <- doEntrezIDAnot(anot.id, paramSet$data.org, paramSet$data.idType)
    fast.write(cbind(EntrezID=anot.id, signif (dataSet$comp.res,5), Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
  } else if (file.exists("annotation.qs")){ # annotation information available
    #msg("[GetSigfeatures] DEBUG: annotation.qs exists - loading...")
    anot.id <- ov_qs_read("annotation.qs");
    #msg("[GetSigfeatures] DEBUG: anot.id loaded: class=", class(anot.id), " length=", length(anot.id))
    #msg("[GetSigfeatures] DEBUG: anot.id is.named=", !is.null(names(anot.id)))
    if(!is.null(names(anot.id))) {
      #msg("[GetSigfeatures] DEBUG: anot.id names head: ", paste(head(names(anot.id), 10), collapse=", "))
    }
    #msg("[GetSigfeatures] DEBUG: anot.id values head: ", paste(head(anot.id, 10), collapse=", "))

    feature.vec <- rownames(dataSet$comp.res);
    gene.anot <- doEntrezIDAnot(feature.vec, paramSet$data.org, paramSet$data.idType);

    #msg("[GetSigfeatures] DEBUG: About to cbind for fast.write...")
    # entrez.vec is already flattened from lines 296-297
    fast.write(signif(dataSet$comp.res,5), row.names=T, file=filename);
    #msg("[GetSigfeatures] DEBUG: fast.write completed.")

    #msg("[GetSigfeatures] DEBUG: About to assign rownames to gene.anot...")
    # Ensure both gene.anot rownames and feature.vec are clean before assignment
    if(!is.null(rownames(gene.anot))) {
      rn_temp <- rownames(gene.anot)
      dim(rn_temp) <- NULL
      #msg("[GetSigfeatures] DEBUG: Cleared existing rownames dimensions")
    }
    # feature.vec should already be flattened from lines 289-290, but ensure it again
    dim(feature.vec) <- NULL
    feature.vec <- as.character(feature.vec)
    #msg("[GetSigfeatures] DEBUG: feature.vec for rownames: class=", class(feature.vec), " length=", length(feature.vec))
    rownames(gene.anot) <- feature.vec;
    #msg("[GetSigfeatures] DEBUG: rownames assignment completed.")
  } else {
    #msg("[GetSigfeatures] DEBUG: No annotation - writing plain results...")
    gene.anot <- NULL;
    fast.write(signif(resTable,5), file=filename);
  }
  #msg("[GetSigfeatures] DEBUG: Annotation section complete.")
  if (is_phospho) {
    rn_comp_res <- rownames(dataSet$comp.res)
    dim(rn_comp_res) <- NULL
    rn_comp_res <- as.character(rn_comp_res)
    phospho_labels <- rn_comp_res

    phospho.map <- try(readDataQs("phospho_symbol_map.qs", paramSet$anal.type, dataSet$name), silent = TRUE)
    if (!inherits(phospho.map, "try-error") && !is.null(phospho.map) && nrow(phospho.map) > 0 &&
        "symbol" %in% colnames(phospho.map)) {
      hit.ids <- intersect(rn_comp_res, rownames(phospho.map))
      if (length(hit.ids) > 0) {
        map.syms <- as.character(phospho.map[hit.ids, "symbol", drop = TRUE])
        for (ii in seq_along(rn_comp_res)) {
          pid <- rn_comp_res[ii]
          if (!(pid %in% hit.ids)) next
          sym <- map.syms[which(hit.ids == pid)[1]]
          if (is.na(sym) || sym == "" || sym == "NA" || sym == pid) next
          # Use symbol directly - it already contains the full display name
          # with isoform and site suffix (e.g., "DOCK10-2_S_12")
          phospho_labels[ii] <- sym
        }
      }
    }
    analSet$comp.features.symbols <- phospho_labels
  } else if(is.null(gene.anot)){
    rn_comp_res <- rownames(dataSet$comp.res)
    dim(rn_comp_res) <- NULL
    rn_comp_res <- as.character(rn_comp_res)
    analSet$comp.features.symbols <- rn_comp_res; # use the id provided
  }else{
    analSet$comp.features.symbols <- gene.anot$symbol;
  }
  dataSet$cls.stat <- cls;
  dataSet$meta.stat <- meta.info;

  # now do protein mapping for network only applicable for annotated
  #msg("[GetSigfeatures] DEBUG: About to extract gene from rownames(resTable)...")
  gene <- rownames(resTable);
  dim(gene) <- NULL
  gene <- as.character(gene)
  #msg("[GetSigfeatures] DEBUG: gene extracted: class=", class(gene), " length=", length(gene))
  
  logFC <- unname(logfc.mat[,1]);
  geneList <- paste(gene, logFC, collapse="\n");
  up <- nrow(resTable[which(logfc.mat[,paramSet$selectedFactorInx]> fc.lvl),])
  down <- nrow(resTable[which(logfc.mat[,paramSet$selectedFactorInx]< -fc.lvl),])
  saveSet(msgSet, "msgSet");
  
  data.norm <- dataSet$data.norm
  colnames(data.norm) <- NULL
  lst <- list(colnames(dataSet$data.norm),data.norm, dataSet$meta.info, dataSet$comp.res, rownames(data.norm), org=paramSet$data.org)
  # OPTIMIZED: Use jsonlite::write_json instead of rjson + sink/cat
  jsonlite::write_json(lst, "ProteoAnalyst_matrix.json", auto_unbox = TRUE, pretty = FALSE);

        if (dataSet$de.method %in% c("deseq2", "edger", "limma", "deqms", "wtt")) {

  significant_gene_table <- list()    # holds one data-frame per comparison

  for (inx in seq_along(dataSet$comp.res.list)) {

    resTable <- dataSet$comp.res.list[[inx]]

    resTable <- resTable[!is.na(resTable[, 1]), ]

    deg.pass <- if (FDR == "true")  resTable$adj.P.Val <= p.lvl
                else                resTable$P.Value   <= p.lvl

    lfc.pass <- abs(resTable[ , "logFC"]) >= fc.lvl

    all_pass <- deg.pass & lfc.pass
    if (!any(all_pass)) {                            # nothing passed
      significant_gene_table[[inx]] <- data.frame() # keep list length constant
      next
    }

    sig <- resTable[all_pass, ]
    sig$GeneID      <- rownames(sig)                # preserve raw ID
    sig$Comparison  <- names(dataSet$comp.res.list)[[inx]]

    # PHOSPHO DATA CHECK: Skip Entrez annotation for phosphosites
    if (is_phospho) {
      sig$Symbol <- sig$GeneID  # Use phosphosite ID as symbol
      sig$Name   <- "Phosphosite"
    } else {
      res.anot <- doEntrezIDAnot(sig$GeneID,
                                 paramSet$data.org,
                                 paramSet$data.idType)
      sig$Symbol <- res.anot$symbol
      sig$Name   <- res.anot$name
    }

    significant_gene_table[[inx]] <- sig
  }

  final_table <- do.call(rbind, significant_gene_table)  # may have duplicates

  output_file <- paste0(dataName, "_logFC_",format(as.numeric(fc.lvl), digits = 2, nsmall = 0, trim = TRUE, scientific = FALSE),
                        "_Significant_features.csv")

  if (nrow(final_table) > 0) {
    write.csv(final_table[ , setdiff(names(final_table), "GeneID")],
              file = output_file, row.names = FALSE)

    all_significant_features <- unique(final_table$GeneID)  # de-duplicate here
    de.Num.total          <- length(all_significant_features)

    #msg("Significant features table exported to: ", output_file)
  } else {
    de.Num.total <- 0
    #msg("No significant features identified to export.")
  }

  ## ---------- bookkeeping -----------------------------------------------
  if (de.Num.total == 0) {
    msgSet$current.msg <- paste(
      msgSet$current.msg,
      "No significant features were identified using the given design and cutoff."
    )
  }
  analSet$sig.gene.count.total <- de.Num.total


    ## 2) Build & write the binary 0/1 incidence table (one row per gene)
    comp_list  <- dataSet$comp.res.list
    comp_names <- names(comp_list)
    use_fdr    <- (FDR == "true") || isTRUE(FDR)

    # collect DE gene sets per comparison
    pass_sets <- vector("list", length(comp_list))
    names(pass_sets) <- comp_names
    for (i in seq_along(comp_list)) {
      rt <- comp_list[[i]]
      if (is.null(rt) || nrow(rt) == 0) { pass_sets[[i]] <- character(0); next }

      pcol <- if (use_fdr && "adj.P.Val" %in% names(rt)) "adj.P.Val" else "P.Value"
      lfc_col <- if ("logFC" %in% names(rt)) "logFC"
                 else if ("log2FoldChange" %in% names(rt)) "log2FoldChange"
                 else NULL

      deg_pass <- is.finite(rt[[pcol]]) & (rt[[pcol]] <= p.lvl)
      lfc_pass <- if (!is.null(lfc_col)) is.finite(rt[[lfc_col]]) & (abs(rt[[lfc_col]]) >= fc.lvl) else TRUE

      pass_sets[[i]] <- rownames(rt)[deg_pass & lfc_pass]
    }

    features <- sort(unique(unlist(pass_sets)))

    # derive a second filename from your existing output_file
    out_bin_file <- sub("_Significant_features\\.csv$", "_DE_binary_matrix.csv", output_file)
    if (identical(out_bin_file, output_file)) {
      out_bin_file <- paste0(tools::file_path_sans_ext(output_file), "_DE_binary_matrix.csv")
    }

    if (length(features) > 0) {
      de_mat <- matrix(0L, nrow = length(features), ncol = length(comp_names),
                       dimnames = list(features, comp_names))
      for (i in seq_along(pass_sets)) {
        if (length(pass_sets[[i]]) > 0) de_mat[pass_sets[[i]], i] <- 1L
      }

      # optional annotation; if it fails, keep Symbol/Name as NA
      # PHOSPHO DATA CHECK: Skip Entrez annotation for phosphosites
      if (is_phospho) {
        annot <- data.frame(symbol = features,  # Use phosphosite ID as symbol
                           name   = rep("Phosphosite", length(features)))
      } else {
        annot <- tryCatch(
          doEntrezIDAnot(features, paramSet$data.org, paramSet$data.idType),
          error = function(e) data.frame(symbol = rep(NA_character_, length(features)),
                                         name   = rep(NA_character_, length(features)))
        )
      }

      de_df <- data.frame(
        GeneID = features,
        Symbol = annot$symbol,
        Name   = annot$name,
        as.data.frame(de_mat, check.names = FALSE),
        check.names = FALSE
      )
      de_df$DE_Count <- rowSums(de_mat)
      de_df$Comparisons <- apply(de_mat, 1, function(z) {
        hits <- which(z == 1L)
        if (length(hits) == 0) "" else paste(comp_names[hits], collapse = ";")
      })

      write.csv(de_df, file = out_bin_file, row.names = FALSE)
      #msg("Binary DE matrix written to: ", out_bin_file)
    } else {
      # write an empty shell (optional)
      write.csv(data.frame(GeneID = character(0)), file = out_bin_file, row.names = FALSE)
      #msg("No features passed any comparison; wrote empty binary matrix to: ", out_bin_file)
    }
}
dataSet$comp.res.list <- lapply(dataSet$comp.res.list, function(tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tbl)

  pcol <- if ( "P.Value" %in% names(tbl)) "P.Value"
          else if ("padj" %in% names(tbl)) "padj"
          else "P.Value"
  pvec <- suppressWarnings(as.numeric(tbl[[pcol]]))

  # tie-break by |logFC|
  if ("logFC" %in% names(tbl)) {
    lfc <- abs(suppressWarnings(as.numeric(tbl$logFC)))
  } else if ("log2FoldChange" %in% names(tbl)) {
    lfc <- abs(suppressWarnings(as.numeric(tbl$log2FoldChange)))
  } else {
    lfc <- rep(0, nrow(tbl))
  }

  sig_idx <- (pvec <= p.lvl) & (lfc >= fc.lvl)

  sig  <- tbl[sig_idx,  , drop = FALSE]
  rest <- tbl[!sig_idx, , drop = FALSE]

  if (nrow(sig)  > 0) sig  <- sig [order(pvec[sig_idx],  -lfc[sig_idx],  na.last = NA),  , drop = FALSE]
  if (nrow(rest) > 0) rest <- rest[order(pvec[!sig_idx],            na.last = TRUE), , drop = FALSE]

  rbind(sig, rest)
})

## Rebuild the combined table with “sig first” and sorted by p
# p column for the combined table
pcol <- if (paramSet$use.fdr && "adj.P.Val" %in% names(resTable)) "adj.P.Val"
        else if ("padj" %in% names(resTable)) "padj"
        else "P.Value"

# ensure the sig block (resTable) itself is ordered
resTable <- resTable[order(as.numeric(resTable[[pcol]]),
                           -abs(suppressWarnings(as.numeric(resTable[[if ("logFC" %in% names(resTable)) "logFC" else 1]]))),
                           na.last = NA),
                     , drop = FALSE]

# the remainder (“other”), ordered too
other <- dataSet$comp.res[!(rownames(dataSet$comp.res) %in% rownames(resTable)), , drop = FALSE]
if (nrow(other) > 0 && pcol %in% names(other)) {
  # pick tie-break column for 'other'
  tie_col <- if ("logFC" %in% names(other)) "logFC" else 1
  other <- other[order(as.numeric(other[[pcol]]),
                       -abs(suppressWarnings(as.numeric(other[[tie_col]]))),
                       na.last = TRUE),
                 , drop = FALSE]
}

dataSet$comp.res <- rbind(resTable, other)

  output_file <- paste0(dataName, "_logFC_",format(as.numeric(fc.lvl), digits = 2, nsmall = 0, trim = TRUE, scientific = FALSE),
                        "_Significant_features.csv")
    write.csv(dataSet$comp.res,
              file = output_file, row.names = FALSE)

  analSet$sig.gene.count <- de.Num;
  saveSet(analSet, "analSet");

  dataSet$pval <- p.lvl;
  dataSet$fc.val <- fc.lvl;
  dataSet$comp.res.filename <- filename;
  res <- RegisterData(dataSet);
  saveSet(paramSet, "paramSet");
  return(c(output_file, de.Num, geneList, total, up, down, non.de.Num));
}
