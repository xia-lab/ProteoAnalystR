# Helper function to suppress ALL graphics device popups on macOS
.suppress_quartz <- function(expr) {
  # Save current state
  old_device <- getOption("device")
  old_bitmapType <- getOption("bitmapType")

  # Override device options to prevent Quartz
  options(device = function(...) pdf(file = NULL))
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    options(bitmapType = "cairo")
  }

  # Close any stray devices
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)

  # Ensure cleanup happens
  on.exit({
    while (dev.cur() > 1) try(dev.off(), silent = TRUE)
    options(device = old_device, bitmapType = old_bitmapType)
  })

  # Execute the expression
  eval(expr, envir = parent.frame())
}

.trim_cem_object_for_save <- function(cem) {
  # Simplified version: only clear slots we know are safe to clear
  # Avoid problematic slots that may not exist in all CEMiTool versions

  # Just return the object as-is - trimming is optional optimization
  # If slot access fails, better to have a larger saved object than crash
  return(cem)
}

my.build.cemi.net <- function(dataName,
                              filter      = TRUE,
                              min_ngen    = 30,
                              cor_method  = "pearson",
                              verbose     = FALSE,
                              classCol    = NULL,
                              auto_impute = TRUE) {   # <-- optional: auto-impute NAs before analysis
  tryCatch({

    ## 1 · load dataset -------------------------------------------------
    dataSet  <- readDataset(dataName)

    # Validate dataSet was loaded
    if (is.null(dataSet)) {
      stop("Failed to load dataset: ", dataName)
    }
    if (is.null(dataSet$data.norm)) {
      stop("Dataset '", dataName, "' has no normalized data (data.norm is NULL)")
    }
    if (is.null(dataSet$meta.info)) {
      stop("Dataset '", dataName, "' has no metadata (meta.info is NULL)")
    }

    expr_mat <- as.data.frame(dataSet$data.norm)    # features × samples

    ## metadata: keep *all* columns, coerce factors -> character
    meta_df <- dataSet$meta.info
    meta_df[] <- lapply(meta_df, \(x) {
      if (is.factor(x)) {
        as.character(x)
      } else {
        x
      }
    })

    ## decide which column is the class
    if (is.null(classCol) || is.na(classCol) || !nzchar(classCol)) {
      classCol <- colnames(meta_df)[1]              # first column by default
    }
    if (!is.character(classCol) || length(classCol) != 1) {
      stop("classCol must be a single character string")
    }
    if (!(classCol %in% colnames(meta_df))) {
      stop("classCol '", classCol, "' not found in meta.info. Available columns: ",
           paste(colnames(meta_df), collapse=", "))
    }

    ## build annotation table (SampleName + all meta)
    annot_df <- data.frame(SampleName = rownames(meta_df),
                           meta_df,
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

    # Ensure SampleName column has no NA values
    if (any(is.na(annot_df$SampleName))) {
      stop("Annotation table has NA values in SampleName column")
    }

    # Check that expr_mat and annot_df have matching samples
    if (!all(colnames(expr_mat) %in% annot_df$SampleName)) {
      stop("Expression matrix contains samples not found in metadata")
    }

    # Check the class column for NA values - this is critical!
    if (classCol %in% colnames(annot_df)) {
      class_values <- annot_df[[classCol]]
      if (any(is.na(class_values))) {
        na_samples <- annot_df$SampleName[is.na(class_values)]
        stop("Class column '", classCol, "' contains NA values for samples: ",
             paste(na_samples, collapse=", "),
             ". Please ensure all samples have a valid class assignment.")
      }
      # Check if class column is character/factor with valid values
      if (length(unique(class_values)) < 2) {
        stop("Class column '", classCol, "' must have at least 2 different groups, found: ",
             length(unique(class_values)))
      }
      msg("Class column '", classCol, "' has ", length(unique(class_values)),
          " groups: ", paste(unique(class_values), collapse=", "))
    }

    ## 2 · run CEMiTool -------------------------------------------------
    suppressPackageStartupMessages({
      library(CEMiTool)
      library(WGCNA)
    })

    # Configure WGCNA options to avoid thread-related issues
    # Disable WGCNA threading to avoid potential NA/logical issues
    WGCNA::disableWGCNAThreads()

    # Force single-threaded BLAS/OpenMP to prevent Rserve hanging during matrix operations
    # This is critical when running under Rserve - multi-threaded BLAS can deadlock
    tryCatch({
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
        msg("Set BLAS/OpenMP to single-threaded mode");
      } else {
        # Alternative: set environment variables
        Sys.setenv(OMP_NUM_THREADS = 1)
        Sys.setenv(OPENBLAS_NUM_THREADS = 1)
        Sys.setenv(MKL_NUM_THREADS = 1)
        msg("Set threading environment variables to 1");
      }
    }, error = function(e) {
      msg("Warning: Could not configure BLAS threading: ", conditionMessage(e));
    })

    # Set WGCNA options to be more tolerant
    options(stringsAsFactors = FALSE)


    # Validate parameters to avoid "missing value where TRUE/FALSE needed" errors

    # Validate filter - ensure it's a proper logical value
    if (is.null(filter) || length(filter) == 0 || is.na(filter)) {
      filter <- TRUE  # default to TRUE if invalid
      msg("Warning: filter was NULL/NA, defaulting to TRUE");
    } else {
      filter <- as.logical(filter)[1]  # coerce to logical, take first element
      if (is.na(filter)) {
        filter <- TRUE
        msg("Warning: filter could not be coerced to logical, defaulting to TRUE");
      }
    }

    # Validate verbose - ensure it's a proper logical value
    if (is.null(verbose) || length(verbose) == 0 || is.na(verbose)) {
      verbose <- FALSE  # default to FALSE if invalid
      msg("Warning: verbose was NULL/NA, defaulting to FALSE");
    } else {
      verbose <- as.logical(verbose)[1]  # coerce to logical, take first element
      if (is.na(verbose)) {
        verbose <- FALSE
        msg("Warning: verbose could not be coerced to logical, defaulting to FALSE");
      }
    }

    # Validate cor_method - ensure it has a valid value before passing to match.arg
    if (is.null(cor_method) || length(cor_method) == 0 || is.na(cor_method) || !nzchar(as.character(cor_method))) {
      cor_method <- "pearson"  # default to pearson if invalid
      msg("Warning: cor_method was NULL/NA/empty, defaulting to 'pearson'");
    } else {
      cor_method <- as.character(cor_method)[1]  # ensure it's character, take first element
      # Ensure it's one of the valid choices
      if (!(cor_method %in% c("pearson", "spearman"))) {
        msg("Warning: cor_method '", cor_method, "' not valid, defaulting to 'pearson'");
        cor_method <- "pearson"
      }
    }

    # Validate min_ngen - ensure it's a positive integer
    if (is.null(min_ngen) || length(min_ngen) == 0 || is.na(min_ngen)) {
      min_ngen <- 30  # default to 30 if invalid
      msg("Warning: min_ngen was NULL/NA, defaulting to 30");
    } else {
      min_ngen <- as.integer(min_ngen)[1]  # coerce to integer, take first element
      if (is.na(min_ngen) || min_ngen < 1) {
        min_ngen <- 30
        msg("Warning: min_ngen was invalid (<1), defaulting to 30");
      }
    }


    # Check for NA values in expression matrix - this is critical!
    if (any(is.na(expr_mat))) {
      na_count <- sum(is.na(expr_mat))
      na_pct <- 100 * na_count / (nrow(expr_mat) * ncol(expr_mat))
      msg("Expression matrix contains ", na_count, " NA values (",
          sprintf("%.2f", na_pct), "%)");

      # Filter genes with too many missing values per gene
      # For co-expression, we need genes with sufficient valid observations
      na_per_gene <- rowSums(is.na(expr_mat)) / ncol(expr_mat)
      genes_too_many_na <- na_per_gene > 0.5  # Gene missing in >50% of samples

      if (sum(genes_too_many_na) > 0) {
        msg("Filtering ", sum(genes_too_many_na), " genes with >50% missing values");
        expr_mat <- expr_mat[!genes_too_many_na, , drop = FALSE]

        # Recalculate overall NA percentage after filtering
        na_count_after <- sum(is.na(expr_mat))
        na_pct_after <- 100 * na_count_after / (nrow(expr_mat) * ncol(expr_mat))
        msg("After filtering: ", na_count_after, " NA values (",
            sprintf("%.2f", na_pct_after), "%)");
        na_pct <- na_pct_after  # Update for subsequent checks
      }

      # Option to automatically impute remaining NAs
      if (auto_impute && na_pct > 0) {
        msg("Auto-imputing ", sum(is.na(expr_mat)), " remaining NA values using minimum per gene...");

        # Simple but effective: replace NAs with minimum value per gene (conservative approach)
        expr_mat_imputed <- t(apply(expr_mat, 1, function(row) {
          if (any(is.na(row))) {
            min_val <- min(row, na.rm = TRUE)
            row[is.na(row)] <- min_val
          }
          row
        }))
        colnames(expr_mat_imputed) <- colnames(expr_mat)
        expr_mat <- as.data.frame(expr_mat_imputed)

        msg("After imputation: 0 NA values (100% complete)");
        na_pct <- 0  # Update for subsequent checks
      }

      # For co-expression analysis, we can handle some NAs using pairwise complete observations
      # But if too many remain and auto_impute is disabled, warn or stop
      if (na_pct > 50) {
        stop("Expression matrix still contains too many NA values (", sprintf("%.2f", na_pct),
             "%) after filtering. This will prevent proper network construction. ",
             "Please impute or filter missing values before co-expression analysis, ",
             "or set auto_impute=TRUE to automatically impute remaining NAs.")
      } else if (na_pct > 20) {
        msg("WARNING: Moderate percentage of NA values (", sprintf("%.2f", na_pct),
            "%). Co-expression network may be less reliable. ",
            "Consider setting auto_impute=TRUE for better results.");
      }
    }

    # Check for infinite or NaN values
    if (any(!is.finite(as.matrix(expr_mat)))) {
      inf_count <- sum(!is.finite(as.matrix(expr_mat)))
      stop("Expression matrix contains ", inf_count, " infinite or NaN values. ",
           "Please clean the data before co-expression analysis.")
    }

    # Check for zero-variance genes (all same value) - these cause issues with correlation
    gene_vars <- apply(expr_mat, 1, var, na.rm = TRUE)
    zero_var_genes <- sum(gene_vars == 0 | is.na(gene_vars))
    if (zero_var_genes > 0) {
      msg("Removing ", zero_var_genes, " genes with zero variance");
      expr_mat <- expr_mat[gene_vars > 0 & !is.na(gene_vars), , drop = FALSE]
    }

    # Final check: ensure we have enough genes
    if (nrow(expr_mat) < min_ngen) {
      stop("After filtering, only ", nrow(expr_mat), " genes remain, ",
           "but min_ngen is set to ", min_ngen, ". ",
           "Please reduce min_ngen or relax filtering.")
    }


    # PRACTICAL LIMIT: Cap at top 5000 features by IQR to ensure reasonable computation time
    # Co-expression analysis is computationally expensive (O(n²) for correlation matrix)
    MAX_FEATURES <- 5000
    if (nrow(expr_mat) > MAX_FEATURES) {
      msg("Dataset contains ", nrow(expr_mat), " features, selecting top ", MAX_FEATURES, " by IQR...");

      # Compute IQR for each feature
      feature_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)

      # Select top features by IQR
      top_features <- order(feature_iqr, decreasing = TRUE)[1:MAX_FEATURES]
      expr_mat <- expr_mat[top_features, , drop = FALSE]

      msg("Selected top ", MAX_FEATURES, " features by IQR for co-expression analysis");
    }


    # Additional validation: check that correlation matrix can be computed
    # This is what CEMiTool will try to do, so test it first
    tryCatch({
      test_cor <- cor(t(expr_mat), method = cor_method, use = "pairwise.complete.obs")

      # Check for NA values in correlation matrix
      if (any(is.na(test_cor))) {
        na_cor_count <- sum(is.na(test_cor))
        warning("Correlation matrix contains ", na_cor_count, " NA values. ",
                "This may cause issues with network construction.")
      }

      # Check for genes with all-NA correlations (these cause the 'a < 0' error)
      all_na_cors <- rowSums(is.na(test_cor)) == ncol(test_cor)
      if (any(all_na_cors)) {
        bad_genes <- rownames(test_cor)[all_na_cors]
        msg("Removing ", length(bad_genes), " genes with undefined correlations");
        expr_mat <- expr_mat[!rownames(expr_mat) %in% bad_genes, , drop = FALSE]
      }

    }, error = function(e) {
      stop("Failed to compute correlation matrix: ", conditionMessage(e),
           ". This suggests the expression data has issues that prevent co-expression analysis.")
    })

    # Final dimension check
    if (nrow(expr_mat) < min_ngen) {
      stop("After all filtering, only ", nrow(expr_mat), " genes remain, ",
           "but min_ngen is set to ", min_ngen, ". ",
           "Cannot proceed with co-expression analysis.")
    }

    # FIX: Suppress Quartz popup on macOS - completely disable plotting during cemitool
    # We'll generate plots separately using the other functions

    # Use callr to run CEMiTool in isolated subprocess to reduce memory bandwidth
    # This prevents memory bloat in the main Rserve process

    # Enable callr for memory efficiency
    USE_CALLR <- TRUE

    if (USE_CALLR && requireNamespace("callr", quietly = TRUE)) {
      msg("Using callr subprocess...");
      # Run cemitool in separate process with stdout/stderr capture
      result <- callr::r(
        func = function(expr_mat, annot_df, filter, min_ngen, cor_method, classCol, verbose) {
          suppressPackageStartupMessages({
            library(CEMiTool)
            library(WGCNA)
          })

          # Disable threading in subprocess too
          WGCNA::disableWGCNAThreads()

          message("[callr subprocess] Starting CEMiTool with ", nrow(expr_mat), " features x ", ncol(expr_mat), " samples")
          message("[callr subprocess] Parameters: filter=", filter, ", min_ngen=", min_ngen, ", cor_method=", cor_method)

          cem <- cemitool(expr              = expr_mat,
                          annot             = annot_df,
                          filter            = filter,
                          min_ngen          = min_ngen,
                          cor_method        = cor_method,
                          class_column      = classCol,
                          verbose           = verbose,
                          plot              = FALSE,
                          plot_diagnostics  = FALSE)

          message("[callr subprocess] CEMiTool completed successfully")
          return(cem)
        },
        args = list(
          expr_mat   = expr_mat,
          annot_df   = annot_df,
          filter     = filter,
          min_ngen   = min_ngen,
          cor_method = cor_method,
          classCol   = classCol,
          verbose    = verbose
        ),
        package = TRUE,  # Use package environment for better isolation
        stdout = "|",    # Capture stdout
        stderr = "|"     # Capture stderr
      )

      # Relay subprocess output to main session
      if (!is.null(attr(result, "stdout")) && nzchar(attr(result, "stdout"))) {
        msg("[callr stdout] ", attr(result, "stdout"));
      }
      if (!is.null(attr(result, "stderr")) && nzchar(attr(result, "stderr"))) {
        msg("[callr stderr] ", attr(result, "stderr"));
      }

      cem <- result
    } else {
      msg("Running in-process (callr disabled or not available)");
      cem <- cemitool(expr              = expr_mat,
                      annot             = annot_df,
                      filter            = filter,
                      min_ngen          = min_ngen,
                      cor_method        = cor_method,
                      class_column      = classCol,
                      verbose           = verbose,
                      plot              = FALSE,           # Disable all plotting
                      plot_diagnostics  = FALSE)
    }


    ## 3 · save & return -----------------------------------------------
    mod <- attr(cem, "module")
    expr_check <- attr(cem, "expression")

    # CRITICAL FIX: Ensure expression matrix is preserved
    # Some CEMiTool versions may not properly set the expression attribute
    if (is.null(attr(cem, "expression"))) {
      attr(cem, "expression") <- expr_mat
    }

    # NOTE: Trimming is DISABLED to prevent losing the expression matrix
    # tryCatch({
    #   cem <- .trim_cem_object_for_save(cem)
    # }, error = function(e) {
    #   msg("Warning: Could not trim CEMiTool object: ", conditionMessage(e))
    # })

    qs::qsave(cem, "cem.qs")

    n_samples <- nrow(cem@sample_annotation)
    n_genes <- nrow(expr_mat)

    if (is.null(mod) || !is.data.frame(mod) || nrow(mod) == 0 || !("modules" %in% colnames(mod))) {
      # No modules found - provide detailed diagnostic information
      msg("No modules found. Samples: ", n_samples, ", Genes: ", n_genes);

      error_msg <- paste0(
        "ERROR: No modules found. CEMiTool could not identify co-expression modules.\n\n",
        "Possible causes:\n",
        "1. Beta parameter selection failed (network not scale-free)\n",
        "2. Not enough genes (current: ", n_genes, ", try with more genes)\n",
        "3. Not enough samples (current: ", n_samples, ", recommend 10+ samples)\n",
        "4. Low variance in gene expression (genes too similar)\n",
        "5. Strong batch effects masking biological signal\n\n",
        "Suggestions:\n",
        "- Reduce filtering (lower min_ngen from ", min_ngen, ")\n",
        "- Set filter=FALSE to skip CEMiTool's internal filtering\n",
        "- Use more genes (include more features from your dataset)\n",
        "- Check if samples cluster by biological groups (PCA plot)\n",
        "- Try different correlation method (spearman vs pearson)"
      )

      return(error_msg)
    } else {
      n_modules <- length(unique(mod$modules))
      msg("Network construction successful: ", n_modules, " modules found");
      return("OK")
    }
  }, error = function(e) {
    # Capture detailed error information
    err_msg <- conditionMessage(e)
    err_trace <- paste(capture.output(traceback()), collapse="\n")

    # Log the error
    message("==== CEMiTool Error Details ====")
    message("Error message: ", err_msg)
    message("Error class: ", class(e))
    message("Call: ", deparse(e$call))
    message("\nTraceback:")
    message(err_trace)
    message("================================")

    # Return a detailed error message
    return(paste0("Error: ", err_msg,
                  " [Location: ", if(!is.null(e$call)) deparse(e$call)[1] else "unknown", "]"))
  })
}


PlotCEMiDendro <- function(mode      = c("sample", "module"),
                           metaClass = "NA",
                           imgName   = "cem_dendro",
                           dpi       = 72,
                           format    = "png") {

  library(Cairo); library(WGCNA)

  # FIX: Suppress Quartz popup on macOS - override device at function start
  old_device <- getOption("device")
  old_bitmapType <- getOption("bitmapType")
  options(device = function(...) pdf(file = NULL))
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    options(bitmapType = "cairo")
  }
  on.exit({
    options(device = old_device, bitmapType = old_bitmapType)
  }, add = TRUE)

  cem <- qs::qread("cem.qs")
  if (!inherits(cem, "CEMiTool"))
    stop("'cem.qs' does not contain a valid CEMiTool object.")

  mode <- match.arg(mode)
  expr <- attr(cem, "expression")       # features × samples
  mod <- attr(cem, "module")

  ## helper ----------------------------------------------------------
  plotDendroColoured <- function(hc, colNamed, label, file, legendPal) {

    leaves <- labels(as.dendrogram(hc))
    if (!all(leaves %in% names(colNamed)))
      stop("Colour vector missing leaves:\n  ",
           paste(setdiff(leaves, names(colNamed)), collapse = ", "))

    colMat <- matrix(colNamed[leaves], nrow = length(leaves), ncol = 1,
                     dimnames = list(leaves, NULL))

    if (dpi == 72) dpi <- 96
    width_in <- 10
    height_in <- 6

    # FIX: Suppress Quartz popup on macOS - close any existing devices first
    while (dev.cur() > 1) dev.off()

    Cairo(file, width = width_in, height = height_in, dpi = dpi,
          bg = "white", type = format, units = "in")

    oldMar <- par("mar"); par(mar = oldMar + c(0, 0, 0, 4))
    plotDendroAndColors(
      dendro       = hc,
      colors       = colMat,
      groupLabels  = label,
      dendroLabels = FALSE,
      addGuide     = TRUE,
      guideHang    = 0.05,
      main         = paste("Module dendrogram"),
      ylab         = "1 − Pearson correlation")

    ## legend
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), xpd = NA)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("topright",
           legend = names(legendPal),
           fill   = legendPal,
           border = NA,
           cex    = 0.8,
           bty    = "n")

    # FIX: Suppress Quartz popup on macOS
    par(oldMar)
    invisible(dev.off())
  }

  # ── MODULE dendrogram ────────────────────────────────────────────
  if (mode == "module") {
    # Detect feature column name (genes vs features)
    feature_col_dendro <- if ("genes" %in% colnames(mod)) "genes" else "features";

    ME  <- moduleEigengenes(t(expr[mod[[feature_col_dendro]], ]),
                            colors = mod$modules)$eigengenes
    hc  <- hclust(as.dist(1 - cor(ME)), method = "average")

    ids   <- colnames(ME)
    pal   <- setNames(rainbow(length(ids)), ids)
    idVec <- setNames(ids, ids)

    file <- sprintf("%s_module_dendro_dpi%d.%s", imgName, dpi, format)
    plotDendroColoured(hc, idVec, "modules", file, pal)
    return(1)
  }

  # ── SAMPLE dendrogram ────────────────────────────────────────────
  sa <- cem@sample_annotation

  ## choose metadata column sensibly
  if (is.na(metaClass) || metaClass == "NA") metaClass <- 2  # 1 is SampleName
  if (is.numeric(metaClass)) {
    if (metaClass < 1 || metaClass > ncol(sa))
      stop("'metaClass' index out of range.")
    if (metaClass == 1)
      stop("metaClass 1 is 'SampleName'; choose a metadata column (>=2).")
    classes <- sa[[metaClass]]
  } else {
    if (!metaClass %in% colnames(sa))
      stop("metaClass '", metaClass, "' not found in sample_annotation.")
    if (metaClass == "SampleName")
      stop("metaClass 'SampleName' is invalid; choose real metadata.")
    classes <- sa[[metaClass]]
  }

  names(classes) <- sa$SampleName

  hc  <- hclust(as.dist(1 - cor(expr)), method = "average")
  pal <- setNames(rainbow(length(unique(classes))), unique(classes))
  colNamed <- setNames(pal[classes], names(classes))

  file <- sprintf("%sdpi%d.%s", imgName, dpi, format)
  plotDendroColoured(hc, colNamed, "sample class", file, pal)
  imgSet <- readSet(imgSet, "imgSet");
  imgSet$coexp_dendrogram <- file;
  saveSet(imgSet, "imgSet");
  return("OK")
}

# =======================================================================
# PlotCEMiTreatmentHeatmap
# -----------------------------------------------------------------------
# factorName : name (character) or index (numeric) of the categorical
#              column in cem@sample_annotation to expand into dummies
# imgName    : file stem for the output
# dpi        : resolution
# format     : "png" | "pdf"
# -----------------------------------------------------------------------
# returns 1 on success, 0 on failure
# =======================================================================
PlotCEMiTreatmentHeatmap <- function(factorName,
                                     imgName = "cem_treatment_heatmap",
                                     dpi     = 96,
                                     format  = c("png", "pdf")) {

  tryCatch({

    library(CEMiTool); library(WGCNA); library(Cairo)

    # FIX: Suppress Quartz popup on macOS - override device at function start
    old_device <- getOption("device")
    old_bitmapType <- getOption("bitmapType")
    options(device = function(...) pdf(file = NULL))
    if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
      options(bitmapType = "cairo")
    }
    on.exit({
      options(device = old_device, bitmapType = old_bitmapType)
    }, add = TRUE)

    cem <- qs::qread("cem.qs")
    stopifnot(inherits(cem, "CEMiTool"))

    sa <- cem@sample_annotation
    msg("Sample annotation has ", nrow(sa), " samples and ", ncol(sa), " columns: ", paste(colnames(sa), collapse=", "));

    ## ── 1 · validate factorName  ----------------------------------
    # Default: use first metadata column (column 2, since column 1 is SampleName)
    if (is.na(factorName) || is.null(factorName) || factorName == "NA" || factorName == "") {
      if (ncol(sa) < 2) {
        stop("Sample annotation has no metadata columns (only SampleName)");
      }
      colLabel <- colnames(sa)[2]
      fac <- sa[[2]]
    } else if (is.numeric(factorName)) {
      if (factorName < 2 || factorName > ncol(sa)) {
        stop("Numeric factorName (", factorName, ") is out of range. Must be between 2 and ", ncol(sa),
             ". Available columns: ", paste(colnames(sa), collapse=", "));
      }
      fac      <- sa[[factorName]]
      colLabel <- colnames(sa)[factorName]
    } else {

      # Validate column exists - try case-insensitive match
      available_cols <- colnames(sa)

      if (factorName %in% available_cols) {
        # Exact match
        fac      <- sa[[factorName]]
        colLabel <- factorName
      } else {
        # Try case-insensitive match
        matched_col <- available_cols[tolower(available_cols) == tolower(factorName)]

        if (length(matched_col) > 0) {
          fac <- sa[[matched_col[1]]]
          colLabel <- matched_col[1]
        } else {
          # No match - provide helpful error
          stop("factorName '", factorName, "' not found in sample annotation (case-insensitive search failed). ",
               "Available columns: ", paste(available_cols, collapse=", "), ". ",
               "Please check that your metadata column name matches exactly.");
        }
      }

      if (colLabel == "SampleName") {
        stop("factorName cannot be 'SampleName'. Available columns: ",
               paste(available_cols[available_cols != "SampleName"], collapse=", "));
      }
    }

    fac <- as.factor(fac)

    ## ── 2 · dummy matrix  -----------------------------------------
    mm <- model.matrix(~ 0 + fac)
    colnames(mm) <- levels(fac)
    rownames(mm) <- sa$SampleName

    ## ── 3 · module eigenfeatures  ------------------------------------
    expr   <- attr(cem, "expression")
    modTbl <- attr(cem, "module")

    if (is.null(modTbl) || !is.data.frame(modTbl) || nrow(modTbl) == 0) {
      stop("Module table is NULL or empty - no modules found in CEMiTool object");
    }


    # CRITICAL FIX: CEMiTool uses different column names in different versions
    # Some use "genes", some use "features"
    feature_col <- if ("genes" %in% colnames(modTbl)) {
      "genes"
    } else if ("features" %in% colnames(modTbl)) {
      "features"
    } else {
      stop("Module table has no 'genes' or 'features' column. Available: ",
           paste(colnames(modTbl), collapse=", "));
    }


    g      <- intersect(rownames(expr), modTbl[[feature_col]])

    if (length(g) == 0) {
      stop("No overlapping features between expression matrix and module table");
    }

    colors <- modTbl$modules[match(g, modTbl[[feature_col]])]
    MEs    <- moduleEigengenes(t(expr[g, , drop = FALSE]), colors)$eigengenes

    mm <- mm[rownames(MEs), , drop = FALSE]

    ## ── 4 · correlations  -----------------------------------------
    corMat <- cor(MEs, mm, use = "p")
    pMat   <- corPvalueStudent(corMat, nSamples = nrow(MEs))
    textMat <- paste0(formatC(corMat, 2), "\n(",
                      formatC(pMat , 1, format = "e"), ")")

    ## ── 5 · device  (min 8 × 6 in)  -------------------------------
    colfun <- colorRampPalette(c("royalblue4", "white", "tomato"))

    outFile <- sprintf("%sdpi%d.%s", imgName, dpi, match.arg(format))
    if(dpi == 72){
        dpi <- 96;
    }
    width_px  <- max(50 + 40 * ncol(corMat), 8 * dpi)
    height_px <- max(50 + 20 * nrow(corMat), 6 * dpi)
    width_in  <- width_px  / dpi
    height_in <- height_px / dpi


    # FIX: Suppress Quartz popup on macOS - close any existing devices first
    while (dev.cur() > 1) dev.off()

    if (tolower(format) == "png") {
      Cairo(file   = outFile,
            width  = width_in,
            height = height_in,
            dpi    = dpi,
            bg     = "white",
            type   = "png",
            units  = "in")
    } else {
      Cairo(file   = outFile,
            width  = width_in,
            height = height_in,
            bg     = "white",
            type   = "pdf")
    }


    ## optional: tighten default margins a little
    oldMar <- par("mar"); on.exit(par(oldMar), add = TRUE)
    par(mar = c(5, 9, 4, 2) + 0.1)

    labeledHeatmap(Matrix          = corMat,
                   xLabels         = colnames(mm),
                   yLabels         = rownames(corMat),
                   ySymbols        = rownames(corMat),
                   colorLabels     = FALSE,
                   colors          = colfun(50),
                   textMatrix      = textMat,
                   setStdMargins   = FALSE,
                   cex.text        = 0.7,
                   zlim            = c(-1, 1),
                   main            = paste("Module ×", colLabel, "Levels"))

    # FIX: Suppress Quartz popup on macOS
    invisible(dev.off())

    imgSet <- readSet(imgSet, "imgSet");
    imgSet$coexp_traitheat <- outFile;
    saveSet(imgSet, "imgSet");

    1

  }, error = function(e) {
    message("ERROR in PlotCEMiTreatmentHeatmap: ", conditionMessage(e));
    0
  })
}

PlotCemiScaleFree <- function(imgName = "coexp_scalefree",
                                     dpi = 72,
                                     format = "png") {
  library(Cairo); library(CEMiTool)

  # FIX: Suppress Quartz popup on macOS - override device at function start
  old_device <- getOption("device")
  old_bitmapType <- getOption("bitmapType")
  options(device = function(...) pdf(file = NULL))
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    options(bitmapType = "cairo")
  }
  on.exit({
    options(device = old_device, bitmapType = old_bitmapType)
  }, add = TRUE)

  cem <- qs::qread("cem.qs")
  stopifnot(inherits(cem, "CEMiTool"))

  # Ensure the plot exists (some versions only populate it after calling plot_beta)
  if (is.null(cem@beta_r2_plot)) {
    plot_beta <- try(getFromNamespace("plot_beta", "CEMiTool"), silent = TRUE)
    if (!inherits(plot_beta, "try-error") && is.function(plot_beta)) {
      cem <- plot_beta(cem)  # fills cem@beta_r2_plot
    }
  }

  # Extract ggplot object from the slot
  g <- NULL
  if (!is.null(cem@beta_r2_plot)) {
    if (is.list(cem@beta_r2_plot) && "beta_r2_plot" %in% names(cem@beta_r2_plot)) {
      g <- cem@beta_r2_plot$beta_r2_plot
    } else {
      g <- cem@beta_r2_plot
    }
  }
  if (is.null(g)) return("Error: beta_r2_plot not available on the CEMiTool object.")

  # Save
  file <- sprintf("%sdpi%d.%s", imgName, dpi, format)
    if (dpi == 72) dpi <- 96
  width_in <- 10
  height_in <- 6

    # FIX: Suppress Quartz popup on macOS - close any existing devices first
    while (dev.cur() > 1) dev.off()

    Cairo(file, width = width_in, height = height_in, dpi = dpi, bg = "white", type = format, units = "in")
  invisible(print(g))    # ggplot draw
  invisible(dev.off())
    imgSet <- readSet(imgSet, "imgSet");
    imgSet$coexp_scalefree <- file;
    saveSet(imgSet, "imgSet");
  return(1);
}
