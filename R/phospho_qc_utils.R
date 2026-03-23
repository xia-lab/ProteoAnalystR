##################################################
## R script for ProteoAnalyst
## Description: Phosphoproteomics-specific QC plots
## Authors:
## Generated for phospho module implementation
###################################################

#' Plot distribution of localization probabilities
#' Shows histogram with Class I (>0.75), II, and III sites
PlotLocalizationProbDist <- function(fileName, imgNm, dpi = 72, format = "png") {

  # Ensure no interactive devices on macOS
  if (Sys.info()['sysname'] == 'Darwin') {
      options(bitmapType = 'cairo')
      options(device = function(...) grDevices::png(..., type = "cairo"))
  }

  dataSet <- readDataset(fileName)
  paramSet <- readSet(paramSet, "paramSet")

  # Check if this is phospho data
  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") {
    #msg("[PlotLocalizationProbDist] Not phospho data, skipping.")
    return("NA")
  }

  # Check if feature.info has localization probability
  if (is.null(dataSet$feature.info) || !("Localization prob" %in% colnames(dataSet$feature.info))) {
    #msg("[PlotLocalizationProbDist] Localization probability not available.")
    return("NA")
  }

  loc_probs <- dataSet$feature.info$`Localization prob`

  # Create plot
  qc.loc.prob.hist(loc_probs, imgNm, dpi, format)

  return("NA")
}

#' Plot breakdown of phosphosite classes
#' Shows STY residue distribution and missing values per residue type
PlotSiteClassBreakdown <- function(fileName, imgNm, dpi = 72, format = "png") {

  # Ensure no interactive devices on macOS
  if (Sys.info()['sysname'] == 'Darwin') {
      options(bitmapType = 'cairo')
      options(device = function(...) grDevices::png(..., type = "cairo"))
  }

  dataSet <- readDataset(fileName)
  paramSet <- readSet(paramSet, "paramSet")

  # Check if this is phospho data
  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") {
    #msg("[PlotSiteClassBreakdown] Not phospho data, skipping.")
    return("NA")
  }

  # Check if feature.info has amino acid information
  if (is.null(dataSet$feature.info) || !("Amino acid" %in% colnames(dataSet$feature.info))) {
    #msg("[PlotSiteClassBreakdown] Amino acid information not available.")
    return("NA")
  }

  amino_acids <- dataSet$feature.info$`Amino acid`
  data_mat <- as.matrix(dataSet$data.norm)
  if (file.exists("data.raw.qs")) {
    raw_mat <- try(readDataQs("data.raw.qs", paramSet$anal.type, fileName), silent = TRUE)
    if (!inherits(raw_mat, "try-error") && !is.null(raw_mat)) {
      data_mat <- as.matrix(raw_mat)
    }
  }

  # Create plot
  qc.site.class.breakdown(amino_acids, data_mat, imgNm, dpi, format)

  return("NA")
}

#' Plot QC panel (quantify) for phosphoproteomics data (disabled)
PlotPhosphoQC <- function(fileName, imgNm, dpi = 96, format = "png", panel = "quantify") {
  return("NA")
}

#' Internal function: Localization probability histogram
#' Runs in quiet mode (backend) without triggering Quartz/GUI
qc.loc.prob.hist <- function(loc_probs, imgNm, dpi = 96, format = "png") {

  # 1. SAFETY: Force non-interactive graphics engine for this function scope
  # This prevents any accidental GUI calls on macOS Rserve
  if (Sys.info()['sysname'] == 'Darwin') {
      old_bitmap <- getOption("bitmapType")
      options(bitmapType = 'cairo')
      on.exit(options(bitmapType = old_bitmap), add = TRUE) # Restore option when done
  }

  require('ggplot2')

  # Ensure no interactive devices are opened
  options(device = function(...) grDevices::png(..., type = "cairo"))

  dpi <- as.numeric(dpi)
  finalFileNm <- paste0(imgNm, "dpi", dpi, ".", format)

  # Remove NA values
  loc_probs <- loc_probs[!is.na(loc_probs)]

  if (length(loc_probs) == 0) {
    #msg("[qc.loc.prob.hist] No localization probability data available.")
    return(NULL)
  }

  # Classify sites
  class_I <- sum(loc_probs >= 0.75)
  class_II <- sum(loc_probs >= 0.5 & loc_probs < 0.75)
  class_III <- sum(loc_probs < 0.5)

  # Create data frame for plotting
  df <- data.frame(
    loc_prob = loc_probs
  )

  # Create histogram with class boundaries
  # The plot object 'p' is created in memory but NOT drawn to screen
  p <- ggplot(df, aes(x = loc_prob)) +
    geom_histogram(binwidth = 0.05, fill = "#3498db", color = "white", alpha = 0.8) +
    geom_vline(xintercept = 0.75, linetype = "dashed", color = "#e74c3c", size = 1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "#f39c12", size = 1) +
    annotate("text", x = 0.875, y = Inf, label = paste0("Class I\n(≥0.75)\n", class_I, " sites"),
             vjust = 1.5, hjust = 0.5, size = 3.5, color = "#2c3e50") +
    annotate("text", x = 0.625, y = Inf, label = paste0("Class II\n(0.5-0.75)\n", class_II, " sites"),
             vjust = 1.5, hjust = 0.5, size = 3.5, color = "#2c3e50") +
    labs(
      title = "Phosphosite Localization Probability Distribution",
      x = "Localization Probability",
      y = "Number of Sites"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )

  # Save plot directly to file using safe devices
  if (format == "png") {
    # type="cairo" avoids the Quartz backend on macOS
    ggsave(finalFileNm, plot = p, dpi = dpi, width = 8, height = 6, bg = "white", type = "cairo")
  } else if (format == "pdf") {
    ggsave(finalFileNm, plot = p, width = 8, height = 6, device = cairo_pdf, bg = "white")
  }

  return(finalFileNm)
}

#' Internal function: Site class (STY) breakdown
# fixed for MacOS
#' Internal function: Site class (STY) breakdown
qc.site.class.breakdown <- function(amino_acids, data_mat, imgNm, dpi = 96, format = "png") {
  # FIX 1: Ensure headless graphics on macOS to prevent Rserve crash
  if (Sys.info()['sysname'] == 'Darwin') {
      options(bitmapType = 'cairo')
  }

  require('ggplot2')
  require('gridExtra')

  # Ensure no interactive devices are opened
  options(device = function(...) grDevices::png(..., type = "cairo"))

  dpi <- as.numeric(dpi)
  # Cleaned up filename generation
  finalFileNm <- paste0(imgNm, "dpi", dpi, ".", format)

  # Count residue types
  aa_counts <- table(amino_acids)

  # Create data frame for pie chart
  df_pie <- data.frame(
    Residue = names(aa_counts),
    Count = as.numeric(aa_counts),
    Percentage = round(100 * as.numeric(aa_counts) / sum(aa_counts), 1)
  )

  df_pie$Label <- paste0(df_pie$Residue, "\n", df_pie$Count, " (", df_pie$Percentage, "%)")

  # Color palette for STY
  colors <- c("S" = "#3498db", "T" = "#e74c3c", "Y" = "#2ecc71")

  # Pie chart
  p1 <- ggplot(df_pie, aes(x = "", y = Count, fill = Residue)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors, labels = df_pie$Label) +
    labs(title = "Phosphosite Residue Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )

  # Calculate missing value percentage per residue type
  missing_per_residue <- sapply(unique(amino_acids), function(aa) {
    aa_rows <- which(amino_acids == aa)
    if (length(aa_rows) == 0) return(0)
    aa_data <- data_mat[aa_rows, , drop = FALSE]
    
    # Counts NA and Zeros
    missing_pct <- 100 * sum(is.na(aa_data) | aa_data == 0) / length(aa_data)
    
    return(missing_pct)
  })

  df_missing <- data.frame(
    Residue = names(missing_per_residue),
    Missing_Pct = as.numeric(missing_per_residue)
  )

  # Bar chart for missing values
  p2 <- ggplot(df_missing, aes(x = Residue, y = Missing_Pct, fill = Residue)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = colors) +
    labs(
      title = "Missing Values by Residue Type",
      x = "Residue",
      y = "Missing Values (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none"
    )

  # FIX 2: Use arrangeGrob instead of grid.arrange
  # grid.arrange() draws to the screen immediately (causing the crash).
  # arrangeGrob() saves the layout to an object WITHOUT drawing.
  combined_plot <- arrangeGrob(p1, p2, ncol = 2)

  # FIX 3: Explicitly define the device type for ggsave
  if (format == "png") {
    # Use type="cairo" to force non-GUI rendering if supported, or rely on global option
    ggsave(finalFileNm, plot = combined_plot, dpi = dpi, width = 12, height = 5, bg = "white", type = "cairo")
  } else if (format == "pdf") {
    ggsave(finalFileNm, plot = combined_plot, width = 12, height = 5, device = cairo_pdf, bg = "white")
  }

  return(finalFileNm)
}


PlotProteinCV <- function(fileName, imgNm, dpi = 72, format = "png") {

  # Ensure no interactive devices on macOS
  if (Sys.info()['sysname'] == 'Darwin') {
      options(bitmapType = 'cairo')
      options(device = function(...) grDevices::png(..., type = "cairo"))
  }

  dataSet <- readDataset(fileName)
  paramSet <- readSet(paramSet, "paramSet")

  # Check if grouping information is available
  if (is.null(dataSet$meta.info) && is.null(dataSet$cls)) {
    return("NA")
  }

  # Determine groups (conditions)
  # Prefer 'cls' if available, otherwise use the first column of meta.info
  if (!is.null(dataSet$cls)) {
    groups <- dataSet$cls
  } else {
    groups <- dataSet$meta.info[, 1]
  }
  groups <- as.factor(groups)

  # Extract data matrix
  # Use raw data for pre-normalization CV; use normalized for post-normalization.
  data_mat <- as.matrix(dataSet$data.norm)
  if (!grepl("norm", imgNm, ignore.case = TRUE) && file.exists("data.raw.qs")) {
    raw_mat <- try(readDataQs("data.raw.qs", paramSet$anal.type, fileName), silent = TRUE)
    if (!inherits(raw_mat, "try-error") && !is.null(raw_mat)) {
      data_mat <- as.matrix(raw_mat)
    }
  }

  # Create plot
  qc.protein.cv.hist(data_mat, groups, imgNm, dpi, format)

  return("NA")
}

#' Internal function: Protein CV Distribution Histogram
qc.protein.cv.hist <- function(data_mat, groups, imgNm, dpi = 96, format = "png") {
  # FIX 1: Ensure headless graphics on macOS to prevent Rserve crash
  if (Sys.info()['sysname'] == 'Darwin') {
      options(bitmapType = 'cairo')
  }

  require('ggplot2')
  require('gridExtra')

  # Ensure no interactive devices are opened
  options(device = function(...) grDevices::png(..., type = "cairo"))

  dpi <- as.numeric(dpi)
  finalFileNm <- paste0(imgNm, "dpi", dpi, ".", format)

  # Initialize data storage
  unique_groups <- unique(groups)
  df_list <- list()

  # Calculate CV for each group
  for (grp in unique_groups) {
    # Get sample indices for this group
    samples_idx <- which(groups == grp)
    
    # Need at least 2 samples to calculate variance/CV
    if (length(samples_idx) < 2) {
      next 
    }

    # Subset data
    sub_data <- data_mat[, samples_idx, drop = FALSE]
    finite_vals <- sub_data[is.finite(sub_data)]
    if (length(finite_vals) > 0 && stats::median(finite_vals, na.rm = TRUE) < 100) {
      sub_data <- 2^sub_data
    }

    # Calculate Mean and SD per protein (row)
    # Using rowMeans/apply is more efficient for matrices
    p_means <- rowMeans(sub_data, na.rm = TRUE)
    p_sds <- apply(sub_data, 1, sd, na.rm = TRUE)

    # Calculate CV (%)
    # Handle division by zero if mean is 0
    p_cvs <- (p_sds / p_means) * 100

    # Filter out infinite/NaN values and zero means
    valid_idx <- which(is.finite(p_cvs) & p_means != 0)

    if (length(valid_idx) > 0) {
      temp_df <- data.frame(
        CV = p_cvs[valid_idx],
        Condition = as.character(grp)
      )
      df_list[[grp]] <- temp_df
    }
  }

  # Combine all groups into one data frame
  if (length(df_list) == 0) {
      finite_vals <- data_mat[is.finite(data_mat)]
      if (length(finite_vals) > 0 && stats::median(finite_vals, na.rm = TRUE) < 100) {
        data_mat <- 2^data_mat
      }
      p_means <- rowMeans(data_mat, na.rm = TRUE)
      p_sds <- apply(data_mat, 1, sd, na.rm = TRUE)
      p_cvs <- (p_sds / p_means) * 100
      valid_idx <- which(is.finite(p_cvs) & p_means != 0)
      if (length(valid_idx) == 0) {
        return("NA")
      }
      df_all <- data.frame(CV = p_cvs[valid_idx], Condition = "All samples")
  } else {
      df_all <- do.call(rbind, df_list)
  }

  # Calculate Median CV for vertical lines
  medians <- aggregate(CV ~ Condition, data = df_all, median)
  medians$Label <- round(medians$CV, 1)

  # Create Plot
  p <- ggplot(df_all, aes(x = CV)) +
    geom_histogram(fill = "#56B4E9", color = "white", bins = 40, alpha = 0.8) +
    geom_vline(data = medians, aes(xintercept = CV), color = "red", linetype = "dashed", size = 0.8) +
    geom_text(data = medians, aes(x = CV, y = Inf, label = Label), 
              vjust = 1.5, hjust = -0.2, color = "red", size = 3.5, fontface = "bold") +
    facet_wrap(~Condition, scales = "fixed") +
    labs(
      title = "Distribution of Protein CV per Condition",
      x = "Coefficient of Variation (%)",
      y = "Frequency"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#e0e0e0")
    )

  # FIX 3: Explicitly define the device type for ggsave
  if (format == "png") {
    ggsave(finalFileNm, plot = p, dpi = dpi, width = 10, height = 6, bg = "white", type = "cairo")
  } else if (format == "pdf") {
    ggsave(finalFileNm, plot = p, width = 10, height = 6, device = cairo_pdf, bg = "white")
  }

  return(finalFileNm)
}
