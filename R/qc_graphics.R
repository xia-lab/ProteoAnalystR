##################################################
## R script for ProteoAnalyst
## Description: functions for quality check boxplot
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

PlotDataBox <- function(fileName, boxplotName, dpi, format){
  dataSet <- readDataset(fileName);
  qc.boxplot(as.matrix(dataSet$data.norm), boxplotName, dpi, format, F, meta = dataSet$meta.info);
  return("NA");
}

PlotQcOverview <- function(fileName, imgNm, dpi, format){
  dataSet <- readDataset(fileName);
  qc.overview.patchwork(as.matrix(dataSet$data.norm), imgNm, dpi, format, meta = dataSet$meta.info);
  return("NA");
}

PlotDataMA <- function(fileName, imgNm, dpi, format, sampleIndex=1){
  dataSet <- readDataset(fileName);
  qc.maplot(as.matrix(dataSet$data.norm), imgNm, dpi, format, F, meta = dataSet$meta.info, sampleIndex = sampleIndex, mode="single");
  return("NA");
}

PlotDataMAInteractive <- function(fileName, imgNm){
  dataSet <- readDataset(fileName);
  qc.maplot.json(as.matrix(dataSet$data.norm), imgNm, meta = dataSet$meta.info);
  return("NA");
}

PlotDataMAOverview <- function(fileName, imgNm, dpi, format){
  dataSet <- readDataset(fileName);
  qc.maplot(as.matrix(dataSet$data.norm), imgNm, dpi, format, F, meta = dataSet$meta.info, mode="overview");
  return("NA");
}

coerce.numeric.matrix <- function(dat) {
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  }

  if (is.null(dat)) {
    return(dat)
  }

  rn <- rownames(dat)
  cn <- colnames(dat)
  dat <- as.matrix(dat)
  suppressWarnings(storage.mode(dat) <- "numeric")
  rownames(dat) <- rn
  colnames(dat) <- cn
  dat
}

get.qc.sample.matrix <- function(dataSet, imgNm = NULL) {
  if (!is.null(imgNm) && grepl("norm", imgNm)) {
    data.for.plot <- dataSet$data.norm;
  } else if (file.exists("data.with.missing.qs")) {
    data.for.plot <- qs::qread("data.with.missing.qs");
  } else {
    data.for.plot <- dataSet$data.norm;
  }

  if (is.data.frame(data.for.plot)) {
    data.for.plot <- as.matrix(data.for.plot);
  }

  suppressWarnings(storage.mode(data.for.plot) <- "numeric");
  data.for.plot
}

build.sample.qc.table <- function(dat, meta = NULL) {
  if (is.data.frame(dat)) dat <- as.matrix(dat)

  sample_names <- colnames(dat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    colnames(dat) <- sample_names
  }

  miss.ind <- is.na(dat) | dat == "NA"
  if (is.numeric(dat)) {
    miss.ind <- is.na(dat) | dat <= 0
  }

  miss.count <- colSums(miss.ind)
  n.features <- nrow(dat)
  nonmiss.pct <- 100 * (n.features - miss.count) / n.features
  num.detected <- n.features - miss.count

  # For mean calculation, exclude missing values (NA and <=0)
  # by converting them to NA first
  dat_valid <- dat
  dat_valid[miss.ind] <- NA

  mean.int <- apply(dat_valid, 2, mean, na.rm = TRUE)

  primary.group <- rep("Group1", length(sample_names))
  if (!is.null(meta) && ncol(meta) >= 1) {
    meta.df <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (!is.null(rownames(meta.df))) {
      matched <- meta.df[sample_names, 1, drop = TRUE]
      matched[is.na(matched) | matched == ""] <- "NA"
      primary.group <- as.character(matched)
    }
  }

  tukey.flag <- function(x, reverse = FALSE) {
    x <- as.numeric(x)
    finite <- is.finite(x)
    out <- rep(FALSE, length(x))
    if (sum(finite) < 4) return(out)
    q1 <- stats::quantile(x[finite], 0.25, na.rm = TRUE, names = FALSE)
    q3 <- stats::quantile(x[finite], 0.75, na.rm = TRUE, names = FALSE)
    iqr <- q3 - q1
    if (!is.finite(iqr) || iqr == 0) return(out)
    lower <- q1 - 1.5 * iqr
    upper <- q3 + 1.5 * iqr
    if (reverse) {
      out[finite] <- x[finite] < lower
    } else {
      out[finite] <- x[finite] < lower | x[finite] > upper
    }
    out
  }

  outlier.flag <- tukey.flag(nonmiss.pct, reverse = TRUE) |
    tukey.flag(num.detected, reverse = TRUE) |
    tukey.flag(mean.int)

  data.frame(
    Sample = sample_names,
    Group = primary.group,
    NonMissingPct = round(as.numeric(nonmiss.pct), 1),
    NumDetected = as.integer(num.detected),
    MeanIntensity = signif(as.numeric(mean.int), 4),
    OutlierFlag = ifelse(outlier.flag, "Review", "OK"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

get.primary.group.values <- function(sample_names, meta = NULL) {
  primary.group <- rep("Group1", length(sample_names))
  if (!is.null(meta) && ncol(meta) >= 1) {
    meta.df <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (!is.null(rownames(meta.df))) {
      matched <- meta.df[sample_names, 1, drop = TRUE]
      matched[is.na(matched) | matched == ""] <- "NA"
      primary.group <- as.character(matched)
    }
  }
  stats::setNames(primary.group, sample_names)
}

get.qc.group.palette <- function(group_values) {
  unique_groups <- unique(as.character(group_values))
  n_groups <- length(unique_groups)
  # Reordered to maximize contrast among the first several categorical groups.
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                 "#D55E00", "#CC79A7", "#F0E442", "#999999")

  if (n_groups <= length(okabe_ito)) {
    group.colors <- okabe_ito[seq_len(n_groups)]
  } else {
    group.colors <- grDevices::colorRampPalette(okabe_ito)(n_groups)
  }

  stats::setNames(as.character(group.colors), unique_groups)
}

get.qc.annotation.palette <- function(values, annotation_index = 1L, use_master_group_palette = FALSE) {
  values <- as.character(values)
  values[is.na(values) | values == ""] <- "NA"

  if (use_master_group_palette) {
    return(get.qc.group.palette(values))
  }

  annotation.palettes <- list(
    c("#1F77B4", "#2CA02C", "#17BECF", "#4DB6AC", "#6BAED6", "#31A354", "#9EDAE5", "#74C476"),
    c("#D62728", "#FF7F0E", "#9467BD", "#E377C2", "#8C564B", "#BCBD22", "#C49C94", "#F7B6D2"),
    c("#7F7F7F", "#BCBD22", "#8C564B", "#AEC7E8", "#FFBB78", "#98DF8A", "#C5B0D5", "#C7C7C7")
  )

  levs <- unique(values)
  base.palette <- annotation.palettes[[((annotation_index - 1L) %% length(annotation.palettes)) + 1L]]
  if (length(levs) <= length(base.palette)) {
    pal <- base.palette[seq_along(levs)]
  } else {
    pal <- grDevices::colorRampPalette(base.palette)(length(levs))
  }

  stats::setNames(as.character(pal), levs)
}

get.sample.dendrogram.order <- function(dat, linkage = "complete") {
  dat <- coerce.numeric.matrix(dat)
  if (is.null(dat) || ncol(dat) < 2) {
    sample_names <- colnames(dat)
    if (is.null(sample_names)) {
      sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    }
    return(sample_names)
  }

  suppressWarnings({
    cor.mat <- stats::cor(dat, use = "pairwise.complete.obs", method = "pearson")
  })
  if (is.null(cor.mat) || all(!is.finite(cor.mat))) {
    sample_names <- colnames(dat)
    if (is.null(sample_names)) {
      sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    }
    return(sample_names)
  }

  cor.mat[!is.finite(cor.mat)] <- 0
  diag(cor.mat) <- 1

  hc <- tryCatch(
    hclust(as.dist(1 - cor.mat), method = linkage),
    error = function(e) NULL
  )
  if (is.null(hc)) {
    sample_names <- colnames(dat)
    if (is.null(sample_names)) {
      sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    }
    return(sample_names)
  }

  colnames(cor.mat)[hc$order]
}

build.dendrogram.plot <- function(dat, group.colors, group.map) {
  dat <- coerce.numeric.matrix(dat)
  if (is.null(dat) || ncol(dat) < 2) {
    sample_names <- colnames(dat)
    if (is.null(sample_names)) {
      sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    }
    return(
      ggplot(data.frame(x = 1, y = 1, label = "Need >= 2 samples"), aes(x, y)) +
        geom_text(aes(label = label), size = 4, color = "#54616c") +
        theme_void()
    )
  }

  suppressWarnings({
    cor.mat <- stats::cor(dat, use = "pairwise.complete.obs", method = "pearson")
  })
  cor.mat[!is.finite(cor.mat)] <- 0
  diag(cor.mat) <- 1

  hc <- hclust(as.dist(1 - cor.mat), method = "complete")
  dend <- as.dendrogram(hc)
  dend_data <- ggdendro::dendro_data(dend, type = "rectangle")

  # Get leaf positions and map to group colors
  label_df <- dend_data$labels
  label_df <- label_df[order(label_df$x), , drop = FALSE]
  label_df$group <- group.map[as.character(label_df$label)]
  label_df$color <- group.colors[label_df$group]

  ggplot() +
    geom_segment(
      data = dend_data$segments,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.4,
      color = "#2f3e46",
      lineend = "round"
    ) +
    # Add colored points at leaf positions (bottom of dendrogram)
    geom_point(
      data = label_df,
      aes(x = x, y = 0),
      color = label_df$color,
      size = 2.5,
      shape = 15  # Square shape
    ) +
    scale_x_continuous(expand = expansion(add = 0.6)) +  # Match discrete scale default
    scale_y_reverse(expand = expansion(mult = c(0.02, 0.05))) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = "Distance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 9, color = "#233647"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "#233647", linewidth = 0.4),
      axis.line.y = element_blank(),
      plot.margin = margin(4, 0, 4, 2)
    )
}

build.violin.overview.plot <- function(dat, sample.order, group.colors, group.map) {
  dat <- coerce.numeric.matrix(dat)
  subgene <- 10000

  if (nrow(dat) > subgene) {
    set.seed(28051968)
    dat <- dat[sample(nrow(dat), subgene), , drop = FALSE]
  }

  sample_names <- colnames(dat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    colnames(dat) <- sample_names
  }

  # Use provided group mapping
  primary.group <- group.map[sample_names]

  df <- data.frame(
    values = as.numeric(dat),
    sample_id = rep(sample_names, each = nrow(dat)),
    Group = rep(primary.group, each = nrow(dat)),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$values), , drop = FALSE]
  df$sample_id <- factor(df$sample_id, levels = sample.order)
  df$Group <- factor(df$Group)

  q_vals <- stats::quantile(df$values, probs = c(0.01, 0.99), na.rm = TRUE)
  x_breaks <- pretty(unname(q_vals), n = 5)
  x_breaks <- x_breaks[x_breaks >= q_vals[[1]] & x_breaks <= q_vals[[2]]]
  if (length(x_breaks) > 1) {
    x_breaks <- x_breaks[-length(x_breaks)]
  }
  violin.df <- do.call(
    rbind,
    lapply(split(df, df$sample_id), function(subdf) {
      iqr_bounds <- stats::quantile(subdf$values, probs = c(0.25, 0.75), na.rm = TRUE)
      subdf[subdf$values >= iqr_bounds[[1]] & subdf$values <= iqr_bounds[[2]], , drop = FALSE]
    })
  )
  point.df <- do.call(
    rbind,
    lapply(split(df, df$sample_id), function(subdf) {
      iqr_bounds <- stats::quantile(subdf$values, probs = c(0.25, 0.75), na.rm = TRUE)
      subdf <- subdf[subdf$values < iqr_bounds[[1]] | subdf$values > iqr_bounds[[2]], , drop = FALSE]
      keep_n <- min(nrow(subdf), 900L)
      if (nrow(subdf) > keep_n) {
        subdf <- subdf[sample.int(nrow(subdf), keep_n), , drop = FALSE]
      }
      subdf
    })
  )

  p <- ggplot(df, aes(x = values, y = sample_id, fill = Group, color = Group)) +
    geom_violin(
      data = violin.df,
      scale = "width",
      trim = TRUE,
      adjust = 0.5,
      alpha = 0.88,
      linewidth = 0.35,
      width = 0.62
    ) +
    geom_point(
      data = point.df,
      aes(x = values, y = sample_id, color = Group),
      inherit.aes = FALSE,
      position = position_jitter(height = 0.085, width = 0),
      size = 0.6,
      alpha = 0.28,
      stroke = 0
    ) +
    geom_boxplot(
      width = 0.095,
      outlier.shape = NA,
      fill = "#fbfbfb",
      color = "#4f4f4f",
      linewidth = 0.3
    ) +
    stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 21,
      size = 1.8,
      stroke = 0.35,
      fill = "#ffffff",
      color = "#1f2d3d"
    ) +
    scale_x_continuous(
      limits = unname(q_vals),
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_fill_manual(values = group.colors) +
    scale_color_manual(values = group.colors) +
    labs(x = "Intensity", y = NULL) +
    guides(fill = guide_legend(title = "Group"), color = "none") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(size = 7.5, hjust = 1, color = "#233647"),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#e7ebee", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "#233647", linewidth = 0.4),
      axis.line.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(4, 2, 4, 0)
    )

  return(p)
}

build.sample.median.plot <- function(dat, sample.order, meta = NULL) {
  df <- build.sample.qc.table(dat, meta = meta)
  df <- df[match(sample.order, df$Sample), , drop = FALSE]
  df$Sample <- factor(df$Sample, levels = sample.order)
  df$Group <- factor(df$Group)

  ggplot(df, aes(x = Sample, y = MedianIntensity, fill = Group)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = get.qc.group.palette(as.character(df$Group))) +
    labs(title = "Sample median intensity", x = NULL, y = "Median intensity") +
    guides(fill = "none") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 8, 10, 8)
    )
}

qc.overview.patchwork <- function(dat, imgNm, dpi = 96, format = "png", meta = NULL) {
  dat <- coerce.numeric.matrix(dat)
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  fullPath <- paste0(fileNm, format, sep = "")

  if (is.null(dat) || ncol(dat) == 0 || nrow(dat) == 0) {
    return("NA")
  }

  if (!requireNamespace("patchwork", quietly = TRUE) ||
      !requireNamespace("ggdendro", quietly = TRUE)) {
    return("NA")
  }
  suppressMessages(require("ggplot2"))

  # Create group mapping (sample name -> group)
  sample_names <- colnames(dat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    colnames(dat) <- sample_names
  }

  group.map <- get.primary.group.values(sample_names, meta = meta)
  group.colors <- get.qc.group.palette(group.map)

  # Get sample order from dendrogram clustering
  sample.order <- get.sample.dendrogram.order(dat, linkage = "complete")

  # Build plots with shared sample order, colors, and group mapping
  dendro.plot <- build.dendrogram.plot(dat, group.colors, group.map)
  violin.plot <- build.violin.overview.plot(dat, sample.order, group.colors, group.map)

  overview.plot <- dendro.plot + violin.plot +
    patchwork::plot_layout(ncol = 2, widths = c(1.1, 1.9), guides = "collect") &
    theme(
      legend.position = "bottom",
      plot.margin = margin(2, 2, 2, 2)
    )

  sample_count <- length(sample.order)
  height_in <- max(6.5, min(16, 2.8 + sample_count * 0.22))
  width_in <- 12  # Reduced from 18 (2 panels instead of 3)

  if (dpi == 72) {
    dpi <- 96
  }

  if (format == "png" && capabilities("cairo")) {
    png(filename = fullPath, width = width_in * dpi, height = height_in * dpi, units = "px", res = dpi, type = "cairo")
  } else {
    Cairo(file = fullPath, width = width_in, height = height_in, unit = "in", dpi = dpi, type = format, bg = "white")
  }

  tryCatch({
    print(overview.plot)
  }, finally = {
    dev.off()
  })

  imgSet <- readSet(imgSet, "imgSet")
  if (grepl("norm", imgNm)) {
    imgSet$qc_norm_overview <- fullPath
  } else {
    imgSet$qc_overview <- fullPath
  }
  saveSet(imgSet, "imgSet")
  return("NA")
}

GetSampleQcTable <- function(dataName){
  dataSet <- readDataset(dataName)
  dat <- get.qc.sample.matrix(dataSet)
  df <- build.sample.qc.table(dat, meta = dataSet$meta.info)
  return(df)
}

# Missing value count per sample
PlotMissingPerSample <- function(fileName, imgNm, dpi, format){
  dataSet <- readDataset(fileName);

  # Determine which dataset to use based on image name:
  # - If imgNm contains "norm", we're AFTER normalization -> use data.norm
  # - Otherwise, we're BEFORE normalization -> use original data with missing values
  data.for.plot <- get.qc.sample.matrix(dataSet, imgNm);

  qc.nonmissing.per.sample(as.matrix(data.for.plot), imgNm, dpi, format, interactive = FALSE, meta = dataSet$meta.info);
  return("NA");
}

# Sample correlation heatmap
PlotSampleCorr <- function(fileName, imgNm, dpi, format, method = "pearson"){
  dataSet <- readDataset(fileName);
  qc.sample.corr.json(dataSet, imgNm, method = method);
  return("NA");
}

# Sample distance dendrogram
PlotSampleDendro <- function(fileName, imgNm, dpi, format, method = "correlation", linkage = "complete"){
  dataSet <- readDataset(fileName);
  qc.sample.dendro(as.matrix(dataSet$data.norm), imgNm, dpi, format, interactive = FALSE, method = method, linkage = linkage, meta = dataSet$meta.info);
  return("NA");
}

qc.boxplot <- function(dat, imgNm, dpi=96, format="png", interactive=F, meta = NULL){
  dpi <- as.numeric(dpi)
  require('ggplot2')
  require('lattice');
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  fullPath <- paste0(fileNm, format, sep="");  subgene <- 10000;

  if(class(dat)[1] == "data.frame"){
    dat <- as.matrix(dat);
  }

  full_sample_names <- colnames(dat)
  if (is.null(full_sample_names)) {
    full_sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    colnames(dat) <- full_sample_names
  }
  full_group_palette <- get.qc.group.palette(get.primary.group.values(full_sample_names, meta = meta))

  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  <- sample(nrow(dat), subgene);
    Mss <- dat[sg,,drop=FALSE];
  } else {
    Mss <- dat;
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  <- sample(ncol(Mss), subsmpl)
    Mss <- Mss[,ss,drop=FALSE]
  } else {
    Mss <- Mss
  }

  sample_names <- colnames(Mss)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_len(ncol(Mss)))
    colnames(Mss) <- sample_names
  }
  sample_order <- sort(sample_names, method = "radix")

  group.map <- get.primary.group.values(sample_names, meta = meta)

  df <- data.frame(
    values = as.numeric(Mss),
    sample_id = factor(rep(sample_names, each = nrow(Mss)), levels = sample_order),
    Group = factor(rep(group.map[sample_names], each = nrow(Mss)))
  )

  q_vals <- quantile(df$values, probs = c(0.01, 0.99), na.rm = TRUE)
  xlower <- unname(q_vals[1])
  xupper <- unname(q_vals[2])
  height <- length(unique(df$sample_id)) *20;
  if(height<450){
    height <- 450
  }
  
  bp <- ggplot(df, aes(sample_id, values, fill = Group)) +
    ylab("Values") +
    xlab("Samples") +
    scale_x_discrete(limits = sample_order, labels = sample_order) +
    ylim(xlower, xupper) +
    stat_boxplot(geom = "errorbar", color="black") +
    geom_boxplot(outlier.size=0.5, outlier.alpha=0.4) +
    theme_bw() + coord_flip() +
    guides(fill = guide_legend(title = "Group"));

  bp <- bp + scale_fill_manual(values = full_group_palette)

  if(interactive){
    library(plotly);
    w <- 800;  # Single panel width
    return(layout(ggplotly(bp), autosize = FALSE, width = w, height = 600));
  }else{
    
    # --- FIX: Safe Device Handling ---
    if(dpi == 72){ dpi <- 96 }
    
    # Use native cairo device if available to avoid Quartz crash
    if (format == "png" && capabilities("cairo")) {
      png(filename = fullPath, width = 600*dpi/72, height = height*dpi/72, units = "px", res = dpi, type = "cairo")
    } else {
      require("Cairo")
      Cairo(file=fullPath, width=600*dpi/72, height=height*dpi/72, unit="px", dpi=dpi, type=format, bg="white");
    }

    tryCatch({
        print(bp);
    }, finally = {
        dev.off();
    })
    # ---------------------------------

    imgSet <- readSet(imgSet, "imgSet");
    if(grepl("norm", imgNm)){
      imgSet$qc_norm_boxplot <- fullPath;
    }else{
      imgSet$qc_boxplot <- fullPath;
    }
    saveSet(imgSet, "imgSet");
    return("NA")
  }
}

qc.nonmissing.per.sample <- function(dat, imgNm, dpi = 96, format = "png",
                                     interactive = FALSE, meta = NULL) {
  dat <- coerce.numeric.matrix(dat)
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  fullPath <- paste0(fileNm, format, sep = "")

  miss.ind <- is.na(dat) | dat == "NA"
  if (is.numeric(dat)) {
    miss.ind <- is.na(dat) | dat <= 0
  }

  miss.count <- colSums(miss.ind)
  n.features <- nrow(dat)
  nonmiss.pct   <- 100 * (n.features - miss.count) / n.features

  group.fac <- NULL
  if (!is.null(meta) && ncol(meta) >= 1) {
    group.fac <- factor(meta[, 1])
    names(group.fac) <- rownames(meta)
  }

  miss.df <- data.frame(
    Sample = names(nonmiss.pct),
    NonMissingPct = as.numeric(nonmiss.pct),
    Group = if (!is.null(group.fac)) as.character(group.fac[colnames(dat)]) else "Group1",
    stringsAsFactors = FALSE
  )
  miss.df$Group <- factor(miss.df$Group)
  miss.df <- miss.df[order(miss.df$NonMissingPct, decreasing = TRUE), , drop = FALSE]
  miss.df$Sample <- factor(miss.df$Sample, levels = miss.df$Sample)

  mean_line <- mean(miss.df$NonMissingPct, na.rm = TRUE)
  bp <- ggplot(miss.df, aes(x = Sample, y = NonMissingPct, color = Group)) +
    geom_point(size = 3) +
    geom_hline(yintercept = mean_line, color = "blue", linetype = "dashed") +
    labs(x = "Sample", y = "Non-missing values (%)", color = "") +
    ylim(0, 100) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    ) +
    scale_color_manual(values = scales::hue_pal()(length(levels(miss.df$Group))))

  if (interactive) {
    library(plotly)
    return(ggplotly(bp))
  }

  # --- FIX: Safe Device Handling ---
  if (dpi == 72) { dpi <- 96 }
  num_samples <- nrow(miss.df)
  width <- ifelse(num_samples < 50, 800, 800 + (num_samples - 50) * 10)
  height <- 600

  if (format == "png" && capabilities("cairo")) {
      png(filename = fullPath, width = width * dpi / 72, height = height * dpi / 72, units = "px", res = dpi, type = "cairo")
  } else {
      Cairo(file = fullPath, width  = width * dpi / 72, height = height * dpi / 72, unit   = "px", dpi    = dpi, type   = format, bg     = "white")
  }

  tryCatch({
      print(bp)
  }, finally = {
      dev.off()
  })
  # ---------------------------------

  imgSet <- readSet(imgSet, "imgSet")
  if (grepl("norm", imgNm)) {
    imgSet$qc_norm_nonmissing <- fullPath
  } else {
    imgSet$qc_nonmissing <- fullPath
  }
  saveSet(imgSet, "imgSet")
  return("NA")
}

# MA plot (sample vs feature mean) for QC
qc.maplot <- function(dat, imgNm, dpi = 96, format = "png", interactive = FALSE, meta = NULL, sampleIndex = 1, mode = "single") {
  dat <- coerce.numeric.matrix(dat)
  dpi <- as.numeric(dpi)
  require("ggplot2")

  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  fullPath <- paste0(fileNm, format, sep = "")

  if (is.null(dat) || nrow(dat) < 2 || ncol(dat) < 1) {
    return("NA")
  }

  dat[!is.finite(dat)] <- NA_real_
  finite_vals <- as.numeric(dat[is.finite(dat)])
  frac_nonpos <- if (length(finite_vals) == 0) 1 else mean(finite_vals <= 0)
  use_logged_scale <- frac_nonpos > 0.2

  # Subsample features for responsiveness on large matrices
  max_features <- 5000
  if (nrow(dat) > max_features) {
    set.seed(28051968)
    keep <- sample(seq_len(nrow(dat)), max_features)
    dat <- dat[keep, , drop = FALSE]
  }

  sample_names <- colnames(dat)
  if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_len(ncol(dat)))
  mode <- ifelse(is.null(mode), "single", as.character(mode)[1])
  mode <- tolower(mode)
  if (!(mode %in% c("single", "overview"))) mode <- "single"

  # QC MA plot: selected sample vs pooled/global reference (feature-wise mean across samples)
  sampleIndex <- suppressWarnings(as.integer(sampleIndex))
  if (is.na(sampleIndex) || sampleIndex < 1) sampleIndex <- 1
  if (sampleIndex > ncol(dat)) sampleIndex <- 1

  row_mean <- rowMeans(dat, na.rm = TRUE)
  row_mean[!is.finite(row_mean)] <- NA_real_
  sample_indices <- sampleIndex
  if (mode == "overview") {
    sample_indices <- seq_len(min(ncol(dat), 30))
  }

  ma_list <- lapply(sample_indices, function(si) {
    x <- dat[, si]
    if (use_logged_scale) {
      a <- (x + row_mean) / 2
      m <- x - row_mean
    } else {
      x[x <= 0] <- NA_real_
      rm2 <- row_mean
      rm2[rm2 <= 0] <- NA_real_
      a <- (log2(x) + log2(rm2)) / 2
      m <- log2(x) - log2(rm2)
    }
    data.frame(
      Panel = sample_names[si],
      A = a,
      M = m,
      stringsAsFactors = FALSE
    )
  })
  ma_df <- do.call(rbind, ma_list)
  ma_df <- ma_df[is.finite(ma_df$A) & is.finite(ma_df$M), , drop = FALSE]

  if (is.null(ma_df) || nrow(ma_df) == 0) {
    return("NA")
  }

  # Calculate per-sample statistics for labeling
  ma_stats <- aggregate(M ~ Panel, data = ma_df, FUN = function(m) {
    c(median = stats::median(m, na.rm = TRUE),
      iqr = stats::IQR(m, na.rm = TRUE))
  })
  ma_stats <- data.frame(
    Panel = ma_stats$Panel,
    Median = ma_stats$M[, "median"],
    IQR = ma_stats$M[, "iqr"]
  )
  ma_stats$Label <- paste0(
    "Median(M): ", format(signif(ma_stats$Median, 4), trim = TRUE),
    "\nIQR(M): ", format(signif(ma_stats$IQR, 4), trim = TRUE)
  )

  pt_alpha <- if (mode == "overview") 0.45 else 0.65
  pt_size  <- if (mode == "overview") 0.7 else 0.9

  p <- ggplot(ma_df, aes(x = A, y = M)) +
    geom_point(alpha = pt_alpha, size = pt_size, na.rm = TRUE, color = "#2c7fb8") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_smooth(method = "loess", se = FALSE, color = "#b30000", linewidth = 0.4, span = 0.7) +
    geom_label(data = ma_stats, aes(x = Inf, y = Inf, label = Label),
               hjust = 1.18, vjust = 1.18, size = 3,
               label.size = 0.2, fill = "white", alpha = 0.9, inherit.aes = FALSE) +
    facet_wrap(~Panel, scales = "free_x") +
    labs(
      x = if (use_logged_scale) "A = average abundance (input scale)" else "A = average log2 intensity",
      y = if (use_logged_scale) "M = sample - pooled reference (input scale)" else "M = log2(sample) - log2(pooled reference)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  n_panels <- length(unique(ma_df$Panel))
  if (mode == "overview") {
    ncol_facets <- if (n_panels <= 9) 3 else if (n_panels <= 16) 4 else 5
    p <- p + facet_wrap(~Panel, scales = "free_x", ncol = ncol_facets)
    fig_w <- max(1000, min(2200, 320 * ncol_facets))
    fig_h <- max(700, min(2600, 250 * ceiling(n_panels / ncol_facets)))
  } else {
    fig_w <- 1000
    fig_h <- 650
  }

  if (interactive) {
    library(plotly)
    return(layout(ggplotly(p), autosize = FALSE, width = fig_w, height = fig_h))
  } else {
    if (dpi == 72) dpi <- 96
    if (format == "png" && capabilities("cairo")) {
      png(filename = fullPath, width = fig_w * dpi/96, height = fig_h * dpi/96, units = "px", res = dpi, type = "cairo")
    } else {
      require("Cairo")
      Cairo(file = fullPath, width = fig_w * dpi/96, height = fig_h * dpi/96, unit = "px", dpi = dpi, type = format, bg = "white")
    }
    tryCatch({
      print(p)
    }, finally = {
      dev.off()
    })

    imgSet <- readSet(imgSet, "imgSet")
    if (grepl("norm", imgNm)) {
      if (grepl("overview", imgNm, ignore.case = TRUE)) {
        imgSet$qc_norm_ma_overview <- fullPath
      } else {
        imgSet$qc_norm_ma <- fullPath
      }
    } else {
      if (grepl("overview", imgNm, ignore.case = TRUE)) {
        imgSet$qc_ma_overview <- fullPath
      } else {
        imgSet$qc_ma <- fullPath
      }
    }
    saveSet(imgSet, "imgSet")
    return("NA")
  }
}

qc.maplot.json <- function(dat, imgNm, meta = NULL) {
  dat <- coerce.numeric.matrix(dat)
  if (is.null(dat) || nrow(dat) < 2 || ncol(dat) < 1) {
    return("NA")
  }

  dat[!is.finite(dat)] <- NA_real_
  finite_vals <- as.numeric(dat[is.finite(dat)])
  frac_nonpos <- if (length(finite_vals) == 0) 1 else mean(finite_vals <= 0)
  use_logged_scale <- frac_nonpos > 0.2

  max_features <- 5000
  if (nrow(dat) > max_features) {
    set.seed(28051968)
    keep <- sample(seq_len(nrow(dat)), max_features)
    dat <- dat[keep, , drop = FALSE]
  }

  sample_names <- colnames(dat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_len(ncol(dat)))
    colnames(dat) <- sample_names
  }

  group_values <- get.primary.group.values(sample_names, meta = meta)
  group_levels <- unique(group_values)
  group_palette <- get.qc.group.palette(group_values)

  row_mean <- rowMeans(dat, na.rm = TRUE)
  row_mean[!is.finite(row_mean)] <- NA_real_

  traces <- list()
  group_a_values <- setNames(vector("list", length(group_levels)), group_levels)
  group_m_values <- setNames(vector("list", length(group_levels)), group_levels)
  for (si in seq_len(ncol(dat))) {
    x <- dat[, si]
    if (use_logged_scale) {
      a <- (x + row_mean) / 2
      m <- x - row_mean
    } else {
      x[x <= 0] <- NA_real_
      rm2 <- row_mean
      rm2[rm2 <= 0] <- NA_real_
      a <- (log2(x) + log2(rm2)) / 2
      m <- log2(x) - log2(rm2)
    }

    keep <- is.finite(a) & is.finite(m)
    if (sum(keep) < 10) {
      next
    }

    ord <- order(a[keep])
    a_sorted <- a[keep][ord]
    m_sorted <- m[keep][ord]
    smooth_df <- tryCatch(
      stats::lowess(a_sorted, m_sorted, f = 0.25, iter = 1),
      error = function(e) NULL
    )
    if (is.null(smooth_df) || length(smooth_df$x) < 2) {
      keep_n <- min(length(a_sorted), 800)
      idx <- unique(round(seq(1, length(a_sorted), length.out = keep_n)))
      smooth_df <- list(x = a_sorted[idx], y = m_sorted[idx])
    }
    if (length(smooth_df$x) < 2) {
      next
    }

    sample_name <- sample_names[si]
    group_name <- group_values[[sample_name]]
    group_a_values[[group_name]] <- c(group_a_values[[group_name]], a_sorted)
    group_m_values[[group_name]] <- c(group_m_values[[group_name]], m_sorted)
    sample_median <- stats::median(m_sorted, na.rm = TRUE)
    sample_iqr <- stats::IQR(m_sorted, na.rm = TRUE)
    hover_text <- rep(
      paste0(
        "Sample: ", sample_name,
        "<br>Group: ", group_name,
        "<br>Median(M): ", signif(sample_median, 4),
        "<br>IQR(M): ", signif(sample_iqr, 4)
      ),
      length(smooth_df$x)
    )

    traces[[length(traces) + 1]] <- list(
      x = unname(smooth_df$x),
      y = unname(smooth_df$y),
      type = "scatter",
      mode = "lines",
      name = sample_name,
      legendgroup = group_name,
      showlegend = FALSE,
      line = list(color = unname(group_palette[[group_name]]), width = 1.5),
      text = hover_text,
      hovertemplate = paste0(
        "%{text}<br>",
        "A: %{x:.3f}<br>",
        "M: %{y:.3f}<extra></extra>"
      )
    )
  }

  if (length(traces) == 0) {
    return("NA")
  }

  group_summary <- list()
  legend_traces <- list()
  for (idx in seq_along(group_levels)) {
    group_name <- group_levels[idx]
    group_a <- as.numeric(group_a_values[[group_name]])
    group_vals <- as.numeric(group_m_values[[group_name]])
    keep_group <- is.finite(group_a) & is.finite(group_vals)
    group_a <- group_a[keep_group]
    group_vals <- group_vals[keep_group]
    if (length(group_vals) < 10) {
      next
    }

    group_median <- format(signif(stats::median(group_vals, na.rm = TRUE), 4), trim = TRUE)
    group_iqr <- format(signif(stats::IQR(group_vals, na.rm = TRUE), 4), trim = TRUE)
    group_summary[[length(group_summary) + 1]] <- list(
      group = group_name,
      sample_count = sum(group_values == group_name),
      median_m = group_median,
      iqr_m = group_iqr,
      color = unname(group_palette[[group_name]])
    )

    legend_traces[[length(legend_traces) + 1]] <- list(
      x = list(NULL),
      y = list(NULL),
      type = "scatter",
      mode = "lines",
      name = group_name,
      legendgroup = group_name,
      showlegend = TRUE,
      line = list(color = unname(group_palette[[group_name]]), width = 3),
      hoverinfo = "skip"
    )
  }

  plot_data <- list(
    data = c(traces, legend_traces),
    groupSummary = group_summary,
    layout = list(
      title = list(text = "MA Trend by Sample"),
      plot_bgcolor = "#FFFFFF",
      paper_bgcolor = "#FFFFFF",
      hovermode = "closest",
      xaxis = list(
        title = if (use_logged_scale) "A = average abundance (input scale)" else "A = average log2 intensity",
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "#000000",
        ticks = "outside"
      ),
      yaxis = list(
        title = if (use_logged_scale) "M = sample - pooled reference (input scale)" else "M = log2(sample) - log2(pooled reference)",
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "#000000",
        ticks = "outside"
      ),
      shapes = list(
        list(
          type = "line",
          xref = "paper",
          x0 = 0,
          x1 = 1,
          yref = "y",
          y0 = 0,
          y1 = 0,
          line = list(color = "#444444", dash = "dash")
        )
      ),
      legend = list(
        title = list(text = "Primary Group"),
        orientation = "v",
        x = 1.02,
        y = 1,
        xanchor = "left",
        yanchor = "top"
      ),
      margin = list(l = 70, r = 180, t = 50, b = 60)
    )
  )

  json_nm <- paste0(imgNm, ".json")
  jsonlite::write_json(plot_data, json_nm, auto_unbox = TRUE, digits = 16, null = "null")
  invisible("NA")
}

qc.sample.dendro <- function(dat, imgNm, dpi = 96, format = "png",
                             interactive = FALSE, method = "correlation",
                             linkage = "complete", meta = NULL) {
  if (is.data.frame(dat)) dat <- as.matrix(dat)
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  fullPath <- paste0(fileNm, format, sep = "")

  if (ncol(dat) < 2) {
    return("NA")
  }

  suppressWarnings({
    cor.mat <- stats::cor(dat, use = "pairwise.complete.obs", method = "pearson")
  })
  if (is.null(cor.mat) || all(!is.finite(cor.mat))) {
    return("NA")
  }
  cor.mat[!is.finite(cor.mat)] <- 0
  diag(cor.mat) <- 1

  hc <- hclust(as.dist(1 - cor.mat), method = linkage)
  ord <- hc$order
  cor.ord <- cor.mat[ord, ord, drop = FALSE]
  sample.nms <- colnames(cor.ord)

  corr.range <- cor.ord[row(cor.ord) != col(cor.ord)]
  corr.range <- corr.range[is.finite(corr.range)]
  if (length(corr.range) == 0) {
    corr.min <- -1
    corr.max <- 1
  } else {
    corr.min <- min(corr.range)
    corr.max <- max(corr.range)
    if (!is.finite(corr.min) || !is.finite(corr.max) || corr.min == corr.max) {
      pad <- max(0.01, abs(corr.min) * 0.05)
      corr.min <- corr.min - pad
      corr.max <- corr.max + pad
    }
  }

  ann.df <- NULL
  ann.cols <- list()
  ann.pals <- list()
  if (!is.null(meta) && ncol(meta) >= 1) {
    ann.df <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (!is.null(rownames(ann.df))) {
      ann.df <- ann.df[sample.nms, , drop = FALSE]
    }
    ann.keep <- setdiff(colnames(ann.df), "SampleName")
    ann.keep <- ann.keep[seq_len(min(2, length(ann.keep)))]
    if (length(ann.keep) > 0) {
      ann.df <- ann.df[, ann.keep, drop = FALSE]
      for (idx in seq_along(ann.keep)) {
        nm <- ann.keep[idx]
        vals <- as.character(ann.df[[nm]])
        vals[is.na(vals) | vals == ""] <- "NA"
        pal <- get.qc.annotation.palette(vals, annotation_index = idx, use_master_group_palette = idx == 1L)
        ann.cols[[nm]] <- pal[vals]
        ann.pals[[nm]] <- pal
      }
    } else {
      ann.df <- NULL
    }
  }

  ann.count <- if (is.null(ann.df)) 0 else ncol(ann.df)
  width_in <- max(7, min(16, 4 + ncol(cor.ord) * 0.28))
  height_in <- max(6, min(16, 4 + nrow(cor.ord) * 0.22 + ann.count * 0.45 + 1.0))
  legend_height <- if (ann.count > 0) 1.4 else 0.9
  
  # --- FIX: Safe Device Handling ---
  if (dpi == 72) dpi <- 96

  if (format == "png" && capabilities("cairo")) {
     png(filename = fullPath, width = width_in * dpi, height = height_in * dpi, units = "px", res = dpi, type = "cairo")
  } else {
     Cairo(file = fullPath, width = width_in, height = height_in, unit = "in", dpi = dpi, type = format, bg = "white")
  }

  tryCatch({
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)

      layout(matrix(seq_len(ann.count + 2), ncol = 1), heights = c(rep(0.45, ann.count), 5.5, legend_height))

      if (ann.count > 0) {
        n <- length(sample.nms)
        for (nm in colnames(ann.df)) {
          par(mar = c(0, 8, 0.2, 2))
          plot.new()
          plot.window(xlim = c(0.5, n + 0.5), ylim = c(0, 1))
          rect(seq_len(n) - 0.5, 0, seq_len(n) + 0.5, 1, col = ann.cols[[nm]], border = NA)
          text(0.3, 0.5, labels = nm, xpd = NA, adj = 1, font = 2, cex = 0.85)
          box()
        }
      }

      n <- ncol(cor.ord)
      par(mar = c(9, 8, 2, 3))
      heat_cols <- grDevices::colorRampPalette(c("#2c7bb6", "#f7f7f7", "#d7191c"))(100)
      image(
        x = seq_len(n),
        y = seq_len(n),
        z = t(cor.ord[n:1, , drop = FALSE]),
        col = heat_cols,
        zlim = c(corr.min, corr.max),
        axes = FALSE,
        xlab = "",
        ylab = "",
        useRaster = TRUE
      )
      axis(1, at = seq_len(n), labels = sample.nms, las = 2, cex.axis = max(0.45, min(0.85, 12 / n)))
      axis(2, at = seq_len(n), labels = rev(sample.nms), las = 2, cex.axis = max(0.45, min(0.85, 12 / n)))
      box()
      title(main = "Sample similarity correlation heatmap", line = 0.2)

      par(mar = c(0, 0, 0, 0))
      plot.new()
      usr <- par("usr")
      x0 <- usr[1] + 0.03 * diff(usr[1:2])
      y0 <- usr[3] + 0.60 * diff(usr[3:4])
      grad_w <- 0.22 * diff(usr[1:2])
      grad_h <- 0.18 * diff(usr[3:4])
      grad_cols <- grDevices::colorRampPalette(c("#2c7bb6", "#f7f7f7", "#d7191c"))(50)
      for (i in seq_along(grad_cols)) {
        rect(x0 + (i - 1) * grad_w / length(grad_cols), y0,
             x0 + i * grad_w / length(grad_cols), y0 + grad_h,
             col = grad_cols[i], border = NA)
      }
      rect(x0, y0, x0 + grad_w, y0 + grad_h, border = "grey40")
      text(x0, y0 - 0.08, labels = format(signif(corr.min, 3), trim = TRUE), adj = c(0, 1), cex = 0.8)
      text(x0 + grad_w / 2, y0 - 0.08, labels = format(signif((corr.min + corr.max) / 2, 3), trim = TRUE), adj = c(0.5, 1), cex = 0.8)
      text(x0 + grad_w, y0 - 0.08, labels = format(signif(corr.max, 3), trim = TRUE), adj = c(1, 1), cex = 0.8)
      text(x0, y0 + grad_h + 0.08, labels = "Pearson correlation", adj = c(0, 0), cex = 0.85, font = 2)

      if (ann.count > 0) {
        legend_x <- x0 + grad_w + 0.08 * diff(usr[1:2])
        legend_y <- usr[4] - 0.05 * diff(usr[3:4])
        for (nm in names(ann.pals)) {
          legend(legend_x, legend_y, legend = names(ann.pals[[nm]]), fill = ann.pals[[nm]],
                 border = NA, title = nm, bty = "n", cex = 0.75, xjust = 0, yjust = 1)
          legend_y <- legend_y - (0.11 + 0.045 * length(ann.pals[[nm]])) * diff(usr[3:4])
        }
      }
  }, finally = {
      dev.off()
  })
  # ---------------------------------

  imgSet <- readSet(imgSet, "imgSet");
  if (grepl("norm", imgNm)) {
    imgSet$qc_norm_dendro <- fullPath;
  } else {
    imgSet$qc_dendro <- fullPath;
  }
  saveSet(imgSet, "imgSet");
  return("NA")
}

#' Sample correlation heatmap (metadata-style aesthetic)
qc.sample.corr <- function(dat, imgNm, dpi=96, format="png", interactive=FALSE, method="pearson",
                           colorGradient = "default") {
  if (is.data.frame(dat)) dat <- as.matrix(dat)
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="")
  imgNm <- paste0(fileNm, format, sep="")

  if (ncol(dat) < 2 || nrow(dat) == 0) {
    #msg("[qc.sample.corr] Not enough samples to compute correlation")
    return("NA")
  }
  cor.mat <- tryCatch(stats::cor(dat, use = "pairwise.complete.obs", method = method), error = function(e) NULL)
  if (is.null(cor.mat)) {
    #msg("[qc.sample.corr] Correlation failed")
    return("NA")
  }
  cor.mat <- round(cor.mat, 3)

  get_upper_tri <- function(mat) { mat[lower.tri(mat)] <- NA; mat }
  melted <- reshape2::melt(get_upper_tri(cor.mat), na.rm = TRUE)
  melted$value <- signif(melted$value, 3)

  # Colour scale similar to PlotMetaCorrHeatmap defaults
  colour_scale <- switch(tolower(colorGradient),
                         gbr     = scale_fill_gradient2(low = "green", mid = "black", high = "red",
                                                        midpoint = 0, limits = c(-1, 1), name = "Correlation"),
                         heat    = scale_fill_gradientn(colours = heat.colors(10)),
                         topo    = scale_fill_gradientn(colours = topo.colors(10)),
                         gray    = scale_fill_gradientn(colours = c("grey90", "grey10")),
                         byr     = scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "RdYlBu"))),
                         viridis = scale_fill_viridis_c(option = "viridis"),
                         plasma  = scale_fill_viridis_c(option = "plasma"),
                         npj     = scale_fill_gradient2(low = "#00A087FF", mid = "white", high = "#E64B35FF",
                                                        midpoint = 0, limits = c(-1, 1), name = "Correlation"),
                         aaas    = scale_fill_gradient2(low = "#4DBBD5FF", mid = "white", high = "#E64B35FF",
                                                        midpoint = 0, limits = c(-1, 1), name = "Correlation"),
                         d3      = scale_fill_gradient2(low = "#2CA02CFF", mid = "white", high = "#FF7F0EFF",
                                                        midpoint = 0, limits = c(-1, 1), name = "Correlation"),
                         scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                              midpoint = 0, limits = c(-1, 1), name = "Correlation"))

  samp.num <- ncol(dat)
  if (samp.num > 25) { width <- 24; textSize <- 3.5 }
  else if (samp.num > 10) { width <- 16; textSize <- 4 }
  else { width <- 10; textSize <- 4 }
  height <- max(10, samp.num * 0.6)

  p <- ggplot(melted, aes(Var2, Var1, fill = value)) +
    geom_tile(colour = "white") +
    scale_y_discrete(position = "right") +
    colour_scale +
    theme_minimal() +
    theme(axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.text.y.right = element_text(),
          legend.direction = "vertical",
          legend.position  = "left") +
    coord_fixed() +
    geom_text(aes(label = value), size = textSize, colour = "black") +
    ggtitle(paste("Sample correlation (", method, ")", sep=""))

  imgSet <- readSet(imgSet, "imgSet");
  if (grepl("norm", imgNm)) {
    imgSet$qc_norm_corr <- imgNm;
  } else {
    imgSet$qc_corr <- imgNm;
  }
  saveSet(imgSet, "imgSet");
  if(dpi == 72){
    dpi <- dpi *1.34
  }
  # Scale dimensions with sample size and convert to inches (Cairo unit="in")
  Cairo(file=imgNm, unit="in", dpi=dpi, width=width, height=height, type=format, bg="white");
  print(p);
  dev.off();
  return("NA")
}

qc.sample.corr.json <- function(dataSet, imgNm, method = "pearson") {
  suppressMessages({
    require(jsonlite)
    require(RColorBrewer)
  })

  dat <- dataSet$data.norm
  if (is.null(dat) || ncol(dat) < 2 || nrow(dat) == 0) {
    return("NA")
  }

  dat <- as.matrix(dat)
  cor.mat <- tryCatch(
    stats::cor(dat, use = "pairwise.complete.obs", method = method),
    error = function(e) NULL
  )
  if (is.null(cor.mat)) {
    return("NA")
  }

  cor.mat[!is.finite(cor.mat)] <- 0
  diag(cor.mat) <- 1

  corr.range <- cor.mat[row(cor.mat) != col(cor.mat)]
  corr.range <- corr.range[is.finite(corr.range)]
  if (length(corr.range) == 0) {
    corr.min <- -1
    corr.max <- 1
  } else {
    corr.min <- min(corr.range)
    corr.max <- max(corr.range)
    if (!is.finite(corr.min) || !is.finite(corr.max) || corr.min == corr.max) {
      pad <- max(0.01, abs(corr.min) * 0.05)
      corr.min <- corr.min - pad
      corr.max <- corr.max + pad
    }
  }

  hc <- tryCatch(
    stats::hclust(stats::as.dist(pmax(0, 1 - cor.mat)), method = "average"),
    error = function(e) NULL
  )
  ord <- if (is.null(hc)) seq_len(ncol(cor.mat)) else hc$order

  cor.mat <- cor.mat[ord, ord, drop = FALSE]
  sample_names <- colnames(cor.mat)

  meta.info <- dataSet$meta.info
  annotation.names <- character(0)
  annotation.df <- NULL
  if (!is.null(meta.info) && ncol(meta.info) > 0) {
    annotation.names <- colnames(meta.info)[seq_len(min(2, ncol(meta.info)))]
    annotation.df <- as.data.frame(meta.info[sample_names, annotation.names, drop = FALSE], stringsAsFactors = FALSE)
  }

  sample_count <- length(sample_names)
  base_height <- max(560, min(1200, sample_count * 24))
  anno_count <- length(annotation.names)
  anno_total <- if (anno_count > 0) min(0.18, anno_count * 0.07 + 0.02) else 0
  main_domain <- c(0, 1 - anno_total)

  heatmap_trace <- list(
    type = "heatmap",
    x = sample_names,
    y = sample_names,
    z = unname(cor.mat),
    zmin = corr.min,
    zmax = corr.max,
    xgap = 1,
    ygap = 1,
    colorscale = list(
      list(0, "#2166ac"),
      list(0.5, "#f7f7f7"),
      list(1, "#b2182b")
    ),
    colorbar = list(
      title = list(text = "Correlation"),
      x = 1.02,
      y = 0.5,
      len = 0.75
    ),
    hovertemplate = paste0(
      "Sample X: %{x}<br>",
      "Sample Y: %{y}<br>",
      "Correlation: %{z:.3f}<extra></extra>"
    )
  )

  traces <- list(heatmap_trace)

  layout <- list(
    title = list(text = paste0("Sample Similarity (", method, " correlation)")),
    xaxis = list(
      domain = c(0, 1),
      side = "bottom",
      tickangle = -45,
      automargin = TRUE,
      showgrid = FALSE,
      zeroline = FALSE
    ),
    yaxis = list(
      domain = main_domain,
      autorange = "reversed",
      automargin = TRUE,
      showgrid = FALSE,
      zeroline = FALSE
    ),
    margin = list(l = 130, r = 120, t = 60, b = 120),
    height = base_height,
    paper_bgcolor = "#ffffff",
    plot_bgcolor = "#ffffff"
  )

  if (!is.null(annotation.df) && ncol(annotation.df) > 0) {
    anno_band <- anno_total / anno_count
    legend_traces <- list()
    for (i in seq_len(ncol(annotation.df))) {
      col_name <- colnames(annotation.df)[i]
      vals <- annotation.df[[i]]
      axis_name <- paste0("yaxis", i + 1)
      trace_axis <- paste0("y", i + 1)
      domain_end <- 1 - (i - 1) * anno_band
      domain_start <- domain_end - anno_band + 0.01

      if (is.numeric(vals)) {
        zvals <- matrix(as.numeric(vals), nrow = 1)
        num.range <- range(as.numeric(vals), na.rm = TRUE)
        if (!all(is.finite(num.range))) {
          num.range <- c(0, 1)
        } else if (num.range[1] == num.range[2]) {
          num.range <- num.range + c(-0.5, 0.5)
        }
        colorscale <- "Viridis"
        custom_vals <- matrix(format(signif(vals, 4), trim = TRUE), nrow = 1)
        zmin.val <- num.range[1]
        zmax.val <- num.range[2]
      } else {
        vals <- as.character(vals)
        vals[is.na(vals) | vals == ""] <- "NA"
        levs <- unique(vals)
        pal <- get.qc.annotation.palette(vals, annotation_index = i, use_master_group_palette = i == 1L)
        idx <- match(vals, levs) - 1
        if (length(levs) == 1) {
          colorscale <- list(list(0, pal[1]), list(1, pal[1]))
        } else {
          denom <- length(levs)
          colorscale <- unlist(lapply(seq_along(levs), function(j) {
            start <- (j - 1) / denom
            end <- j / denom
            list(list(start, pal[j]), list(end, pal[j]))
          }), recursive = FALSE)
        }
        zvals <- matrix(idx, nrow = 1)
        custom_vals <- matrix(vals, nrow = 1)
        zmin.val <- 0
        zmax.val <- max(length(levs) - 1, 1)

        legend_traces[[length(legend_traces) + 1]] <- lapply(seq_along(levs), function(j) {
          list(
            type = "scatter",
            mode = "markers",
            x = list(NULL),
            y = list(NULL),
            name = paste0(col_name, ": ", levs[j]),
            marker = list(color = pal[j], size = 10, symbol = "square"),
            hoverinfo = "skip",
            showlegend = TRUE
          )
        })
      }

      traces[[length(traces) + 1]] <- list(
        type = "heatmap",
        x = sample_names,
        y = list(col_name),
        z = zvals,
        xaxis = "x",
        yaxis = trace_axis,
        colorscale = colorscale,
        zmin = zmin.val,
        zmax = zmax.val,
        showscale = FALSE,
        xgap = 1,
        ygap = 1,
        customdata = custom_vals,
        hovertemplate = paste0(col_name, ": %{customdata}<br>Sample: %{x}<extra></extra>")
      )

      layout[[axis_name]] <- list(
        domain = c(domain_start, domain_end),
        automargin = TRUE,
        showgrid = FALSE,
        zeroline = FALSE,
        ticks = "",
        showticklabels = TRUE
      )

      if (i < ncol(annotation.df)) {
        layout[[paste0("xaxis", i + 1)]] <- list(
          domain = c(0, 1),
          anchor = trace_axis,
          matches = "x",
          showticklabels = FALSE,
          showgrid = FALSE,
          zeroline = FALSE
        )
      }

      if (length(legend_traces) > 0) {
        for (legend_group in legend_traces) {
          for (legend_trace in legend_group) {
            traces[[length(traces) + 1]] <- legend_trace
          }
        }
        legend_traces <- list()
      }
    }

    layout$legend <- list(
      orientation = "v",
      x = 1.14,
      y = 1,
      xanchor = "left",
      yanchor = "top",
      font = list(size = 11)
    )
  }

  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE,
    digits = 16,
    null = "null"
  )

  invisible("NA")
}


PlotDataDensity <- function(fileName, imgNm, dpi,format){
  dataSet <- readDataset(fileName);
  res <- qc.density(dataSet, imgNm, dpi, format, FALSE);
  return(res);
}

qc.density<- function(dataSet, imgNm="abc", dpi=96, format, interactive){
  require("ggplot2")
  dat <- dataSet$data.norm
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  imgNm <- paste0(fileNm, format, sep="");
  dpi <- as.numeric(dpi)

  ######check.names=F is important for names not following conventions (i.e. with -, starts with number)
  df <- data.frame(dataSet$data.norm, stringsAsFactors = FALSE, check.names = FALSE)
  df <- stack(df)
  sampleNms <-colnames(dataSet$data.norm)
  if(length(dataSet$meta.info) == 2){

    # OPTIMIZED: Single merge instead of sequential merges to eliminate intermediate copies
    Factor1 <- dataSet$meta.info[,1]
    Factor2 <- dataSet$meta.info[,2]
    factorNm1 <- colnames(dataSet$meta.info)[1]
    factorNm2 <- colnames(dataSet$meta.info)[2]

    # Build combined metadata once
    conv <- data.frame(
      ind = sampleNms,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    conv[[factorNm1]] <- Factor1
    conv[[factorNm2]] <- Factor2

    # Single merge operation
    df1 <- merge(df, conv, by="ind")
    df2 <- reshape::melt(df1, measure.vars=c(factorNm1,factorNm2))
    colnames(df2)[4] <- "Conditions"
    
    g = ggplot(df2, aes(x=values)) + 
      geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) + 
      facet_grid(. ~ variable) +
      theme_bw()
    
    width <- 12
    height <- 6
  }else{
    Conditions <- dataSet$meta.info[,1];
    conv <- data.frame(ind=sampleNms, Conditions=Conditions, stringsAsFactors = FALSE, check.names = FALSE)
    df1 <- merge(df, conv, by="ind")
    
    g = ggplot(df1, aes(x=values)) + 
      geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) +
      theme_bw()
    
    width <- 8
    height <- 6
  }
  
  if(interactive){
    library(plotly);
    m <- list(
      l = 50,
      r = 50,
      b = 20,
      t = 20,
      pad = 0.5
    )
    if(length(dataSet$meta.info) == 2){
      w=1000;
    }else{
      w=800;
    }
    ggp_build <- layout(ggplotly(g), autosize = FALSE, width = w, height = 600, margin = m)
    return(ggp_build);
  }else{
    imgSet <- readSet(imgSet, "imgSet");
    imgSet$qc.density_norm <- imgNm;
    saveSet(imgSet, "imgSet");
  if(dpi == 72){
  dpi <- dpi *1.34
  }
    Cairo(file=imgNm, width=width, height=height, type=format, bg="white", dpi=dpi, unit="in");
    print(g)
    dev.off();
    return("NA")
  }
}


PlotLibSizeView<-function(fileName, imgNm,dpi=96, format="png"){
  library("ggrepel");
  require("ggplot2");

  dataSet <- readDataset(fileName);
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  imgNm <- paste0(fileNm, format, sep="");
  dpi <- as.numeric(dpi);

  data.anot <- .get.annotated.data();
  data_bef<-data.matrix(data.anot);
  
  # Use NA-tolerant column sums; raw proteomics matrices often contain missing values
  smpl.sums <- colSums(data_bef, na.rm = TRUE);
  
  names(smpl.sums) <- sampleNms <- colnames(data_bef);
  df <- data.frame(count=smpl.sums,ind=colnames(data_bef))
  
  if(length(dataSet$meta.info) == 2){
    Factor1 <- as.vector(dataSet$meta.info[,1])
    factor1Nm <- colnames(dataSet$meta.info)[1]
    conv <- data.frame(ind=sampleNms, Factor1=Factor1)
    colnames(conv) <- c("ind", factor1Nm)
    df1 <- merge(df, conv, by="ind")
    Factor2 <- as.vector(dataSet$meta.info[,2])
    factor2Nm <- colnames(dataSet$meta.info)[2]
    conv <- data.frame(ind=sampleNms, Factor2=Factor2)
    colnames(conv) <- c("ind", factor2Nm)
    df1 <- merge(df1, conv, by="ind")
    df2 <- reshape::melt(df1, measure.vars=c(factor1Nm,factor2Nm))
    colnames(df2)[4] <- "Conditions"
    if(length(df2$ind)>20){

      g <- ggplot(df2, aes(x = Conditions, y = count, fill=Conditions, label= ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',position = position_dodge(), dotsize=0.7) + 
        geom_text() + 
        ylab("Sum") + 
        facet_grid(. ~ variable) +
        theme_bw()

      plotData <- ggplot_build(g)
      g$layers[[2]] = NULL;
    }else{

      g <- ggplot(df2, aes(x = Conditions, y = count, fill=Conditions, label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) + 
        geom_text_repel(force=5) + 
        ylab("Sum") + 
        facet_grid(. ~ variable) +
        theme_bw()

      plotData <- ggplot_build(g)
    }
    width <- 12
    height <- 6
    
  }else{
    Conditions= as.character(dataSet$meta.info[,1]);
    conv <- data.frame(ind=sampleNms, Conditions=Conditions)
    df1 <- merge(df, conv, by="ind")
    if(length(df1$ind)>20){
      g <- ggplot(df1, aes(x = Conditions, y = count, fill=Conditions, label= ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) + 
        xlab("Sum") + 
        geom_text() +
        theme_bw()

      plotData <- ggplot_build(g)
      g$layers[[2]] = NULL;
    }else{

      g <- ggplot(df1, aes(x = Conditions, y = count, label=ind, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), dotsize=0.7) + 
        geom_text_repel(force=5) + 
        xlab("Sum") +
        theme_bw()

      plotData <- ggplot_build(g)
      
    }
    width <- 8
    height <- 6
    
  }
  if(dpi == 72){
  dpi <- dpi *1.34
  }
  Cairo(file=imgNm, width=width, height=height, type=format, bg="white", unit="in", dpi=dpi);
  print(g);
  dev.off();
  str <- "NA"

  imgSet <- readSet(imgSet, "imgSet");
  imgSet$libsize <- imgNm;
    saveSet(imgSet, "imgSet");

  return(str);
}

PlotDataMeanStd <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName);
  if(grepl("_norm", imgName)){
    res <- qc.meanstd(dataSet$data.norm, imgName, dpi, format);
  }else{
    data.anot <- .get.annotated.data();
    res <- qc.meanstd(data.anot, imgName, dpi, format);
  }
  return(res);
}

qc.meanstd <- function(dat, imgNm, dpi=96, format="png"){
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="");
  fullPath <- paste0(fileNm, format, sep="");
  
  # --- FIX: Safe Device Handling ---
  if(dpi == 72){ dpi <- 96 }
  
  if (format == "png" && capabilities("cairo")) {
      png(filename = fullPath, width=8, height=6, type="cairo", units="in", res=dpi)
  } else {
      Cairo(file=fullPath, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  }
  
  tryCatch({
      # Call meanSdPlot with plot=FALSE to prevent it from printing to the wrong device
      res <- meanSdPlot(dat, ranks=FALSE, plot=FALSE) 
      print(res$gg)
  }, finally = {
      dev.off();
  })
  # ---------------------------------

  imgSet <- readSet(imgSet, "imgSet");
  imgSet$qc.meanstd <- fullPath;
  saveSet(imgSet, "imgSet");

  return("NA");
}

meanSdPlot <- function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
                       ylab = "sd", pch, plot = TRUE, bins = 50, ...) {
  
  stopifnot(is.logical(ranks), length(ranks) == 1, !is.na(ranks))
  
  n <- nrow(x)
  if (n == 0L) {
    warning("In 'meanSdPlot': input matrix 'x' has 0 rows. There is nothing to be done.")
    return()
  }

  px   = rowMeans(x, na.rm = TRUE)
  py   = sqrt(rowV(x, mean = px, na.rm = TRUE))
  rpx  = rank(px, na.last = FALSE, ties.method = "random")
  
  dm        = 0.025
  midpoints = seq(dm, 1-dm, by = dm)
  within    = function(x, x1, x2) { (x >= x1) & (x <= x2) }
  mediwind  = function(mp) median(py[within(rpx/n, mp - 2*dm, mp + 2*dm)], na.rm = TRUE)
  rq.sds    = sapply(midpoints, mediwind)
  
  res = if(ranks) {
    list(rank = midpoints*n, sd = rq.sds, px = rpx, py = py)
  } else {
    list(quantile = quantile(px, probs = midpoints, na.rm = TRUE), sd = rq.sds, px = px, py = py)
  }
  
  fmt = function() function(x) format(round(x, 0), nsmall = 0L, scientific = FALSE)
  
  res$gg = ggplot(data.frame(px = res$px, py = res$py),
                  aes_string(x = "px", y = "py")) + 
    xlab(xlab) + ylab(ylab) +
    geom_hex(bins = bins, ...) +
    scale_fill_gradient(name = "count", trans = "log", labels = fmt()) + 
    geom_line(aes_string(x = "x", y = "y"),
              data = data.frame(x = res[[1]], y = res$sd), color = "red") +
    theme_bw();
  
  # Only print if explicitly requested (qc.meanstd calls this with FALSE)
  if (plot) print(res$gg)
  
  return(invisible(res))
}

rowV = function(x, mean, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<1]  = NA
  if(missing(mean))
    mean=rowMeans(x, ...)
  return(rowSums(sqr(x-mean), ...)/(n-1))
}



PlotDataPCA <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName);

  if(grepl("_norm", imgName)){
    qc.pcaplot(dataSet, dataSet$data.norm, imgName, dpi, format, F);
    if (paramSet$oneDataAnalType == "dose") {
      qc.pcaplot.outliers.json(dataSet, dataSet$data.norm, imgName);
    } else {
      qc.pcaplot.json(dataSet, dataSet$data.norm, imgName);
    }
  } else {
    data.anot <- .get.annotated.data();
    qc.pcaplot(dataSet, data.anot, imgName, dpi, format, F);
    if (paramSet$oneDataAnalType == "dose") {
      qc.pcaplot.outliers.json(dataSet, data.anot, imgName);
    } else {
      qc.pcaplot.json(dataSet, data.anot, imgName);
    }
  }
  return("NA");
}


qc.pcaplot <- function(dataSet, x, imgNm, dpi=96, format="png", interactive=FALSE) {
  dpi <- as.numeric(dpi)
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep="")
  fullPath <- paste0(fileNm, format, sep="")
  
  require('lattice')
  require('ggplot2')
  require('reshape')
  require('see')
  require('ggrepel')

  if (!is.matrix(x)) x <- try(as.matrix(x), silent = TRUE)
  if (inherits(x, "try-error") || is.null(x)) return("NA")

  x[!is.finite(x)] <- NA
  clean <- na.omit(x)
  if (nrow(clean) < 2 || ncol(clean) < 2) return("NA")

  pca <- prcomp(t(clean))
  imp.pca <- summary(pca)$importance
  if (ncol(pca$x) < 2) return("NA")
  
  xlabel <- paste0("PC1"," (", 100 * round(imp.pca[2,][1], 3), "%)")
  ylabel <- paste0("PC2"," (", 100 * round(imp.pca[2,][2], 3), "%)")
  pca.res <- as.data.frame(pca$x)[, 1:2]
  
  # Increase xlim and ylim for text label
  xlim <- GetExtendRange(pca.res$PC1)
  ylim <- GetExtendRange(pca.res$PC2)
  
  common_smpls <- intersect(rownames(pca.res), rownames(dataSet$meta.info))
  pca.res <- pca.res[common_smpls, , drop = FALSE]
  dataSet$meta.info <- dataSet$meta.info[common_smpls, , drop = FALSE]

  # Always use only the primary (first) metadata column for coloring
  Factor <- dataSet$meta.info[, 1]
  legend_title <- colnames(dataSet$meta.info)[1]
  pca.res$Conditions <- Factor
  pcafig <- ggplot(pca.res, aes(x = PC1, y = PC2, color = Conditions)) +
      geom_point(size = 3, alpha = 0.5) + xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) +
      theme_bw() + labs(color = legend_title)

  width <- 10; height <- 6
  if(grepl("norm", imgNm)){
      n_conditions <- length(unique(pca.res$Conditions))

      if (paramSet$oneDataAnalType == "dose") {
        pal <- colorRampPalette(c("#2196F3", "#DE690D"))
        col.pal <- pal(n_conditions)
        pcafig <- pcafig + scale_fill_manual(values = col.pal) + scale_color_manual(values = col.pal)
      } else {
        color_pal <- get.qc.group.palette(pca.res$Conditions)
        pcafig <- pcafig + scale_fill_manual(values = color_pal) + scale_color_manual(values = color_pal)
      }
  }

  analSet <- readSet(analSet, "analSet")
  num_pcs <- min(3, ncol(pca$x))
  analSet$pca <- list(x = pca$x[, 1:num_pcs, drop = FALSE], xlabel = xlabel, ylabel = ylabel)
  saveSet(analSet, "analSet");

  if (interactive) {
    library(plotly)
    w <- 800  # Single panel width (no longer creating dual panels)
    return(layout(ggplotly(pcafig), autosize=FALSE, width=w, height=600))
  } else {
    
    # --- FIX: Safe Device Handling ---
    if(dpi == 72){ dpi <- 96 }

    if (format == "png" && capabilities("cairo")) {
      png(filename = fullPath, width=width, height=height, type="cairo", units="in", res=dpi)
    } else {
      Cairo(file = fullPath, width=width, height=height, type=format, bg="white", unit="in", dpi=dpi)
    }
    
    tryCatch({
        print(pcafig)
    }, finally = {
        dev.off()
    })
    # ---------------------------------
    return("NA")
  }
}

GetPcaOutliers <- function(){
    paramSet <- readSet(paramSet, "paramSet")
    return(paramSet$pca.outliers);
}

PlotDataNcov5 <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName)
  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }

  ncov5_df <- dataSet$summary_df[, c("Sample", "HighCoverageGeneCount")]
  
  ## Compute outlier limits
  Q1  <- quantile(ncov5_df$HighCoverageGeneCount, 0.25)
  Q3  <- quantile(ncov5_df$HighCoverageGeneCount, 0.75)
  IQRv <- IQR(ncov5_df$HighCoverageGeneCount)
  lower <- Q1 - 3 * IQRv
  upper <- Q3 + 3 * IQRv

  ncov5_df$Status <- ifelse(ncov5_df$HighCoverageGeneCount < lower |
                            ncov5_df$HighCoverageGeneCount > upper,
                            "Outlier", "Normal")

  qc.ncov5.plot(ncov5_df, imgName, lower, upper, dpi, format);
  qc.ncov5plot.json(ncov5_df, imgName, lower, upper);
  return("NA")
}

qc.ncov5.plot <- function(ncov5_df,
                          imgNm = "NCov5_plot",
                          lower,
                          upper,
                          dpi = 72,
                          format = "png",
                          interactive = FALSE) {
  require(ggplot2)
  require(ggrepel)
  require(Cairo)
  
  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be a positive number.")

  g <- ggplot(ncov5_df, aes(x = "", y = HighCoverageGeneCount)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80") +
    geom_jitter(aes(color = Status), width = 0.25, height = 0) +
    geom_hline(yintercept = c(lower, upper),
               linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(ncov5_df, Status == "Outlier"),
                    aes(label = Sample), nudge_x = 0.35, size = 3) +
    scale_color_manual(values = c(Normal = "grey40", Outlier = "red"),
                       name = "Sample status") +
    theme_minimal(base_size = 11) +
    labs(x = NULL,
         y = "features with > 5 uniquely mapped reads") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  width  <- 8
  height <- 6
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  imgNm  <- paste0(fileNm, format)

  if (interactive) {
    require(plotly)
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = imgNm, width = width, height = height,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return("NA")
  }
}


PlotDataNsig <- function(fileName, imgName, dpi, format){
  dataSet <- readDataset(fileName)
  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }

  nsig_df <- dataSet$summary_df[, c("Sample", "NSig80")]


  ## identify outliers (± 3×IQR)
  Q1  <- quantile(nsig_df$NSig80, 0.25)
  Q3  <- quantile(nsig_df$NSig80, 0.75)
  IQRv <- IQR(nsig_df$NSig80)
  lower <- Q1 - 3 * IQRv
  upper <- Q3 + 3 * IQRv

  nsig_df$outlier <- ifelse(nsig_df$NSig80 < lower | nsig_df$NSig80 > upper,
                            "Outlier", "Normal")

  qc.nsig.plot(nsig_df, imgName, lower, upper, dpi, format)
  qc.nsigplot.json(nsig_df, imgName, lower, upper); 
  return("NA")
}

qc.nsig.plot <- function(nsig_df,
                         imgNm = "NSig80_plot",
                         lower,
                         upper,
                         dpi = 72,
                         format = "png",
                         interactive = FALSE) {
  require("ggplot2")
  require("Cairo")
  require("ggrepel")

  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be a positive number.")

  g <- ggplot(nsig_df, aes(x = "", y = NSig80)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80") +
    geom_jitter(aes(color = outlier), width = 0.25, height = 0) +
    scale_color_manual(values = c("Normal" = "grey40", "Outlier" = "red")) +
    geom_text_repel(data = subset(nsig_df, outlier == "Outlier"),
                    aes(label = Sample), nudge_x = 0.35, size = 3) +
    geom_hline(yintercept = c(lower, upper), linetype = "dashed",
               color = "blue") +
    theme_minimal(base_size = 11) +
    labs(x = NULL,
         y = "NSig80 (features reaching 80 % of signal)",
         color = "Sample Status") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  width  <- 8
  height <- 6
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  imgNm  <- paste0(fileNm, format)

  if (interactive) {
    require("plotly")
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = imgNm, width = width, height = height,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return("NA")
  }
}

PlotDataDendrogram <- function(fileName, imgName, threshold, dpi, format){
  dataSet <- readDataset(fileName)

  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }
  dendro_df <- dataSet$summary_df[, c("Sample", "Dendrogram_Distance")]
  finite_rows <- is.finite(dendro_df$Dendrogram_Distance)
  if (!any(finite_rows)) {
    msgSet <- readSet(msgSet, "msgSet")
    msgSet$current.msg <- "Dendrogram plot skipped: no finite pair-wise distances available."
    saveSet(msgSet, "msgSet")
    return("NA")
  }
  dendro_df <- dendro_df[finite_rows, , drop = FALSE]
  dendro_df$Status <- ifelse(dendro_df$Dendrogram_Distance > threshold, "Outlier", "Normal")

  ## Decide label set
  out_idx <- which(dendro_df$Status == "Outlier")
  label_idx <- if (length(out_idx) <= 20) {
    out_idx
  } else {
    out_idx[order(dendro_df$Dendrogram_Distance[out_idx], decreasing = TRUE)[1:20]]
  }
  dendro_df$LabelMe <- FALSE
  dendro_df$LabelMe[label_idx] <- TRUE


  qc.dendrogram.plot(dendro_df, threshold, imgName, dpi, format)
  qc.dendrogram.json(dendro_df, imgName);
  return("NA")
}

qc.dendrogram.plot <- function(dendro_df,
                               threshold = 0.1,
                               imgNm = "Dendrogram_plot",
                               dpi = 72,
                               format = "png",
                               interactive = FALSE) {
  require(ggplot2)
  require(ggrepel)
  require(Cairo)

  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be positive.")

  set.seed(1)  # For reproducible jitter
  dendro_df$xj <- jitter(rep(1, nrow(dendro_df)), amount = 0.25)

  g <- ggplot(dendro_df, aes(x = xj, y = Dendrogram_Distance)) +
    geom_boxplot(aes(x = 1), outlier.shape = NA,
                 width = 0.4, fill = "grey80") +
    geom_point(aes(color = Status), size = 2.2) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
    geom_text_repel(data = dendro_df[dendro_df$LabelMe, ],
                    aes(label = Sample),
                    max.overlaps = Inf,
                    box.padding = 0.35,
                    segment.size = 0.2,
                    size = 4.2) +
    scale_color_manual(values = c(Normal = "grey40", Outlier = "red"),
                       name = "Sample status") +
    theme_minimal(base_size = 12) +
    labs(x = NULL, y = "Max pair-wise distance (1 − Pearson ρ)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  outFile <- paste0(imgNm, "dpi", dpi, ".", format)

  if (interactive) {
    require(plotly)
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = outFile, width = 8, height = 6,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return(outFile)
  }
}



GetSummaryTable <- function(dataName){
  dataSet <- readDataset(dataName)
  df <- dataSet$summary_df;
  df_rounded <- df;
  df_rounded[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], signif, digits = 3)
  return(df_rounded);
}

calculate_gini <- function(x) {
  n <- length(x)
  sorted_x <- sort(x)
  index <- 1:n
  gini <- (2 * sum(index * sorted_x) / sum(sorted_x)) - (n + 1)
  return(gini / n)
}


ComputePERMANOVA <- function(pc1, pc2, cls, numPermutations = 999) {
  # Combine PC1 and PC2 scores into a matrix
  pc.mat <- cbind(pc1, pc2)
  
  # Calculate PERMANOVA significance
  res <- .calculateDistSig(pc.mat, cls)
  
  # Extract the main results
  resTab <- res[[1]][1, ]
  
  # Format and create the PERMANOVA summary statistics
  stat.info <- paste("[PERMANOVA] F-value: ", signif(resTab$F, 5),
                     "; R-squared: ", signif(resTab$R2, 5),
                     "; p-value (based on ", numPermutations, " permutations): ",
                     signif(resTab$Pr, 5), sep = "")
  
  # Create a named vector for the statistics
  stat.info.vec <- c(F_value = signif(resTab$F, 5), 
                     R_squared = signif(resTab$R2, 5), 
                     p_value = signif(resTab$Pr, 5))
  names(stat.info.vec) <- c("F-value", "R-squared", "p-value");

  # Extract pairwise PERMANOVA results if available
  pair.res <- res[[2]]
  
  # Return the results as a list
  list(
    stat.info = stat.info,
    stat.info.vec = stat.info.vec,
    pair.res = pair.res
  )
}

# use a PERMANOVA to partition the euclidean distance by groups based on current score plot:
.calculateDistSig <- function(pc.mat, grp){

    data.dist <- dist(as.matrix(pc.mat), method = 'euclidean');
    res <- vegan::adonis2(formula = data.dist ~ grp);

    # pairwise for multi-grp
    if(length(levels(grp)) > 2){
      pair.res <- .permanova_pairwise(x = data.dist, grp);
      rownames(pair.res) <- pair.res$pairs;
      pair.res$pairs <- NULL;
      pair.res <- signif(pair.res,5);
      fast.write.csv(pair.res, file="pca_pairwise_permanova.csv");
    }else{
      pair.res <- NULL;
    }

    return(list(res, pair.res));
}

###adopted from ecole package https://rdrr.io/github/phytomosaic/ecole/
.permanova_pairwise <- function(x,
                                 grp,
                                 permutations = 999,
                                 method = 'bray',
                                 padj = 'fdr', ...) {
  f     <- grp
  if (!all(table(f) > 1)) warning('factor has singletons! perhaps lump them?')
  co    <- combn(unique(as.character(f)),2)
  nco   <- NCOL(co)
  out   <- data.frame(matrix(NA, nrow=nco, ncol=5))
  dimnames(out)[[2]] <- c('pairs', 'SumOfSqs', 'F.Model', 'R2', 'pval')
  if (!inherits(x, 'dist')) {
    D <- vegan::vegdist(x, method=method)
  } else {
    D <- x
  }
  #cat('Now performing', nco, 'pairwise comparisons. Percent progress:\n')
  for(j in 1:nco) {
    cat(round(j/nco*100,0),'...  ')
    ij  <- f %in% c(co[1,j],co[2,j])
    Dij <- as.dist(as.matrix(D)[ij,ij])
    fij <- data.frame(fij = f[ij])
    a   <- vegan::adonis2(Dij ~ fij, data=fij, permutations = permutations, ...);
    out[j,1] <- paste(co[1,j], 'vs', co[2,j])
    out[j,2] <- a$SumOfSqs[1]
    out[j,3] <- a$F[1]
    out[j,4] <- a$R2[1]
    out[j,5] <- a$`Pr(>F)`[1]
  }
  #cat('\n')
  out$p.adj <- p.adjust(out$pval, method=padj)
  out$SumOfSqs <-NULL
  #attr(out, 'p.adjust.method') <- padj
  #cat('\np-adjust method:', padj, '\n\n');
  return(out)
}

qc.pcaplot.json <- function(dataSet, x, imgNm) {
  jsonFile <- paste0(imgNm, ".json")

  # libs (only jsonlite is really needed now)
  suppressMessages({
    require(jsonlite)
  })

  # ---------- Load PCA & metadata ----------
  analSet <- readSet(analSet, "analSet")
  pca     <- analSet$pca
  xlabel  <- pca$xlabel
  ylabel  <- pca$ylabel

  # Safety check: ensure at least 2 PCs exist
  if (ncol(pca$x) < 2) {
    stop("Insufficient principal components for PCA plot (need >= 2, found ", ncol(pca$x), ")")
  }

  pca.res <- as.data.frame(pca$x)[, 1:2, drop = FALSE]
  colnames(pca.res) <- c("PC1", "PC2")

  match_indices <- match(rownames(dataSet$meta.info), rownames(pca.res))

  pca.res <- pca.res[match_indices, , drop = FALSE]

  # metadata1 → color (group)
  pca.res$group  <- as.character(dataSet$meta.info[[1]])
  pca.res$sample <- rownames(pca.res)

  # ---------- Detect 2nd metadata for shapes ----------
  doShape <- FALSE
  shape.levels <- character(0)
  shape.map <- NULL

  if (ncol(dataSet$meta.info) >= 2) {
    second <- dataSet$meta.info[[2]]
    # treat non-numeric as discrete for shapes
    isDisc  <- !is.numeric(second)
    levs    <- unique(as.character(second))
    if (isDisc && length(levs) <= 8) {
      doShape <- TRUE
      pca.res$shape <- as.character(second)
      symbols <- c("circle","square","diamond",
                   "cross","x","triangle-up",
                   "triangle-down","star")
      shape.map    <- stats::setNames(symbols[seq_along(levs)], levs)
      shape.levels <- levs
    }
  }

  # ---------- Color mapping (dose-aware like qc.pcaplot) ----------
  paramSet    <- readSet(paramSet, "paramSet")
  unique_grps <- unique(pca.res$group)

  if (grepl("norm", imgNm) &&
      !is.null(paramSet$oneDataAnalType) &&
      paramSet$oneDataAnalType == "dose") {
    pal <- grDevices::colorRampPalette(c("#2196F3", "#DE690D"))(length(unique_grps))
  } else {
    pal <- unname(get.qc.group.palette(unique_grps))
  }
  col.map <- stats::setNames(pal, unique_grps)

  # ---------- Build traces ----------
  traces <- list()

  if (doShape) {
    # one trace per (group × shape)
    combos <- unique(pca.res[, c("group","shape")])
    for (i in seq_len(nrow(combos))) {
      g  <- combos$group[i]
      sh <- combos$shape[i]

      df <- pca.res[pca.res$group == g & pca.res$shape == sh, , drop = FALSE]
      if (nrow(df) == 0) {
        next
      }

      mkr <- list(
        color = unname(col.map[g]),
        size  = 8,
        line  = list(color = "white", width = 0.5),
        symbol = unname(shape.map[[sh]]) # scalar symbol per trace
      )

      traces[[length(traces) + 1]] <- list(
        x            = df$PC1,
        y            = df$PC2,
        type         = "scatter",
        mode         = if (nrow(df) > 20) "markers" else "markers+text",
        name         = paste0(g, " • ", sh),
        legendgroup  = g,          # groups align in legend
        marker       = mkr,
        text         = if (nrow(df) <= 20) df$sample else NULL,
        hoverinfo    = "text",
        textposition = "top center"
      )
    }
  } else {
    # one trace per color group (no shapes)
    for (g in unique_grps) {
      df <- pca.res[pca.res$group == g, , drop = FALSE]
      if (nrow(df) == 0) {
        next
      }

      mkr <- list(
        color = unname(col.map[g]),
        size  = 8,
        line  = list(color = "white", width = 0.5)
      )

      traces[[length(traces) + 1]] <- list(
        x            = df$PC1,
        y            = df$PC2,
        type         = "scatter",
        mode         = if (nrow(df) > 20) "markers" else "markers+text",
        name         = g,
        legendgroup  = g,
        marker       = mkr,
        text         = if (nrow(df) <= 20) df$sample else NULL,
        hoverinfo    = "text",
        textposition = "top center"
      )
    }
  }

  # ---------- Layout ----------
  layout <- list(
    title = "",
    xaxis = list(title = xlabel, zeroline = FALSE),
    yaxis = list(title = ylabel, zeroline = FALSE),
    legend = list(
      orientation = "v",
      x           = 1.02,
      y           = 1,
      xanchor     = "left",
      yanchor     = "top"
    )
  )

  # ---------- Dump JSON ----------
  plot_data <- list(data = traces, layout = layout)
  json.obj  <- jsonlite::toJSON(plot_data, auto_unbox = TRUE, null = "null", digits = NA)
  writeLines(json.obj, jsonFile)

  return("NA")
}


PlotDataGini <- function(fileName, imgName, threshold, dpi, format){
  dataSet <- readDataset(fileName)
  if (is.null(dataSet$summary_df)) {
    stop("summary_df not found in dataSet. Please run SummarizeQC first.")
  }
  
  # Select Gini data
  gini_df <- dataSet$summary_df[, c("Sample", "Gini")]
  gini_df$Status <- ifelse(gini_df$Gini > threshold, "Outlier", "Normal")
  
  ## Plot
  qc.gini.plot(gini_df, imgName, threshold, dpi, format)
  qc.giniplot.json(gini_df, imgName);
  return("NA")
}

qc.gini.plot <- function(gini_df,
                         imgNm   = "Gini_plot",
                         threshold = 0.95,
                         dpi     = 72,
                         format  = "png",
                         interactive = FALSE) {
  require(ggplot2)
  require(ggrepel)
  require(Cairo)
  
  dpi <- as.numeric(dpi)
  if (dpi <= 0) stop("DPI must be a positive number.")
  
  g <- ggplot(gini_df, aes(x = "", y = Gini)) +
    geom_boxplot(outlier.shape = NA, fill = "grey80") +
    geom_jitter(aes(color = Status), width = 0.25, height = 0) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(gini_df, Status == "Outlier"),
                    aes(label = Sample), nudge_x = 0.35, size = 3) +
    scale_color_manual(values = c(Normal = "grey40", Outlier = "red"),
                       name = "Sample status") +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = "Gini coefficient") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  width  <- 8
  height <- 6
  fileNm <- paste(imgNm, "dpi", dpi, ".", sep = "")
  imgNm  <- paste0(fileNm, format)
  
  if (interactive) {
    require(plotly)
    m <- list(l = 50, r = 50, b = 20, t = 20, pad = 0.5)
    return(layout(plotly::ggplotly(g),
                  autosize = FALSE, width = 1000, height = 600, margin = m))
  } else {
    if (dpi == 72) dpi <- dpi * 1.34
    Cairo(file = imgNm, width = width, height = height,
          type = format, bg = "white", dpi = dpi, unit = "in")
    print(g)
    dev.off()
    return("NA")
  }
}

SummarizeQC <- function(fileName, imgNameBase, threshold = 0.1) {
  # save.image("summarize.RData");
  dataSet <- readDataset(fileName)

  summary_df <- data.frame(
    Sample = character(),
    HighCoverageGeneCount = numeric(),
    NSig80 = numeric(),
    Gini = numeric(),
    Dendrogram_Distance = numeric(),
    Outlier_HighCoverageGeneCount = numeric(),
    Outlier_NSig80 = numeric(),
    Outlier_Gini = numeric(),
    Outlier_Dendrogram = numeric(),
    stringsAsFactors = FALSE
  )

  if (grepl("norm", imgNameBase)) {
    data <- dataSet$data.norm
  } else {
    data <- .get.annotated.data()
  }

  ## --- basic per-sample metrics ---
  HighCoverageGeneCount <- colSums(data > 5, na.rm = TRUE)
  ncov5_df <- data.frame(
    Sample = names(HighCoverageGeneCount),
    HighCoverageGeneCount = as.numeric(HighCoverageGeneCount),
    stringsAsFactors = FALSE
  )

  NSig80 <- apply(data, 2, function(col) {
    col[is.na(col)] <- 0
    s <- sum(col)
    if (s <= 0) return(0)
    sum(cumsum(sort(col, decreasing = TRUE)) <= 0.8 * s)
  })
  nsig_df <- data.frame(
    Sample = names(NSig80),
    NSig80 = as.numeric(NSig80),
    stringsAsFactors = FALSE
  )

  gini_scores <- apply(data, 2, calculate_gini)
  gini_df <- data.frame(
    Sample = colnames(data),
    Gini = as.numeric(gini_scores),
    stringsAsFactors = FALSE
  )

  ## --- correlation distance (0..1) & within-group mean distance ---
  pearson_corr <- suppressWarnings(cor(data, method = "pearson", use = "pairwise.complete.obs"))
  # keep diagonal sane and replace NA correlations
  if (!is.null(pearson_corr)) {
    diag(pearson_corr) <- 1
    pearson_corr[is.na(pearson_corr)] <- 0
  }
  distance_matrix <- as.dist((1 - pearson_corr) / 2)  # 0..1
  dist_mat <- as.matrix(distance_matrix)

  group_info <- dataSet$meta.info[, 1]
  names(group_info) <- rownames(dataSet$meta.info)

  mean_distances <- sapply(colnames(data), function(sample) {
    sample_group <- group_info[sample]
    same_group_samples <- names(group_info)[group_info == sample_group]
    same_group_samples <- setdiff(same_group_samples, sample)
    # ensure they exist in the distance matrix
    same_group_samples <- intersect(same_group_samples, colnames(data))
    if (length(same_group_samples) == 0 || is.null(dist_mat)) return(NA_real_)
    mean(dist_mat[sample, same_group_samples], na.rm = TRUE)
  })

  dendrogram_df <- data.frame(
    Sample = names(mean_distances),
    Dendrogram_Distance = as.numeric(mean_distances),
    stringsAsFactors = FALSE
  )

  ## --- Merge first (align by Sample) ---
  summary_df <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE),
                       list(ncov5_df, nsig_df, gini_df, dendrogram_df))

  ## --- Outlier calls on the merged columns ---
  # NSig80 IQR rule (3*IQR)
  Q1_nsig <- quantile(summary_df$NSig80, 0.25, na.rm = TRUE)
  Q3_nsig <- quantile(summary_df$NSig80, 0.75, na.rm = TRUE)
  IQR_nsig <- IQR(summary_df$NSig80, na.rm = TRUE)
  summary_df$Outlier_NSig80 <- as.integer(
    summary_df$NSig80 < (Q1_nsig - 3 * IQR_nsig) |
      summary_df$NSig80 > (Q3_nsig + 3 * IQR_nsig)
  )

  # HighCoverageGeneCount IQR rule (3*IQR)
  Q1_cov <- quantile(summary_df$HighCoverageGeneCount, 0.25, na.rm = TRUE)
  Q3_cov <- quantile(summary_df$HighCoverageGeneCount, 0.75, na.rm = TRUE)
  IQR_cov <- IQR(summary_df$HighCoverageGeneCount, na.rm = TRUE)
  summary_df$Outlier_HighCoverageGeneCount <- as.integer(
    summary_df$HighCoverageGeneCount < (Q1_cov - 3 * IQR_cov) |
      summary_df$HighCoverageGeneCount > (Q3_cov + 3 * IQR_cov)
  )

  # Gini hard cutoff
  summary_df$Outlier_Gini <- as.integer(summary_df$Gini > 0.95)

  # Dendrogram distance threshold (0..1 scale)
  summary_df$Outlier_Dendrogram <- as.integer(summary_df$Dendrogram_Distance > threshold)

  ## --- finalize & register ---
  # Optional: stable ordering by Sample
  summary_df <- summary_df[order(summary_df$Sample), , drop = FALSE]

  dataSet$summary_df <- summary_df
  RegisterData(dataSet)

  return(1)
}


# -------------------------------------------------------------------------
#  qc.giniplot.json()
#  ------------------------------------------------------------------------
#  gini_df      data.frame with columns: Sample, Gini, Status
#  imgNm        stem for the JSON file (".json" is appended automatically)
#  threshold    horizontal dashed reference line
#  jitter.w     half-width of horizontal jitter (0–0.5 recommended)
# -------------------------------------------------------------------------
qc.giniplot.json <- function(gini_df,
                             imgNm     = "Gini_plot",
                             threshold = 0.95,
                             jitter.w  = 0.45) {

  stopifnot(all(c("Sample", "Gini", "Status") %in% names(gini_df)))

  ## 1 · Tukey fences & statistical-outlier flag -------------------------
  finite_gini <- is.finite(gini_df$Gini)
  if (!any(finite_gini)) {
    warning("qc.giniplot.json: no finite Gini values; skipping plot")
    return(invisible("NA"))
  }
  stats      <- boxplot.stats(gini_df$Gini[finite_gini], coef = 1.5)$stats
  q1         <- stats[2]; q3 <- stats[4]; iqr <- q3 - q1
  lowFence   <- q1 - 1.5 * iqr
  highFence  <- q3 + 1.5 * iqr
  gini_df$stat_out <- with(gini_df, Gini < lowFence | Gini > highFence)
  gini_df$stat_out[is.na(gini_df$stat_out)] <- FALSE
  non_out <- sum(!gini_df$stat_out)

  ## 2 · Semantic palette (Normal / Outlier) -----------------------------
  status_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## 3 · Traces ----------------------------------------------------------
  # 3a · Box built from in-fence values only
  tr_box <- list(
    x              = rep(0, non_out),
    y              = I(gini_df$Gini[!gini_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b · Invisible all-points trace (for autoscale)
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(gini_df), -jitter.w, jitter.w)),
    y          = I(gini_df$Gini),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c · Visible points (semantic colour, stat-outlines)
  show_labels <- nrow(gini_df) <= 20
  set.seed(2)
  points_trace <- list(
    x    = I(runif(nrow(gini_df), -jitter.w, jitter.w)),
    y    = I(gini_df$Gini),
    type = "scatter",
    mode = if (show_labels) "markers+text" else "markers",
    text = if (show_labels) gini_df$Sample else "",
    textposition = "right",
    name = "Samples",
    hoverinfo = "text",
    hovertext = paste0(
      "Sample: ", gini_df$Sample,
      "<br>Gini: ", signif(gini_df$Gini, 3),
      "<br>Status: ", gini_df$Status
    ),
    marker = list(
      color = status_cols[gini_df$Status],
      size  = 8,
      line  = list(
        color = ifelse(gini_df$stat_out, "black", "rgba(0,0,0,0)"),
        width = ifelse(gini_df$stat_out, 1, 0)
      )
    ),
    showlegend = FALSE
  )

  traces <- list(tr_box, tr_all, points_trace)

  ## 4 · Layout ----------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF",
    paper_bgcolor = "#FFFFFF",
    xaxis = list(
      title          = "",
      range          = c(-jitter.w - 0.1, jitter.w + 0.1),
      zeroline       = FALSE,
      showticklabels = FALSE,
      showline       = TRUE,
      linecolor      = "#000000"
    ),
    yaxis = list(
      title     = list(text = "Gini coefficient"),
      zeroline  = FALSE,
      ticks     = "outside",
      showline  = TRUE,
      linecolor = "#000000",
      showgrid  = TRUE,
      gridcolor = "rgba(200,200,200,0.4)"
    ),
    shapes = list(list(
      type  = "line",
      xref  = "paper", x0 = 0, x1 = 1,
      yref  = "y",     y0 = threshold, y1 = threshold,
      line  = list(color = "#0026FF", dash = "dot")
    )),
    legend = list(
      title       = list(text = "Sample Status"),
      orientation = "v",
      x = 1.02, y = 1,
      xanchor = "left", yanchor = "top"
    ),
    margin = list(l = 60, r = 110, t = 20, b = 40)
  )

  ## 5 · Write JSON ------------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )

  invisible("NA")
}


# -------------------------------------------------------------------------
#  dendro_df   data.frame with columns:
#              Sample, Dendrogram_Distance, Status (Normal / Outlier),
#              LabelMe (TRUE/FALSE → label on plot)
#  imgNm       stem for JSON file ("<imgNm>.json")
#  threshold   horizontal dashed cut-off
#  jitter.w    half-width for horizontal jitter of points
# -------------------------------------------------------------------------
qc.dendrogram.json <- function(dendro_df,
                               imgNm     = "Dendrogram_plot",
                               threshold = 0.10,
                               jitter.w  = 0.45) {

  stopifnot(all(c("Sample", "Dendrogram_Distance", "Status", "LabelMe") %in% names(dendro_df)))

  ## ── 1 · Tukey fences and statistical-outlier flag -------------------
  finite_dendro <- is.finite(dendro_df$Dendrogram_Distance)
  if (!any(finite_dendro)) {
    warning("qc.dendrogram.json: no finite dendrogram distances; skipping plot")
    return(invisible("NA"))
  }
  stats     <- boxplot.stats(dendro_df$Dendrogram_Distance[finite_dendro], coef = 1.5)$stats
  q1        <- stats[2];  q3 <- stats[4];  iqr <- q3 - q1
  lowFence  <- q1 - 1.5 * iqr
  highFence <- q3 + 1.5 * iqr
  dendro_df$stat_out <- with(dendro_df,
                             Dendrogram_Distance < lowFence |
                             Dendrogram_Distance > highFence)
  dendro_df$stat_out[is.na(dendro_df$stat_out)] <- FALSE
  non_out <- sum(!dendro_df$stat_out)

  ## ── 2 · Semantic palette -------------------------------------------
  status_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## ── 3 · Traces ------------------------------------------------------
  # 3a · Box from in-fence points
  tr_box <- list(
    x              = rep(0, non_out),
    y              = I(dendro_df$Dendrogram_Distance[!dendro_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b · Invisible all-points scatter (autoscale helper)
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(dendro_df), -jitter.w, jitter.w)),
    y          = I(dendro_df$Dendrogram_Distance),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c · Visible points (semantic colouring, stat outlines)
  set.seed(2)
  points_trace <- list(
    x    = I(runif(nrow(dendro_df), -jitter.w, jitter.w)),
    y    = I(dendro_df$Dendrogram_Distance),
    type = "scatter",
    mode = "markers+text",
    text = ifelse(dendro_df$LabelMe, dendro_df$Sample, ""),
    textposition = "right",
    name = "Samples",
    hoverinfo = "text",
    hovertext = paste0(
      "Sample: ", dendro_df$Sample,
      "<br>Distance: ", signif(dendro_df$Dendrogram_Distance, 3),
      "<br>Status: ", dendro_df$Status
    ),
    marker = list(
      color = status_cols[dendro_df$Status],
      size  = 8,
      line  = list(
        color = ifelse(dendro_df$stat_out, "black", "rgba(0,0,0,0)"),
        width = ifelse(dendro_df$stat_out, 1, 0)
      )
    ),
    showlegend = FALSE
  )

  traces <- list(tr_box, tr_all, points_trace)

  ## ── 4 · Layout ------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF",
    paper_bgcolor = "#FFFFFF",
    xaxis = list(
      title          = "",
      range          = c(-jitter.w - 0.1, jitter.w + 0.1),
      zeroline       = FALSE,
      showticklabels = FALSE,
      showline       = TRUE,
      linecolor      = "#000000"
    ),
    yaxis = list(
      title     = list(text = "Max pair-wise distance (1 \u2212 Pearson \u03c1)"),
      zeroline  = FALSE,
      ticks     = "outside",
      showline  = TRUE,
      linecolor = "#000000",
      showgrid  = TRUE,
      gridcolor = "rgba(200,200,200,0.4)"
    ),
    shapes = list(list(
      type  = "line",
      xref  = "paper", x0 = 0, x1 = 1,
      yref  = "y",     y0 = threshold, y1 = threshold,
      line  = list(color = "#0026FF", dash = "dot")
    )),
    legend = list(
      title       = list(text = "Sample Status"),
      orientation = "v",
      x = 1.02, y = 1,
      xanchor = "left", yanchor = "top"
    ),
    margin = list(l = 70, r = 110, t = 20, b = 40)
  )

  ## ── 5 · Write JSON ---------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )
  invisible("NA")
}

qc.ncov5plot.json <- function(ncov5_df,
                              imgNm    = "NCov5_plot",
                              lower,
                              upper,
                              jitter.w = 0.45) {

  stopifnot(all(c("Sample", "HighCoverageGeneCount", "Status") %in% names(ncov5_df)),
            is.numeric(lower), length(lower) == 1,
            is.numeric(upper), length(upper) == 1)

  ## ── 1 · Tukey fences & statistical-outlier flag --------------------
  finite_cov <- is.finite(ncov5_df$HighCoverageGeneCount)
  if (!any(finite_cov)) {
    warning("qc.ncov5plot.json: no finite coverage counts; skipping plot")
    return(invisible("NA"))
  }
  stats      <- boxplot.stats(ncov5_df$HighCoverageGeneCount[finite_cov], coef = 1.5)$stats
  q1         <- stats[2]; q3 <- stats[4]; iqr <- q3 - q1
  lowFence   <- q1 - 1.5 * iqr
  highFence  <- q3 + 1.5 * iqr
  ncov5_df$stat_out <- with(ncov5_df,
                             HighCoverageGeneCount < lowFence |
                             HighCoverageGeneCount > highFence)
  ncov5_df$stat_out[is.na(ncov5_df$stat_out)] <- FALSE
  non_out <- sum(!ncov5_df$stat_out)

  ## ── 2 · palettes ----------------------------------------------------
  stat_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## ── 3 · Traces ------------------------------------------------------
  # 3a · box (only in-fence points)
  tr_box <- list(
    x              = rep(0, non_out),
    y              = I(ncov5_df$HighCoverageGeneCount[!ncov5_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b · invisible all-points scatter (for autoscale)
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(ncov5_df), -jitter.w, jitter.w)),
    y          = I(ncov5_df$HighCoverageGeneCount),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c · labelled points (semantic status colouring, stat-outlines)
  # 3c · visible points (one trace for all samples)
set.seed(2)
points_trace <- list(
  x    = I(runif(nrow(ncov5_df), -jitter.w, jitter.w)),
  y    = I(ncov5_df$HighCoverageGeneCount),
  type = "scatter",
  mode = if (any(ncov5_df$stat_out)) "markers+text" else "markers",
  text = ifelse(ncov5_df$stat_out, ncov5_df$Sample, ""),
  textposition = "right",
  name = "Samples",
  hoverinfo = "text",
  hovertext = paste0(
    "Sample: ", ncov5_df$Sample,
    "<br>Count: ", ncov5_df$HighCoverageGeneCount,
    "<br>Status: ", ncov5_df$Status
  ),
  marker = list(
    color = stat_cols[ncov5_df$Status],           # vector OK
    size  = 8,
    line  = list(
      color = ifelse(ncov5_df$stat_out, "black", "rgba(0,0,0,0)"),
      width = ifelse(ncov5_df$stat_out, 1, 0)
    )
  ),
  showlegend = FALSE
)

  traces <- list(tr_box, tr_all, points_trace)

  ## ── 4 · Layout ------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF",
    paper_bgcolor = "#FFFFFF",
    xaxis = list(
      title          = "",
      range          = c(-jitter.w - 0.1, jitter.w + 0.1),
      zeroline       = FALSE,
      showticklabels = FALSE,
      showline       = TRUE,
      linecolor      = "#000000"
    ),
    yaxis = list(
      title     = list(text = "features with > 5 uniquely mapped reads"),
      zeroline  = FALSE,
      ticks     = "outside",
      showline  = TRUE,
      linecolor = "#000000",
      showgrid  = TRUE,
      gridcolor = "rgba(200,200,200,0.4)"
    ),
    shapes = list(
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=lower, y1=lower,
           line=list(color="#0026FF", dash="dot")),
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=upper, y1=upper,
           line=list(color="#0026FF", dash="dot"))
    ),
    legend = list(
      title       = list(text = "Sample Status"),
      orientation = "v",
      x = 1.02, y = 1,
      xanchor = "left", yanchor = "top"
    ),
    margin = list(l = 70, r = 110, t = 20, b = 40)
  )

  ## ── 5 · Write JSON ---------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )

  invisible("NA")
}

qc.nsigplot.json <- function(nsig_df,
                             imgNm    = "NSig80_plot",
                             lower,
                             upper,
                             jitter.w = 0.45) {

  stopifnot(all(c("Sample", "NSig80", "outlier") %in% names(nsig_df)))

  ## ------------------------------------------------------------------
  ## 1 · Compute Tukey fences and flag statistical outliers
  ## ------------------------------------------------------------------
  finite_nsig <- is.finite(nsig_df$NSig80)
  if (!any(finite_nsig)) {
    warning("qc.nsigplot.json: no finite NSig80 values; skipping plot")
    return(invisible("NA"))
  }
  stats       <- boxplot.stats(nsig_df$NSig80[finite_nsig], coef = 1.5)$stats
  q1          <- stats[2]; q3 <- stats[4]; iqr <- q3 - q1
  lowFence    <- q1 - 1.5 * iqr
  highFence   <- q3 + 1.5 * iqr
  nsig_df$stat_out <- with(nsig_df, NSig80 < lowFence | NSig80 > highFence)
  nsig_df$stat_out[is.na(nsig_df$stat_out)] <- FALSE
  non_out <- sum(!nsig_df$stat_out)

  ## ------------------------------------------------------------------
  ## 2 · Palette for semantic status (Normal / Outlier)
  ## ------------------------------------------------------------------
  status_cols <- c(Normal = "#666666", Outlier = "#E41A1C")

  ## ------------------------------------------------------------------
  ## 3 · Plotly traces
  ## ------------------------------------------------------------------
  # 3a ─ Box: only in-fence points
  tr_box <- list(
    x              = rep(0, non_out),
    y              = I(nsig_df$NSig80[!nsig_df$stat_out]),
    quartilemethod = "linear",
    type           = "box",
    width          = 0.8,
    name           = "",
    boxpoints      = FALSE,
    fillcolor      = "rgba(200,200,200,0.6)",
    line           = list(color = "#000000"),
    hoverinfo      = "skip",
    showlegend     = FALSE
  )

  # 3b ─ Invisible “all” trace for autoscaling
  set.seed(1)
  tr_all <- list(
    x          = I(runif(nrow(nsig_df), -jitter.w, jitter.w)),
    y          = I(nsig_df$NSig80),
    type       = "scatter",
    mode       = "markers",
    marker     = list(color = "rgba(0,0,0,0)", size = 0),
    hoverinfo  = "skip",
    showlegend = FALSE
  )

  # 3c ─ Points (semantic status colouring, stat-outliers outlined)
  set.seed(2)
  points_trace <- list(
    x    = I(runif(nrow(nsig_df), -jitter.w, jitter.w)),
    y    = I(nsig_df$NSig80),
    type = "scatter",
    mode = "markers+text",
    text = ifelse(nsig_df$stat_out, nsig_df$Sample, ""),
    textposition = "right",
    name = "Samples",
    hoverinfo = "text",
    hovertext = paste0(
      "Sample: ", nsig_df$Sample,
      "<br>NSig80: ", nsig_df$NSig80,
      "<br>Status: ", nsig_df$outlier
    ),
    marker = list(
      color = status_cols[nsig_df$outlier],
      size  = 8,
      line  = list(
        color = ifelse(nsig_df$stat_out, "black", "rgba(0,0,0,0)"),
        width = ifelse(nsig_df$stat_out, 1, 0)
      )
    ),
    showlegend = FALSE
  )

  traces <- list(tr_box, tr_all, points_trace)

  ## ------------------------------------------------------------------
  ## 4 · Layout (unchanged)
  ## ------------------------------------------------------------------
  layout <- list(
    plot_bgcolor  = "#FFFFFF", paper_bgcolor = "#FFFFFF",
    xaxis = list(title="", range=c(-jitter.w-0.1, jitter.w+0.1),
                 zeroline=FALSE, showticklabels=FALSE,
                 showline=TRUE, linecolor="#000"),
    yaxis = list(title=list(text="NSig80 (features reaching 80% of signal)"),
                 zeroline=FALSE, ticks="outside", showline=TRUE,
                 linecolor="#000", showgrid=TRUE,
                 gridcolor="rgba(200,200,200,0.4)"),
    shapes = list(
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=lower, y1=lower,
           line=list(color="#0026FF", dash="dot")),
      list(type="line", xref="paper", x0=0, x1=1,
           yref="y", y0=upper, y1=upper,
           line=list(color="#0026FF", dash="dot"))
    ),
    legend = list(title=list(text="Sample Status"),
                  orientation="v", x=1.02, y=1,
                  xanchor="left", yanchor="top"),
    margin = list(l=70, r=110, t=20, b=40)
  )

  ## ------------------------------------------------------------------
  ## 5 · Write JSON
  ## ------------------------------------------------------------------
  jsonlite::write_json(
    list(data = traces, layout = layout),
    paste0(imgNm, ".json"),
    auto_unbox = TRUE, digits = 16
  )
}
qc.pcaplot.outliers.json <- function(dataSet, x, imgNm,
                                     uniq_map_col = "uniq_map",
                                     min_per_dose = 2,
                                     min_vehicles = 3) {
  jsonFile <- paste0(imgNm, ".json")
  csvFile  <- paste0(imgNm, "_outliers.csv")

  require(plotly)
  require(rjson)

  # ----- load PCA & align to meta -----
  analSet <- readSet(analSet, "analSet")
  pca     <- analSet$pca
  xlabel  <- pca$xlabel
  ylabel  <- pca$ylabel

  #msg("[PCA outliers] pca$x dimensions: ", nrow(pca$x), " x ", ncol(pca$x))

  # Safety check: ensure at least 2 PCs exist
  if (ncol(pca$x) < 2) {
    #msg("[PCA outliers] ERROR: Need at least 2 principal components, found ", ncol(pca$x))
    stop("Insufficient principal components for PCA outlier plot (need >= 2, found ", ncol(pca$x), ")")
  }

  pca.res <- as.data.frame(pca$x)[, 1:2, drop = FALSE]
  colnames(pca.res) <- c("PC1","PC2")
  pca.res <- pca.res[match(rownames(dataSet$meta.info), rownames(pca.res)), , drop = FALSE]
  pca.res$sample_id <- rownames(pca.res)

  meta <- dataSet$meta.info
  stopifnot(nrow(meta) == nrow(pca.res))
  pca.res$group <- as.character(meta[[1]])

  doShape <- FALSE; shape.map <- NULL; shape.levels <- NULL
  if (ncol(meta) >= 2) {
    second <- meta[[2]]
    isDisc <- !is.numeric(second)
    levs   <- unique(as.character(second))
    if (isDisc && length(levs) <= 6) {
      doShape <- TRUE
      pca.res$shape <- as.character(second)
      symbols <- c("circle","square","diamond","cross","x","triangle-up","triangle-down","star")
      shape.map <- setNames(symbols[seq_along(levs)], levs)
      shape.levels <- levs
    }
  }

  nR <- nrow(meta)
  pca.res$dose <- if ("dose" %in% colnames(meta)) as.character(meta[["dose"]]) else rep(NA_character_, nR)
  pca.res$is_vehicle <- if ("is_vehicle" %in% colnames(meta)) {
    as.logical(as.character(meta[["is_vehicle"]]))
  } else rep(FALSE, nR)
  pca.res$uniq_map <- if (uniq_map_col %in% colnames(meta)) as.numeric(meta[[uniq_map_col]]) else rep(NA_real_, nR)

  # ----- color mapping (your palettes) -----
  paramSet <- readSet(paramSet, "paramSet")
  unique_grps <- unique(pca.res$group)
  if (grepl("norm", imgNm) && !is.null(paramSet$oneDataAnalType) && paramSet$oneDataAnalType == "dose") {
    pal <- colorRampPalette(c("#2196F3", "#DE690D"))(length(unique_grps))
  } else {
    okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
    pal <- rep(okabe, length.out = length(unique_grps))
  }
  col.map <- setNames(pal, unique_grps)

  axis_class <- function(vals_named) {
    labs <- names(vals_named); vals <- as.numeric(vals_named)
    n <- length(vals)
    cls <- rep("none", n)

    med <- median(vals, na.rm = TRUE)
    md  <- mad(vals, constant = 1, na.rm = TRUE)
    if (is.na(md) || md == 0) {
      q <- quantile(vals, probs = c(0.10, 0.90), na.rm = TRUE, names = FALSE)
      core_idx <- vals >= q[1] & vals <= q[2]
    } else {
      core_idx <- abs(vals - med) <= 3 * md
    }
    if (sum(core_idx, na.rm = TRUE) < 3) {
      ord <- order(vals); k1 <- max(1, floor(0.10 * n)); k2 <- min(n, ceiling(0.90 * n))
      core_idx <- FALSE; core_idx[ord[k1:k2]] <- TRUE
    }

    core_vals <- vals[core_idx]
    span_core <- max(core_vals, na.rm = TRUE) - min(core_vals, na.rm = TRUE)
    if (!is.finite(span_core) || span_core == 0) span_core <- diff(range(vals, na.rm = TRUE))

    for (k in seq_len(n)) {
      xi <- vals[k]
      if (core_idx[k]) { cls[k] <- "none"; next }
      sep_core <- min(abs(xi - core_vals))
      if (sep_core > 2 * span_core)      cls[k] <- "strong"
      else if (sep_core > 1 * span_core) cls[k] <- "moderate"
      else                                cls[k] <- "none"
    }
    setNames(cls, labs)
  }

  euclid_med <- function(df) {
    cx <- median(df$PC1, na.rm = TRUE); cy <- median(df$PC2, na.rm = TRUE)
    sqrt((df$PC1 - cx)^2 + (df$PC2 - cy)^2)
  }
  within_dose_far <- function(df) {
    if (all(is.na(df$dose))) return(rep(FALSE, nrow(df)))
    unlist(
      by(df, df$dose, function(dd) {
        if (nrow(dd) < 2) return(rep(NA, nrow(dd)))
        cx <- median(dd$PC1); cy <- median(dd$PC2)
        d  <- sqrt((dd$PC1-cx)^2 + (dd$PC2-cy)^2)
        thr <- median(d, na.rm = TRUE) * 2
        d > thr
      }),
      use.names = FALSE
    )
  }

  .safe_chr <- function(x) { if (length(x) == 0 || is.na(x) || x == "NA") "" else as.character(x) }
  .safe_reason <- function(x) { if (length(x) == 0 || is.na(x) || x == "NA" || x == "") "" else as.character(x) }
  .make_customdata <- function(subdf) {
    n <- nrow(subdf)
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- list(
        dose       = .safe_chr(subdf$dose[i]),
        is_vehicle = .safe_chr(subdf$is_vehicle[i]),
        reason     = .safe_reason(subdf$reason[i])
      )
    }
    out
  }

  # trace appending guards
  .is_single_trace <- function(x) is.list(x) && !is.null(x$type) && !is.null(x$x) && !is.null(x$y)
  .append_trace <- function(dst, tr) {
    if (is.null(tr) || !is.list(tr)) return(dst)
    if (.is_single_trace(tr)) { dst[[length(dst) + 1]] <- tr; return(dst) }
    for (k in seq_along(tr)) {
      tk <- tr[[k]]
      if (.is_single_trace(tk)) dst[[length(dst) + 1]] <- tk
    }
    dst
  }

  df <- pca.res
  ax1 <- axis_class(setNames(df$PC1, df$sample_id))
  ax2 <- axis_class(setNames(df$PC2, df$sample_id))
  df$ax_PC1 <- ax1[df$sample_id]
  df$ax_PC2 <- ax2[df$sample_id]
  df$axis_class <- ifelse(df$ax_PC1 == "strong" | df$ax_PC2 == "strong", "strong",
                          ifelse(df$ax_PC1 == "moderate" | df$ax_PC2 == "moderate", "moderate", "none"))
  df$moderate_both_axes <- (df$ax_PC1 == "moderate" & df$ax_PC2 == "moderate")

  D  <- euclid_med(df)
  df$far_euclid <- D > (2 * median(D, na.rm = TRUE))
  df$far_repl  <- within_dose_far(df)

  df$reason  <- NA_character_
  df$exclude <- df$axis_class == "strong"
  df$reason[df$exclude] <- "Strong axis separation vs. core (>2× core span)"

  if (!all(is.na(df$dose))) {
    strong_rows <- which(df$exclude)
    if (length(strong_rows) > 1) {
      dups <- duplicated(df$dose[strong_rows]) | duplicated(df$dose[strong_rows], fromLast = TRUE)
      if (any(dups, na.rm = TRUE)) {
        df$exclude[strong_rows] <- FALSE
        df$reason[strong_rows]  <- NA_character_
      }
    }
  }

  m_idx <- which(!df$exclude & df$axis_class == "moderate")
  for (i in m_idx) {
    reasons <- character(0)
    if (isTRUE(df$moderate_both_axes[i]) && isTRUE(df$far_euclid[i]))
      reasons <- c(reasons, "Moderate on both axes with large Euclidean distance")
    if (isTRUE(df$far_repl[i]))
      reasons <- c(reasons, "Far from replicate/similar dose cluster")
    if (isTRUE(df$worse_qc[i]))
      reasons <- c(reasons, "Lower sequencing quality")
    if (length(reasons)) {
      df$exclude[i] <- TRUE
      df$reason[i]  <- paste(reasons, collapse = "; ")
    }
  }

  if (!all(is.na(df$dose))) {
    kept_by_dose <- tapply(!df$exclude & !df$is_vehicle, df$dose, sum)
    drop_doses <- names(kept_by_dose[!is.na(kept_by_dose) & kept_by_dose < min_per_dose])
    if (length(drop_doses)) {
      hit <- which(df$dose %in% drop_doses & !df$is_vehicle)
      df$exclude[hit] <- TRUE
      df$reason[hit]  <- ifelse(is.na(df$reason[hit]),
                                "Dose dropped (<2 samples after QC/outlier)",
                                paste(df$reason[hit], "Dose dropped (<2 samples)", sep="; "))
    }
  }

  veh_kept <- sum(!df$exclude & df$is_vehicle, na.rm = TRUE)
  vehicle_note <- if (any("is_vehicle" == colnames(meta)) && veh_kept < min_vehicles)
    sprintf("Warning: only %d vehicle samples kept (< %d).", veh_kept, min_vehicles) else NULL

  status_lab <- ifelse(df$exclude, "Outlier",
                       ifelse(df$axis_class == "moderate", "Moderate", "Kept"))
  df$.__status__ <- status_lab

  status_styles <- list(
    Kept     = list(line = list(color = "white", width = 0.5), size = 8,  opacity = 0.9),
    Moderate = list(line = list(color = "orange", width = 2),  size = 9,  opacity = 1.0),
    Outlier = list(line = list(color = "red",    width = 3),  size = 10, opacity = 1.0)
  )

  mk_trace <- function(subdf, name, color) {
    st <- as.character(unique(subdf$.__status__))[1]
    base_marker <- c(list(color = as.character(color)), status_styles[[st]])
    mode_val <- if (st == "Kept") "markers" else "markers+text"
    text_val <- if (st == "Kept") NULL else subdf$sample_id
    leg_name <- if (st == "Kept") name else paste0(name, " • ", st)

    if (doShape && "shape" %in% colnames(subdf)) {
      spl <- split(subdf, subdf$shape)
      out <- vector("list", length(spl)); i <- 0L
      for (sh in names(spl)) {
        i <- i + 1L
        ss <- spl[[sh]]
        mkr <- base_marker; mkr$symbol <- unname(shape.map[sh])
        shaped_name <- if (st == "Kept") paste0(name, " • ", sh) else paste0(name, " • ", sh, " • ", st)
        out[[i]] <- list(
          x = ss$PC1, y = ss$PC2, type = "scatter",
          mode = mode_val,
          name = shaped_name,
          showlegend = TRUE,
          marker = mkr,
          text = if (st == "Kept") NULL else ss$sample_id,
          customdata = .make_customdata(ss),  # list-of-objects
          hoverinfo = "text",
          textposition = "top center"
        )
      }
      return(out)
    }

    list(
      x = subdf$PC1, y = subdf$PC2, type = "scatter",
      mode = mode_val,
      name = leg_name,
      showlegend = TRUE,
      marker = base_marker,
      text = text_val,
      customdata = .make_customdata(subdf),  # list-of-objects
      hoverinfo = "text",
      textposition = "top center"
    )
  }

  traces <- list()
  for (g in unique_grps) {
    gdf <- df[df$group == g, , drop = FALSE]
    for (st in c("Kept","Moderate","Outlier")) {
      sdf <- gdf[gdf$.__status__ == st, , drop = FALSE]
      if (nrow(sdf) == 0) next
      tr <- mk_trace(sdf, g, col.map[[g]])
      traces <- .append_trace(traces, tr)
    }
  }

  if (doShape) {
    for (sh in shape.levels) {
      traces[[length(traces) + 1]] <- list(
        x = c(NA), y = c(NA), type = "scatter", mode = "markers",
        name = paste0("Shape: ", sh), showlegend = TRUE,
        marker = list(symbol = shape.map[[sh]], color = "black", size = 8)
      )
    }
  }

  subtitle <- if (!is.null(vehicle_note)) vehicle_note else ""
  layout <- list(
    title = "",
    xaxis = list(title = xlabel),
    yaxis = list(title = ylabel),
    legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left", yanchor = "top"),
    `shape.map` = shape.map,          # keep for JS legend helpers if you use them
    meta2Name = if (doShape) names(meta)[2] else NULL,
    annotations = if (nzchar(subtitle)) list(list(
      x = 0, y = 1.08, xref = "paper", yref = "paper",
      xanchor = "left", yanchor = "bottom",
      text = subtitle, showarrow = FALSE, font = list(size = 12)
    )) else NULL
  )

  plot_data <- list(data = traces, layout = layout)
  json.obj  <- toJSON(plot_data)
  sink(jsonFile); cat(json.obj); sink()

  out_cols <- c("sample_id","group","dose","is_vehicle","uniq_map",
                "PC1","PC2","ax_PC1","ax_PC2","axis_class",
                "moderate_both_axes","far_euclid","far_repl",
                "worse_qc","exclude","reason","__.__status__")
  keep_cols <- intersect(out_cols, colnames(df))
  out_tab   <- df[, keep_cols, drop = FALSE]
  colnames(out_tab)[colnames(out_tab) == ".__status__"] <- "status"
  utils::write.csv(out_tab, file = csvFile, row.names = FALSE)

  return("NA")
}
