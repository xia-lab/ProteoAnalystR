##################################################
## R script for ProteoAnalyst
## Description: plot violin plot for individual gene across metadata
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################
# given a gene id, plot its expression profile as violin plot
NormalizeFeatureId <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- as.character(x[1])
  x <- trimws(x)
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  x <- strsplit(x, ";", fixed = TRUE)[[1]][1]
  x <- trimws(x)
  if (grepl("\\|", x)) {
    parts <- strsplit(x, "\\|")[[1]]
    if (length(parts) >= 2 && nzchar(parts[2])) {
      x <- parts[2]
    }
  }
  x <- sub("_.*$", "", x)
  x <- sub("-\\d+$", "", x)
  x
}

ResolveFeatureRowId <- function(mat, feature.id) {
  if (is.null(mat) || is.null(rownames(mat))) return(feature.id)
  rn <- rownames(mat)
  if (feature.id %in% rn) return(feature.id)

  feature.norm <- NormalizeFeatureId(feature.id)
  if (!is.na(feature.norm)) {
    rn.norm <- vapply(rn, NormalizeFeatureId, character(1))
    idx <- which(rn.norm == feature.norm)
    if (length(idx) == 1) {
      return(rn[idx[1]])
    }
  }
  feature.id
}

GetGroupPalette <- function(groups, paletteOpt = "default") {
  groups <- as.factor(groups)
  paletteOpt <- tolower(as.character(paletteOpt)[1])
  if (identical(paletteOpt, "okabeito")) {
    n <- length(levels(groups))
    okabe_ito_base <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                        "#0072B2", "#D55E00", "#CC79A7", "#999999")
    if (n <= 8) {
      return(setNames(okabe_ito_base[seq_len(n)], levels(groups)))
    }
    cols <- grDevices::colorRampPalette(okabe_ito_base)(n)
    return(setNames(cols, levels(groups)))
  }
  if (identical(paletteOpt, "tol_bright")) {
    n <- length(levels(groups))
    tol_bright_base <- c("#4477AA", "#EE6677", "#228833", "#CCBB44",
                         "#66CCEE", "#AA3377", "#BBBBBB")
    if (n <= length(tol_bright_base)) {
      return(setNames(tol_bright_base[seq_len(n)], levels(groups)))
    }
    cols <- grDevices::colorRampPalette(tol_bright_base)(n)
    return(setNames(cols, levels(groups)))
  }
  cols <- unique(GetColorSchema(groups))
  names(cols) <- levels(groups)
  cols
}

PlotProteinPeptideOverview <- function(dataName = "", imageName = "", protein.id = "", format = "png", dpi = 96, paletteOpt = "default", plotType = "boxplot", maxPeptides = 35) {
  require(ggplot2)
  require(Cairo)

  if (is.null(protein.id) || !nzchar(protein.id)) {
    msg("[PlotProteinPeptideOverview] empty protein id")
    return(0)
  }

  plotType <- tolower(as.character(plotType)[1])
  if (!(plotType %in% c("violin", "boxplot"))) {
    plotType <- "boxplot"
  }
  use.violin <- plotType == "violin"

  dataSet <- readDataset(dataName)
  if (is.null(dataSet) || is.null(dataSet$data.norm)) {
    msg("[PlotProteinPeptideOverview] dataset or protein matrix is missing")
    return(0)
  }

  data.norm <- dataSet$data.norm
  if (length(dataSet$rmidx) > 0) {
    data.norm <- data.norm[, -dataSet$rmidx, drop = FALSE]
  }

  protein.row.id <- ResolveFeatureRowId(data.norm, protein.id)
  if (!(protein.row.id %in% rownames(data.norm))) {
    msg("[PlotProteinPeptideOverview] protein not found in protein matrix: ", protein.id)
    return(0)
  }

  peptide.file <- if (file.exists("peptide_level_data.qs")) {
    "peptide_level_data.qs"
  } else if (file.exists("data.stat.qs")) {
    "data.stat.qs"
  } else {
    ""
  }
  if (!nzchar(peptide.file) || !file.exists("peptide_to_protein_map.qs")) {
    msg("[PlotProteinPeptideOverview] peptide matrix or peptide-to-protein map is missing")
    return(0)
  }

  pep.mat <- ov_qs_read(peptide.file)
  pep.cache <- if (exists(".paGetPeptideFeatureCache", mode = "function")) {
    try(.paGetPeptideFeatureCache(dataName), silent = TRUE)
  } else {
    NULL
  }
  pep.map <- if (!inherits(pep.cache, "try-error") && !is.null(pep.cache) && !is.null(pep.cache$pep.map)) {
    pep.cache$pep.map
  } else {
    ov_qs_read("peptide_to_protein_map.qs")
  }
  if (is.null(pep.mat) || is.null(pep.map) || nrow(pep.map) == 0) {
    msg("[PlotProteinPeptideOverview] peptide inputs are empty")
    return(0)
  }

  pep.col <- if ("Peptide" %in% colnames(pep.map)) "Peptide" else colnames(pep.map)[1]
  prot.col <- if ("Protein" %in% colnames(pep.map)) "Protein" else colnames(pep.map)[2]
  prot.norm <- NormalizeFeatureId(protein.id)
  map.prot.norm <- if ("Protein.norm" %in% colnames(pep.map)) {
    as.character(pep.map$Protein.norm)
  } else {
    vapply(pep.map[[prot.col]], NormalizeFeatureId, character(1))
  }
  map.prot.raw <- trimws(as.character(pep.map[[prot.col]]))
  peps <- unique(as.character(pep.map[[pep.col]][map.prot.raw == protein.id | map.prot.norm == prot.norm]))
  peps <- peps[!is.na(peps) & nzchar(peps)]
  peps <- peps[peps %in% rownames(pep.mat)]
  if (length(peps) == 0) {
    msg("[PlotProteinPeptideOverview] no mapped peptides found for protein: ", protein.id)
  }

  peptide.count <- length(peps)
  if (peptide.count > maxPeptides) {
    peps <- peps[seq_len(maxPeptides)]
  }

  common.samples <- intersect(colnames(data.norm), colnames(pep.mat))
  if (length(common.samples) == 0) {
    common.samples <- colnames(data.norm)
  }
  protein.samples <- intersect(common.samples, colnames(data.norm))
  peptide.samples <- intersect(common.samples, colnames(pep.mat))

  meta <- dataSet$meta.info
  if (!is.null(meta) && nrow(meta) > 0) {
    group.col <- dataSet$analysisVar
    if (is.numeric(group.col) && length(group.col) > 0 && group.col[1] >= 1 && group.col[1] <= ncol(meta)) {
      group.col <- colnames(meta)[group.col[1]]
    }
    if (is.null(group.col) || !(group.col %in% colnames(meta))) {
      group.col <- colnames(meta)[1]
    }
    meta <- meta[rownames(meta) %in% common.samples, , drop = FALSE]
    sample.groups <- as.character(meta[match(common.samples, rownames(meta)), group.col])
  } else {
    sample.groups <- rep("All samples", length(common.samples))
  }
  sample.groups[is.na(sample.groups) | !nzchar(sample.groups)] <- "Unknown"
  group.map <- setNames(sample.groups, common.samples)

  make_label <- function(type, id) {
    id <- as.character(id)
    ifelse(nchar(id) > 64, paste0(substr(id, 1, 61), "..."), id)
  }

  protein.display <- protein.row.id
  protein.symbol <- tryCatch({
    paramSet <- readSet(paramSet, "paramSet")
    id.type <- tolower(paramSet$data.idType)
    org <- paramSet$data.org
    mapped <- if (id.type == "uniprot") {
      doUniprot2SymbolMapping(protein.row.id, org, paramSet$data.idType)
    } else {
      doEntrez2SymbolMapping(protein.row.id, org, paramSet$data.idType)
    }
    if (is.null(mapped) || is.na(mapped) || !nzchar(as.character(mapped)[1])) protein.row.id else as.character(mapped)[1]
  }, error = function(e) protein.row.id)
  protein.label <- make_label("Protein", protein.symbol)
  plot.df <- data.frame(
    Type = "Protein",
    Feature = protein.label,
    Sample = protein.samples,
    Group = group.map[protein.samples],
    Value = as.numeric(data.norm[protein.row.id, protein.samples]),
    stringsAsFactors = FALSE
  )

  if (length(peps) > 0) {
    pep.df <- do.call(rbind, lapply(peps, function(pep.id) {
      data.frame(
        Type = "Peptide",
        Feature = make_label("Peptide", pep.id),
        Sample = peptide.samples,
        Group = group.map[peptide.samples],
        Value = as.numeric(pep.mat[pep.id, peptide.samples]),
        stringsAsFactors = FALSE
      )
    }))
    plot.df <- rbind(plot.df, pep.df)
  }

  plot.df <- plot.df[!is.na(plot.df$Value) & !is.na(plot.df$Group), , drop = FALSE]
  if (nrow(plot.df) == 0) {
    msg("[PlotProteinPeptideOverview] no plottable values for protein: ", protein.id)
    return(0)
  }

  feature.levels <- unique(plot.df$Feature)
  plot.df$Feature <- factor(plot.df$Feature, levels = feature.levels)
  plot.df$Group <- factor(plot.df$Group)
  col <- GetGroupPalette(plot.df$Group, paletteOpt)
  make_axis_expr <- function(lbl, is.protein = FALSE) {
    lbl <- as.character(lbl)
    lbl <- gsub("\\\\", "\\\\\\\\", lbl)
    lbl <- gsub("\"", "\\\\\"", lbl)
    if (is.protein) {
      return(paste0("bold(\"", lbl, "\")"))
    }
    paste0("\"", lbl, "\"")
  }
  feature.labels <- setNames(
    vapply(as.character(feature.levels), function(lbl) make_axis_expr(lbl, identical(lbl, protein.label)), character(1)),
    as.character(feature.levels)
  )
  group.levels <- levels(plot.df$Group)
  group.count <- length(group.levels)
  dodge.width <- if (group.count > 1) 0.72 else 0

  n.feat <- length(feature.levels)
  band.df <- data.frame(
    xmin = c(-Inf, seq_len(n.feat - 1) + 0.5),
    xmax = c(seq_len(n.feat - 1) + 0.5, Inf),
    ymin = -Inf,
    ymax = Inf,
    Fill = rep(c("#d9d9d9", "#ececec"), length.out = n.feat),
    stringsAsFactors = FALSE
  )

  plot.geom <- if (use.violin) {
    geom_violin(aes(color = Group),
                trim = FALSE,
                position = position_dodge(width = dodge.width),
                width = if (group.count > 1) 0.66 else 0.72,
                show.legend = TRUE)
  } else {
    geom_boxplot(aes(color = Group),
                 outlier.shape = NA,
                 position = position_dodge(width = dodge.width),
                 width = if (group.count > 1) 0.56 else 0.62,
                 show.legend = TRUE)
  }

  myplot <- ggplot(plot.df, aes(x = Feature, y = Value, fill = Group)) +
    geom_rect(data = band.df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = band.df$Fill,
              alpha = 0.7) +
    plot.geom +
    geom_point(color = "black",
               position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0,
                                                dodge.width = dodge.width),
               size = 1.25,
               alpha = 0.65,
               show.legend = FALSE) +
    stat_summary(fun = mean, colour = "yellow", geom = "point",
                 shape = 18, size = 2.4,
                 position = position_dodge(width = dodge.width),
                 show.legend = FALSE) +
    scale_fill_manual(values = col) +
    scale_color_manual(values = col) +
    scale_x_discrete(labels = parse(text = feature.labels)) +
    theme_gray(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.ticks.x = element_line(colour = "#555555"),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "white", size = 0.35),
      panel.background = element_rect(fill = "#e5e5e5", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_blank()
    ) +
    ylab("Expression")

  imgName <- paste(imageName, "dpi", dpi, ".", format, sep = "")
  plot.width <- min(20, max(9.2, 4.2 + 0.55 * length(feature.levels)))
  Cairo(file = imgName, width = plot.width, height = 5.8, unit = "in",
        dpi = dpi, bg = "white", type = format)
  invisible(print(myplot))
  invisible(dev.off())
  return(0)
}

PlotSelectedGene <-function(dataName="",imageName="", gene.id="", type="notvolcano", format="png", dpi=96, fc = T, plotType = "boxplot", dataType = "default", paletteOpt = "default"){

  require(see)
  require(ggplot2)
  require(lattice)


  dataType <- tolower(dataType)
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  dataSet <- readDataset(dataName);
  anal.type <- paramSet$anal.type;

  # Protein list analysis doesn't have conditions/groups, so boxplot/violin doesn't make sense
  if(anal.type == "proteinlist"){
    msg("[PlotSelectedGene] Skipping plot for proteinlist - no group comparisons available")
    return(0)
  }

  plotType <- tolower(plotType)
  if (!(plotType %in% c("violin", "boxplot"))) {
    plotType <- "boxplot"
  }
  use.violin <- plotType == "violin"
  imgName <- paste(imageName,"dpi",dpi,".",format,sep="");
  
  if(length(dataSet$rmidx)>0){
    data.norm <- dataSet$data.norm[,-dataSet$rmidx]  
  }else{
    data.norm <- dataSet$data.norm;
  }
  if(anal.type %in% c("onedata", "proteinlist")){
    # Look up symbol directly by gene ID (not positional index, which can desync)
    cmpdNm <- tryCatch({
      id.type <- tolower(paramSet$data.idType)
      if (id.type == "uniprot") {
        doUniprot2SymbolMapping(gene.id, paramSet$data.org, paramSet$data.idType)
      } else {
        doEntrez2SymbolMapping(gene.id, paramSet$data.org, paramSet$data.idType)
      }
    }, error = function(e) gene.id)
    if (dataType == "peptide" || is.null(cmpdNm) || is.na(cmpdNm) || length(cmpdNm) == 0 || cmpdNm == "") {
      cmpdNm <- gene.id
    }
    if(type== "volcano"){
      cmpdNm <- "";
    }
    if (length(dataSet$sec.cls) <= 1) {
      gene.row.id <- gene.id
      if (dataType != "peptide") {
        gene.row.id <- ResolveFeatureRowId(data.norm, gene.id)
        if (!(gene.row.id %in% rownames(data.norm))) {
          msg("[R DEBUG] PlotSelectedGene: feature ID not found in protein matrix: ", gene.id)
          return(0)
        }
      }

      # Collect classes & expression ---------------------------------------
      if (dataType == "peptide") {
        if (!file.exists("peptide_level_data.qs")) {
          return(0)
        }
        pep.mat <- ov_qs_read("peptide_level_data.qs")
        pep.row.id <- ResolveFeatureRowId(pep.mat, gene.id)
        if (!(pep.row.id %in% rownames(pep.mat))) {
          msg("[R DEBUG] PlotSelectedGene: feature ID not found in peptide matrix: ", gene.id)
          return(0)
        }
        if (!is.null(dataSet$comp.type) &&
            length(dataSet$comp.type) > 0 &&
            dataSet$comp.type == "custom") {
          grp.nms <- dataSet$grp.nms
          if (dataSet$cont.inx[dataSet$analysisVar] ||
              any(grepl("(^[0-9]+).*", as.character(dataSet$cls)))) {
            grp.nms <- sub(paste0(dataSet$analysisVar, "_"), "", grp.nms)
          }
          keep <- dataSet$cls %in% grp.nms
          cls  <- droplevels(dataSet$cls[keep])
          expr <- as.numeric(pep.mat[pep.row.id, keep])
        } else {
          expr <- as.numeric(pep.mat[pep.row.id, ])
          meta <- dataSet$meta.info[
                    rownames(dataSet$meta.info) %in% colnames(pep.mat),
                    , drop = FALSE]
          cls  <- droplevels(
                    meta[match(rownames(meta), colnames(pep.mat)),
                         dataSet$analysisVar])
        }
      } else if (!is.null(dataSet$comp.type) &&
                 length(dataSet$comp.type) > 0 &&
                 dataSet$comp.type == "custom") {
        grp.nms <- dataSet$grp.nms
        if (dataSet$cont.inx[dataSet$analysisVar] ||
            any(grepl("(^[0-9]+).*", as.character(dataSet$cls)))) {
          grp.nms <- sub(paste0(dataSet$analysisVar, "_"), "", grp.nms)
        }
        keep <- dataSet$cls %in% grp.nms
        cls  <- droplevels(dataSet$cls[keep])
        expr <- unlist(data.norm[gene.row.id, keep])
      } else {
        expr <- unlist(data.norm[gene.row.id, ])
        meta <- dataSet$meta.info[
                  rownames(dataSet$meta.info) %in% colnames(data.norm),
                  , drop = FALSE]
        cls  <- droplevels(
                  meta[match(rownames(meta), colnames(data.norm)),
                       dataSet$analysisVar])
      }

      # Optional fold-change ------------------------------------------------
      is.disc <- dataSet$disc.inx[dataSet$analysisVar]
      if (is.disc && fc) {
        baseline <- levels(cls)[1]                                     # control level
        baseAvg  <- mean(expr[cls == baseline], na.rm = TRUE)          # log2 mean

        expr    <- expr - baseAvg
        ylabStr <- "Relative expression (vs. baseline mean)"
      } else {
        ylabStr  <- "Expression"
      }

      # Build plot ----------------------------------------------------------
      # Validate data before creating data frame
      if(is.null(expr) || length(expr) == 0){
        print(paste("Error: expr is NULL or empty for gene.id:", gene.id))
        print(paste("expr class:", class(expr), "length:", length(expr)))
        return(0)
      }

      if(is.null(cls) || length(cls) == 0){
        print(paste("Error: cls is NULL or empty for gene.id:", gene.id))
        print(paste("cls class:", class(cls), "length:", length(cls)))
        return(0)
      }

      # Ensure expr is a proper vector (handle matrix row extraction)
      if(is.matrix(expr) || is.data.frame(expr)){
        print(paste("Warning: expr is", class(expr), "converting to vector"))
        expr <- as.vector(expr)
      }

      # Ensure expr is numeric
      expr <- as.numeric(expr)

      # Ensure cls is a vector
      if(is.matrix(cls) || is.data.frame(cls)){
        print(paste("Warning: cls is", class(cls), "converting to vector"))
        cls <- as.vector(cls)
      }

      # Ensure cls is factor
      cls <- as.factor(cls)

      if(length(expr) != length(cls)){
        print(paste("Error: Length mismatch - expr:", length(expr), "cls:", length(cls)))
        return(0)
      }

      # Remove any NA values
      valid.idx <- !is.na(expr) & !is.na(cls)
      if(sum(valid.idx) == 0){
        print(paste("Error: All values are NA for gene.id:", gene.id))
        return(0)
      }

      expr <- expr[valid.idx]
      cls <- cls[valid.idx]

      # Create data frame with explicit column types
      df <- data.frame(value = as.numeric(expr), name = as.character(cls), stringsAsFactors = FALSE)
      col <- GetGroupPalette(df$name, paletteOpt)

      # Validate data frame
      if(!is.data.frame(df) || nrow(df) == 0 || ncol(df) != 2){
        print(paste("Error: Invalid data frame for gene.id:", gene.id))
        print(paste("dim(df):", paste(dim(df), collapse=" x ")))
        print(paste("class(df):", class(df)))
        return(0)
      }

      plot.geom <- if (use.violin) {
        geom_violin(trim = FALSE, aes(color = name), show.legend = FALSE)
      } else {
        geom_boxplot(aes(color = name), outlier.shape = NA, show.legend = FALSE)
      }

      myplot <- ggplot(df, aes(x = name, y = value, fill = name)) +
             plot.geom +
             geom_jitter(width = 0.05, height = 0, color = "black", show.legend = FALSE)       +
             stat_summary(fun = mean, colour = "yellow", geom = "point",
                          shape = 18, size = 3, show.legend = FALSE)          +
             scale_fill_manual(values = col) + scale_color_manual(values = col) +
             theme_bw(base_size = 10) +
             theme(axis.title.x   = element_blank(),
                   legend.position = "none",
                   panel.grid      = element_blank(),
                   plot.title      = element_text(size = 11, hjust = 0.5))    +
             ylab(ylabStr) +
             ggtitle(cmpdNm)

      # Save ----------------------------------------------------------------
      # FIX: Suppress Quartz popup on macOS by using invisible() and checking device
      Cairo(file = imgName, width = 5, height = 5, unit = "in",
            dpi  = dpi, bg   = "white", type = format)
      invisible(print(myplot))
      invisible(dev.off())
      return;
    }else{
      gene.row.id <- gene.id
      if (dataType == "peptide") {
        if (!file.exists("peptide_level_data.qs")) {
          return(0)
        }
        pep.mat <- ov_qs_read("peptide_level_data.qs")
        pep.row.id <- ResolveFeatureRowId(pep.mat, gene.id)
        if (!(pep.row.id %in% rownames(pep.mat))) {
          msg("[R DEBUG] PlotSelectedGene: feature ID not found in peptide matrix: ", gene.id)
          return(0)
        }
        data.norm <- pep.mat
        gene.row.id <- pep.row.id
      } else {
        gene.row.id <- ResolveFeatureRowId(data.norm, gene.id)
        if (!(gene.row.id %in% rownames(data.norm))) {
          msg("[R DEBUG] PlotSelectedGene: feature ID not found in protein matrix: ", gene.id)
          return(0)
        }
      }
      out.fac <- dataSet$sec.cls
      in.fac <- dataSet$fst.cls
      xlab <- colnames(dataSet$meta.info[,1]);
      col <- GetGroupPalette(in.fac, paletteOpt);
      
      img.num <- length(levels(out.fac));
      row.num <- ceiling(img.num/2)
      
      if(row.num == 6){
        layout <- c(row.num, 1);
        h=320;
        w=160*row.num;
      }else{
        rn <- round(row.num/2);
        layout <- c(rn, 2);
        h=500;
        w=160*rn;
      }
      
      if (length(out.fac) == 0 || length(levels(out.fac)) == 0) {
        print("Error: sec.cls is empty; cannot build multi-class plot")
        return(0)
      }

      p_all <- list()
      
      sample_ids <- colnames(data.norm)
      for(lv in levels(out.fac)){
        inx <- out.fac == lv
        expr_vals <- as.numeric(data.norm[gene.row.id, inx])
        cls_vals <- in.fac[inx]
        sample_vals <- sample_ids[inx]
        df.orig <- data.frame(
          facA = lv,
          value = expr_vals,
          sample = sample_vals,
          group = cls_vals,
          stringsAsFactors = FALSE
        )
        p_all[[lv]] <- df.orig
      }
      # FIX: Suppress Quartz popup on macOS
      Cairo(file = imgName, width=5, height=5, type=format, bg="white", dpi=dpi,unit="in");

      alldata <- do.call(rbind, p_all)
      alldata$facA <- factor(as.character(alldata$facA), levels=levels(out.fac))
      plot.geom <- if (use.violin) {
        geom_violin(aes(x = group, color = group), trim = FALSE, show.legend = FALSE)
      } else {
        geom_boxplot(aes(x = group, color = group), outlier.shape = NA, show.legend = FALSE)
      }

      p.time <- ggplot2::ggplot(alldata, aes(x=group, y=value, fill=group, label=sample)) +
        plot.geom +
        geom_jitter(aes(x = group), height = 0, width = 0.05, color = "black", show.legend = FALSE) +
        facet_wrap(~facA, nrow = row.num) +
        theme(axis.title.x = element_blank(), legend.position = "none", axis.text.x = element_text(angle=90, hjust=1),
              plot.title = element_text(size = 11, hjust=0.5, face = "bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        scale_fill_manual(values = col) +
        scale_color_manual(values = col) +
        ggtitle(cmpdNm) +
        ylab("Expression") +
        theme_bw()
      myplot <- p.time + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
      invisible(print(myplot))
      invisible(dev.off())
      return;
    }

  }else{ # metadata
    mdata.all <- paramSet$mdata.all;
    inmex.meta <- ov_qs_read("inmex_meta.qs");
    if(inmex.meta$id.type == "entrez"){
      cmpdNm <- inmex.meta$gene.symbls[gene.id];
    }else{
      cmpdNm <- gene.id;
    }
    num <- length(mdata.all);
    gene.row.id.meta <- ResolveFeatureRowId(inmex.meta$plot.data, gene.id)
    if (!(gene.row.id.meta %in% rownames(inmex.meta$plot.data))) {
      msg("[R DEBUG] PlotSelectedGene: feature ID not found in metadata matrix: ", gene.id)
      return(0)
    }
    # calculate width based on the dateset number
    if(num == 1){
      # FIX: Suppress Quartz popup on macOS
      Cairo(file = imgName, width=5, height=5, type=format, bg="white", dpi=dpi);
      col <- GetGroupPalette(as.character(inmex.meta$cls.lbl), paletteOpt);
      expr_vals <- as.numeric(inmex.meta$plot.data[gene.row.id.meta,])
      cls_vals <- as.character(inmex.meta$cls.lbl)
      df.norm <- data.frame(value=expr_vals, name = cls_vals, stringsAsFactors = FALSE)
      plot.geom <- if (use.violin) {
        geom_violin(trim = FALSE, aes(color = name), show.legend = FALSE)
      } else {
        geom_boxplot(aes(color = name), outlier.shape = NA, show.legend = FALSE)
      }

      p.norm <- ggplot2::ggplot(df.norm, aes(x=name, y=value, fill=name))
      p.norm <- p.norm + plot.geom + geom_jitter(height = 0, width = 0.05, color = "black", show.legend = FALSE)  + theme_bw()
      p.norm <- p.norm + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
      p.norm <- p.norm + stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
      p.norm <- p.norm + scale_fill_manual(values=col) +
        scale_color_manual(values=col) +
        ggtitle(cmpdNm) + theme(axis.text.x = element_text(angle=90, hjust=1))
      p.norm <- p.norm + theme(plot.title = element_text(size = 11, hjust=0.5))
      p.norm <- p.norm + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
      myplot <- p.norm + theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.5, "cm"), axis.text = element_text(size=10))
      invisible(print(myplot))
      invisible(dev.off())
      return;
    }else{
      # calculate layout
      h=420;
      if(num <= 3){
        w=630;
      }else{
        w=num*170
      }
      
      width = w*dpi/72
      height = h*dpi/72
      #print(w)
      #print(h)
      # FIX: Suppress Quartz popup on macOS
      Cairo(file = imgName, width=width, height=height, type=format, bg="white", dpi=dpi);
      data.lbl <- as.character(inmex.meta$data.lbl);
      data.lbl <- substr(data.lbl, 0, nchar(data.lbl)-4);
      
      # get counts in each data, same order as a levels
      counts <- table(data.lbl);
      # back to factor 
      data.lbl <- factor(data.lbl);
      
      # get new lbls to cut potential long names, and add sample numbers
      nlbls <- data.lbl;
      levels(nlbls) <- abbreviate(levels(nlbls),9);
      nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")");
      # update labels
      data.lbl <- factor(data.lbl, labels=nlbls);
      # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
      p_all <- list();
      
      out.fac <- data.lbl;
      in.fac <- inmex.meta$cls.lbl;
      xlab <- colnames(dataSet$meta.info[,1]);
      
      col <- GetGroupPalette(as.character(inmex.meta$cls.lbl), paletteOpt);   
      
      for(lv in levels(out.fac)){
        inx <- out.fac == lv;
        df.orig <- data.frame(facA = lv, value = inmex.meta$plot.data[gene.row.id.meta, inx], name = in.fac[inx])
        p_all[[lv]] <- df.orig
      }
      
      alldata <- do.call(rbind, p_all)
      alldata$Dataset <- factor(as.character(alldata$facA), levels=levels(out.fac))
      colnames(alldata) <- c("Dataset", "value", "Factor")
      
      # Assuming alldata, col, and cmpdNm are already defined
      plot.geom <- if (use.violin) {
        geom_violin(trim = FALSE, aes(color = Factor))
      } else {
        geom_boxplot(aes(color = Factor), outlier.shape = NA)
      }

      p.time <- ggplot(alldata, aes(x = Factor, y = value, fill=Factor)) +
        theme_bw() +
        plot.geom +
        geom_jitter(height = 0, width = 0.05, color = "black", show.legend = FALSE) +
        scale_fill_manual(values = col) +
        theme(axis.title.x = element_blank(),   # Remove x-axis title
              axis.text.x = element_blank(),     # Remove x-axis tick text
              plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        ggtitle(cmpdNm) +
        ylab("Expression")
      
      # Adding facet_wrap to create a separate panel for each Dataset
      p.time <- p.time + facet_wrap(~ Dataset, scales = "free_x")

      myplot <- p.time + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
      invisible(print(myplot))
      invisible(dev.off())
      return;
    }
  }
}

UpdateMultifacPlot <-function(dataName="",imgName, gene.id, boxmeta,format="png", dpi=96, paletteOpt="default", plotType="boxplot"){
  
  require(ggplot2);
  require(see);
  require(lattice);
  
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  dataSet <- readDataset(dataName);
  anal.type <- paramSet$anal.type;
  imgName <- paste(imgName,"dpi",dpi,".",format,sep="");
  plotType <- tolower(as.character(plotType)[1])
  if (!(plotType %in% c("violin", "boxplot"))) {
    plotType <- "boxplot"
  }
  meta <- dataSet$meta.info[dataSet$meta.info[,boxmeta]!="NA",boxmeta,drop=F];
  cls <- droplevels(meta[,boxmeta]);
  data.norm <- dataSet$data.norm[,colnames(dataSet$data.norm) %in% rownames(meta)];
  
  if(anal.type == "onedata"){
    cmpdNm <- tryCatch({
      id.type <- tolower(paramSet$data.idType)
      if (id.type == "uniprot") {
        doUniprot2SymbolMapping(gene.id, paramSet$data.org, paramSet$data.idType)
      } else {
        doEntrez2SymbolMapping(gene.id, paramSet$data.org, paramSet$data.idType)
      }
    }, error = function(e) gene.id)
    if (is.null(cmpdNm) || is.na(cmpdNm) || cmpdNm == "") cmpdNm <- gene.id

    # FIX: Suppress Quartz popup on macOS
    Cairo(file = imgName,  width=320*dpi/72, height=380*dpi/72, type=format, dpi=dpi, bg="white");
    dat <- data.norm
    
    df.norm <- data.frame(value=dat[gene.id,], name = cls);
    if(dataSet$disc.inx[boxmeta]){
      col <- GetGroupPalette(df.norm$name, paletteOpt)
      plot.geom <- if (identical(plotType, "boxplot")) {
        geom_boxplot(aes(color = name), outlier.shape = NA, width = 0.55, show.legend = FALSE)
      } else {
        geom_violin(trim = FALSE, aes(color = name), show.legend = FALSE)
      }
      p.norm <- ggplot2::ggplot(df.norm, aes(x = name, y = value, fill = name)) +
        plot.geom +
        geom_jitter(height = 0, width = 0.05, color = "black", show.legend = FALSE) +
        theme_bw()+
        theme(legend.position = "none") +  xlab(boxmeta) +
        stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE) +
        scale_fill_manual(values = col) + 
        scale_color_manual(values = col) + 
        ggtitle(cmpdNm) + 
        theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(size = 11, hjust=0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    }else{
      df.norm$name <- as.numeric(as.character(df.norm$name ))
      p.norm <- ggplot2::ggplot(df.norm, aes(x=name, y=value))+
        geom_point(size=2) + theme_bw()  + geom_smooth(method=lm,se=T)+
        xlab(boxmeta) +
        theme(axis.text.x = element_text(angle=90, hjust=1)) + guides(size="none")+
        ggtitle(cmpdNm) + theme(plot.title = element_text(size = 11, hjust=0.5, face = "bold")) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    }
    myplot <- p.norm + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
    invisible(print(myplot))
    invisible(dev.off())
    return;

  }else{ # metadata
    mdata.all <- paramSet$mdata.all;
    inmex.meta <- ov_qs_read("inmex_meta.qs");
    if(inmex.meta$id.type == "entrez"){
      cmpdNm <- inmex.meta$gene.symbls[gene.id];
    }else{
      cmpdNm <- gene.id;
    }
    num <- sum(mdata.all == 1);
    # calculate width based on the dateset number
    if(num == 1){
      # FIX: Suppress Quartz popup on macOS
      Cairo(file = imgName, width=280*dpi/72, height=320*dpi/72, type=format, dpi=dpi, bg="white");
      
      col <- GetGroupPalette(as.character(inmex.meta$cls.lbl), paletteOpt);
      df.norm <- data.frame(value=inmex.meta$plot.data[gene.id,], name = as.character(inmex.meta$cls.lbl))
      plot.geom <- if (identical(plotType, "boxplot")) {
        geom_boxplot(aes(color = name), outlier.shape = NA, show.legend = FALSE)
      } else {
        geom_violin(trim = FALSE, aes(color = name), show.legend = FALSE)
      }
      p.norm <- ggplot2::ggplot(df.norm, aes(x=name, y=value, fill=name))  
      p.norm <- p.norm + plot.geom + geom_jitter(height = 0, width = 0.05, color = "black", show.legend = FALSE)  + theme_bw()
      p.norm <- p.norm + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
      p.norm <- p.norm + stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
      p.norm <- p.norm + scale_fill_manual(values=col) + 
        scale_color_manual(values=col) +
        ggtitle(cmpdNm) + theme(axis.text.x = element_text(angle=90, hjust=1))
      p.norm <- p.norm + theme(plot.title = element_text(size = 11, hjust=0.5))
      p.norm <- p.norm + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
      myplot <- p.norm + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
      invisible(print(myplot))
      invisible(dev.off())
      return;
    }
  }
}

PlotSelectedGeneRaw <- function(gene.id="",imgName="", format="png", dpi=96) {

  library(ggplot2)
  library(Cairo)

  dataSet <- readDataset("dataset_rawmodule")
  dat <- dataSet$data.raw
  cls <- dataSet$meta.info[,1]  # Extract metadata info
  
  col <- unique(GetColorSchema(cls))   
  
  df.raw <- data.frame(value = unlist(dat[which(rownames(dat) == gene.id), ]), name = cls)
  imgName <- paste(imgName,"dpi",dpi,".",format,sep="");

  p <- ggplot2::ggplot(df.raw, aes(x = name, y = value, fill = name)) +
          geom_violin(trim = FALSE, aes(color = name), show.legend = FALSE) + 
          geom_jitter(height = 0, width = 0.05, color = "black", show.legend = FALSE) +
          theme(legend.position = "none") +  
          xlab("Metadata") +
          stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE) +
          ggtitle(gene.id) + 
          theme(axis.title.x = element_blank(), plot.title = element_text(size = 11, hjust = 0.5), 
                panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
          theme_bw()

  # Generate and save the plot using Cairo
  # FIX: Suppress Quartz popup on macOS
  Cairo(file = imgName, width = 420 * dpi / 72, height = 560 * dpi / 72, type = format, dpi = dpi, bg = "white")
  invisible(print(p))
  invisible(dev.off())
}
