

### Example for HRiD / HRiP estimation ###


## ------------------------------------------------------------------------
## 1. Import phased VCF data
## ------------------------------------------------------------------------

## This block concatenates autosomes (1–21) and X chromosome into a single data.frame.
## Adjust the chromosome range if your dataset has a different naming scheme.
## In this example, chromosome X is stored in vcf as 22, as there are 21 autosomes.

vcf <- do.call(rbind, lapply(1:22, function(i) {
  if (i == 22) {
    file_name <- "FinalReport_afterQC_chrX.vcf"
  } else {
    file_name <- paste0("FinalReport_afterQC_autosomes_chr", i, ".vcf")
  }
  file_data <- read.table(file_name,
                          skip = 5, 
                          header = TRUE, 
                          comment.char = "")   # Prevent '#' from being ignored
  # Clean up column names
  colnames(file_data) <- sub("^#", "", colnames(file_data))  # Remove leading '#'
  if ("X.CHROM" %in% colnames(file_data)) {
    colnames(file_data)[colnames(file_data) == "X.CHROM"] <- "CHROM"
  }
  file_data
}))

## ensure CHROM is numeric (map X->22 beforehand if needed)


## ------------------------------------------------------------------------
## 2. Load HaplotypeRichness.Estimation() function
## ------------------------------------------------------------------------
## Make sure the path points to the function file inside the repo
source("R/HaplotypeRichness_Estimation.R")


## ------------------------------------------------------------------------
## 3. Run tests with all different modes (4)
## ------------------------------------------------------------------------

## Each call below demonstrates a different combination of approach and windowing mode.
## Results are stored in objects Test1–Test4.

# Sliding windows, bp-based windows in estimation and slide
Test1 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = FALSE,
                                      approach = "Bp.based")

# Sliding windows, SNP-based windows in estimation and slide
Test2 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = FALSE,
                                      approach = "SNP.based")

# Per-SNP (consecutive), bp-based windows in estimation
Test3 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = TRUE,
                                      approach = "Bp.based")

# Per-SNP (consecutive), SNP-based windows in estimation
Test4 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = TRUE,
                                      approach = "SNP.based")


## ------------------------------------------------------------------------
## 4. Add Mb positions for visualization
## ------------------------------------------------------------------------

## Sliding modes: use window centers; consecutive modes: use anchor SNP position

if (all(c("start.position","end.position") %in% names(Test1))) {
  Test1$Center_Mb <- ((Test1$start.position + Test1$end.position) / 2) / 1e6
}
if (all(c("start.position","end.position") %in% names(Test2))) {
  Test2$Center_Mb <- ((Test2$start.position + Test2$end.position) / 2) / 1e6
}
if ("position" %in% names(Test3)) {
  Test3$Mb <- Test3$position / 1e6
}
if ("position" %in% names(Test4)) {
  Test4$Mb <- Test4$position / 1e6
}


## ------------------------------------------------------------------------
## 5. HRiD outlier detection and signatures (consecutiveSNP = TRUE)
## ------------------------------------------------------------------------

## Goal:
##  - identify "solo outliers": significant anchors that are NOT in a consecutive run
##  - identify "consecutive outliers - true signatures": runs of >= min_run_len significant anchors
##
## Significance is defined on HRiD_LogPvalue using a user-defined threshold 'thr'.
## NOTE: Choose 'thr' via simpleM or FDR in your study design.

thr <- 5            # example threshold on HRiD_LogPvalue
min_run_len <- 3    # minimal length for consecutive run

## Helper to annotate significant runs per chromosome
## Adds only is_sig, run_id, is_consecutive_sig
annotate_runs <- function(df, chr_col = "CHROM",
                          stat_col = "HRiD_LogPvalue",
                          thr = 5, min_run_len = 3) {
  stopifnot(all(c(chr_col, stat_col) %in% names(df)))
  # sort within chromosome
  if ("position" %in% names(df)) {
    df <- df[order(df[[chr_col]], df[["position"]]), ]
  } else if ("start.Index" %in% names(df)) {
    df <- df[order(df[[chr_col]], df[["start.Index"]]), ]
  } else {
    df <- df[order(df[[chr_col]]), ]
  }
  df$is_sig <- !is.na(df[[stat_col]]) & (df[[stat_col]] >= thr)
  df$run_id <- NA_integer_
  run_counter <- 0L
  for (chr in unique(df[[chr_col]])) {
    idx <- which(df[[chr_col]] == chr)
    if (!any(df$is_sig[idx], na.rm = TRUE)) next
    r <- rle(df$is_sig[idx])
    ends <- cumsum(r$lengths)
    starts <- c(1, head(ends, -1) + 1)
    for (k in seq_along(r$lengths)) {
      if (r$values[k]) {
        seg <- idx[starts[k]:ends[k]]
        if (length(seg) >= min_run_len) {
          run_counter <- run_counter + 1L
          df$run_id[seg] <- run_counter
        }
      }
    }
  }
  df$is_consecutive_sig <- !is.na(df$run_id)
  return(df)
}

## Apply to Test3 and Test4 (per-SNP modes)
if (all(c("CHROM","HRiD_LogPvalue") %in% names(Test3))) {
  Test3_annot <- annotate_runs(Test3, chr_col = "CHROM",
                               stat_col = "HRiD_LogPvalue",
                               thr = thr, min_run_len = min_run_len)
  
  # requested outputs
  Test3_outliers        <- subset(Test3_annot, is_sig)                      # all significant anchors
  Test3_signatures      <- subset(Test3_annot, is_consecutive_sig)          # consecutive runs
  Test3_con_outliers    <- Test3_signatures$SNPname                         # vector of SNP names 
  Test3_solo_signatures <- subset(Test3_annot, is_sig & !is_consecutive_sig)# solo outliers (significant but not in a run)
  Test3_solo_outliers   <- Test3_solo_signatures$SNPname                    # vector of SNP names
  
  # ordering for readability
  if ("position" %in% names(Test3_outliers)) {
    Test3_outliers        <- Test3_outliers[order(Test3_outliers$CHROM, Test3_outliers$position), ]
    Test3_signatures      <- Test3_signatures[order(Test3_signatures$CHROM, Test3_signatures$run_id, Test3_signatures$position), ]
    Test3_solo_signatures <- Test3_solo_signatures[order(Test3_solo_signatures$CHROM, Test3_solo_signatures$position), ]
  }
}

if (all(c("CHROM","HRiD_LogPvalue") %in% names(Test4))) {
  Test4_annot <- annotate_runs(Test4, chr_col = "CHROM",
                               stat_col = "HRiD_LogPvalue",
                               thr = thr, min_run_len = min_run_len)
  
  Test4_outliers        <- subset(Test4_annot, is_sig)                     
  Test4_signatures      <- subset(Test4_annot, is_consecutive_sig)
  Test4_con_outliers    <- Test4_signatures$SNPname
  Test4_solo_signatures <- subset(Test4_annot, is_sig & !is_consecutive_sig)
  Test4_solo_outliers   <- Test4_solo_signatures$SNPname
  
  if ("position" %in% names(Test4_outliers)) {
    Test4_outliers        <- Test4_outliers[order(Test4_outliers$CHROM, Test4_outliers$position), ]
    Test4_signatures      <- Test4_signatures[order(Test4_signatures$CHROM, Test4_signatures$run_id, Test4_signatures$position), ]
    Test4_solo_signatures <- Test4_solo_signatures[order(Test4_solo_signatures$CHROM, Test4_solo_signatures$position), ]
  }
}


## ------------------------------------------------------------------------
## 6) HRiD outlier detection and signatures (consecutiveSNP = FALSE)
## ------------------------------------------------------------------------

## Significance is defined on HRiD_LogPvalue using a user-defined threshold 'thr'.
## NOTE: Choose 'thr' via simpleM or FDR in your study design.
## All significant windows are considered as outliers/signatures

Test1_signatures <- Test1[Test1$HRiD_LogPvalue >= thr & !is.na(Test1$HRiD_LogPvalue), ]
Test1_outliers   <- Test1_signatures$Window

Test2_signatures <- Test2[Test2$HRiD_LogPvalue >= thr & !is.na(Test2$HRiD_LogPvalue), ]
Test2_outliers   <- Test2_signatures$Window


## ------------------------------------------------------------------------
## 7. Quick inspection
## ------------------------------------------------------------------------

head(Test1); head(Test2); head(Test3); head(Test4)

## Inspect counts of detected signatures (consecutiveSNP=FALSE)
cat("\nTest1 (Bp.based, sliding):\n")
cat("  Outlier windows (>=",thr,"): ", nrow(Test1_signatures), "\n", sep = '')
cat("\nTest2 (SNP.based, sliding):\n")
cat("  Outlier windows (>=",thr,"): ", nrow(Test2_signatures), "\n", sep = '')

## Inspect counts of detected signals (consecutiveSNP=TRUE)
if (exists("Test3_outliers")) {
  cat("\nTest3 (Bp.based, consecutiveSNP=TRUE):\n", sep = '')
  cat("  All outliers (>=",thr,"): ", nrow(Test3_outliers), "\n", sep = '')
  cat("  Solo outliers:       ", nrow(Test3_solo_signatures), "\n", sep = '')
  cat("  Consecutive signatures:", length(unique(Test3_signatures$run_id)), "runs\n")
}
if (exists("Test4_outliers")) {
  cat("\nTest4 (SNP.based, consecutiveSNP=TRUE):\n", sep = '')
  cat("  All outliers (>=",thr,"): ", nrow(Test4_outliers), "\n", sep = '')
  cat("  Solo outliers:       ", nrow(Test4_solo_signatures), "\n", sep = '')
  cat("  Consecutive signatures:", length(unique(Test4_signatures$run_id)), "runs\n")
}


## ------------------------------------------------------------------------
## 8) Visualization helpers
## ------------------------------------------------------------------------

# Light wrapper adapted from qqman's manhattan()
manhatan1 <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                       col = c("ivory4", "ivory3"),
                       chrlabs = NULL,
                       suggestiveline = -log10(1e-05),
                       genomewideline = -log10(5e-08),
                       highlight1 = NULL, highlight2 = NULL,
                       logp = TRUE, annotatePval = NULL,
                       annotateTop = TRUE, ...) {
  
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))  stop(paste("Column", bp,  "not found!"))
  if (!(p %in% names(x)))   stop(paste("Column", p,   "not found!"))
  if (!(snp %in% names(x))) warning("No SNP column found. OK unless you're trying to highlight.")
  
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric."))
  if (!is.numeric(x[[bp]]))  stop(paste(bp,  "column should be numeric."))
  if (!is.numeric(x[[p]]))   stop(paste(p,   "column should be numeric."))
  
  if (!is.null(x[[snp]])) {
    d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]],
                    pos = NA, index = NA, SNP = x[[snp]],
                    stringsAsFactors = FALSE)
  } else {
    d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]],
                    pos = NA, index = NA)
  }
  
  d <- d[order(d$CHR, d$BP), ]
  
  d$logp <- if (logp) -log10(d$P) else d$P
  
  d$index <- rep.int(seq_along(unique(d$CHR)),
                     times = tapply(if (!is.null(d$SNP)) d$SNP else d$BP, d$CHR, length))
  nchr <- length(unique(d$CHR))
  
  if (nchr == 1) {
    d$pos <- d$BP
    xlabel <- paste("Chromosome", unique(d$CHR), "position")
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos <- d[d$index == i, ]$BP
      } else {
        lastbase <- lastbase + max(d[d$index == (i - 1), "BP"])
        d[d$index == i, "BP"]  <- d[d$index == i, "BP"] - min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] <- d[d$index == i, "BP"] + lastbase
      }
    }
    ticks  <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel <- "Chromosome"
    labs   <- unique(d$CHR)
  }
  
  # robust ylim (finite only)
  ylim_max <- suppressWarnings(max(d$logp[is.finite(d$logp)], na.rm = TRUE))
  if (!is.finite(ylim_max)) ylim_max <- 1
  
  xmax <- ceiling(max(d$pos) * 1.03)
  xmin <- floor(max(d$pos) * -0.03)
  
  def_args <- list(
    xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", las = 1, pch = 20,
    xlim = c(xmin, xmax), ylim = c(0, ceiling(ylim_max)),
    xlab = xlabel, ylab = if (logp) expression(-log[10](italic(P))) else "Statistic"
  )
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(unique(d$CHR))) {
        labs <- chrlabs
      } else {
        warning("Number of chromosome labels != number of chromosomes.")
      }
    } else {
      warning("chrlabs must be a character vector.")
    }
  }
  
  if (nchr == 1) axis(1, ...) else axis(1, at = ticks, labels = labs, ...)
  
  col <- rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  } else {
    icol <- 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], col = col[icol], pch = 20, ...)
      icol <- icol + 1
    }
  }
  
  if (isTRUE(suggestiveline))  abline(h = suggestiveline, col = "blue")
  if (isTRUE(genomewideline) || is.numeric(genomewideline)) abline(h = genomewideline, col = "magenta", lwd = 2)
  
  if (!is.null(highlight1)) {
    d.h1 <- d[!is.null(d$SNP) & d$SNP %in% highlight1, , drop = FALSE]
    if (nrow(d.h1) > 0) with(d.h1, points(pos, logp, col = "lightsalmon", pch = 20, ...))
  }
  if (!is.null(highlight2)) {
    d.h2 <- d[!is.null(d$SNP) & d$SNP %in% highlight2, , drop = FALSE]
    if (nrow(d.h2) > 0) with(d.h2, points(pos, logp, col = "dodgerblue", pch = 20, ...))
  }
}


## ------------------------------------------------------------------------
## 9) Manhattan plots
## ------------------------------------------------------------------------

new.names <- c(1:21,'X')
par(mfrow = c(2, 2)) 

# Sliding, bp-based (Test1)
manhatan1(Test1, bp= "Center_Mb", chr= "CHROM", p= "HRiD_Pvalue", snp="Window",highlight2 = Test1_outliers, logp=TRUE,ylim=c(0,max(Test1$HRiD_LogPvalue, na.rm = TRUE)+2), xlab = 'Chromosome',suggestiveline = FALSE, genomewideline = 5,chrlabs = new.names,cex.axis = 0.8, las = 2)

# Sliding, SNP-based (Test2)
manhatan1(Test2, bp= "Center_Mb", chr= "CHROM", p= "HRiD_Pvalue", snp="Window",highlight2 = Test2_outliers, logp=TRUE,ylim=c(0,max(Test2$HRiD_LogPvalue, na.rm = TRUE)+2), xlab = 'Chromosome',suggestiveline = FALSE, genomewideline = 5,chrlabs = new.names,cex.axis = 0.8, las = 2)

# Consecutive, bp-based (Test3)
manhatan1 (Test3, bp= "Mb", chr= "CHROM", p= "HRiD_Pvalue", snp="SNPname",highlight1 = Test3_solo_outliers,highlight2 = Test3_con_outliers, logp=TRUE,ylim=c(0,max(Test3$HRiD_LogPvalue, na.rm = TRUE)+2), xlab = 'Chromosome',suggestiveline = FALSE, genomewideline = 5,chrlabs = new.names,cex.axis = 0.8, las = 2)

# Consecutive, SNP-based (Test4)
manhatan1 (Test4, bp= "Mb", chr= "CHROM", p= "HRiD_Pvalue", snp="SNPname",highlight1 = Test4_solo_outliers,highlight2 = Test4_con_outliers, logp=TRUE,ylim=c(0,max(Test4$HRiD_LogPvalue, na.rm = TRUE)+2), xlab = 'Chromosome',suggestiveline = FALSE, genomewideline = 5,chrlabs = new.names,cex.axis = 0.8, las = 2)


## ------------------------------------------------------------------------
## 10) (Optional) Environment cleanup — keep only key objects
## ------------------------------------------------------------------------

rm(list = setdiff(ls(), c(
  "vcf", "HaplotypeRichness.Estimation","manhatan1", 
  "thr", "min_run_len", "new.names",
  "Test1","Test2","Test3","Test4",
  "Test1_signatures","Test2_signatures","Test1_outliers","Test2_outliers",
  "Test3_outliers","Test3_signatures","Test3_con_outliers","Test3_solo_signatures","Test3_solo_outliers",
  "Test4_outliers","Test4_signatures","Test4_con_outliers","Test4_solo_signatures","Test4_solo_outliers"
)))


###########################################################################



