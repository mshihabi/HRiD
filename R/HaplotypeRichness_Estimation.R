HaplotypeRichness.Estimation <- function(
    # --- Input data ---
  vcf               = NULL,                       # data.frame with standard VCF columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) + sample genotype columns
  vcf.file.names    = NULL,                       # character vector of VCF file paths to read and row-bind; use either 'vcf' OR 'vcf.file.names'
  
  # --- Estimation mode & size ---
  approach          = c("Bp.based","SNP.based"),  # estimation mode: windows used in estimation are defined by base-pair length ("Bp.based") or by SNP count ("SNP.based")
  consecutiveSNP    = FALSE,                      # TRUE: compute estimation centered at each SNP; FALSE: compute estimation for each window using slide windows by 'slide'
  start             = 0,                          # Bp.based: starting bp coordinate for the first window, usually leave at 0
                                                  # SNP.based: starting SNP index (offset), usually leave at 0
  end               = NULL,                       # If omitted, defaults depend on 'approach':
                                                  #   Bp.based  -> 1,000,000 bp (1 Mb)
                                                  #   SNP.based -> 30 SNPs
                                                  # If both 'start' and 'end' are set, window length is (end - start)
  slide             = end/2,                      # Step when consecutiveSNP = FALSE; units follow 'approach' (bp in Bp.based, SNPs in SNP.based)
                                                  # If omitted, it defaults to half of the chosen 'end' (i.e., end/2) after 'end' is set
                                                  # How to set start/end/slide:
                                                  #   Bp.based  : set 'end' to desired bp window (e.g., 1 Mb = 1000000), keep 'start' at 0 (starts at first position), choose 'slide' in bp
                                                  #   SNP.based : set 'end' to desired SNP count (e.g., 30), keep 'start' at 0 (starts at first SNP), choose 'slide' in SNPs
                                                  #   When consecutiveSNP = TRUE, 'slide' is ignored
  
  # --- Parameters of quality / constraints ---
  minSNP            = 1,                          # If omitted, defaults depend on 'approach':
                                                  #   Bp.based  -> 15 (absolute SNP count per bp window)
                                                  #   SNP.based -> 2/3 (fraction of the SNP window that must remain after splitting by 'maxGap')
                                                  # Bp.based  : to disable the filter, set 1
                                                  # SNP.based : interpreted as FRACTION (0–1); e.g., 0.05 with end=200 => need ≥10 SNPs in the kept block
                                                  #             1.0 requires the full window (no large gaps after splitting); >1 makes all windows fail
                                                  #             to disable the filter here, set 0 (not recommended)
  maxGap            = 1000000,                    # SNP.based only: maximum allowed distance (bp) between adjacent SNPs inside a SNP window
                                                  #                 If any gap > maxGap, split into contiguous blocks and keep a single block
                                                  #                 (the block containing the center SNP when consecutiveSNP = TRUE; otherwise the longest block)
                                                  #                 HR is computed on that kept block if 'minSNP' (fractional) is satisfied
  
  # --- Haploid support ---
  HaploidExistance  = FALSE,                      # set TRUE if some chromosomes are haploid in one or both sexes (e.g., X/Y systems)
  HemizygosityProportion = 0.95,                  # share of identical alleles (h1 == h2) on a haploid chromosome above which a sample is flagged as hemizygous
  
  # --- Normalization & p-values ---
  WithNormalisation = TRUE,                       # also compute Z-scores and p-values (HRiD for decrease; HRiP for increase)
  NormalisationByChr= FALSE                       # TRUE: normalize within each chromosome; FALSE: normalize across all chromosomes together
) {
  # Validate/select windowing approach
  approach <- match.arg(approach)
  
  # --- Mode-specific defaults (only if the user did not supply these args) ---
  if (missing(end) || is.null(end)) {
    # 1 Mb for Bp.based; 30 SNPs for SNP.based
    end <- if (approach == "Bp.based") 1000000 else 30
  }
  if (missing(minSNP) || is.null(minSNP)) {
    # 15 for Bp.based (absolute); 2/3 for SNP.based (fraction)
    minSNP <- if (approach == "Bp.based") 15 else 2/3
  }
  if (missing(slide) || is.null(slide)) {
    # keep intended behavior: slide = end/2 once 'end' is set
    slide <- end/2
  }
  
  # Helper — minimal SNP count required in a window (minSNP interpreted as a fraction if < 1, used only in SNP.based mode)
  needed_snp_count <- function(minSNP, win_len) {
    if (win_len <= 0) return(Inf)
    ceiling(minSNP * win_len)
  }
  
  # Helper — enforce maxGap, and if violated salvage the window by keeping a single contiguous block, with optional anchoring
  # pos_vec: vector of genomic positions for the whole chromosome table
  # idx: integer indices (rows) that form the tentative window
  # anchor_idx: if provided, choose the block that contains this global row index, otherwise choose the longest block
  salvage_by_gap <- function(pos_vec, idx, maxGap, min_needed, anchor_idx = NULL) {
    if (length(idx) == 0) return(list(valid=FALSE, use_idx=integer(0)))
    if (length(idx) == 1) {
      return(list(valid = (1 >= min_needed), use_idx = idx))
    }
    pos <- pos_vec[idx]
    gaps <- diff(pos)
    # All gaps within the limit — use the full window if it meets min_needed
    if (all(gaps <= maxGap, na.rm = TRUE)) {
      return(list(valid = (length(idx) >= min_needed), use_idx = idx))
    }
    # Split at large gaps
    brk <- c(0, which(gaps > maxGap), length(pos))
    rel_runs <- lapply(seq_len(length(brk) - 1), function(i) (brk[i] + 1):brk[i+1])
    runs <- lapply(rel_runs, function(rr) idx[rr])
    # Choose the block to keep
    choose_run <- function() {
      if (!is.null(anchor_idx)) {
        hit <- vapply(runs, function(r) (anchor_idx >= r[1] && anchor_idx <= r[length(r)]), logical(1))
        if (any(hit)) return(runs[[which(hit)[1]]])
      }
      # Fallback to the longest block
      lens <- vapply(runs, length, integer(1))
      runs[[which.max(lens)]]
    }
    chosen <- choose_run()
    list(valid = (length(chosen) >= min_needed), use_idx = chosen)
  }
  
  # Input check — at least one of vcf or vcf.file.names must be provided
  if (is.null(vcf) & is.null(vcf.file.names)) {
    stop("Either a VCF data.frame (vcf) or a vector of VCF file names (vcf.file.names) must be provided.\nPlease ensure that the file names in 'vcf.file.names' are listed in the correct order — e.g., chromosome X should appear after autosomes.")
  }
  
  # Read and row-bind VCF files when provided via vcf.file.names
  if (!is.null(vcf.file.names)) {
    vcf <- tryCatch({
      do.call(rbind, lapply(vcf.file.names, function(file_path) {
        file_data <- tryCatch({
          all_lines <- readLines(file_path)
          header_line_index <- grep("^#CHROM", all_lines)
          if (length(header_line_index) == 0) {
            stop(paste0("Header line starting with '#CHROM' not found in file: ", file_path))
          }
          file_data <- read.table(
            file_path,
            skip = header_line_index - 1,
            header = TRUE,
            sep = "\t",
            quote = "",
            comment.char = "#",
            check.names = FALSE
          )
          colnames(file_data) <- sub("^#", "", colnames(file_data))
          if ("X.CHROM" %in% colnames(file_data)) {
            colnames(file_data)[colnames(file_data) == "X.CHROM"] <- "CHROM"
          }
          return(file_data)
        }, error = function(e) {
          stop(paste0("Error while processing file: ", file_path, "\n", e$message))
        })
      }))
    }, error = function(e) {
      stop("Error while reading VCF file(s): format is not suitable.\nPlease ensure your VCF files follow the expected structure.\nSee example input file format (e.g., ExampleImportVCF.txt).")
    })
  }
  
  # Validate VCF data.frame when provided directly via vcf
  if (!is.null(vcf)) {
    required_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    actual_columns <- colnames(vcf)[1:9]
    if (!all(required_columns == actual_columns)) {
      stop("Error: VCF object is not in the expected format.\nPlease ensure the VCF contains standard columns in order:\nCHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT.\nSee example input file format (e.g., ExampleImportVCF.txt).")
    }
  }
  
  vcf_org <- vcf
  
  # Optional — detect haploid chromosomes and process them separately if present
  if (HaploidExistance == TRUE) {
    haploid.chromosomes <- NULL
    total_chromosomes <- length(unique(vcf$CHROM))
    message("Haploid chromosome detection is enabled.")
    message("Total loaded chromosomes: ", total_chromosomes)
    message("Please specify which chromosome(s) are haploid in at least one sex.")
    message("Enter chromosome numbers only (e.g. 23), if more of them are present, separate them with commas (e.g., 23, 24):")
    user_input <- readline("Haploid chromosome(s): ")
    haploid.chromosomes <- as.numeric(unlist(strsplit(user_input, ",")))
    if (any(is.na(haploid.chromosomes))) {
      stop("Invalid input: only numeric chromosome identifiers are allowed.")
    }
    message("Haploid chromosome(s) set to: ", paste(haploid.chromosomes, collapse = ", "))
    
    vcf_haplo <- vcf[vcf$CHROM %in% haploid.chromosomes, ]
    vcf       <- vcf[!(vcf$CHROM %in% haploid.chromosomes), ]
    
    message(Sys.time(), " | Detection of hemizygous IDs...")
    
    # Store per-chromosome phased haplotypes for haploids
    haplo_fazes <- setNames(vector("list", length(haploid.chromosomes)),
                            as.character(haploid.chromosomes))
    
    # Track hemizygous individuals and their chromosomes for reporting
    hemizygous.IDs <- character(0)
    hemizygous.Chr <- numeric(0)
    
    if (nrow(vcf_haplo) > 0) {
      for (f in 10:ncol(vcf_haplo)) {
        sample_genotypes <- vcf_haplo[, c("CHROM", "POS", colnames(vcf_haplo)[f]), drop = FALSE]
        first_non_na <- sample_genotypes[which(!is.na(sample_genotypes[, 3]))[1], 3]
        if (grepl("/", first_non_na, fixed = TRUE)) {
          separator <- "/"
        } else if (grepl("|", first_non_na, fixed = TRUE)) {
          separator <- "|"
        } else {
          stop(paste("Unknown genotype separator in sample:", colnames(vcf_haplo)[f]))
        }
        alleles <- do.call("rbind",
                           strsplit(as.character(sample_genotypes[, 3]), separator, fixed = TRUE))
        colnames(alleles) <- c("h1", "h2")
        sample_name <- colnames(vcf_haplo)[f]
        
        sample_tbl <- cbind(sample_genotypes[, c("CHROM", "POS")], alleles, stringsAsFactors = FALSE)
        for (chr in haploid.chromosomes) {
          h_sample <- sample_tbl[sample_tbl$CHROM == chr, , drop = FALSE]
          prop_identical <- mean(h_sample$h1 == h_sample$h2, na.rm = TRUE)
          key <- as.character(chr)
          if (is.null(haplo_fazes[[key]])) {
            haplo_fazes[[key]] <- data.frame(POS = h_sample$POS, stringsAsFactors = FALSE)
          } else {
            if (!identical(haplo_fazes[[key]]$POS, h_sample$POS)) {
              stop(sprintf("POS mismatch detected for chromosome %s across samples;", key),
                   " ensure VCF rows are consistent/sorted before phasing.")
            }
          }
          if (prop_identical > HemizygosityProportion) {
            haplo_fazes[[key]][[paste0(sample_name, "_h1")]] <- h_sample$h1
            hemizygous.IDs <- c(hemizygous.IDs, sample_name)
            hemizygous.Chr <- c(hemizygous.Chr, chr)
          } else {
            haplo_fazes[[key]][[paste0(sample_name, "_h1")]] <- h_sample$h1
            haplo_fazes[[key]][[paste0(sample_name, "_h2")]] <- h_sample$h2
          }
        }
      }
    }
    
    if (length(hemizygous.IDs) > 0) {
      for (ind in seq_along(hemizygous.IDs)) {
        message("Hemizygosity (>", HemizygosityProportion, ") detected in individual ",
                hemizygous.IDs[ind], " on chromosome: ", hemizygous.Chr[ind])
      }
    } else {
      message("No hemizygosity detected above threshold (>", HemizygosityProportion, ").")
    }
    
    assign("haplo_fazes", haplo_fazes, inherits = TRUE)
    message("If you have a Pseudo-Autosomal Region (PAR), please include it as a separate chromosome or exclude it before analysis.")
    
    message(Sys.time(), " | HR computation on haploid chromosomes...")
    
    HR_haplo <- NULL
    for (chr in haploid.chromosomes) {
      key <- as.character(chr)
      H <- haplo_fazes[[key]]
      H$RowIndex <- seq_len(nrow(H))
      CHROM_h <- start.Index_h <- end.Index_h <- start.position_h <- end.position_h <- na_h <- nhh <- HR_value_h <- NULL
      anchor_pos_h_chr <- NULL
      
      if (!isTRUE(consecutiveSNP)) {
        # Bp.based sliding by bp distance, non-consecutive, haploid
        if (approach == "Bp.based") {
          startt <- start
          endd   <- end
          if (endd <= H$POS[1]) {
            startt <- H$POS[1]
            endd   <- startt + end
          }
          while (endd <= H$POS[nrow(H)]) {
            if (endd != H$POS[nrow(H)] & (endd + slide) > H$POS[nrow(H)]) {
              endd <- H$POS[nrow(H)]
            }
            idx <- which(H$POS >= startt & H$POS <= endd)
            if (length(idx) >= minSNP) {
              CHROM_h <- c(CHROM_h, chr)
              subH <- H[idx, , drop = FALSE]
              index_pos <- c(subH$RowIndex[1], subH$POS[1], subH$RowIndex[nrow(subH)], subH$POS[nrow(subH)])
              hap_cols <- setdiff(colnames(subH), c("POS", "RowIndex"))
              k <- as.matrix(subH[, hap_cols, drop = FALSE])
              if (!is.character(k)) storage.mode(k) <- "character"
              keys <- apply(k, 2, paste0, collapse = "\r")
              grp  <- match(keys, unique(keys))
              cnt  <- tabulate(grp, nbins = max(grp))
              freq <- cnt / length(keys)
              eff  <- 1 / sum(freq^2)
              na_current <- length(cnt)
              start.position_h <- c(start.position_h, index_pos[2])
              end.position_h   <- c(end.position_h,   index_pos[4])
              start.Index_h    <- c(start.Index_h,    index_pos[1])
              end.Index_h      <- c(end.Index_h,      index_pos[3])
              nhh  <- c(nhh, eff)
              na_h <- c(na_h, na_current)
            } else if (length(idx) > 0) {
              subH <- H[idx, , drop = FALSE]
              index_pos <- c(subH$RowIndex[1], subH$POS[1], subH$RowIndex[nrow(subH)], subH$POS[nrow(subH)])
              CHROM_h <- c(CHROM_h, chr)
              start.position_h <- c(start.position_h, index_pos[2])
              end.position_h   <- c(end.position_h,   index_pos[4])
              start.Index_h    <- c(start.Index_h,    index_pos[1])
              end.Index_h      <- c(end.Index_h,      index_pos[3])
              na_h <- c(na_h, NA)
              nhh  <- c(nhh, NA)
            }
            startt <- startt + slide
            endd   <- endd + slide
          }
          m <- length(nhh)
          HR_value_h <- rep(NA_real_, m)
          if (m >= 2) {
            HR_value_h[1] <- nhh[2] / nhh[1]
            HR_value_h[m] <- nhh[m-1] / nhh[m]
          }
          if (m > 2) {
            mid <- 2:(m-1)
            HR_value_h[mid] <- (nhh[mid-1] + nhh[mid+1]) / (2 * nhh[mid])
          }
        } else {
          # SNP.based, non-consecutive, haploid — index windows with maxGap salvage
          win_len <- if ((end - start) > 0) (end - start) else end
          if (win_len <= 0) stop("In SNP.based mode (haploid, non-consecutive) window length (end-start) must be > 0.")
          stepS <- max(1, as.integer(slide))
          min_needed <- needed_snp_count(minSNP, win_len)
          
          for (iL in seq.int(1 + start, nrow(H), by = stepS)) {
            iR <- min(iL + win_len - 1, nrow(H))
            idx <- iL:iR
            
            sv <- salvage_by_gap(H$POS, idx, maxGap, min_needed, anchor_idx = NULL)
            CHROM_h <- c(CHROM_h, chr)
            if (sv$valid) {
              subH <- H[sv$use_idx, , drop = FALSE]
              index_pos <- c(subH$RowIndex[1], subH$POS[1], subH$RowIndex[nrow(subH)], subH$POS[nrow(subH)])
              hap_cols <- setdiff(colnames(subH), c("POS","RowIndex"))
              M <- as.matrix(subH[, hap_cols, drop = FALSE]); if (!is.character(M)) storage.mode(M) <- "character"
              keys <- apply(M, 2, paste0, collapse = "\r"); grp <- match(keys, unique(keys))
              cnt <- tabulate(grp, nbins = max(grp)); freq <- cnt / length(keys); eff <- 1 / sum(freq^2)
              start.Index_h <- c(start.Index_h, index_pos[1]); end.Index_h <- c(end.Index_h, index_pos[3])
              start.position_h <- c(start.position_h, index_pos[2]); end.position_h <- c(end.position_h, index_pos[4])
              na_h <- c(na_h, length(cnt)); nhh <- c(nhh, eff)
            } else {
              # Keep attempted bounds even if the window was not salvageable
              start.Index_h <- c(start.Index_h, H$RowIndex[idx[1]])
              end.Index_h   <- c(end.Index_h,   H$RowIndex[idx[length(idx)]])
              start.position_h <- c(start.position_h, H$POS[idx[1]])
              end.position_h   <- c(end.position_h,   H$POS[idx[length(idx)]])
              na_h <- c(na_h, NA); nhh <- c(nhh, NA)
            }
          }
          m <- length(nhh)
          HR_value_h <- rep(NA_real_, m)
          if (m >= 2) { HR_value_h[1] <- nhh[2] / nhh[1]; HR_value_h[m] <- nhh[m-1] / nhh[m] }
          if (m > 2) { mid <- 2:(m-1); HR_value_h[mid] <- (nhh[mid-1] + nhh[mid+1]) / (2 * nhh[mid]) }
        }
        
      } else {
        # consecutiveSNP = TRUE, SNP-wise on haploids
        first_pos <- H$POS[1]
        last_pos  <- H$POS[nrow(H)]
        half_end  <- end / 2
        
        if (approach == "Bp.based") {
          # Bp.based, consecutive = TRUE, haploid
          compute_hap_stats_h <- function(L, R) {
            L0 <- max(L, first_pos); R0 <- min(R, last_pos)
            idx <- which(H$POS >= L0 & H$POS <= R0)
            if (length(idx) < minSNP) {
              sIdx <- if (length(idx) > 0) H$RowIndex[idx[1]] else NA
              eIdx <- if (length(idx) > 0) H$RowIndex[idx[length(idx)]] else NA
              return(list(valid=FALSE, na=NA, nh=NA,
                          startIndex=sIdx, endIndex=eIdx,
                          startPos=L0, endPos=R0))
            }
            subH <- H[idx, , drop = FALSE]
            index_pos <- c(subH$RowIndex[1], subH$POS[1], subH$RowIndex[nrow(subH)], subH$POS[nrow(subH)])
            hap_cols <- setdiff(colnames(subH), c("POS", "RowIndex"))
            M <- as.matrix(subH[, hap_cols, drop = FALSE])
            if (!is.character(M)) storage.mode(M) <- "character"
            keys <- apply(M, 2, paste0, collapse = "\r")
            grp  <- match(keys, unique(keys))
            cnt  <- tabulate(grp, nbins = max(grp))
            freq <- cnt / length(keys)
            list(valid=TRUE,
                 na=length(cnt),
                 nh=1 / sum(freq^2),
                 startIndex=index_pos[1], endIndex=index_pos[3],
                 startPos=index_pos[2], endPos=index_pos[4])
          }
          
          i0 <- which(H$POS >= max(start, first_pos))[1]
          if (!is.na(i0)) {
            for (i in i0:nrow(H)) {
              p <- H$POS[i]
              left_edge  <- (p - end) < first_pos
              right_edge <- (p + end) > last_pos
              
              if (!left_edge && !right_edge) {
                st_center <- compute_hap_stats_h(p - half_end, p + half_end)
                st_left   <- compute_hap_stats_h(p - end, p)
                st_right  <- compute_hap_stats_h(p, p + end)
                if (st_center$valid) {
                  CHROM_h <- c(CHROM_h, chr)
                  start.Index_h    <- c(start.Index_h,    st_center$startIndex)
                  end.Index_h      <- c(end.Index_h,      st_center$endIndex)
                  start.position_h <- c(start.position_h, st_center$startPos)
                  end.position_h   <- c(end.position_h,   st_center$endPos)
                  na_h <- c(na_h, st_center$na)
                  nhh  <- c(nhh,  st_center$nh)
                  HR_value_h <- c(HR_value_h, if (st_left$valid && st_right$valid) (st_left$nh + st_right$nh)/(2*st_center$nh) else NA)
                } else {
                  CHROM_h <- c(CHROM_h, chr)
                  start.Index_h    <- c(start.Index_h,    st_center$startIndex)
                  end.Index_h      <- c(end.Index_h,      st_center$endIndex)
                  start.position_h <- c(start.position_h, st_center$startPos)
                  end.position_h   <- c(end.position_h,   st_center$endPos)
                  na_h <- c(na_h, NA); nhh <- c(nhh, NA); HR_value_h <- c(HR_value_h, NA)
                }
                anchor_pos_h_chr <- c(anchor_pos_h_chr, p)
                
              } else if (left_edge) {
                st_near <- compute_hap_stats_h(p, p + end)
                st_far  <- compute_hap_stats_h(p + half_end, p + 3*half_end)
                if (st_near$valid) {
                  CHROM_h <- c(CHROM_h, chr)
                  start.Index_h    <- c(start.Index_h,    st_near$startIndex)
                  end.Index_h      <- c(end.Index_h,      st_near$endIndex)
                  start.position_h <- c(start.position_h, st_near$startPos)
                  end.position_h   <- c(end.position_h,   st_near$endPos)
                  na_h <- c(na_h, st_near$na)
                  nhh  <- c(nhh,  st_near$nh)
                  HR_value_h <- c(HR_value_h, if (st_far$valid) st_far$nh / st_near$nh else NA)
                } else {
                  CHROM_h <- c(CHROM_h, chr)
                  start.Index_h    <- c(start.Index_h,    st_near$startIndex)
                  end.Index_h      <- c(end.Index_h,      st_near$endIndex)
                  start.position_h <- c(start.position_h, st_near$startPos)
                  end.position_h   <- c(end.position_h,   st_near$endPos)
                  na_h <- c(na_h, NA); nhh <- c(nhh, NA); HR_value_h <- c(HR_value_h, NA)
                }
                anchor_pos_h_chr <- c(anchor_pos_h_chr, p)
                
              } else { # right edge
                st_near <- compute_hap_stats_h(p - end, p)
                st_far  <- compute_hap_stats_h(p - 3*half_end, p - half_end)
                if (st_near$valid) {
                  CHROM_h <- c(CHROM_h, chr)
                  start.Index_h    <- c(start.Index_h,    st_near$startIndex)
                  end.Index_h      <- c(end.Index_h,      st_near$endIndex)
                  start.position_h <- c(start.position_h, st_near$startPos)
                  end.position_h   <- c(end.position_h,   st_near$endPos)
                  na_h <- c(na_h, st_near$na)
                  nhh  <- c(nhh,  st_near$nh)
                  HR_value_h <- c(HR_value_h, if (st_far$valid) st_far$nh / st_near$nh else NA)
                } else {
                  CHROM_h <- c(CHROM_h, chr)
                  start.Index_h    <- c(start.Index_h,    st_near$startIndex)
                  end.Index_h      <- c(end.Index_h,      st_near$endIndex)
                  start.position_h <- c(start.position_h, st_near$startPos)
                  end.position_h   <- c(end.position_h,   st_near$endPos)
                  na_h <- c(na_h, NA); nhh <- c(nhh, NA); HR_value_h <- c(HR_value_h, NA)
                }
                anchor_pos_h_chr <- c(anchor_pos_h_chr, p)
              }
            }
          }
        } else {
          # SNP.based, consecutive = TRUE, haploid — index windows with maxGap salvage and anchoring
          win_len <- if ((end - start) > 0) (end - start) else end
          if (win_len <= 0) stop("In SNP.based mode (haploid, consecutiveSNP=TRUE) window length (end-start) must be > 0.")
          halfW <- floor(win_len/2)
          min_needed <- needed_snp_count(minSNP, win_len)
          
          for (i in 1:nrow(H)) {
            left_ok  <- (i - win_len + 1) >= 1
            right_ok <- (i + win_len - 1) <= nrow(H)
            
            if (left_ok && right_ok) {
              c_start <- i - floor((win_len - 1)/2)
              c_end   <- c_start + win_len - 1
              idx_c <- c_start:c_end
              idx_l <- (i - win_len + 1):i
              idx_r <- i:(i + win_len - 1)
              
              # Center, left, and right windows are salvaged with anchoring at i
              sC <- salvage_by_gap(H$POS, idx_c, maxGap, min_needed, anchor_idx = i)
              sL <- salvage_by_gap(H$POS, idx_l, maxGap, min_needed, anchor_idx = i)
              sR <- salvage_by_gap(H$POS, idx_r, maxGap, min_needed, anchor_idx = i)
              
              CHROM_h <- c(CHROM_h, chr)
              if (sC$valid) {
                subC <- H[sC$use_idx, , drop=FALSE]
                hap_cols <- setdiff(colnames(subC), c("POS","RowIndex"))
                MC <- as.matrix(subC[, hap_cols, drop=FALSE]); if (!is.character(MC)) storage.mode(MC) <- "character"
                keysC <- apply(MC, 2, paste0, collapse="\r"); grpC <- match(keysC, unique(keysC))
                cntC <- tabulate(grpC, nbins=max(grpC)); nhC <- 1/sum((cntC/length(keysC))^2)
                start.Index_h <- c(start.Index_h, subC$RowIndex[1]); end.Index_h <- c(end.Index_h, subC$RowIndex[nrow(subC)])
                start.position_h <- c(start.position_h, subC$POS[1]); end.position_h <- c(end.position_h, subC$POS[nrow(subC)])
                na_h <- c(na_h, length(cntC)); nhh <- c(nhh, nhC)
                
                if (sL$valid && sR$valid) {
                  subL <- H[sL$use_idx, , drop=FALSE]; subR <- H[sR$use_idx, , drop=FALSE]
                  ML <- as.matrix(subL[, hap_cols, drop=FALSE]); MR <- as.matrix(subR[, hap_cols, drop=FALSE])
                  if (!is.character(ML)) storage.mode(ML) <- "character"
                  if (!is.character(MR)) storage.mode(MR) <- "character"
                  keysL <- apply(ML, 2, paste0, collapse="\r"); grpL <- match(keysL, unique(keysL))
                  keysR <- apply(MR, 2, paste0, collapse="\r"); grpR <- match(keysR, unique(keysR))
                  cntL <- tabulate(grpL, nbins=max(grpL)); cntR <- tabulate(grpR, nbins=max(grpR))
                  nhL <- 1/sum((cntL/length(keysL))^2); nhR <- 1/sum((cntR/length(keysR))^2)
                  HR_value_h <- c(HR_value_h, (nhL + nhR)/(2*nhC))
                } else {
                  HR_value_h <- c(HR_value_h, NA)
                }
              } else {
                # Keep attempted center bounds when not salvageable
                start.Index_h    <- c(start.Index_h,    H$RowIndex[idx_c[1]])
                end.Index_h      <- c(end.Index_h,      H$RowIndex[idx_c[length(idx_c)]])
                start.position_h <- c(start.position_h, H$POS[idx_c[1]])
                end.position_h   <- c(end.position_h,   H$POS[idx_c[length(idx_c)]])
                na_h <- c(na_h, NA); nhh <- c(nhh, NA); HR_value_h <- c(HR_value_h, NA)
              }
              anchor_pos_h_chr <- c(anchor_pos_h_chr, H$POS[i])
              
            } else if (!left_ok) {
              idx_n <- i:(min(i + win_len - 1, nrow(H)))          # Near window anchored at i
              idx_f <- (i + halfW):(min(i + halfW + win_len - 1, nrow(H)))  # Far window without anchoring
              
              sN <- salvage_by_gap(H$POS, idx_n, maxGap, min_needed, anchor_idx = i)
              sF <- salvage_by_gap(H$POS, idx_f, maxGap, min_needed, anchor_idx = NULL)
              
              CHROM_h <- c(CHROM_h, chr)
              if (sN$valid) {
                subN <- H[sN$use_idx, , drop=FALSE]
                hap_cols <- setdiff(colnames(subN), c("POS","RowIndex"))
                MN <- as.matrix(subN[, hap_cols, drop=FALSE]); if (!is.character(MN)) storage.mode(MN) <- "character"
                keysN <- apply(MN, 2, paste0, collapse="\r"); grpN <- match(keysN, unique(keysN))
                cntN <- tabulate(grpN, nbins=max(grpN)); nhN <- 1/sum((cntN/length(keysN))^2)
                start.Index_h <- c(start.Index_h, subN$RowIndex[1]); end.Index_h <- c(end.Index_h, subN$RowIndex[nrow(subN)])
                start.position_h <- c(start.position_h, subN$POS[1]); end.position_h <- c(end.position_h, subN$POS[nrow(subN)])
                na_h <- c(na_h, length(cntN)); nhh <- c(nhh, nhN)
                
                HR_value_h <- c(HR_value_h, if (sF$valid) {
                  subF <- H[sF$use_idx, , drop=FALSE]
                  MF <- as.matrix(subF[, hap_cols, drop=FALSE]); if (!is.character(MF)) storage.mode(MF) <- "character"
                  keysF <- apply(MF, 2, paste0, collapse="\r"); grpF <- match(keysF, unique(keysF))
                  cntF <- tabulate(grpF, nbins=max(grpF)); nhF <- 1/sum((cntF/length(keysF))^2)
                  nhF / nhN
                } else NA)
              } else {
                start.Index_h    <- c(start.Index_h,    H$RowIndex[idx_n[1]])
                end.Index_h      <- c(end.Index_h,      H$RowIndex[idx_n[length(idx_n)]])
                start.position_h <- c(start.position_h, H$POS[idx_n[1]])
                end.position_h   <- c(end.position_h,   H$POS[idx_n[length(idx_n)]])
                na_h <- c(na_h, NA); nhh <- c(nhh, NA); HR_value_h <- c(HR_value_h, NA)
              }
              anchor_pos_h_chr <- c(anchor_pos_h_chr, H$POS[i])
              
            } else {
              idx_n <- (max(1, i - win_len + 1)):i                      # Near window anchored at i
              idx_f <- (max(1, i - halfW - win_len + 1)):(max(1, i - halfW))  # Far window without anchoring
              
              sN <- salvage_by_gap(H$POS, idx_n, maxGap, min_needed, anchor_idx = i)
              sF <- salvage_by_gap(H$POS, idx_f, maxGap, min_needed, anchor_idx = NULL)
              
              CHROM_h <- c(CHROM_h, chr)
              if (sN$valid) {
                subN <- H[sN$use_idx, , drop=FALSE]
                hap_cols <- setdiff(colnames(subN), c("POS","RowIndex"))
                MN <- as.matrix(subN[, hap_cols, drop=FALSE]); if (!is.character(MN)) storage.mode(MN) <- "character"
                keysN <- apply(MN, 2, paste0, collapse="\r"); grpN <- match(keysN, unique(keysN))
                cntN <- tabulate(grpN, nbins=max(grpN)); nhN <- 1/sum((cntN/length(keysN))^2)
                start.Index_h <- c(start.Index_h, subN$RowIndex[1]); end.Index_h <- c(end.Index_h, subN$RowIndex[nrow(subN)])
                start.position_h <- c(start.position_h, subN$POS[1]); end.position_h <- c(end.position_h, subN$POS[nrow(subN)])
                na_h <- c(na_h, length(cntN)); nhh <- c(nhh, nhN)
                
                HR_value_h <- c(HR_value_h, if (sF$valid) {
                  subF <- H[sF$use_idx, , drop=FALSE]
                  MF <- as.matrix(subF[, hap_cols, drop=FALSE]); if (!is.character(MF)) storage.mode(MF) <- "character"
                  keysF <- apply(MF, 2, paste0, collapse="\r"); grpF <- match(keysF, unique(keysF))
                  cntF <- tabulate(grpF, nbins=max(grpF)); nhF <- 1/sum((cntF/length(keysF))^2)
                  nhF / nhN
                } else NA)
              } else {
                start.Index_h    <- c(start.Index_h,    H$RowIndex[idx_n[1]])
                end.Index_h      <- c(end.Index_h,      H$RowIndex[idx_n[length(idx_n)]])
                start.position_h <- c(start.position_h, H$POS[idx_n[1]])
                end.position_h   <- c(end.position_h,   H$POS[idx_n[length(idx_n)]])
                na_h <- c(na_h, NA); nhh <- c(nhh, NA); HR_value_h <- c(HR_value_h, NA)
              }
              anchor_pos_h_chr <- c(anchor_pos_h_chr, H$POS[i])
            }
          }
        }
      }
      
      if (length(nhh) > 0 || length(CHROM_h) > 0) {
        HR_chr <- as.data.frame(cbind(
          CHROM         = CHROM_h,
          start.Index   = start.Index_h,
          start.position= start.position_h,
          end.Index     = end.Index_h,
          end.position  = end.position_h,
          na            = na_h,
          nh            = nhh,
          HR_value      = HR_value_h
        ))
        HR_chr <- as.data.frame(lapply(HR_chr, as.numeric))
        if (isTRUE(consecutiveSNP)) {
          HR_chr$anchor.pos <- as.numeric(anchor_pos_h_chr)
        }
        HR_haplo <- rbind(HR_haplo, HR_chr)
      }
    }
    
    assign("HR_haplo", HR_haplo, inherits = TRUE)
    message(Sys.time(), " | Haploid chromosomes HR computed.")
  }
  
  message(Sys.time(), " | Autosomal data preparation...")
  
  # Initialize storage objects
  fazes <- CHROM <- na <- start.Index <- end.Index <- start.position <- end.position <- nhh <- nh <- HR_value <- HRdata <- NULL
  anchor_pos <- NULL   # Used only when consecutiveSNP = TRUE on autosomes
  
  # Parse genotype columns into phased alleles h1 and h2 for each individual
  for (f in 10:ncol(vcf)) {
    sample_genotypes <- as.character(vcf[, f])
    first_non_na <- sample_genotypes[which(!is.na(sample_genotypes))[1]]
    if (grepl("/", first_non_na, fixed = TRUE)) {
      separator <- "/"
    } else if (grepl("|", first_non_na, fixed = TRUE)) {
      separator <- "|"
    } else {
      stop(paste("Unknown genotype separator in sample:", colnames(vcf)[f]))
    }
    alleles <- do.call("rbind", strsplit(sample_genotypes, separator, fixed = TRUE))
    fazes <- cbind(fazes, alleles[, 1], alleles[, 2])
  }
  
  vcf <- cbind(vcf[, 1:9], fazes)
  
  # Identify chromosomes present in the VCF
  Chr <- as.numeric(unique(vcf$CHROM))
  vcf <- cbind(1:nrow(vcf), vcf)
  
  message(Sys.time(), " | Autosomal data prepared.")
  message(Sys.time(), " | HR computation on autosomes...")
  
  pb <- txtProgressBar(min = 1, max = length(Chr), style = 3)
  
  for (s in 1:length(Chr)) {
    setTxtProgressBar(pb, s)
    
    vcf_chr <- subset(vcf, vcf$CHROM == Chr[s])
    pos_vec <- vcf_chr$POS
    last_pos <- pos_vec[length(pos_vec)]
    
    if (!isTRUE(consecutiveSNP)) {
      # Sliding by bp distance
      if (approach == "Bp.based") {
        # Bp.based, non-consecutive, autosomes
        startt <- start
        endd <- end
        if (endd <= vcf_chr[1, 3]) {
          startt <- vcf_chr[1, 3]
          endd <- startt + endd
        }
        while (endd <= vcf_chr[nrow(vcf_chr), 3]) {
          if (endd != vcf_chr[nrow(vcf_chr), 3] & (endd + slide) > vcf_chr[nrow(vcf_chr), 3]) {
            endd <- vcf_chr[nrow(vcf_chr), 3]
          }
          k <- vcf_chr[which(vcf_chr$POS >= startt & vcf_chr$POS <= endd), c(-2, -4:-10)]
          if (nrow(k) >= minSNP) {
            CHROM <- c(CHROM, vcf_chr[1, 2])
            index_pos <- k[c(1, nrow(k)), 1:2]
            M <- as.matrix(k[, -1:-2, drop = FALSE])
            if (!is.character(M)) storage.mode(M) <- "character"
            keys <- apply(M, 2, paste0, collapse = "\r")
            grp  <- match(keys, unique(keys))
            cnt  <- tabulate(grp, nbins = max(grp))
            freq <- cnt / length(keys)
            eff  <- 1 / sum(freq^2)
            na_current <- length(cnt)
            start.position <- c(start.position, index_pos[1, 2])
            end.position   <- c(end.position,   index_pos[2, 2])
            start.Index    <- c(start.Index,    index_pos[1, 1])
            end.Index      <- c(end.Index,      index_pos[2, 1])
            nhh <- c(nhh, eff)
            na  <- c(na,  na_current)
          } else if (nrow(k) > 0) {
            index_pos <- k[c(1, nrow(k)), 1:2]
            CHROM <- c(CHROM, vcf_chr[1, 2])
            start.position <- c(start.position, index_pos[1, 2])
            end.position   <- c(end.position,   index_pos[2, 2])
            start.Index    <- c(start.Index,    index_pos[1, 1])
            end.Index      <- c(end.Index,      index_pos[2, 1])
            na  <- c(na,  NA)
            nhh <- c(nhh, NA)
          }
          startt <- startt + slide
          endd   <- endd + slide
        }
        m <- length(nhh)
        HR_value_chr <- rep(NA_real_, m)
        if (m >= 2) {
          HR_value_chr[1] <- nhh[2] / nhh[1]
          HR_value_chr[m] <- nhh[m-1] / nhh[m]
        }
        if (m > 2) {
          mid <- 2:(m-1)
          HR_value_chr[mid] <- (nhh[mid-1] + nhh[mid+1]) / (2 * nhh[mid])
        }
        HR_value <- c(HR_value, HR_value_chr)
        nh  <- c(nh, nhh)
        nhh <- NULL
      } else {
        # SNP.based, non-consecutive, autosomes — index windows with maxGap salvage
        win_len <- if ((end - start) > 0) (end - start) else end
        if (win_len <= 0) stop("In SNP.based mode (autosomes, non-consecutive) window length (end-start) must be > 0.")
        stepS <- max(1, as.integer(slide))
        min_needed <- needed_snp_count(minSNP, win_len)
        
        tmp_nh <- tmp_na <- tmp_si <- tmp_ei <- tmp_sp <- tmp_ep <- NULL
        for (iL in seq.int(1 + start, nrow(vcf_chr), by = stepS)) {
          iR <- min(iL + win_len - 1, nrow(vcf_chr))
          idx <- iL:iR
          
          sv <- salvage_by_gap(vcf_chr$POS, idx, maxGap, min_needed, anchor_idx = NULL)
          CHROM <- c(CHROM, vcf_chr[1, 2])
          
          if (sv$valid) {
            k  <- vcf_chr[sv$use_idx, c(-2, -4:-10)]
            ip <- k[c(1, nrow(k)), 1:2]
            M  <- as.matrix(k[, -1:-2, drop = FALSE]); if (!is.character(M)) storage.mode(M) <- "character"
            keys <- apply(M, 2, paste0, collapse = "\r"); grp <- match(keys, unique(keys))
            cnt  <- tabulate(grp, nbins = max(grp)); eff <- 1/sum((cnt/ncol(M))^2)
            tmp_si <- c(tmp_si, ip[1,1]); tmp_ei <- c(tmp_ei, ip[2,1])
            tmp_sp <- c(tmp_sp, ip[1,2]); tmp_ep <- c(tmp_ep, ip[2,2])
            tmp_na <- c(tmp_na, length(cnt)); tmp_nh <- c(tmp_nh, eff)
          } else {
            k0 <- vcf_chr[idx, c(-2, -4:-10)]
            ip0 <- k0[c(1, nrow(k0)), 1:2]
            tmp_si <- c(tmp_si, ip0[1,1]); tmp_ei <- c(tmp_ei, ip0[2,1])
            tmp_sp <- c(tmp_sp, ip0[1,2]); tmp_ep <- c(tmp_ep, ip0[2,2])
            tmp_na <- c(tmp_na, NA);       tmp_nh <- c(tmp_nh, NA)
          }
        }
        start.Index <- c(start.Index, tmp_si); end.Index <- c(end.Index, tmp_ei)
        start.position <- c(start.position, tmp_sp); end.position <- c(end.position, tmp_ep)
        na <- c(na, tmp_na); nh <- c(nh, tmp_nh)
        m <- length(tmp_nh); HR_value_chr <- rep(NA_real_, m)
        if (m >= 2) { HR_value_chr[1] <- tmp_nh[2] / tmp_nh[1]; HR_value_chr[m] <- tmp_nh[m-1] / tmp_nh[m] }
        if (m > 2) { mid <- 2:(m-1); HR_value_chr[mid] <- (tmp_nh[mid-1] + tmp_nh[mid+1])/(2*tmp_nh[mid]) }
        HR_value <- c(HR_value, HR_value_chr)
      }
      
    } else {
      # consecutiveSNP = TRUE, SNP-wise on autosomes
      first_pos <- pos_vec[1]
      last_pos  <- pos_vec[length(pos_vec)]
      half_end  <- end / 2
      
      if (approach == "Bp.based") {
        # Bp.based, consecutive = TRUE, autosomes
        compute_hap_stats <- function(L, R) {
          L0 <- max(L, first_pos); R0 <- min(R, last_pos)
          idx <- which(vcf_chr$POS >= L0 & vcf_chr$POS <= R0)
          if (length(idx) < minSNP) {
            sIdx <- if (length(idx) > 0) vcf_chr[idx[1], 1] else NA
            eIdx <- if (length(idx) > 0) vcf_chr[idx[length(idx)], 1] else NA
            return(list(valid=FALSE, na=NA, nh=NA,
                        startIndex=sIdx, endIndex=eIdx,
                        startPos=L0, endPos=R0))
          }
          k  <- vcf_chr[idx, c(-2, -4:-10)]
          ip <- k[c(1, nrow(k)), 1:2]
          M  <- as.matrix(k[, -1:-2, drop = FALSE])
          if (!is.character(M)) storage.mode(M) <- "character"
          keys <- apply(M, 2, paste0, collapse = "\r")
          grp  <- match(keys, unique(keys))
          cnt  <- tabulate(grp, nbins = max(grp))
          freq <- cnt / length(keys)
          list(valid=TRUE,
               na=length(cnt),
               nh=1 / sum(freq^2),
               startIndex=ip[1, 1], endIndex=ip[2, 1],
               startPos=ip[1, 2],  endPos=ip[2, 2])
        }
        
        i0 <- which(pos_vec >= max(start, first_pos))[1]
        if (!is.na(i0)) {
          for (i in i0:nrow(vcf_chr)) {
            p <- vcf_chr$POS[i]
            left_edge  <- (p - end) < first_pos
            right_edge <- (p + end) > last_pos
            
            if (!left_edge && !right_edge) {
              st_center <- compute_hap_stats(p - half_end, p + half_end)
              st_left   <- compute_hap_stats(p - end, p)
              st_right  <- compute_hap_stats(p, p + end)
              if (st_center$valid) {
                CHROM <- c(CHROM, vcf_chr[1, 2])
                start.Index    <- c(start.Index,    st_center$startIndex)
                end.Index      <- c(end.Index,      st_center$endIndex)
                start.position <- c(start.position, st_center$startPos)
                end.position   <- c(end.position,   st_center$endPos)
                na  <- c(na,  st_center$na)
                nhh <- c(nhh, st_center$nh)
                HR_value <- c(HR_value,
                              if (st_left$valid && st_right$valid)
                                (st_left$nh + st_right$nh) / (2 * st_center$nh)
                              else
                                NA)
              } else {
                CHROM <- c(CHROM, vcf_chr[1, 2])
                start.Index    <- c(start.Index, st_center$startIndex)
                end.Index      <- c(end.Index,   st_center$endIndex)
                start.position <- c(start.position, st_center$startPos)
                end.position   <- c(end.position,   st_center$endPos)
                na <- c(na, NA); nhh <- c(nhh, NA); HR_value <- c(HR_value, NA)
              }
              anchor_pos <- c(anchor_pos, p)
              
            } else if (left_edge) {
              st_near <- compute_hap_stats(p, p + end)
              st_far  <- compute_hap_stats(p + half_end, p + 3*half_end)
              if (st_near$valid) {
                CHROM <- c(CHROM, vcf_chr[1, 2])
                start.Index    <- c(start.Index,    st_near$startIndex)
                end.Index      <- c(end.Index,      st_near$endIndex)
                start.position <- c(start.position, st_near$startPos)
                end.position   <- c(end.position,   st_near$endPos)
                na  <- c(na,  st_near$na)
                nhh <- c(nhh, st_near$nh)
                HR_value <- c(HR_value, if (st_far$valid) st_far$nh / st_near$nh else NA)
              } else {
                CHROM <- c(CHROM, vcf_chr[1, 2])
                start.Index    <- c(start.Index, st_near$startIndex)
                end.Index      <- c(end.Index,   st_near$endIndex)
                start.position <- c(start.position, st_near$startPos)
                end.position   <- c(end.position,   st_near$endPos)
                na <- c(na, NA); nhh <- c(nhh, NA); HR_value <- c(HR_value, NA)
              }
              anchor_pos <- c(anchor_pos, p)
              
            } else { # right edge
              st_near <- compute_hap_stats(p - end, p)
              st_far  <- compute_hap_stats(p - 3*half_end, p - half_end)
              if (st_near$valid) {
                CHROM <- c(CHROM, vcf_chr[1, 2])
                start.Index    <- c(start.Index,    st_near$startIndex)
                end.Index      <- c(end.Index,      st_near$endIndex)
                start.position <- c(start.position, st_near$startPos)
                end.position   <- c(end.position,   st_near$endPos)
                na  <- c(na,  st_near$na)
                nhh <- c(nhh, st_near$nh)
                HR_value <- c(HR_value, if (st_far$valid) st_far$nh / st_near$nh else NA)
              } else {
                CHROM <- c(CHROM, vcf_chr[1, 2])
                start.Index    <- c(start.Index, st_near$startIndex)
                end.Index      <- c(end.Index,   st_near$endIndex)
                start.position <- c(start.position, st_near$startPos)
                end.position   <- c(end.position,   st_near$endPos)
                na <- c(na, NA); nhh <- c(nhh, NA); HR_value <- c(HR_value, NA)
              }
              anchor_pos <- c(anchor_pos, p)
            }
          }
        }
      } else {
        # SNP.based, consecutive = TRUE, autosomes — index windows with maxGap salvage and anchoring
        win_len <- if ((end - start) > 0) (end - start) else end
        if (win_len <= 0) stop("In SNP.based mode (autosomes, consecutiveSNP=TRUE) window length (end-start) must be > 0.")
        halfW <- floor(win_len/2)
        min_needed <- needed_snp_count(minSNP, win_len)
        
        for (i in 1:nrow(vcf_chr)) {
          left_ok  <- (i - win_len + 1) >= 1
          right_ok <- (i + win_len - 1) <= nrow(vcf_chr)
          
          if (left_ok && right_ok) {
            c_start <- i - floor((win_len - 1)/2)
            c_end   <- c_start + win_len - 1
            idx_c <- c_start:c_end
            idx_l <- (i - win_len + 1):i
            idx_r <- i:(i + win_len - 1)
            
            sC <- salvage_by_gap(vcf_chr$POS, idx_c, maxGap, min_needed, anchor_idx = i)
            sL <- salvage_by_gap(vcf_chr$POS, idx_l, maxGap, min_needed, anchor_idx = i)
            sR <- salvage_by_gap(vcf_chr$POS, idx_r, maxGap, min_needed, anchor_idx = i)
            
            CHROM <- c(CHROM, vcf_chr[1, 2])
            # Record center bounds — salvaged if valid, otherwise attempted bounds
            if (sC$valid) {
              kC <- vcf_chr[sC$use_idx, c(-2, -4:-10)]
              ipC <- kC[c(1, nrow(kC)), 1:2]
              start.Index    <- c(start.Index, ipC[1,1]); end.Index <- c(end.Index, ipC[2,1])
              start.position <- c(start.position, ipC[1,2]); end.position <- c(end.position, ipC[2,2])
              MC <- as.matrix(kC[, -1:-2, drop=FALSE]); if (!is.character(MC)) storage.mode(MC) <- "character"
              keysC <- apply(MC, 2, paste0, collapse="\r"); grpC <- match(keysC, unique(keysC))
              cntC <- tabulate(grpC, nbins=max(grpC)); nhC <- 1/sum((cntC/ncol(MC))^2)
              na <- c(na, length(cntC)); nhh <- c(nhh, nhC)
              
              if (sL$valid && sR$valid) {
                kL <- vcf_chr[sL$use_idx, c(-2, -4:-10)]
                kR <- vcf_chr[sR$use_idx, c(-2, -4:-10)]
                ML <- as.matrix(kL[, -1:-2, drop=FALSE]); MR <- as.matrix(kR[, -1:-2, drop=FALSE])
                if (!is.character(ML)) storage.mode(ML) <- "character"
                if (!is.character(MR)) storage.mode(MR) <- "character"
                keysL <- apply(ML, 2, paste0, collapse="\r"); grpL <- match(keysL, unique(keysL))
                keysR <- apply(MR, 2, paste0, collapse="\r"); grpR <- match(keysR, unique(keysR))
                cntL <- tabulate(grpL, nbins=max(grpL)); cntR <- tabulate(grpR, nbins=max(grpR))
                nhL <- 1/sum((cntL/ncol(ML))^2); nhR <- 1/sum((cntR/ncol(MR))^2)
                HR_value <- c(HR_value, (nhL + nhR)/(2*nhC))
              } else {
                HR_value <- c(HR_value, NA)
              }
            } else {
              k0 <- vcf_chr[idx_c, c(-2, -4:-10)]
              ip0 <- k0[c(1, nrow(k0)), 1:2]
              start.Index    <- c(start.Index, ip0[1,1]); end.Index <- c(end.Index, ip0[2,1])
              start.position <- c(start.position, ip0[1,2]); end.position <- c(end.position, ip0[2,2])
              na <- c(na, NA); nhh <- c(nhh, NA); HR_value <- c(HR_value, NA)
            }
            anchor_pos <- c(anchor_pos, vcf_chr$POS[i])
            
          } else if (!left_ok) {
            idx_n <- i:(min(i + win_len - 1, nrow(vcf_chr)))                 # Near window anchored at i
            idx_f <- (i + halfW):(min(i + halfW + win_len - 1, nrow(vcf_chr))) # Far window without anchoring
            
            sN <- salvage_by_gap(vcf_chr$POS, idx_n, maxGap, min_needed, anchor_idx = i)
            sF <- salvage_by_gap(vcf_chr$POS, idx_f, maxGap, min_needed, anchor_idx = NULL)
            
            CHROM <- c(CHROM, vcf_chr[1, 2])
            if (sN$valid) {
              kN <- vcf_chr[sN$use_idx, c(-2, -4:-10)]
              ipN <- kN[c(1, nrow(kN)), 1:2]
              start.Index    <- c(start.Index, ipN[1,1]); end.Index <- c(end.Index, ipN[2,1])
              start.position <- c(start.position, ipN[1,2]); end.position <- c(end.position, ipN[2,2])
              MN <- as.matrix(kN[, -1:-2, drop=FALSE]); if (!is.character(MN)) storage.mode(MN) <- "character"
              keysN <- apply(MN, 2, paste0, collapse="\r"); grpN <- match(keysN, unique(keysN))
              cntN <- tabulate(grpN, nbins=max(grpN)); nhN <- 1/sum((cntN/ncol(MN))^2)
              na <- c(na, length(cntN)); nhh <- c(nhh, nhN)
              
              HR_value <- c(HR_value, if (sF$valid) {
                kF <- vcf_chr[sF$use_idx, c(-2, -4:-10)]
                MF <- as.matrix(kF[, -1:-2, drop=FALSE]); if (!is.character(MF)) storage.mode(MF) <- "character"
                keysF <- apply(MF, 2, paste0, collapse="\r"); grpF <- match(keysF, unique(keysF))
                cntF <- tabulate(grpF, nbins=max(grpF)); nhF <- 1/sum((cntF/ncol(MF))^2)
                nhF / nhN
              } else NA)
            } else {
              k0 <- vcf_chr[idx_n, c(-2, -4:-10)]
              ip0 <- k0[c(1, nrow(k0)), 1:2]
              start.Index    <- c(start.Index, ip0[1,1]); end.Index <- c(end.Index, ip0[2,1])
              start.position <- c(start.position, ip0[1,2]); end.position <- c(end.position, ip0[2,2])
              na <- c(na, NA); nhh <- c(nhh, NA); HR_value <- c(HR_value, NA)
            }
            anchor_pos <- c(anchor_pos, vcf_chr$POS[i])
            
          } else {
            idx_n <- (max(1, i - win_len + 1)):i                                # Near window anchored at i
            idx_f <- (max(1, i - halfW - win_len + 1)):(max(1, i - halfW))      # Far window without anchoring
            
            sN <- salvage_by_gap(vcf_chr$POS, idx_n, maxGap, min_needed, anchor_idx = i)
            sF <- salvage_by_gap(vcf_chr$POS, idx_f, maxGap, min_needed, anchor_idx = NULL)
            
            CHROM <- c(CHROM, vcf_chr[1, 2])
            if (sN$valid) {
              kN <- vcf_chr[sN$use_idx, c(-2, -4:-10)]
              ipN <- kN[c(1, nrow(kN)), 1:2]
              start.Index    <- c(start.Index, ipN[1,1]); end.Index <- c(end.Index, ipN[2,1])
              start.position <- c(start.position, ipN[1,2]); end.position <- c(end.position, ipN[2,2])
              MN <- as.matrix(kN[, -1:-2, drop=FALSE]); if (!is.character(MN)) storage.mode(MN) <- "character"
              keysN <- apply(MN, 2, paste0, collapse="\r"); grpN <- match(keysN, unique(keysN))
              cntN <- tabulate(grpN, nbins=max(grpN)); nhN <- 1/sum((cntN/ncol(MN))^2)
              na <- c(na, length(cntN)); nhh <- c(nhh, nhN)
              
              HR_value <- c(HR_value, if (sF$valid) {
                kF <- vcf_chr[sF$use_idx, c(-2, -4:-10)]
                MF <- as.matrix(kF[, -1:-2, drop=FALSE]); if (!is.character(MF)) storage.mode(MF) <- "character"
                keysF <- apply(MF, 2, paste0, collapse="\r"); grpF <- match(keysF, unique(keysF))
                cntF <- tabulate(grpF, nbins=max(grpF)); nhF <- 1/sum((cntF/ncol(MF))^2)
                nhF / nhN
              } else NA)
            } else {
              k0 <- vcf_chr[idx_n, c(-2, -4:-10)]
              ip0 <- k0[c(1, nrow(k0)), 1:2]
              start.Index    <- c(start.Index, ip0[1,1]); end.Index <- c(end.Index, ip0[2,1])
              start.position <- c(start.position, ip0[1,2]); end.position <- c(end.position, ip0[2,2])
              na <- c(na, NA); nhh <- c(nhh, NA); HR_value <- c(HR_value, NA)
            }
            anchor_pos <- c(anchor_pos, vcf_chr$POS[i])
          }
        }
      }
      nh  <- c(nh, nhh)
      nhh <- NULL
    }
  }
  
  if (!isTRUE(consecutiveSNP)) {
    HR <- as.data.frame(cbind(CHROM, start.Index, start.position, end.Index, end.position, na, nh, HR_value))
    HR <- as.data.frame(lapply(HR, as.numeric))
  } else {
    HR <- as.data.frame(cbind(CHROM, start.Index, start.position, end.Index, end.position, na, nh, HR_value))
    HR <- as.data.frame(lapply(HR, as.numeric))
    HR$anchor.pos <- as.numeric(anchor_pos)
  }
  
  close(pb)
  
  # Optional — merge results if haploid chromosomes were included
  if (HaploidExistance == TRUE) {
    if (isTRUE(consecutiveSNP)) {
      # HR_haplo already contains anchor.pos when SNP-wise
      HR <- as.data.frame(rbind(HR, HR_haplo))
    } else {
      HR <- as.data.frame(rbind(HR, HR_haplo))
    }
    
    key_map    <- paste(as.character(vcf_org$CHROM), as.numeric(vcf_org$POS), sep = ":")
    start_keys <- paste(as.character(HR$CHROM),      as.numeric(HR$start.position), sep = ":")
    end_keys   <- paste(as.character(HR$CHROM),      as.numeric(HR$end.position),   sep = ":")
    
    HR$start.Index <- match(start_keys, key_map)
    HR$end.Index   <- match(end_keys,   key_map)
    
    if (anyNA(HR$start.Index) || anyNA(HR$end.Index)) {
      warning(sum(is.na(HR$start.Index)) + sum(is.na(HR$end.Index)),
              " window boundary(ies) could not be mapped to .vcf file (CHROM:POS mismatch).")
    }
    
    HR <- HR[order(HR$start.Index, HR$end.Index), , drop = FALSE]
    rownames(HR) <- NULL
    Chr <- unique(vcf_org$CHROM)
  }
  
  if (!isTRUE(consecutiveSNP)) {
    Window <- paste0("w", seq_len(nrow(HR)))
    HR <- cbind(Window = Window, HR)
  } else {
    # Map SNPname (ID) using CHROM and anchor.pos for both autosomes and haploids
    key_map <- paste(as.character(vcf_org$CHROM), as.numeric(vcf_org$POS), sep=":")
    id_map  <- as.character(vcf_org$ID)
    anchor_keys <- paste(as.character(HR$CHROM), as.numeric(HR$anchor.pos), sep=":")
    SNPname <- id_map[match(anchor_keys, key_map)]
    HR$SNPname  <- SNPname
    HR$position <- HR$anchor.pos
    HR$anchor.pos <- NULL
    HR <- HR[, c("SNPname", "position", setdiff(names(HR), c("SNPname","position")))]
  }
  
  # Optional — normalize and compute p-values
  if (WithNormalisation == TRUE) {
    if (isTRUE(NormalisationByChr)) {
      HRdata <- NULL
      for (s in 1:length(Chr)) {
        HR_Chr <- HR[HR$CHROM == Chr[s], ]
        mean_val <- mean(HR_Chr$HR_value, na.rm = TRUE)
        sd_val   <- sd(HR_Chr$HR_value,   na.rm = TRUE)
        if (is.na(sd_val) || sd_val == 0) {
          HR_Chr$Z_value <- NA
        } else {
          HR_Chr$Z_value <- (HR_Chr$HR_value - mean_val) / sd_val
        }
        HR_Chr$HRiD_Pvalue    <- ifelse(is.na(HR_Chr$Z_value), NA, pnorm(-(HR_Chr$Z_value)))
        HR_Chr$HRiD_LogPvalue <- ifelse(is.na(HR_Chr$HRiD_Pvalue), NA, -log10(HR_Chr$HRiD_Pvalue))
        HR_Chr$HRiP_Pvalue    <- ifelse(is.na(HR_Chr$Z_value), NA, pnorm(HR_Chr$Z_value))
        HR_Chr$HRiP_LogPvalue <- ifelse(is.na(HR_Chr$HRiP_Pvalue), NA, -log10(HR_Chr$HRiP_Pvalue))
        HRdata <- rbind(HRdata, HR_Chr)
      }
      message(Sys.time(), " | Analysis completed successfully.")
      return(HRdata)
    } else {
      mean_val <- mean(HR$HR_value, na.rm = TRUE)
      sd_val   <- sd(HR$HR_value,   na.rm = TRUE)
      if (is.na(sd_val) || sd_val == 0) {
        HR$Z_value <- NA
      } else {
        HR$Z_value <- (HR$HR_value - mean_val) / sd_val
      }
      HR$HRiD_Pvalue    <- ifelse(is.na(HR$Z_value), NA, pnorm(-(HR$Z_value)))
      HR$HRiD_LogPvalue <- ifelse(is.na(HR$HRiD_Pvalue), NA, -log10(HR$HRiD_Pvalue))
      HR$HRiP_Pvalue    <- ifelse(is.na(HR$Z_value), NA, pnorm(HR$Z_value))
      HR$HRiP_LogPvalue <- ifelse(is.na(HR$HRiP_Pvalue), NA, -log10(HR$HRiP_Pvalue))
      message(Sys.time(), " | Analysis completed successfully.")
      return(HR)
    }
  }
  
  message(Sys.time(), " | Analysis completed successfully.")
  return(HR)
}
