

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


## ------------------------------------------------------------------------
## 2. Load HaplotypeRichness.Estimation() function
## ------------------------------------------------------------------------
## Make sure the path points to the function file inside the repo
source("R/HaplotypeRichness_Estimation.R")


## ------------------------------------------------------------------------
## 3. Run tests with different modes
## ------------------------------------------------------------------------

## Each call below demonstrates a different combination of approach and windowing mode.
## Results are stored in objects Test1–Test4.

# Sliding windows, bp-based
Test1 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = FALSE,
                                      approach = "Bp.based")

# Sliding windows, SNP-based
Test2 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = FALSE,
                                      approach = "SNP.based")

# Consecutive-SNP mode, bp-based
Test3 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = TRUE,
                                      approach = "Bp.based")

# Consecutive-SNP mode, SNP-based
Test4 <- HaplotypeRichness.Estimation(vcf = vcf,
                                      HaploidExistance = TRUE,
                                      consecutiveSNP = TRUE,
                                      approach = "SNP.based")




