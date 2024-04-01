
load("reference_data/hapmap3_snp.RDATA")

# Load packages bigsnpr and bigstatsr
library(bigsnpr)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("genotype_data/data.rds")
obj.bigSNP$genotypes <- snp_fastImputeSimple(Gna = obj.bigSNP$genotypes,method = "mean2",ncores = 10)


# Get aliases for useful slot
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

NCORES <- nb_cores()

# Read external summary statistics
sumstats <- bigreadr::fread2("GWAS_data/stroke_Any.std")
str(sumstats)

library(dplyr)
#change all the lowercase alleles into uppercase alleles
sumstats$a1 <- toupper(sumstats$a1)
sumstats$a2 <- toupper(sumstats$a2)

names(sumstats)[names(sumstats) == "a1"] <- "a0"
names(sumstats)[names(sumstats) == "a2"] <- "a1"
names(sumstats)[names(sumstats) == "b"] <- "beta"
str(sumstats)
# sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
# sumstats$n_case <- sumstats$n_control <- NULL
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map)
# To convert physical positions (in bp) to genetic positions (in cM), use
# POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

# To avoid downloading "large" files, `POS2` has been precomputed here
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/tmp_dir/")
corr <- NULL
ld <- NULL

tmp <- tempfile(tmpdir = "/tmp_dir/")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(df_beta$chr == chr)
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
      G,
      ind.col = ind.chr2,
      ncores = 11,
      infos.pos = POS2[ind.chr2],
      size = 3 / 1000
    )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.frame(obj.bigSNP$fam)
# Rename fam order
fam.order <- fam.order %>% 
  rename("FID" = "family.ID",
         "IID" = "sample.ID")
df_beta$n_eff <- df_beta$ncol
df_beta$beta_se <- df_beta$se
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref
ldsc_h2_est <- ldsc[["h2"]]
set.seed(1)  # to get the same result every time
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = ldsc_h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = NCORES,
  # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)
str(multi_auto, max.level = 1)
str(multi_auto[[1]], max.level = 1)
library(ggplot2)
auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
# `range` should be between 0 and 2
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

if(is.null(obj.bigSNP)){
  obj.bigSNP <- snp_attach("genotyp_data/data.rds")
}
genotype <- obj.bigSNP$genotypes
ind.test <- 1:nrow(genotype)

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

#test if the extracted data complies the results

rows_to_extract <- c(1, 2, 3, 4, 5)  # Example: Extracting the first three rows

library(bigstatsr)

# Assuming 'G' is your FBM from a bigSNP object
G <- obj.bigSNP$genotypes

# Path for the new FBM file
new_fbm_path <- tempfile("new_fbm", fileext = ".bk")

# Copying selected rows (and all columns) into a new FBM
G_subset <- big_copy(G, ind.row = rows_to_extract, ind.col = 1:ncol(G), 
                     backingfile = new_fbm_path)

ind.test <- 1:5  # Adjusting to use all rows in the subset
####most important the upper codes just used to test

pred_auto <- big_prodVec(G_subset, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])


beta_auto
pred_auto

PRS <- data.frame("IID" = fam.order$IID,
           "prs" = pred_auto)

# Assuming df is your data frame
write.csv(PRS, "PRS.csv", row.names = FALSE)