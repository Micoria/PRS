# the scripts of writing the function of calculating PRS
PRSAnalysis <- function(reference_data_path, summary_statistics_path, gene_data_path,outputDir, nCores = 10) {
  # Ensure necessary libraries are loaded
  if (!requireNamespace("bigsnpr", quietly = TRUE)) stop("Please install the 'bigsnpr' package.")
  if (!requireNamespace("bigstatsr", quietly = TRUE)) stop("Please install the 'bigstatsr' package.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the 'dplyr' package.")
  #if you don't need graphs to tell the results, then you don't need this
  #if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install the 'ggplot2' package.")
  
  #load reference data
  #rm(list = ls()) remove all the listed docu
  
  load(reference_data_path)
  
  # Read from bed/bim/fam, it generates .bk and .rds files.
  gene_data <- snp_readBed(gene_data_path)
  
  
  # Load genotype data - reference data
  obj.bigSNP <- snp_attach(gene_data) 
  # Perform imputation and other preprocessing steps as required
  # (You might need to adjust these steps based on your data)
  obj.bigSNP$genotypes <- snp_fastImputeSimple(Gna = obj.bigSNP$genotypes, method = "mean2", ncores = nCores)
  
  # Get aliases for useful slot
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  NCORES <- nb_cores()
  
  # Load and preprocess summary statistics
  summary_statistics <- bigreadr::fread2(summary_statistics_path)
  summary_statistics <- summary_statistics %>%
    mutate(a1 = toupper(a1), a2 = toupper(a2)) %>% rename(a0 = a1, a1 = a2, beta = b)
  # Continue with summary statistics preprocessing...
  # summary_statistics$n_eff <- 4 / (1 / summary_statistics$n_case + 1 / summary_statistics$n_control)
  # summary_statistics$n_case <- summary_statistics$n_control <- NULL
  #in this case, the number of the objects is very huge, so which can treat them identically.
  summary_statistics$n_eff <- summary_statistics$N
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  df_beta <- snp_match(summary_statistics, map)
  # To convert physical positions (in bp) to genetic positions (in cM), use
  # POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)
  
  # To avoid downloading "large" files, `POS2` has been precomputed here
  POS2 <- snp_asGeneticPos(CHR, POS, dir = outputDir)
  corr <- NULL
  ld <- NULL
  
  tmp <- tempfile(tmpdir = outputDir)
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
  #if you want to see the graph of multi_auto
  #library(ggplot2)
  #auto <- multi_auto[[1]]  # first chain
  #plot_grid(
    #qplot(y = auto$path_p_est) + 
      #theme_bigstatsr() + 
      #geom_hline(yintercept = auto$p_est, col = "blue") +
      #scale_y_log10() +
      #labs(y = "p"),
    #qplot(y = auto$path_h2_est) + 
      #theme_bigstatsr() + 
      #geom_hline(yintercept = auto$h2_est, col = "blue") +
      #labs(y = "h2"),
    #ncol = 1, align = "hv"
  #)
  # `range` should be between 0 and 2
  (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
  (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))
  
  if(is.null(obj.bigSNP)){
    obj.bigSNP <- snp_attach(gene_data_path) # this should be rds. data
  }
  genotype <- obj.bigSNP$genotypes
  ind.test <- 1:nrow(genotype)
  
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
  beta_auto
  pred_auto
  return(list(pred_auto = pred_auto, beta_auto = beta_auto))
}

pred_auto <- results[["pred_auto"]]
beta_auto <- results[["beta_auto"]]
  
  
  
  # Perform genetic analysis (LD score regression, LDpred2, etc.)
  # This will include matching, computing LD matrices, and so forth
  # Note: Detailed steps should be included here based on your original code
  
  # Compute Polygenic Risk Scores
  # Adjust this part based on how you calculate PRS in your workflow
  
  # Save or return the results
  # You might save files to `outputDir` or return data frames/lists directly

# Example of how to call your function
results <- PRSAnalysis(
  reference_data_path = "/home/shizh/research/tmp/prs/output/hapmap3_snp/hapmap3_snp.RDATA",
  summary_statistics_path = "/home/shizh/research//GIGASTROKE2022/stroke_Any.std",
  gene_data_path = "/home/shizh/research/PRS/rawData/IMPROVE_after_qc.endotype.bed",
  outputDir = "/home/shizh/research/tmp/tmp_dir/"
)

