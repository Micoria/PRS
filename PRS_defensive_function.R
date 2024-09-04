rm(list=ls())
defensive_PRS <- function(path_reference,path_GWAS,path_genotype,path_temp){
  #load the required packages
  require(dplyr)
  require(bigreadr)
  require(bigsnpr)
  require(data.table)
  require(stringr)
  require(R.utils)
  require(runonce)
  require(bigstatsr)
  # load hapmap 3+ in the environment
  load(path_reference)
  # Check if the specific object 'hapmap3' exists in the environment
  if (exists("info")) {
  } else {
    # If 'hapmap3' doesn't exist, download and load the alternative dataset
    hapmap3 <- readRDS(runonce::download_file(
      "https://figshare.com/ndownloader/files/37802721",
      dir = "tmp-data", fname = "map_hm3_plus.rds"))
  }
  #change the pattern and the names of the datasets
  variable_patterns <- list(
    chr = c("chrom", "chromosome", "chromo","Chr"),
    rsid = c("rs", "rsID", "RS"),
    pos = c("position", "Position"),
    a0 = c("a1", "referentalleles","allele1"),
    a1 = c("a2", "riskallelels", "effectalleles","allele2"),
    beta = c("b"),
    beta_se = c("se"),
    n_eff = c("N", "ncol")
  )
  rename_variables <- function(df, patterns) {
    for (standard_name in names(patterns)) {
      pattern_list <- patterns[[standard_name]]
      for (pattern in pattern_list) {
        if (any(grepl(pattern, names(df), ignore.case = TRUE))) {
          df <- df %>%
            rename_with(~ standard_name, matches(pattern, ignore.case = TRUE))
        }
      }
    }
    return(df)
  }
  pattern <- function(df, patterns) {
    # Apply the renaming function
    df <- rename_variables(df, patterns)
    # Check if all standard variables are present
    missing_vars <- setdiff(names(patterns), names(df))
    if (length(missing_vars) > 0) {
      stop(paste("Missing variables after renaming:", paste(missing_vars, collapse = ", ")))
    }
    return(df)
  }
  
  # Load GWAS data
  sumstats <- fread2(path_GWAS)
  sumstats <- pattern(sumstats, variable_patterns)
  sumstats$a0 <- toupper(sumstats$a0)
  sumstats$a1 <- toupper(sumstats$a1)
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
  #load genotype data
  process_genotype_files <- function(file_path) {
    #get the extension name of the document
    file_extension <- tools::file_ext(file_path)
    if (file_extension == "rds") {
      obj.bigSNP <- snp_attach(file_path)
    } else if (file_extension %in% c("bed", "bim", "fam")) {
      snp_readBed(file_path)
      # hypotheses they have the same path and same name
      rds_file_path <- sub("\\.bed$", ".rds", file_path)
      obj.bigSNP <- snp_attach(rds_file_path)
    } else {
      stop("unsuportted file")
    }
    return(obj.bigSNP)
  }
  obj.bigSNP <-process_genotype_files(path_genotype)
  ncores = 11
  available_cores <- parallel::detectCores()
  if (ncores > available_cores) {
    stop(paste("Requested", ncores, "cores, but only", available_cores, "are available."))
  }
  
  # Apply a function over all columns of the FBM to check for any NA values
  has_na <- big_apply(
    X = obj.bigSNP$genotypes,
    a.FUN = function(X, ind) {
      anyNA(X[, ind])
    },
    ind = cols_along(obj.bigSNP$genotypes),  # Apply across all columns
    ncores = 1  # Adjust the number of cores depending on your system
  )
  
  # Check if any column contains NA
  if (any(has_na)) {
    print("The genotypes dataset contains NA values.")
  } else {
    print("The genotypes dataset does not contain NA values.")
  }
  obj.bigSNP$genotypes <- snp_fastImputeSimple(Gna = obj.bigSNP$genotypes,method = "mean2",ncores = 11)
  # Get aliases for useful slots
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  (NCORES <- nb_cores())
  
  withCallingHandlers(
    expr = { 
      map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
    }, 
    error = function(e) {
      print("The names are not followed the requirments, check the variables' names.")
    },
    warning = function(e) {
      print("something does not look right")
    },
    message = function(e) {
      print("You should make sure results are correct.")
    }
  )
  gc()
  withCallingHandlers(
    expr = {df_beta <- snp_match(sumstats, map)
    }, 
    error = function(e) {
      print("Might use different genomics，check the format of rsid in GWAS data")
    },
    warning = function(e) {
      print("something does not look right")
    },
    message = function(e) {
      print("You should make sure results are correct.")
    }
  )
  genotype <- obj.bigSNP$genotypes
  ind.row <- rows_along(genotype)
  maf <- snp_MAF(genotype, ind.row = ind.row, ind.col = df_beta$"_NUM_ID_", ncores = NCORES)
  maf_thr <- 0.05 # threshold privefl likes to use
  df_beta <- df_beta[maf > maf_thr, ]
  
  # We assume the fam order is the same across different chromosomes
  fam.order <- as.data.frame(obj.bigSNP$fam)
  tryCatch(
    expr = {
      fam.order <- fam.order %>% 
        rename("FID" = "family.ID",
               "IID" = "sample.ID")
    }, 
    error = function(e) {
      print("names are not matched")
    },
    warning = function(e) {
      print("something does not look right")
    },
    message = function(e) {
      print("You should make sure results are correct.")
    }
  )
  rm(info)
  rm(has_na)
  rm(map)
  rm(obj.bigSNP)
  rm(sumstats)
  gc()
  POS2 <- snp_asGeneticPos(CHR, POS, dir =path_temp)
  corr <- NULL
  ld <- NULL
  on.exit(file.remove(paste0(path_temp, ".sbk")), add = TRUE)
  gc()
  
  sbk_file <- paste0(path_temp, ".sbk")
  cat("Attempting to delete file:", sbk_file, "\n")
  
  # Check if the file exists
  if (file.exists(sbk_file)) {
    cat("File exists. Attempting to delete...\n")
    success <- file.remove(sbk_file)
    if (success) {
      cat("File successfully deleted.\n")
    } else {
      cat("Failed to delete the file.\n")
    }
  } else {
    cat("File does not exist.\n")
  }
  # Defensive loop to process each chromosome
  for (chr in 1:22) {
    
    # Check if df_beta contains the required columns
    if (!all(c("chr", "_NUM_ID_") %in% names(df_beta))) {
      stop("The data frame 'df_beta' must contain 'chr' and '_NUM_ID_' columns.")
    }
    
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(df_beta$chr == chr)
    
    if (length(ind.chr) == 0) {
      warning(paste("No SNPs found for chromosome", chr, ". Skipping this chromosome."))
      next
    }
    
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    
    # Ensure POS2 has enough positions
    if (any(ind.chr2 > length(POS2))) {
      stop("Indexing error: Some indices in 'ind.chr2' exceed the length of 'POS2'.")
    }
    
    # Calculate the LD
    corr0 <- tryCatch({
      snp_cor(
        G,
        ind.col = ind.chr2,
        ncores = 11,
        infos.pos = POS2[ind.chr2],
        size = 3 / 1000
      )
    }, error = function(e) {
      stop(paste("Error in calculating LD for chromosome", chr, ":", e$message))
    })
    
    if (chr == 1) {
      ld <- tryCatch({
        Matrix::colSums(corr0^2)
      }, error = function(e) {
        stop("Error in calculating 'ld' for the first chromosome: ", e$message)
      })
      
      corr <- tryCatch({
        as_SFBM(corr0, path_temp)
      }, error = function(e) {
        stop("Error in converting 'corr0' to SFBM for the first chromosome: ", e$message)
      })
      
    } else {
      ld <- tryCatch({
        c(ld, Matrix::colSums(corr0^2))
      }, error = function(e) {
        stop("Error in appending 'ld' for chromosome", chr, ":", e$message)
      })
      
      tryCatch({
        corr$add_columns(corr0, nrow(corr))
      }, error = function(e) {
        stop("Error in adding columns to 'corr' for chromosome", chr, ":", e$message)
      })
    }
  }

  ldsc <- snp_ldsc( ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
  ldsc_h2_est <- ldsc[["h2"]]
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = ldsc_h2_est)
  coef_shrink <- 0.95  
  set.seed(1)  # to get the same result every time
  # takes less than 2 min with 4 cores
  multi_auto <- snp_ldpred2_auto(
    corr, df_beta, h2_init = ldsc_h2_est,
    vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = NCORES,
    # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
    allow_jump_sign = FALSE, shrink_corr = coef_shrink)
  
  (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
  (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))
  index_for_prs_calculation <- 1:nrow(G)
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  pred_auto <- big_prodVec(G, 
                           beta_auto, 
                           ind.row = index_for_prs_calculation, 
                           ind.col = df_beta[["_NUM_ID_"]])
  PRS <- data.frame("IID" = fam.order$IID,
                    "prs" = pred_auto)
  return(PRS)
}
rm(PRS)
PRS <- defensive_PRS(
  path_reference = "/home/XXXXX.RDATA",
  path_GWAS = "/home/XXXXXXX.std",
  path_genotype = "/home/XXXXXXX.rds",
  path_temp = "/home/XXXXX"
)

fwrite(PRS, file = "/home/XXXXXX.tsv",
       sep = "\t",quote = F)
