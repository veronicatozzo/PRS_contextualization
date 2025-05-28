suppressPackageStartupMessages({
    library(pROC)
    library(optparse)
    library(readr)
    library(dplyr)
    library(broom)
})

parser <- OptionParser()
parser <- add_option(parser,"--path_to_pgs_file", action = "store", type = "character", default = NULL)
parser <- add_option(parser,"--path_to_phenotype_file", action = "store", type = "character", default = NULL)
parser <- add_option(parser,"--output_file", action = "store", type = "character", default = NULL)
parser <- parse_args(parser)

nagelkerke <- function(model1, model0) {
  ll_null <- logLik(model0)
  ll_model <- logLik(model1)
  n <- length(model1$fitted.values)
  dev_null <- -2 * as.numeric(ll_null)
  dev_model <- -2 * as.numeric(ll_model)
  r2_n <- 1 - (dev_model / dev_null)
  r2_n_adj <- r2_n / (1 - exp(-dev_null / (n / 2)))
  return(r2_n_adj)
}

pgs <- read_tsv(parser$path_to_pgs_file)
data <- read_tsv(parser$path_to_phenotype_file)
data_pgs <- inner_join(data, pgs, by = "ID")

analyze_group <- function(data, ancestry_label, pgs_type) {
  trait <- strsplit(pgs_type, "_")[[1]][1]
  z_col <- paste0(pgs_type, "_Z")
  
  # Standardize the PGS
  data[[z_col]] <- scale(data[[pgs_type]], center = TRUE, scale = TRUE)

  if (trait == "BC") {
    formula_str <- paste0("BC ~ ", z_col, " + BC_age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  } else if (trait == "CHD") {
    formula_str <- paste0("CHD ~ ", z_col, " + CHD_age + self_identified_sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  } else {
    stop(paste("Unsupported trait:", trait))
  }

  model <- glm(as.formula(formula_str), data = data, family = binomial(link = "logit"), na.action = na.exclude)
  
  coef_data <- tidy(model) %>% filter(term == z_col)
  OR <- exp(coef_data$estimate)
  lower_CI <- exp(confint(model)[z_col, 1])
  upper_CI <- exp(confint(model)[z_col, 2])
  p_value <- coef_data$p.value

  return(data.frame(
    pgs_type = pgs_type,
    ancestry = ancestry_label,
    OR = OR,
    lower_CI = lower_CI,
    upper_CI = upper_CI,
    p_value = p_value
  ))
}


ancestry_data <- list(
  EUR = data_pgs %>% filter(ancestry == "EUR"),
  AFR = data_pgs %>% filter(ancestry == "AFR"),
  AMR = data_pgs %>% filter(ancestry == "AMR"),
  EAS = data_pgs %>% filter(ancestry == "EAS"),
  SAS = data_pgs %>% filter(ancestry == "SAS")
)

pgs_types <- c("BC_eMERGE", "BC_PGS000507", "CHD_eMERGE", "CHD_PGS003725")

results <- list()
for (pgs in pgs_types) {
  for (ancestry in names(ancestry_data)) {
    df <- ancestry_data[[ancestry]]
    if (nrow(df) > 0 && pgs %in% names(df)) {
      tryCatch({
        result <- analyze_group(df, ancestry, pgs)
        results <- append(results, list(result))
      }, error = function(e) {
        message(sprintf("Failed for %s in %s: %s", pgs, ancestry, e$message))
      })
    }
  }
}
results_df <- do.call(rbind, results)
write.table(results_df, parser$output_file, sep = '\t', row.names = FALSE)


