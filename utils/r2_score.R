suppressPackageStartupMessages({
    # library(rcompanion)
    library(pROC)
    library(optparse)
})

parser <- OptionParser()
parser <- add_option(parser,"--path_to_pgs_file", action = "store", type = "character", default = NULL)
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


all_data = c()
for (PGS in c('BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725')){
    print(PGS)
    data = read.csv(paste0(parser$path_to_pgs_file, "/", PGS, ".csv"))
    data$self_identified_sex = as.factor(data$self_identified_sex)
    data$all = 1
    formula = "phenotype~age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
    if (substring(PGS, 1, 2)!='BC'){
        formula = paste0(formula, '+ bmi')
    }
                
    for (prs_type in c('PRS_PC', 'PRS_PC_sex', 'PRS_PC_BMI', 'PRS_PC_age', 'PRS_PC_all')){
            for (context in c('all', 'Age_binned', 'BMI_binned', 'self_identified_sex', 'ancestry')){
                for (class in unique(data[[context]])){
                    if (context=='nan'){
                        next
                    }
                    if (context !='self_identified_sex' & substring(PGS, 1, 2)!='BC'){
                        
                        formula_ = paste0(formula, '+self_identified_sex')
                     }
                    else{
                        formula_=formula
                        }
                    data_context = data[data[[context]]==class,]
                    model1 <- glm(paste0(formula_, '+', prs_type),
                                 data = data_context,family='binomial')
                    model0 <-  glm(formula_,
                                 data = data_context,family='binomial')
                    r2 <- nagelkerke(model1, model0)
                    roc_context <- roc(data_context$phenotype, predict(model1))
                    roc_context0 <- roc(data_context$phenotype, predict(model0))
                    
                    df1 <- data.frame(PGS=PGS, 
                                      prs_type=prs_type,
                                      context=context,
                                      context_group=class,
                                      r2=r2, 
                                      aux=auc(roc_context)-auc(roc_context0))
                    all_data = rbind(all_data, df1)
               }
            }
       }
        
 }
write.table(all_data, parser$output_file)
 