# model number passed in from command line (number refers to a line in the "df_enum_formulas.RDS" file)
# example - line 50 from df_enum_formulas.RDS: 

# lambda_form   gamma_form  omega_form  detec_form  rownumber
# ~NorthCove    ~1          ~Year       ~NorthCove  50

args <- commandArgs(TRUE)

# counter indicates which model number
counter <- as.integer(args[1])

# install/load libraries tidyverse and unmarked
if(!require("tidyverse")) install.packages("tidyverse", repos='http://cran.us.r-project.org', dependencies=TRUE); library(tidyverse)
if(!require("unmarked")) install.packages("unmarked", repos='http://cran.us.r-project.org', dependencies=TRUE); library(unmarked)

# define AIC function
uAIC <- function(pcountOpen_model) {
  s <- pcountOpen_model
  L <- unmarked::logLik(pcountOpen_model)
  k <- length(unmarked::coef(s))
  return(2*k-2*L)
}

# define BIC function
uBIC <- function(pcountOpen_model) {
  s <- pcountOpen_model
  L <- unmarked::logLik(pcountOpen_model)
  k <- length(unmarked::coef(s))
  #  n <- nrow(s@data@y)*ncol(s@data@y) # site and time replicates
  n <- nrow(s@data@y) # site only replicates 
  return(log(n)*k-2*L)
}

# load the count data
chick_dat <- read.csv("chicks_F1toF6_1990to2006.csv")
yeardat   <- data.frame(Year=as.factor(rep(1990:2006,6)))

# counts
Y <- chick_dat[,2:18]
# site covariates
X <- data.frame(Funnel=as.factor(chick_dat[,1]), NorthCove=as.factor(chick_dat[,19]))
# temporal covariates
Z <- yeardat                        

uframeO <- unmarkedFramePCO(y=Y, siteCovs=X, yearlySiteCovs=Z, numPrimary = 17)

# define upper bound on summations
defK <- 300

df <- data.frame(matrix(ncol=11, nrow=1))
colnames(df) <- c("AIC",
                  "BIC",
                  "lambda_f", 
                  "gamma_f", 
                  "omega_f", 
                  "p_f", 
                  "loglikelihood",
                  "filename",
                  "warnings",
                  "optim_converged",
                  "optim_message")

enum_formulas <- readRDS(file = "df_enum_formulas.RDS")

model_output <- tryCatch(
  {
    # fit model using unmarked::pcountOpen
    mod <- pcountOpen(as.formula(enum_formulas[counter,]$lambda_form), 
                      as.formula(enum_formulas[counter,]$gamma_form), 
                      as.formula(enum_formulas[counter,]$omega_form), 
                      as.formula(enum_formulas[counter,]$detec_form), 
                      uframeO, 
                      K=defK,
                      control=list(maxit=1000, trace=TRUE, REPORT=1),
                      se = TRUE)
    # save output model
    saveRDS(mod, file = paste0("./output/mod_murr_",counter,".RDS"))
    
    mod_out <- list(AIC=uAIC(mod),
                    BIC=uBIC(mod),
                    lambda_f=mod@formlist[1], 
                    gamma_f=mod@formlist[2], 
                    omega_f=mod@formlist[3], 
                    p_f=mod@formlist[4], 
                    loglikelihood=unmarked::logLik(mod),
                    filename=paste0("mod_murr_",counter,".RDS"),
                    warnings=list(mod@opt$message),
                    optim_converged=list(mod@opt$convergence),
                    optim_message=list(mod@opt$message)
    )
    mod_out
    
  },
  error = function(e){
    mod_out <- list(AIC=NA,
                    BIC=NA,
                    lambda_f=NA, 
                    gamma_f=NA, 
                    omega_f=NA, 
                    p_f=NA, 
                    loglikelihood=NA,
                    filename=paste0("mod_murr_",counter,".RDS"),
                    warnings=list(e),
                    optim_converged=NA,
                    optim_message=NA
    )
    mod_out
  }
)

# save output data
df <- model_output
saveRDS(df, file = paste0("./output/df_mod_murr_",counter,".RDS"))

