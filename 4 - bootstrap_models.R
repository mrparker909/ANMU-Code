# input bootstrap iteration from command line
args <- commandArgs(TRUE)

# counts which bootstrap iteration we are currently on
counter <- as.integer(args[1])

# our best model:
current_model <- 54

# install/load libraries
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

## GENERATE POPULATION ##
# note that the covariates are the fitted parameters from our best fitted model

# site covariates
lambda_cov <- c(rep(159.49,times=4), rep(262.17,times=2))
pdet_cov   <- c(rep(1.00,times=4),rep(0.63,times=2))

# time covariates
omega_cov  <- c(0.64, 1.00, 0.88, 1.00, 0.77, 1.00, 0.88, 0.94, 0.95, 1.00, 0.94, 0.97, 0.90, 0.87, 1.00, 0.76)
gamma_cov  <- c(0.00,27.11, 9.97, 3.67, 0.00,27.12, 0.00, 0.00,18.17, 9.97, 0.00, 6.69, 3.67, 0.00, 0.00,22.20)

sites <- 6
times <- 17

# bootstrap population generating function
gen_pop <- function(num_sites,num_times,l_cov,p_cov,o_cov,g_cov) {
  U    = num_sites
  T    = num_times
  
  # generate initial population for each site from poisson distribution
  Ntemp <- c(rep(rpois(n=U,lambda = l_cov),times=T))
  Ni <- matrix(data=Ntemp, nrow = U, ncol = T)
  
  # extrapolate population for T>1 using binomial survival and poisson recruitment
  for(i in 2:T) {
    Ni[,i] <- rbinom(n = U, size = Ni[,i-1], prob = o_cov[i-1]) + rpois(n = U, lambda = g_cov[i-1])
  }
  
  nit <- Ni
  
  # perform binomial thinning to convert from true population to observed counts 
  for(site in 1:U) {
    nit[site,] <- rbinom(n = T, size = Ni[site,], prob = p_cov[site])
  }
  
  return(list(Nit=Ni, nit=nit))
}

# genY is one set of bootstrap counts
genY <- gen_pop(num_sites = sites, num_times = times, l_cov = lambda_cov, p_cov = pdet_cov, o_cov = omega_cov, g_cov = gamma_cov)


## FIT NMIXTURE MODEL TO BOOTSTRAP SAMPLE ##
chick_dat <- genY$nit
yeardat <- data.frame(Year=as.factor(rep(1990:2006,6)))

# counts
Y <- chick_dat
# site covariates
X <- data.frame(Funnel=as.factor(1:6), NorthCove=as.factor(c(1,1,1,1,0,0)))
# temporal covariates
Z <- yeardat                        

uframeO <- unmarkedFramePCO(y=Y, siteCovs=X, yearlySiteCovs=Z, numPrimary = 17)


defK <- max(300,ceiling(2*max(genY$nit)))

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
    mod <- pcountOpen(as.formula(enum_formulas[current_model,]$lambda_form), 
                      as.formula(enum_formulas[current_model,]$gamma_form), 
                      as.formula(enum_formulas[current_model,]$omega_form), 
                      as.formula(enum_formulas[current_model,]$detec_form), 
                      uframeO, 
                      K=defK,
                      control=list(maxit=1000, trace=TRUE, REPORT=1),
                      se = TRUE,
                      starts=c(5.569,
                               -0.497,
                               -18.1,
                               21.4,
                               20.4,
                               19.4,
                               -14.1,
                               21.4,
                               -34.1,
                               -20.3,
                               21.0,
                               20.4,
                               -11.9,
                               20.0,
                               19.4,
                               -11.0,
                               -10.7,
                               21.2,
                               0.596,
                               16.488,
                               1.371,
                               26.598,
                               0.599,
                               23.709,
                               1.412,
                               2.203,
                               2.289,
                               44.336,
                               2.168,
                               2.945,
                               1.592,
                               1.273,
                               65.946,
                               0.572,
                               0.537,
                               0.513))
    saveRDS(mod, file = paste0("./output/bootstrap/mod_boot_",counter,".RDS"))
    
    mod_out <- list(AIC=uAIC(mod),
                    BIC=uBIC(mod),
                    lambda_f=mod@formlist[1], 
                    gamma_f=mod@formlist[2], 
                    omega_f=mod@formlist[3], 
                    p_f=mod@formlist[4], 
                    loglikelihood=unmarked::logLik(mod),
                    filename=paste0("mod_boot_",counter,".RDS"),
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
                    filename=paste0("mod_boot_",counter,".RDS"),
                    warnings=list(e),
                    optim_converged=NA,
                    optim_message=NA
    )
    mod_out
  }
)


df <- model_output
saveRDS(df, file = paste0("./output/bootstrap/df_mod_boot_",counter,".RDS"))

