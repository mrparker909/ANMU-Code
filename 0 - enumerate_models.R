library(tidyverse)
lambda_opts <- c("1","NorthCove")
gamma_opts <- c("1","NorthCove","Year")
omega_opts <- c("1","NorthCove","Year")
detec_opts <- c("1","NorthCove","Year")

param_list <- list(lambda_opts, gamma_opts, omega_opts, detec_opts)

make_formula <- function(x) {
  suffix <- "~"
  paste0(suffix,paste(x,collapse="+"))
}

enum_temp <- lapply(param_list, function(w) data.frame(f=unlist(lapply(1:length(w), function(x) combn(w, x, FUN=make_formula)))))
enum <- lapply(enum_temp, function(x) x %>% filter(!grepl("~1+", f, fixed=TRUE)))

# enum contains a list of each covariate combination for each parameter, now we need 
# to make a list of all combinations, 1 choice for each paramter

enum_formulas_temp <- expand.grid(lambda_form=enum[[1]]$f, 
                                  gamma_form=enum[[2]]$f, 
                                  omega_form=enum[[3]]$f, 
                                  detec_form=enum[[4]]$f,
                                  stringsAsFactors = FALSE)
enum_formulas <- enum_formulas_temp %>% mutate(lambda_form=as.character(lambda_form),
                                               gamma_form=as.character(gamma_form),
                                               omega_form=as.character(omega_form),
                                               detec_form=as.character(detec_form),
                                               rownumber=1:n())
saveRDS(object = enum_formulas, file = "./df_enum_formulas.RDS")

