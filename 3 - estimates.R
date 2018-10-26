# The purpose of this script is to generate breeding pair abundance estimates using productivity and area expansion

library(unmarked)
library(tidyverse)

our_abundance_estimates <- function( productivity1995, productivity2006 ) {
  
  FUNNEL_AREA <- 34189.9+17906.9
  
  # read in the previously calculated ranef
  r <- readRDS("./output/mod_murr_54_ranef.RDS")
  a <- bup(r)
  
  # group estimates by North and Cabin cove
  chickEstimateN <- round(colSums(a[1:4,]))
  chickEstimateC <- round(colSums(a[5:6,]))
  
  # 95% CIs
  chickConfN <- colSums(confint(r, level = 0.95)[1:4,,])
  chickConfC <- colSums(confint(r, level = 0.95)[5:6,,])
  
  df_chick <- data.frame(Abundance_Estimate=c(chickEstimateN, chickEstimateC), 
                         Year=1990:2006, 
                         NorthCove=rep(c("North Cove","Cabin Cove"),each=17),
                         Area=rep(c(34189.9,17906.9),each=17),
                         Lower2.5=c(chickConfN[1,],chickConfC[1,]), 
                         Upper97.5=c(chickConfN[2,],chickConfC[2,]))
  
  df_productivity <- data.frame(Year=c(1995, 2006), Prod=c(productivity1995, productivity2006), ColonyArea=c(137600,125500))
  
  # breeding pair abundance estimates
  pop_est <- df_chick %>% filter(Year %in% c(1995,2006)) %>%
    left_join(df_productivity, by="Year") %>%
    mutate(
      Abundance = round(Abundance_Estimate/Prod/Area*(Area/FUNNEL_AREA)*ColonyArea),
      Lower2.5 = round(Lower2.5/Prod/Area*(Area/FUNNEL_AREA)*ColonyArea),
      Upper97.5 = round(Upper97.5/Prod/Area*(Area/FUNNEL_AREA)*ColonyArea)
    ) %>%
    select(Year, Abundance, Lower2.5, Upper97.5) %>%
    group_by(Year) %>% 
    summarise_all(.funs = sum)  
  return(pop_est)
}

# calculate abundance estimates given productivity for 1995 and 2006
our_abundance_estimates(1.49, 1.1)
our_abundance_estimates(1.54, 1.54)


