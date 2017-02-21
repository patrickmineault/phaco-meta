# Let's verify something in the data:
# We only get both absolute values of the IOP and changes in the IOP when we have full follow through
source("read-data.R")
df <- read.data(agg.arms=FALSE, impute.change=FALSE)

df %>% filter(!is.na(df$SixMoAbsIOPChangeMean) & !is.na(df$SixMoIOPMean)) %>% dplyr::select(PreOpEyes, SixMoEyes)
df %>% filter(!is.na(df$OneYAbsIOPChangeMean) & !is.na(df$OneYIOPMean)) %>% dplyr::select(PreOpEyes, OneYEyes)
df %>% filter(!is.na(df$LastPeriodIOPMean) & !is.na(df$LastPeriodAbsIOPChangeMean)) %>% dplyr::select(study.name, subtype, PreOpEyes, LastPeriodEyes, PreOpIOPMean, LastPeriodIOPMean, LastPeriodAbsIOPChangeMean)

# That's true for all studies except 3: 
#                       study.name subtype PreOpEyes LastPeriodEyes PreOpIOPMean LastPeriodIOPMean LastPeriodAbsIOPChangeMean
# 4  Leelachaikul et. al (2005)     OAG        58             54     16.50000          15.20000                  -1.600000
# 5     Mathalone et. al (2005)     OAG        58             24     16.30000          15.10000                  -1.900000
# 6        Shoji et. al  (2007)     OAG        35             20     16.70000          15.60000                  -1.000000

# In the Mathalone study, the changes are for those that are not lost to follow up.
# In the Shoji study as well. 
# Interpretation: changes are measured within the surviving group. 

# What about those that cite only changes?
df %>% filter(!is.na(df$SixMoAbsIOPChangeMean) & is.na(df$SixMoIOPMean)) %>% dplyr::select(study.name, PreOpEyes, SixMoEyes)
df %>% filter(!is.na(df$OneYAbsIOPChangeMean) & is.na(df$OneYIOPMean)) %>% dplyr::select(study.name, PreOpEyes, OneYEyes)

# Relative changes: within those that are not lost to follow up.
# Absolutes: within those that are not lost to follow up.

# Q: how many cases do we have where we have one but not the other?
sum((!is.na(df$SixMoAbsIOPChangeMean) & is.na(df$SixMoIOPMean)) |
    (!is.na(df$OneYAbsIOPChangeMean) & is.na(df$OneYIOPMean)) |
    (!is.na(df$LastPeriodAbsIOPChangeMean) & is.na(df$LastPeriodIOPMean)))
# 7 cases where we have change but not mean.

sum((is.na(df$SixMoAbsIOPChangeMean) & !is.na(df$SixMoIOPMean)) |
    (is.na(df$OneYAbsIOPChangeMean) & !is.na(df$OneYIOPMean)) |
    (is.na(df$LastPeriodAbsIOPChangeMean) & !is.na(df$LastPeriodIOPMean)))
# 44 cases where we have mean but not change

sum((!is.na(df$SixMoAbsIOPChangeMean)) |
      (!is.na(df$OneYAbsIOPChangeMean)) |
      (!is.na(df$LastPeriodAbsIOPChangeMean)))
# 22 cases where we have change.

sum((is.na(df$SixMoAbsIOPChangeMean) & !is.na(df$SixMoIOPMean)) |
      (is.na(df$OneYAbsIOPChangeMean) & !is.na(df$OneYIOPMean)) |
      (is.na(df$LastPeriodAbsIOPChangeMean) & !is.na(df$LastPeriodIOPMean)))

a <- mvrnorm(n = 100, mu = c(0, 0), Sigma = matrix(c(1, .5, .5, 1), nrow=2))
cor(a)

sd(a[,1])
sd(a[,2])
sd(a[,1] - a[,2])
b <- rbinom(n=100, 1, .5)

sd(a[b==1,1])
sd(a[b==1,2])
sd(a[b==1,2] - a[b==1,1])

# If we wipe out some portion of the population completely at random, the correlation structure
# and the sds are preserved. 
estimate.corr <- function(a, b, c) {
  return(- (c ** 2 - a ** 2 - b ** 2) / (2 * a * b))
}
 
cs <- c(with(df, estimate.corr(PreOpIOPStdDev, SixMoIOPStdDev, SixMoAbsIOPChangeStdDev)),
        with(df, estimate.corr(PreOpIOPStdDev, OneYIOPStdDev, OneYAbsIOPChangeStdDev)),
        with(df, estimate.corr(PreOpIOPStdDev, LastPeriodIOPStdDev, LastPeriodAbsIOPChangeStdDev)))
cs <- cs[!is.na(cs)]
median(cs)
mean(cs)

# Therefore, impute changes under a MCAR assumption, with 
# StdDevChange := sqrt(StdDevBefore ** 2 + StdDevAfter ** 2 - rho * StdDevBefore * StdDevAfter)
# Where rho is a conservative .35.

# Simulate missing information.

# library(MASS)
# a <- mvrnorm(n = 100, mu = c(0, 0), Sigma = matrix(c(1, .5, .5, 1), nrow=2))
#plot(a[,1], a[,2])
#cor(a)

#sd(a[,1])
#sd(a[,2])
#sd(a[,1] - a[,2])

#
aa <- c(rep(1, 6), rep(2, 31), rep(3, 7))
mean(aa)
sd(aa)

aa <- c(rep(0, 38), rep(1, 6))
mean(aa)
sd(aa)

