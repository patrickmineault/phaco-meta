library(testthat)
library(dplyr)

agg.arms <- function(df, which.arms) {
  # Aggregate different arms according to severity
  df_ <- df[which.arms, ]
  
  stopifnot(nrow(df_) == 2) # Currently only support aggregation over 2 arms.
  row <- df_[1, ]# Use the convention that whatever is not specifically aggregated comes from the first row.
  
  # Weighted means.
  wm <- list()
  wm[['PreOpEyes']] <- c('AgeMean', 'RxPreOpMean', 'PreOpIOPMean', 'VisualAcuityPreOpMean', 'Male', 'Female', 'OAG', 'ACG', 'NTG', 'PXG')
  wm[['SixMoEyes']] <- c('SixMoIOPMean', 'SixMoAbsIOPChangeMean')
  wm[['OneYEyes']] <- c('OneYIOPMean', 'OneYAbsIOPChangeMean')
  wm[['LastPeriodEyes']] <- c('LastPeriodIOPMean', 'LastPeriodAbsIOPChangeMean', 'RxPostOpMean', 'VisualAcuityPostOpMean')
  
  # Weighted standard deviations.
  ws <- list()
  ws[['PreOpEyes']] <- c('AgeStdDev', 'RxPreOpStdDev', 'PreOpIOPStdDev', 'VisualAcuityPreOpStdDev')
  ws[['SixMoEyes']] <- c('SixMoIOPStdDev', 'SixMoAbsIOPChangeStdDev')
  ws[['OneYEyes']] <- c('OneYIOPStdDev', 'OneYAbsIOPChangeStdDev')
  ws[['LastPeriodEyes']] <- c('LastPeriodIOPStdDev', 'LastPeriodAbsIOPChangeStdDev', 'RxPostOpStdDev', 'VisualAcuityPostOpStdDev')
  
  for(i in names(wm)) {
    # Weighted mean: easy.
    row[wm[[i]]] <- sapply(df_[,wm[[i]]], function(x) {sum(x * df_[,i] / sum(df_[,i])) })
    
    # From: http://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
    r <- sapply(1:length(ws[[i]]), function(j) {
      sqrt(
        sum(df_[,i] / sum(df_[,i]) * df_[,ws[[i]][j]] ** 2) +
        sum(df_[,i] / sum(df_[,i]) * df_[,wm[[i]][j]] ** 2) -
        sum(df_[,i] / sum(df_[,i]) * df_[,wm[[i]][j]]) ** 2
      )
      })
    row[ws[[i]]] <- r
    row[i] = sum(df_[,i])
  }
  
  df <- df[!which.arms,]
  df <- rbind(df, row)
  return(df)
}

fill.eyes <- function(df) {
  # Fill in number of eyes.
  bad.info <- (((!is.na(df$SixMoAbsIOPChangeMean) | !is.na(df$SixMoIOPMean)) & is.na(df$SixMoEyes)) | 
               ((!is.na(df$OneYAbsIOPChangeMean) | !is.na(df$OneYIOPMean)) & is.na(df$OneYEyes)) | 
               ((!is.na(df$LastPeriodIOPMean) | !is.na(df$LastPeriodAbsIOPChangeMean)) & is.na(df$LastPeriodEyes)))
  
  # Some eyes being bad and others not is the worst.
  saveable <- bad.info & (!is.na(df$SixMoEyes) | !is.na(df$OneYEyes) | !is.na(df$LastPeriodEyes))
  if(sum(saveable) > 0) {
    cat("This study only has some eye number information missing, fix it:\n")
    cat(df[saveable,]$study.name)
  }
  
  df[, 'prospective'] <- !bad.info
  df[bad.info, 'SixMoEyes'] <- df[bad.info, 'PreOpEyes']
  df[bad.info, 'OneYEyes'] <- df[bad.info, 'PreOpEyes']
  df[bad.info, 'LastPeriodEyes'] <- df[bad.info, 'PreOpEyes']
  return(df)
}

fill.absolutes <- function(df) {
  missing.abs <- (!is.na(df$SixMoAbsIOPChangeMean) & is.na(df$SixMoIOPMean)) | 
                 (!is.na(df$OneYAbsIOPChangeMean) & is.na(df$OneYIOPMean)) |
                 (!is.na(df$LastPeriodAbsIOPChangeMean) & is.na(df$LastPeriodIOPMean))
  df[bad.info, 'SixMoEyes'] <- df[bad.info, 'PreOpEyes']
  df[bad.info, 'OneYEyes'] <- df[bad.info, 'PreOpEyes']
  df[bad.info, 'LastPeriodEyes'] <- df[bad.info, 'PreOpEyes']
  
  return(df)
}

read.data <- function(agg.arms) {
  # Read data from the phaco meta-analysis and fix encodings
  df <- read.csv("phaco.csv", na.strings='-')
  
  df <- df %>% rename(SixMoEyes = Eyes6mo,
                      SixMoIOPMean = MeanIOP6mo,
                      SixMoIOPStdDev = MeanIOPsd6mo,
                      SixMoAbsIOPChangeMean = IOPchangemean6m,
                      SixMoAbsIOPChangeStdDev = IOPchangesd6mo,
                      OneYEyes = Eyes12mo,
                      OneYIOPMean = MeanIOP12mo,
                      OneYIOPStdDev = MeanIOPsd12mo,
                      OneYAbsIOPChangeMean = IOPchangemean12m,
                      OneYAbsIOPChangeStdDev = IOPchangeSD12mo,
                      LastPeriodAbsIOPChangeStdDev = LastPeriodAbsIOPChangeStd,
                      LastPeriodEyes = LastPeriodofEyes
  )
  
  # TODO(Patrick): Aggregate the Euswas study, because it's all the same treatment - 
  # just different severity.
  # Similarly for Lee et. al 2016.
  # Aggregate the two PXG Mierzejewski studies
  # Hayashi et al.: Only
  # patients who were followed for at least 12 months were
  # included in the analyses.
  
  df <- df %>% mutate(subtype = as.factor(
    ifelse(acuteangleclosure == 'Y', 'acute', 
           ifelse(MIGsYorN == 'Y', 'MIGS',
                  ifelse(OAG > 50, 'OAG',
                         ifelse(ACG > 50, 'ACG',
                                ifelse(PXG > 50, 'PXG', NA)))))))
  
  df <- df %>% mutate(study.name = paste0(Author, ' (', Year, ')', ifelse(WashOut == 'Y', '*', '')))
  df <- df %>% dplyr::arrange(Year, study.name)
  
  # Fill in number of eyes. As soon as one has missing missingness information, the study will be tagged as retrospective.
  # TODO(Patrick): read pro - retro information directly from CSV once available.
  df <- fill.eyes(df)
  
  # Fill in absolutes.
  df <- fill.absolutes(df)

  # Aggregate some study arms which correspond to different severities
  if(agg.arms) {
    nrow.before <- nrow(df)
    df <- agg.arms(df, regexpr('Mierzejewski', df$study.name) > 0 & df$subtype == 'PXG') 
    df <- agg.arms(df, regexpr('Lee', df$study.name) > 0 & df$Year == '2016')
    df <- agg.arms(df, regexpr('Euswas', df$study.name) > 0)
    stopifnot(nrow.before - nrow(df) == 3)
  }

  # TODO(Patrick): Right now, Vold is classified as prospective, because we don't have all the eye numbers. 
  # That's inaccurate. Resolve that by going through another pass in the data.
  # stopifnot(df[regexpr('Vold', df$study.name) > 0, 'prospective'])
  return(df)
}

# Unit test agg.arms.
df <- read.data(FALSE)
expect_equal(nrow(df) - 1, nrow())
tgt <- regexpr('Euswas', df$study.name) > 0

df[tgt,'PreOpEyes'] <- 10
df[tgt,'PreOpIOPMean'] <- 0
df[tgt,'PreOpIOPStdDev'] <- 1
df_ <- agg.arms(df, tgt)
tgt_ <- regexpr('Euswas', df_$study.name) > 0
row <- df_[tgt_,]
expect_equal(row$PreOpEyes, 20)
expect_equal(row$PreOpIOPMean, 0)
expect_equal(row$PreOpIOPStdDev, 1)

df[tgt,'PreOpIOPMean'] <- c(0, 10)
df[tgt,'PreOpIOPStdDev'] <- 0.001
df_ <- agg.arms(df, tgt)
tgt_ <- regexpr('Euswas', df_$study.name) > 0
row <- df_[tgt_,]
expect_equal(row$PreOpEyes, 20)
expect_equal(row$PreOpIOPMean, 5, tolerance=1e-5)
expect_equal(row$PreOpIOPStdDev, 5, tolerance=1e-6)

