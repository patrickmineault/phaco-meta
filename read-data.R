read.data <- function() {
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
  df <- df %>% filter(Author != "Mierzejewski et. al")  # Exclude Mierzejewski studies, they're abstracts and weirdly coded.
  df <- df %>% dplyr::arrange(Year, study.name)
  df
}