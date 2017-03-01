library(testthat)
library(dplyr)

agg.arms <- function(df, which.arms) {
  # Aggregate different arms of the same study.
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
  # Fill in number of eyes at each time period where information is missing.
  bad.info <- (((!is.na(df$SixMoAbsIOPChangeMean) | !is.na(df$SixMoIOPMean)) & is.na(df$SixMoEyes)) | 
               ((!is.na(df$OneYAbsIOPChangeMean) | !is.na(df$OneYIOPMean)) & is.na(df$OneYEyes)) | 
               ((!is.na(df$LastPeriodIOPMean) | !is.na(df$LastPeriodAbsIOPChangeMean)) & is.na(df$LastPeriodEyes)))
  
  # Some eyes being bad and others not is the worst.
  saveable <- bad.info & (!is.na(df$SixMoEyes) | !is.na(df$OneYEyes) | !is.na(df$LastPeriodEyes))
  if(sum(saveable) > 0) {
    message("This study only has some eye number information missing, fix it:")
    message(df[saveable,]$study.name)
  }
  
  # Compare whether a study has bad eye info data to prospectiveness
  bad.studies <- df[(df$StudyType == 'Prospective') & bad.info]$study.name
  if(!is.null(bad.studies)) {
    message("These prospective studies are missing eye number information:")
    message(bad.studies)
  }

  df.weird <- df %>% filter(StudyType == 'Retrospective', 
                (SixMoEyes < PreOpEyes) | (OneYEyes < PreOpEyes) | (LastPeriodEyes < PreOpEyes))
  weird <- unique(df.weird$study.name)
  if(length(weird) > 0) {
    message("These retrospective studies are losing eyes per period - not impossible, but unusual:")
    message(paste(weird, '\n'))
  }
  
  df.weird <- df %>% filter((SixMoEyes > PreOpEyes) | (OneYEyes > PreOpEyes) | (LastPeriodEyes > PreOpEyes))
  weird <- unique(df.weird$study.name)
  if(length(weird) > 0) {
    message("These retrospective studies are gaining eyes as the study goes")
    message(paste(weird, '\n'))
  }
  
  df[bad.info, 'SixMoEyes'] <- df[bad.info, 'PreOpEyes']
  df[bad.info, 'OneYEyes'] <- df[bad.info, 'PreOpEyes']
  df[bad.info, 'LastPeriodEyes'] <- df[bad.info, 'PreOpEyes']
  
  # Set LastPeriodEyes to whatever the the latest reading is.
  bad.info <- is.na(df$LastPeriodEyes)
  df[bad.info, 'LastPeriodEyes'] <- df[bad.info, 'OneYEyes']
  bad.info <- is.na(df$LastPeriodEyes)
  df[bad.info, 'LastPeriodEyes'] <- df[bad.info, 'SixMoEyes']
  
  
  return(df)
}

fill.change <- function(df, rho) {
  # Fill in information of absolute IOP numbers and relative numbers.
  missing.abs <- (is.na(df$SixMoAbsIOPChangeMean) & !is.na(df$SixMoIOPMean)) | 
                 (is.na(df$OneYAbsIOPChangeMean) & !is.na(df$OneYIOPMean)) |
                 (is.na(df$LastPeriodAbsIOPChangeMean) & !is.na(df$LastPeriodIOPMean))
  
  df_ <- df[missing.abs,]
  df_ <- df_ %>% mutate(SixMoAbsIOPChangeMean = SixMoIOPMean - PreOpIOPMean,
                 OneYAbsIOPChangeMean = OneYIOPMean - PreOpIOPMean,
                 LastPeriodAbsIOPChangeMean = LastPeriodIOPMean - PreOpIOPMean,
                 SixMoAbsIOPChangeStdDev = sqrt(PreOpIOPStdDev ** 2 + PreOpIOPStdDev ** 2 - 2 * rho * SixMoIOPStdDev * PreOpIOPStdDev),
                 OneYAbsIOPChangeStdDev = sqrt(PreOpIOPStdDev ** 2 + OneYIOPStdDev ** 2 - 2 * rho * OneYIOPStdDev * PreOpIOPStdDev),
                 LastPeriodAbsIOPChangeStdDev = sqrt(PreOpIOPStdDev ** 2 + 
                                                       LastPeriodIOPStdDev ** 2 - 
                                                       2 * rho * LastPeriodIOPStdDev * PreOpIOPStdDev))
  df[missing.abs,] <- df_
  
  df <- df %>% mutate(RxChangeMean = RxPostOpMean - RxPreOpMean,
                      RxChangeStdDev = sqrt(RxPostOpStdDev ** 2 + RxPreOpStdDev ** 2))
  return(df)
}

read.data <- function(agg.arms=TRUE, impute.change=TRUE, impute.eyes=TRUE) {
  # Read data from the phaco meta-analysis and fix encodings
  fileName <- 'phaco.csv'
  read.file <- readChar(fileName, file.info(fileName)$size)
  # We've had problems before with \003 popping up in the stream out of nowhere.
  stopifnot(length(grep('\003', read.file)) == 0)
  
  df <- read.csv(fileName, na.strings='-')
  colnames(df) <- sub('X','', gsub('\\.','', colnames(df)))
  
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
  
  df <- df %>% mutate(subtype = as.factor(
    ifelse(acuteangleclosure == 'Y', 'acute', 
           ifelse(MIGsYorN == 'Y' | OAG > 50, 'OAG',
                         ifelse(ACG > 50, 'ACG',
                                ifelse(PXG > 50, 'PXG', NA))))))
  
  # See investigate-washout.Rmd for an explanation of these labels and a 
  # verification that the labels are correct.
  df <- df %>% mutate(washout.type = 
                        ifelse(is.na(Washoutbaseline) & is.na(WashoutIOP), 'None',
                          ifelse(Washoutbaseline == PreOpIOPMean & is.na(WashoutIOP), 'Pre',
                            ifelse(!is.na(Washoutbaseline) & !is.na(WashoutIOP), 'Both', 'Partial'))))
  
  washout.labels <- c(None = '',
                         Pre = '**',
                         Both = '*',
                         Partial = '')
  
  # To fix Pfeiffer and any other similar studies.
  switcheroo <- with(df, washout.type == 'Both' & Washoutbaseline > PreOpIOPMean)
  expect_equal(sum(switcheroo), 2)
  df[switcheroo, 'PreOpIOPMean'] <- df[switcheroo, 'Washoutbaseline']
  df[switcheroo, 'PreOpIOPStdDev'] <- df[switcheroo, 'BaselinewashoutSD']
  
  df <- df %>% mutate(study.name = paste0(Author, ' (', Year, ')', washout.labels[washout.type]))
  
  # Fix a weird bug that happened recently.
  if(any(df$OneYAbsIOPChangeStdDev < 0, na.rm=TRUE)) {
    message("Fix the negative SD for this study:")
    message(paste(unique(df[!is.na(df$OneYAbsIOPChangeStdDev) &  df$OneYAbsIOPChangeStdDev < 0,]$study.name), '\n'))
    df <- df %>% mutate(OneYAbsIOPChangeStdDev = abs(OneYAbsIOPChangeStdDev))
  }
  
  # Fill in number of eyes. As soon as one has missing missingness information, the study will be tagged as retrospective.
  # TODO(Patrick): read pro - retro information directly from CSV once available.
  if(impute.eyes) {
    df <- fill.eyes(df)
  }

  if(impute.change) {
    df <- fill.change(df, rho = 0.35)
  }

  # Aggregate some study arms which correspond to different severities
  if(agg.arms) {
    nrow.before <- nrow(df)
    df <- agg.arms(df, regexpr('Mierzejewski', df$study.name) > 0 & df$subtype == 'PXG') 
    df <- agg.arms(df, regexpr('Lee', df$study.name) > 0 & df$Year == '2016')
    df <- agg.arms(df, regexpr('Euswas', df$study.name) > 0)
    stopifnot(nrow.before - nrow(df) == 3)
  }
  
  df <- df %>% dplyr::arrange(Year, study.name)

  # TODO(Patrick): Right now, Vold is classified as prospective, because we don't have all the eye numbers. 
  # That's inaccurate. Resolve that by going through another pass in the data.
  # stopifnot(df[regexpr('Vold', df$study.name) > 0, 'prospective'])
  return(df)
}

filter.data <- function(df, level="all") {
  if(level == 'all') {
    return(df)
  } else if(level == 'prospective') {
    return(df %>% filter(StudyType == 'Prospective'))
  } else if(level == 'prospective-nowashout') {
    return(df %>% filter(StudyType == 'Prospective', washout.type %in% c('None', 'Partial')))
  } else if(level == 'nowashout') {
    return(df %>% filter(washout.type %in% c('None', 'Partial')))
  } else if(level == 'low-low') {
    return(df %>% filter(StudyType == 'Prospective', 
                         washout.type %in% c('None', 'Partial'), 
                         LastPeriodEyes >= .85 * PreOpEyes))
  } else {
    stop("Unknown filtering level.")
  }
}

# Unit test agg.arms, without polluting the global namespace.
test.fun <- function() {
  df <- read.data(FALSE)
  expect_equal(nrow(df) - 1, nrow(agg.arms(df, regexpr('Euswas', df$study.name) > 0)))
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
  
  df <- read.data()
  expect_equal(nrow(df %>% filter(washout.type != 'None', 
                                  (Washoutbaseline < PreOpIOPMean) | 
                                    (WashoutIOP < OneYIOPMean))), 0)
}

test.fun()


# TODO: each of these should conditions should not have any studies associated with them.
# df <- read.data()
# df[!is.na(df$SixMoAbsIOPChangeMean) & is.na(df$SixMoAbsIOPChangeStdDev),]$study.name
# df[!is.na(df$OneYAbsIOPChangeMean) & is.na(df$OneYAbsIOPChangeStdDev),]$study.name
# df[!is.na(df$LastPeriodAbsIOPChangeMean) & is.na(df$LastPeriodAbsIOPChangeStdDev),]$study.name
