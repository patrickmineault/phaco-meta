library(testthat)
library(dplyr)
library(stringr)

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
  
  # Set LastPeriodEyes to whatever the latest reading is.
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

read.data <- function(agg.arms=TRUE, 
                      impute.change=TRUE, 
                      impute.eyes=TRUE, 
                      drop.migs=TRUE,
                      fill.last=TRUE,
                      exclude.short=TRUE,
                      no.processing=FALSE) {
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
  
  df <- df %>% mutate(JournalVolumePage =  
                ifelse(regexpr("O,", JournalVolumePage) == 1, 
                       paste0("Ophthalmol,", substr(JournalVolumePage, 3, 100)),
                       as.character(JournalVolumePage)),
                TimeofLastPostOp = as.character(TimeofLastPostOp))
  
  df <- df %>% mutate(subtype = as.factor(
    ifelse(acuteangleclosure == 'Y', 'AACG', 
           ifelse(MIGsYorN == 'Y' | OAG > 50, 'OAG',
                         ifelse(ACG > 50, 'ACG',
                                ifelse(PXG > 50, 'PXG', NA))))))
  
  # For Jacobi (2002), the error bars for change in IOP imply an r > 1. 
  # Fix this by stripping out that data.
  jacobi <- with(df, Year == '2002' & regexpr('Jacobi', Author))
  expect_equal(sum(jacobi), 1)
  df[jacobi, c('LastPeriodAbsIOPChangeMean', 'LastPeriodAbsIOPChangeStdDev')] <- NA
  
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
  
  # To fix Pfeiffer and any other similar studies with partial washout.
  switcheroo <- with(df, washout.type == 'Both' & Washoutbaseline > PreOpIOPMean)
  expect_equal(sum(switcheroo), 2)
  df[switcheroo, 'PreOpIOPMean'] <- df[switcheroo, 'Washoutbaseline']
  df[switcheroo, 'PreOpIOPStdDev'] <- df[switcheroo, 'BaselinewashoutSD']
  
  # Compile study names with dagger and such.
  df <- df %>% 
    mutate(hayashi.symbol = ifelse(regexpr("Hayashi", df$Author) > 0, "†", ""),
           study.name = paste0(Author, ' (', Year, ')', washout.labels[washout.type], hayashi.symbol))
  
  if(no.processing) {
    return(df)
  }
  
  
  # Prevent a weird bug about negative SDs from reappearing.
  stopifnot(any(df$OneYAbsIOPChangeStdDev < 0, na.rm=TRUE) == FALSE)
  
  # Fill in number of eyes. 
  if(impute.eyes) {
    df <- fill.eyes(df)
  }

  if(impute.change) {
      df <- fill.change(df, rho = 0.28)
  }

  # Aggregate some study arms which correspond to different severities
  if(agg.arms) {
    nrow.before <- nrow(df)
    df <- agg.arms(df, regexpr('Mierzejewski', df$study.name) > 0 & df$subtype == 'PXG') 
    df <- agg.arms(df, regexpr('Lee', df$study.name) > 0 & df$Year == '2016')
    df <- agg.arms(df, regexpr('Euswas', df$study.name) > 0)
    stopifnot(nrow.before - nrow(df) == 3)
  }
  
  # Fill in information about missing time information.
  missing.time <- (is.na(df$LastPeriodIOPMean) & is.na(df$LastPeriodAbsIOPChangeMean)) & (!is.na(df$OneYAbsIOPChangeMean) | !is.na(df$OneYIOPMean))
  df[missing.time,]$TimeofLastPostOp <- "12 mo"
  
  missing.time <- is.na(df$TimeofLastPostOp)
  df[missing.time,]$TimeofLastPostOp <- "6 mo"
  
  # Fix this value which doesn't have a space.
  df[df$TimeofLastPostOp == '24mo',]$TimeofLastPostOp <- "24 mo"
  
  last.space <- str_locate(df$TimeofLastPostOp, ' ')[,1]
  
  expect_equal(sum(is.na(last.space)), 0)
  df$StudyLength <- as.numeric(substr(df$TimeofLastPostOp, 1, last.space))
  expect_equal(sum(is.na(df$StudyLength)), 0)
  expect_gt(sum(df$StudyLength == 6), 0)
  expect_gt(sum(df$StudyLength == 12), 0)
  expect_equal(sum(df$StudyLength == 7.2), 1)
  
  # Clean up time of last follow up. Anything less than a year is put into the
  # 6 month bin. This should only affect 3 studies (4 arms) that have between 6 and 
  # 12 month follow-up, and one MIGS which is weirdly coded.
  low.followup <- (df$StudyLength < 12) & (!is.na(df$LastPeriodIOPMean) | !is.na(df$LastPeriodAbsIOPChangeMean))
  expect_true(sum(low.followup)>= 4 && sum(low.followup)<= 5)
  
  df[low.followup,] <- df[low.followup,] %>% mutate(
    SixMoEyes = LastPeriodEyes,
    SixMoIOPMean = LastPeriodIOPMean,
    SixMoIOPStdDev = LastPeriodIOPStdDev,
    SixMoAbsIOPChangeMean = LastPeriodAbsIOPChangeMean,
    SixMoAbsIOPChangeStdDev = LastPeriodAbsIOPChangeStdDev,
    LastPeriodEyes = NA,
    LastPeriodIOPMean = NA,
    LastPeriodIOPStdDev = NA,
    LastPeriodAbsIOPChangeStdDev = NA,
    LastPeriodAbsIOPChangeMean = NA
  )
  
  if(fill.last) {
    # Where not otherwise specified, fill in the data for last follow-up with 
    # the 12 month data.
    missing.data <- is.na(df$LastPeriodAbsIOPChangeMean)
    
    # Fill in last follow-up data.
    df[missing.data,]$LastPeriodIOPMean <- df[missing.data,]$OneYIOPMean
    df[missing.data,]$LastPeriodIOPStdDev <- df[missing.data,]$OneYIOPStdDev
    df[missing.data,]$LastPeriodAbsIOPChangeMean <- df[missing.data,]$OneYAbsIOPChangeMean
    df[missing.data,]$LastPeriodAbsIOPChangeStdDev <- df[missing.data,]$OneYAbsIOPChangeStdDev
    df[missing.data,]$LastPeriodEyes <- df[missing.data,]$OneYEyes
    
  }
  
  if(exclude.short) {
    # Drop experiments where there's < 12 month follow-up.
    df <- df[df$StudyLength >= 12,]
  }
  
  df <- df %>% dplyr::arrange(Year, study.name)
  
  if(drop.migs) {
    # Originally, the paper included the comparison of different MIGS 
    # treatments, but we didn't the space and time to do it properly, so rather
    # than do it sloppily we decided to leave it for a new paper.
    df <- df %>% filter(MIGsYorN == 'N')
  }
  
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

read.column.names <- function() {
  fileName <- 'phaco.csv'
  df <- read.csv(fileName, na.strings='-', check.names = FALSE)
  return(colnames(df))
}

# Unit test agg.arms, without polluting the global namespace.
test.fun <- function() {
  df <- read.data(FALSE, exclude.short=FALSE)
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
  # Check that every IOP measurement has its own SD.
  expect_equal(nrow(df[!is.na(df$SixMoAbsIOPChangeMean) & is.na(df$SixMoAbsIOPChangeStdDev),]), 0)
  expect_equal(nrow(df[!is.na(df$OneYAbsIOPChangeMean) & is.na(df$OneYAbsIOPChangeStdDev),]), 0)
  expect_equal(nrow(df[!is.na(df$LastPeriodAbsIOPChangeMean) & is.na(df$LastPeriodAbsIOPChangeStdDev),]), 0)
}

test.fun()
