#import libraries
# library(data.table)
# library(plyr)
# library(dplyr)
#
# #set working directory
# setwd("/Users/sonjaperkovic/Dropbox/SSI_R_package")
#
# #simulate data ####
# infoSearch <- data.table(participant=rep(c(1:3), times = c(100,120,150)))
# sim <- nrow(infoSearch)
# infoSearch$trial <- ifelse(infoSearch$participant == 1, rep(c(1:5), times = c(20,15,25,18,22)),
#                        ifelse(infoSearch$participant == 2, rep(c(1:5), times = c(30,25,20,25,20)),
#                               ifelse(infoSearch$participant == 3, rep(c(1:5), times = c(30,35,25,32,28)), 0)))
# infoSearch$alternative <- sample(1:3, sim, T)
# infoSearch$attribute <- sample(c("a","b","c"), sim, T)
# infoSearch <- infoSearch[order(participant, trial),]
#
# infoSearchRan <- NULL
# infoSearchRan$participant <- infoSearch$participant
# infoSearchRan$trial <- infoSearch$trial
# infoSearchRan <- as.data.table(infoSearchRan)

# infoSearch <- as.data.table(read.csv("gabor_clean.csv", header = T, sep = ","))
#
# infoSearchRan <- as.data.table(read.csv("gabor_clean.csv", header = T, sep = ","))
# infoSearchRan[, c("alternative", "attribute") := NULL]

#identify alternative-wise patterns ####

altwise <- function(df, participant, trial, alternative, attribute) {
  #create counter variable that assigns a new number whenever the alternative changes
  df <- setDT(df)[, counter:= rleid(trial, alternative)]
  #concatenate attribute values into strings
  #this is done for each alternative, within each trial and for each participant
  df1 <- df[,list(string <- paste(attribute, collapse = ""),
                                  participant = unique(participant),
                                  trial = unique(trial)), by = counter]
  #delete counter varibale
  df1[, "counter" := NULL]
  #function which keeps unique letters in a string and orders them alphabetically
  relaxedFreqOrder <- function(i){
    paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
  }
  #applying the function
  df1$string <- lapply(df1$V1, relaxedFreqOrder)
  #deleting unnecessary variable
  df1[, "V1" := NULL]
  #changing class
  df1$string <- as.character(df1$string)
  #creating counter variable which identifies identical subsequent strings
  df1 <- setDT(df1)[, counter:= rleid(string, trial)]
  #variable that identifies identical subsequent counter variables
  df1$equalCounter <- ifelse(df1$counter == lag(df1$counter, n = 1L)
                            | df1$counter == lead(df1$counter, n = 1L), 1, 0)
  #subset the data
  df1 <- subset(df1, equalCounter != 0)
  #keep string of length at least two letters
  df1 <- subset(df1, nchar(as.character(string)) >= 2)
  #concatenate strings into patterns using the counter variable
  df2 <- df1[,list(string <- paste(string, collapse = ""),
                                          participant = unique(participant),
                                          trial = unique(trial)), by = counter]

  #delete counter variable
  df2[, "counter" := NULL]
  #count how many times each unique pattern occurs
  df3 <- as.data.table(with(df2, table(V1, trial, participant)))
  #subest the data
  df3 <- subset(df3, N != 0)
  #rename column
  setnames(df3, "V1", "pattern")
  df3
}

#test = altwise(infoSearch, "participant", "trial", "alternative", "attribute")

#alternative-wise pattern simulation ####
#the same procedure as for altwise function, just on a randomized data set
altwiseSim <- function(df, participant, trial, num_alt, num_att) {
  #making sure we sample the total number of elements from the original data set
  sim <- nrow(df)
  #creating new alternative variable (random)
  df$alternative <- sample(1:num_alt, sim, T)
  #creating new attribute variable (random)
  attset <- letters[1:num_att]
  df$attribute <- sample(attset, sim, T)
  df <- setDT(df)[, counter:= rleid(trial, alternative)]
  df1 <- df[,list(string <- paste(attribute, collapse = ""),
                  participant = unique(participant),
                  trial = unique(trial)), by = counter]
  df1[, "counter" := NULL]
  relaxedFreqOrder <- function(i){
    paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
  }
  df1$string <- lapply(df1$V1, relaxedFreqOrder)
  df1[, "V1" := NULL]
  df1$string <- as.character(df1$string)
  df1 <- setDT(df1)[, counter:= rleid(string, trial)]
  df1$equalCounter <- ifelse(df1$counter == lag(df1$counter, n = 1L)
                             | df1$counter == lead(df1$counter, n = 1L), 1, 0)
  df1 <- subset(df1, equalCounter != 0)
  df1 <- subset(df1, nchar(as.character(string)) >= 2)
  df2 <- df1[,list(string <- paste(string, collapse = ""),
                   participant = unique(participant),
                   trial = unique(trial)), by = counter]

  df2[, "counter" := NULL]
  df3 <- as.data.table(with(df2, table(V1, trial, participant)))
  df3 <- subset(df3, N != 0)
  colnames(df3)[c(1,4)] <- c("pattern", "N_sim")
  df3
}

#test1 = altwiseSim(infoSearchRan, "participant", "trial", 3, 3)

#replicate 'altwiseSim' function 10000 times ####

altwiseSimRep <- function(df, participant, trial, num_alt, num_att, iter) {
  do.call(rbind, lapply(iter, function(i) altwiseSim(df, participant, trial, num_alt, num_att)))
}

#test2 = altwiseSimRep(infoSearchRan,"participant", "trial", 3, 3, 10)

#calculate probabilities and probability complements for alternative-wise patterns ####

probAltwise <- function(df, df1, iter) {
 #creating a function which will compare original data set with randomized data sets
 altwiseProb <- function(i){
 sum(df1$pattern == df$pattern[i]
     & df1$participant == df$participant[i]
     & df1$trial == df$trial[i]
     & df1$N_sim >= df$N[i])
}
 #applying the function
 df$N_sim <- sapply(1:nrow(df), altwiseProb)
 #calculating probability
 df$probability <- df$N_sim / iter
 #calculating probability complement
 df$prob_complement <- 1 - df$probability
 #calculating pattern leght
 df$patt_length <- nchar(df$pattern)
 df
 #saving file
 write.csv(file="altwisePatterns.csv", x=df)
 df
}

#test3 = probAltwise(test, test2, 10)

#identify attribute-wise patterns ####

attwise <- function(df, participant, trial, alternative, attribute) {
  #deleting dwells (subsequent fixations within the same AOI)
  df$attribute_clean <- ifelse(df[[alternative]] == shift(df[[alternative]], 1L)
                              & df[[attribute]] == shift(df[[attribute]], 1L), 1, 0)
  df <- subset(df, attribute_clean != 1 | is.na(attribute_clean))
  df[, "attribute_clean" := NULL]
  #creating counter variable which identifies fixation change for attributes (fixating on a new attribute)
  df <- setDT(df)[, counterAttwise:= rleid(trial, attribute)]
  df1 <- df[,list(string <- paste(attribute, collapse = ""),
                                  participant = unique(participant),
                                  trial = unique(trial)), by = counterAttwise]
  #subseting patterns of length 3 or less
  df1 <- subset(df1, nchar(as.character(V1)) >= 4)
  #deleting unnecessary variable
  df1[, "counterAttwise" := NULL]
  #counting occurrences of unique patterns
  df2 <- as.data.table(with(df1, table(V1, trial, participant)))
  df2 <- subset(df2, N != 0)
  #renaming column
  setnames(df2, "V1", "pattern")
  df2
}

#test4 = attwise(infoSearch, "participant", "trial", "alternative", "attribute")

#attribute-wise pattern simulation ####

#the same procesure as for attwise function, just for randomized data set
attwiseSim <- function(df, participant, trial, num_alt, num_att) {
  #making sure we sample the total number of elements from the original data set
  sim <- nrow(df)
  #creating new alternative variable (random)
  df$alternative <- sample(1:num_alt, sim, T)
  #creating new attribute variable (random)
  attset <- letters[1:num_att]
  df$attribute <- sample(attset, sim, T)
  df$attribute_clean <- ifelse(df$alternative == shift(df$alternative, 1L)
                               & df$attribute == shift(df$attribute, 1L), 1, 0)
  df <- subset(df, attribute_clean != 1 | is.na(attribute_clean))
  df[, "attribute_clean" := NULL]
  df <- setDT(df)[, counter:= rleid(trial, attribute)]
  df1 <- df[,list(string <- paste(attribute, collapse = ""),
                  participant = unique(participant),
                  trial = unique(trial)), by = counter]
  df1 <- subset(df1, nchar(as.character(V1)) >= 4)
  df1[, "counter" := NULL]
  df2 <- as.data.table(with(df1, table(V1, trial, participant)))
  df2 <- subset(df2, N != 0)
  colnames(df2)[c(1,4)] <- c("pattern", "N_sim")
  df2
}

#test5 = attwiseSim(infoSearchRan, "participant", "trial", 3, 3)

#replicate 'attwiseSim' function 10000 times ####

attwiseSimRep <- function(df, participant, trial, num_alt, num_att, iter) {
  do.call(rbind, lapply(iter, function(i) attwiseSim(df, participant, trial, num_alt, num_att)))
}

#test6 = attwiseSimRep(infoSearchRan, "participant", "trial", 3, 3, 10)

#calculate probabilities and probability complements for attribute-wise patterns####

#creating a function which will compare original data set with randomized data sets
probAttwise <- function(df, df1, iter) {
attwiseProb <- function(i){
  sum(df1$pattern == df$pattern[i]
      & df1$participant == df$participant[i]
      & df1$trial == df$trial[i]
      & df1$N_sim >= df$N[i])
}
 #applying the function
 df$N_sim <- sapply(1:nrow(df), attwiseProb)
 #calculating probability
 df$probability <- df$N_sim / iter
 #calculating probability complement
 df$prob_complement <- 1 - df$probability
 #calculating pattern lenght
 df$patt_length <- nchar(df$pattern)
 df
 #saving file
 write.csv(file="attwisePatterns.csv", x=df)
 df
}

#test7 = probAttwise(test4, test6, 10)

#calculate SSI ####

computeSSI <- function(df, df1, df3, participant, trial, alternative, attribute) {
  #calculating string length for each trial
  df$attribute_clean <- ifelse(df[[alternative]] == shift(df[[alternative]], 1L)
                                & df[[attribute]] == shift(df[[attribute]], 1L), 1, 0)
  df <- subset(df, attribute_clean != 1 | is.na(attribute_clean))
  df[, "attribute_clean" := NULL]
  stringLength <- ddply(df, .(participant, trial), function(df) length(df$attribute))
  #calculating numerator for alternative-wise patterns
  df1$numerator <- df1$N * df1$patt_length * df1$prob_complement
  df2 <- ddply(df1,.(participant, trial), summarize, altwise_sum = sum(numerator))
  df2 <- as.data.table(df2)
  #changing class
  df2$participant <- as.numeric(df2$participant)
  df2$trial <- as.numeric(df2$trial)
  #ordering values within variables
  df2 <- df2[order(participant, trial),]
  #merging in trial string leghts
  df2 <- merge(df2, stringLength, by = c("participant", "trial"), all = T)
  df2[is.na(df2)] <- 0
  #calculating numerator for attribute-wise patterns
  df3$numerator <- df3$N * df3$patt_length * df3$prob_complement
  df4 <- ddply(df3,.(participant, trial), summarize, attwise_sum = sum(numerator))
  df4 <- as.data.table(df4)
  #changing class
  df4$participant <- as.numeric(df4$participant)
  df4$trial <- as.numeric(df4$trial)
  #ordering values within variables
  df4 <- df4[order(participant, trial),]
  #merging in trial string lengths
  df4 <- merge(df4, stringLength, by = c("participant", "trial"), all = T)
  df4[is.na(df4)] <- 0
  #merging data sets with alternative and attribute-wise pattern numerators
  df5 <- merge(df2, df4, by = c("participant", "trial"), all = T)
  #deleting variable
  df5$V1.x <- NULL
  #renaming column
  setnames(df5, "V1.y", "string_length")
  #applying SSI equation
  df5$SSI <- (df5$altwise_sum + df5$attwise_sum) / df5$string_length
  df5
  #saving file
  write.csv(file="SSI.csv", x=df)
  df5
}

#test8 = computeSSI(infoSearch, test3, test7, "participant", "trial", "alternative", "attribute")

##################################################
#a function that wraps all previous functions ####

#identify alternative-wise patterns ####
computeSSIfast <- function(df, dfRan, participant, trial, alternative, attribute, num_alt, num_att, iter) {

  altwise <- function(df, participant, trial, alternative, attribute) {
  #create counter variable that assigns a new number whenever the alternative changes
  df <- setDT(df)[, counter:= rleid(trial, alternative)]
  #concatenate attribute values into strings
  #this is done for each alternative, within each trial and for each participant
  df1 <- df[,list(string <- paste(attribute, collapse = ""),
                  participant = unique(participant),
                  trial = unique(trial)), by = counter]
  #delete counter varibale
  df1[, "counter" := NULL]
  #function which keeps unique letters in a string and orders them alphabetically
  relaxedFreqOrder <- function(i){
    paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
  }
  #applying the function
  df1$string <- lapply(df1$V1, relaxedFreqOrder)
  #deleting unnecessary variable
  df1[, "V1" := NULL]
  #changing class
  df1$string <- as.character(df1$string)
  #creating counter variable which identifies identical subsequent strings
  df1 <- setDT(df1)[, counter:= rleid(string, trial)]
  #variable that identifies identical subsequent counter variables
  df1$equalCounter <- ifelse(df1$counter == lag(df1$counter, n = 1L)
                             | df1$counter == lead(df1$counter, n = 1L), 1, 0)
  #subset the data
  df1 <- subset(df1, equalCounter != 0)
  #keep string of length at least two letters
  df1 <- subset(df1, nchar(as.character(string)) >= 2)
  #concatenate strings into patterns using the counter variable
  df2 <- df1[,list(string <- paste(string, collapse = ""),
                   participant = unique(participant),
                   trial = unique(trial)), by = counter]

  #delete counter variable
  df2[, "counter" := NULL]
  #count how many times each unique pattern occurs
  df3 <- as.data.table(with(df2, table(V1, trial, participant)))
  #subest the data
  df3 <- subset(df3, N != 0)
  #rename column
  setnames(df3, "V1", "pattern")
  df3
}

test = altwise(df, "participant", "trial", "alternative", "attribute")
test
#test0 = computeSSIfast(infoSearch, infoSearchRan, "participant", "trial", "alternative", "attribute", 3, 3)

#alternative-wise pattern simulation ####
#the same procedure as for altwise function, just on a randomized data set
altwiseSim <- function(dfRan, participant, trial, num_alt, num_att) {
  #making sure we sample the total number of elements from the original data set
  sim <- nrow(dfRan)
  #creating new alternative variable (random)
  dfRan$alternative <- sample(1:num_alt, sim, T)
  #creating new attribute variable (random)
  attset <- letters[1:num_att]
  dfRan$attribute <- sample(attset, sim, T)
  dfRan <- setDT(dfRan)[, counter:= rleid(trial, alternative)]
  dfRan1 <- dfRan[,list(string <- paste(attribute, collapse = ""),
                  participant = unique(participant),
                  trial = unique(trial)), by = counter]
  dfRan1[, "counter" := NULL]
  relaxedFreqOrder <- function(i){
    paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
  }
  dfRan1$string <- lapply(dfRan1$V1, relaxedFreqOrder)
  dfRan1[, "V1" := NULL]
  dfRan1$string <- as.character(dfRan1$string)
  dfRan1 <- setDT(dfRan1)[, counter:= rleid(string, trial)]
  dfRan1$equalCounter <- ifelse(dfRan1$counter == lag(dfRan1$counter, n = 1L)
                             | dfRan1$counter == lead(dfRan1$counter, n = 1L), 1, 0)
  dfRan1 <- subset(dfRan1, equalCounter != 0)
  dfRan1 <- subset(dfRan1, nchar(as.character(string)) >= 2)
  dfRan2 <- dfRan1[,list(string <- paste(string, collapse = ""),
                   participant = unique(participant),
                   trial = unique(trial)), by = counter]

  dfRan2[, "counter" := NULL]
  dfRan3 <- as.data.table(with(dfRan2, table(V1, trial, participant)))
  dfRan3 <- subset(dfRan3, N != 0)
  colnames(dfRan3)[c(1,4)] <- c("pattern", "N_sim")
  dfRan3
}

test1 = altwiseSim(dfRan, "participant", "trial", num_alt, num_att)

#test01 = computeSSIfast(infoSearch, infoSearchRan, "participant", "trial", "alternative", "attribute", 3, 3)

#replicate 'altwiseSim' function 10000 times ####

altwiseSimRep <- function(dfRan, participant, trial, num_alt, num_att, iter) {
  do.call(rbind, lapply(iter, function(i) altwiseSim(dfRan, participant, trial, num_alt, num_att)))
}

test2 = altwiseSimRep(dfRan,"participant", "trial", num_alt, num_att, iter)
test2

#test = computeSSIfast(infoSearch, infoSearchRan, "participant", "trial", "alternative", "attribute", 3, 3)

#calculate probabilities and probability complements for alternative-wise patterns ####

probAltwise <- function(df, df1, iter) {
  #creating a function which will compare original data set with randomized data sets
  altwiseProb <- function(i){
    sum(df1$pattern == df$pattern[i]
        & df1$participant == df$participant[i]
        & df1$trial == df$trial[i]
        & df1$N_sim >= df$N[i])
  }
  #applying the function
  df$N_sim <- sapply(1:nrow(df), altwiseProb)
  #calculating probability
  df$probability <- df$N_sim / iter
  #calculating probability complement
  df$prob_complement <- 1 - df$probability
  #calculating pattern length
  df$patt_length <- nchar(df$pattern)
  df
}

test3 = probAltwise(test, test2, iter)
test3

#test_new = computeSSIfast(infoSearch, infoSearchRan, "participant", "trial", "alternative", "attribute", 4, 4)

#identify attribute-wise patterns ####

attwise <- function(df, participant, trial, alternative, attribute) {
  #deleting dwells (subsequent fixations within the same AOI)
  df$attribute_clean <- ifelse(df[[alternative]] == shift(df[[alternative]], 1L)
                               & df[[attribute]] == shift(df[[attribute]], 1L), 1, 0)
  df <- subset(df, attribute_clean != 1 | is.na(attribute_clean))
  df[, "attribute_clean" := NULL]
  #creating counter variable which identifies fixation change for attributes (fixating on a new attribute)
  df <- setDT(df)[, counterAttwise:= rleid(trial, attribute)]
  df1 <- df[,list(string <- paste(attribute, collapse = ""),
                  participant = unique(participant),
                  trial = unique(trial)), by = counterAttwise]
  #subseting patterns of length 3 or less
  df1 <- subset(df1, nchar(as.character(V1)) >= 4)
  #deleting unnecessary variable
  df1[, "counterAttwise" := NULL]
  #counting occurrences of unique patterns
  df2 <- as.data.table(with(df1, table(V1, trial, participant)))
  df2 <- subset(df2, N != 0)
  #renaming column
  setnames(df2, "V1", "pattern")
  df2
}

test4 = attwise(df, "participant", "trial", "alternative", "attribute")

#attribute-wise pattern simulation ####

#the same procesure as for attwise function, just for randomized data set
attwiseSim <- function(dfRan, participant, trial, num_alt, num_att) {
  #making sure we sample the total number of elements from the original data set
  sim <- nrow(df)
  #creating new alternative variable (random)
  dfRan$alternative <- sample(1:num_alt, sim, T)
  #creating new attribute variable (random)
  attset <- letters[1:num_att]
  dfRan$attribute <- sample(attset, sim, T)
  dfRan$attribute_clean <- ifelse(dfRan$alternative == shift(dfRan$alternative, 1L)
                               & dfRan$attribute == shift(dfRan$attribute, 1L), 1, 0)
  dfRan <- subset(dfRan, attribute_clean != 1 | is.na(attribute_clean))
  dfRan[, "attribute_clean" := NULL]
  dfRan <- setDT(dfRan)[, counter:= rleid(trial, attribute)]
  dfRan1 <- dfRan[,list(string <- paste(attribute, collapse = ""),
                  participant = unique(participant),
                  trial = unique(trial)), by = counter]
  dfRan1 <- subset(dfRan1, nchar(as.character(V1)) >= 4)
  dfRan1[, "counter" := NULL]
  dfRan2 <- as.data.table(with(dfRan1, table(V1, trial, participant)))
  dfRan2 <- subset(dfRan2, N != 0)
  colnames(dfRan2)[c(1,4)] <- c("pattern", "N_sim")
  dfRan2
}

#test5 = attwiseSim(infoSearchRan, "participant", "trial", num_alt, num_att)

#replicate 'attwiseSim' function 10000 times ####

attwiseSimRep <- function(dfRan, participant, trial, num_alt, num_att, iter) {
  do.call(rbind, lapply(iter, function(i) attwiseSim(dfRan, participant, trial, num_alt, num_att)))
}

test6 = attwiseSimRep(dfRan, "participant", "trial", num_alt, num_att, iter)
test6

#calculate probabilities and probability complements for attribute-wise patterns####

#creating a function which will compare original data set with randomized data sets
probAttwise <- function(df, df1, iter) {
  attwiseProb <- function(i){
    sum(df1$pattern == df$pattern[i]
        & df1$participant == df$participant[i]
        & df1$trial == df$trial[i]
        & df1$N_sim >= df$N[i])
  }
  #applying the function
  df$N_sim <- sapply(1:nrow(df), attwiseProb)
  #calculating probability
  df$probability <- df$N_sim / iter
  #calculating probability complement
  df$prob_complement <- 1 - df$probability
  #calculating pattern lenght
  df$patt_length <- nchar(df$pattern)
  df
}

test7 = probAttwise(test4, test6, iter)
test7

#test_newest <- computeSSIfast(infoSearch, infoSearchRan, "participant", "trial", "alternative", "attribute", 3, 3)

#calculate SSI ####

computeSSI <- function(df, df1, df3, participant, trial, alternative, attribute) {
  #calculating string length for each trial
  df$attribute_clean <- ifelse(df[[alternative]] == shift(df[[alternative]], 1L)
                               & df[[attribute]] == shift(df[[attribute]], 1L), 1, 0)
  df <- subset(df, attribute_clean != 1 | is.na(attribute_clean))
  df[, "attribute_clean" := NULL]
  stringLength <- ddply(df, .(participant, trial), function(df) length(df$attribute))
  #calculating numerator for alternative-wise patterns
  df1$numerator <- df1$N * df1$patt_length * df1$prob_complement
  df2 <- ddply(df1,.(participant, trial), summarize, altwise_sum = sum(numerator))
  df2 <- as.data.table(df2)
  #changing class
  df2$participant <- as.numeric(df2$participant)
  df2$trial <- as.numeric(df2$trial)
  #ordering values within variables
  df2 <- df2[order(participant, trial),]
  #merging in trial string leghts
  df2 <- merge(df2, stringLength, by = c("participant", "trial"), all = T)
  df2[is.na(df2)] <- 0
  #calculating numerator for attribute-wise patterns
  df3$numerator <- df3$N * df3$patt_length * df3$prob_complement
  df4 <- ddply(df3,.(participant, trial), summarize, attwise_sum = sum(numerator))
  df4 <- as.data.table(df4)
  #changing class
  df4$participant <- as.numeric(df4$participant)
  df4$trial <- as.numeric(df4$trial)
  #ordering values within variables
  df4 <- df4[order(participant, trial),]
  #merging in trial string lengths
  df4 <- merge(df4, stringLength, by = c("participant", "trial"), all = T)
  df4[is.na(df4)] <- 0
  #merging data sets with alternative and attribute-wise pattern numerators
  df5 <- merge(df2, df4, by = c("participant", "trial"), all = T)
  #deleting variable
  df5$V1.x <- NULL
  #renaming column
  setnames(df5, "V1.y", "string_length")
  #applying SSI equation
  df5$SSI <- (df5$altwise_sum + df5$attwise_sum) / df5$string_length
  df5
}

test8 = computeSSI(df, test3, test7, "participant", "trial", "alternative", "attribute")
test8
}

#test_final = computeSSIfast(infoSearch, infoSearchRan, "participant", "trial", "alternative", "attribute", 3, 3, 100)




