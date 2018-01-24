#' A function that computes Systematicity of Search Index (SSI)
#'
#' Integrates all the necessary functions and computes the SSI. The number of alternatives and attributes, simulations
#' as well as the minimum threshold for pattern length can be adjusted.
#'
#' @param df Object of class data frame
#' @param dfRan The same object used for creating random data.
#' @param participant Identifies each unique subject.
#' @param trial Identifies each unique trial.
#' @param alternative Represents column with eye fixations to different alternatives.
#' @param attribute Represents column with eye fixations to different attributes.
#' @param num_alt Number of alternatives in the experiment.
#' @param num_att Number of attributes in the experiment.
#' @param threshold Sets the threshold for pattern length to two or four.
#' @param iter Number of simulation iterations.
#'
#' @return A data set with SSI values for all identified alternative- and attribute-wise patterns.
#'
#' @author Sonja Perkovic, \email{bnsp@leeds.ac.uk}
#' @keywords SSI
#'
#' @examples
#' #IMPORTANT! Variables representing participants, trials, alternatives and attributes in your data should have the names
#' that match the ones in the example below!
#'
#' dataSet <- data.frame(participant = rep(c(1:50), each = 400),
#'                       trial = rep(c(1:200), each = 100),
#'                       alternative = sample(1:4, 20000, TRUE),
#'                       attribute = sample(c("a","b","c","d"), 20000, TRUE))
#'
#' SSI <- computeSSI(dataSet, dataSet, "participant", "trial", "alternative", "attribute", 4, 4, 4, 10)
#'
#' @export
#'

#a function that wraps all necessary functions and computes SSI

#identify alternative- and attribute-wise patterns ####

computeSSI = function(df, dfRan, participant, trial, alternative, attribute, num_alt, num_att, threshold, iter) {

  patterns = function(df, participant, trial, alternative, attribute, threshold) {

    #change class
    df = as.data.table(df)

    #delete dwells (subsequent fixations within the same AOI)
    df$attributeClean = ifelse(df$trial == shift(df$trial, 1L)
                             & df$alternative == shift(df$alternative, 1L)
                             & df$attribute == shift(df$attribute, 1L), 1, 0)
    df = subset(df, attributeClean != 1 | is.na(attributeClean))
    df$attributeClean = NULL

    #identify alternative- and attribute-wise transitions
    n = nrow(df)
    df$transAlt = c(0L, df$trial[-n] == df$trial[-1] & df$alternative[-n] == df$alternative[-1])
    df$transAtt = c(0L, df$trial[-n] == df$trial[-1] & df$attribute[-n] == df$attribute[-1])

    #create counter variable which will be used for concatenating fixations into substrings
    df$transDiff = c(0L, df$transAtt[-n] != df$transAtt[-1] & df$transAlt[-n] != df$transAlt[-1])
    df$diff1 <- c(NA, diff(df$transDiff))

    df$counter = c(1L, df$trial[-n] != df$trial[-1]
                     | c(df$trial[-n] == df$trial[-1] &
                         df$alternative[-n] != df$alternative[-1] &
                         df$attribute[-n] != df$attribute[-1])

                     | c(df$trial[-n] == df$trial[-1] &
                         df$alternative[-n] != df$alternative[-1] &
                         df$attribute[-n] == df$attribute[-1] &
                         df$transAtt[-n] != df$transAtt[-1] &
                         df$transAlt[-n] != df$transAlt[-1] &
                         df$transDiff[-n] != df$transDiff[-1] &
                         df$diff1[-n] != df$diff1[-1])

                     | c(df$trial[-n] == df$trial[-1] &
                         df$alternative[-n] != df$alternative[-1] &
                         df$attribute[-n] == df$attribute[-1] &
                         df$transAtt[-n] != df$transAtt[-1] &
                         df$transAlt[-n] != df$transAlt[-1] &
                         df$transDiff[-n] == df$transDiff[-1] &
                         df$diff1[-n] == df$diff1[-1])

                     | c(df$trial[-n] == df$trial[-1] &
                         df$alternative[-n] == df$alternative[-1] &
                         df$attribute[-n] != df$attribute[-1] &
                         df$transAtt[-n] != df$transAtt[-1] &
                         df$transAlt[-n] != df$transAlt[-1] &
                         df$transDiff[-n] != df$transDiff[-1]))
    df$counter = cumsum(df$counter)

    #create alternative- and attribute-wise substrings based on counter variable
    df1 = df[, list(string = paste(attribute, collapse = ""),
                                   participant = unique(participant),
                                   trial = unique(trial)), by = counter]

    #function that identifies rows with identical elements only
    identElements = function(i){
      length(unique(unlist(strsplit(i, "")))) == 1
    }

    #apply 'identElements' function
    df1$identElements = lapply(df1$string, identElements)

    #subset data with alternative-wise substrings only
    df2 = subset(df1, identElements == "FALSE")

    #function that keeps unique elements and sorts them alphabetically
    relaxedFreqOrder = function(i){
      paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
    }

    #apply 'relaxedFreqOrder' function
    df2$string = lapply(df2$string, relaxedFreqOrder)

    #change class
    df1$string = as.character(df1$string)
    df1$identElements = as.numeric(df1$identElements)
    df2$string = as.character(df2$string)
    df2$identElements = as.numeric(df2$identElements)

    #merge in formatted data set with altwise substrings
    df1 = merge(df1, df2, by = c("counter", "participant", "trial", "identElements"), all.x = T)
    df1$string.y = ifelse(is.na(df1$string.y), df1$string.x, df1$string.y)
    df1$string.x = NULL
    setnames(df1, "string.y", "string")

    #create counter variable for alternative-wise substrings based on string variable within each trial
    df1 = setDT(df1)[, countEqualSubstrings := rleid(string, trial)]

    #if threshold = 4, then compile next two lines, otherwise skip to line 111
    if (threshold == 4) {

    #create variable that assigns 1 to subsequent equal counter variable values
    df1$equalCounter = ifelse(df1$countEqualSubstrings == shift(df1$countEqualSubstrings, 1L) |
                              df1$countEqualSubstrings == shift(df1$countEqualSubstrings, 1L, type = "lead"), 1, 0)

    #subset substrings
    df1 = df1[df1$equalCounter != 0 | df1$identElements != 0]
    }

    #combine substrings into patterns using counter variable
    #i.e., all substrings with equal count should be collapsed into one pattern
    df2 = df1[,list(pattern = paste(string, collapse = ""),
                              participant = unique(participant),
                              trial = unique(trial)), by = countEqualSubstrings]

    #keep substrings of minimum length four
    df2 = subset(df2, nchar(as.character(pattern)) >= threshold)

    #calculate frequencies for each pattern within each trial and participant
    df3 = setDT(df2)[, list(pattFreq = .N), by = c('pattern', 'trial', 'participant')]
    df3
  }

  test = patterns(df, "participant", "trial", "alternative", "attribute", threshold)
  test

  #alternative- and attribute-wise pattern simulation ####
  #the same procedure as for patterns function, just on a random data set

  patternsSim = function(dfRan, participant, trial, num_alt, num_att, threshold) {

    #change class
    dfRan = as.data.table(dfRan)

    #number of rows corresponding to number of fixations in original data set (total string length)
    sim = nrow(dfRan)

    #create new alternative variable (random)
    dfRan$alternative = sample(1:num_alt, sim, T)

    #create new attribute variable (random)
    attset = letters[1:num_att]
    dfRan$attribute = sample(attset, sim, T)

    #delete dwells (subsequent fixations within the same AOI)
    dfRan$attributeClean = ifelse(dfRan$trial == shift(dfRan$trial, 1L)
                                & dfRan$alternative == shift(dfRan$alternative, 1L)
                                & dfRan$attribute == shift(dfRan$attribute, 1L), 1, 0)
    dfRan = subset(dfRan, attributeClean != 1 | is.na(attributeClean))
    dfRan$attributeClean = NULL

    #identify alternative- and attribute-wise transitions
    n = nrow(dfRan)
    dfRan$transAlt = c(0L, dfRan$trial[-n] == dfRan$trial[-1] & dfRan$alternative[-n] == dfRan$alternative[-1])
    dfRan$transAtt = c(0L, dfRan$trial[-n] == dfRan$trial[-1] & dfRan$attribute[-n] == dfRan$attribute[-1])

    #create counter variable which will be used for concatenating fixations into substrings
    dfRan$transDiff = c(0L, dfRan$transAtt[-n] != dfRan$transAtt[-1] & dfRan$transAlt[-n] != dfRan$transAlt[-1])
    dfRan$diff1 <- c(NA, diff(dfRan$transDiff))

    dfRan$counter = c(1L, dfRan$trial[-n] != dfRan$trial[-1]
                   | c(dfRan$trial[-n] == dfRan$trial[-1] &
                         dfRan$alternative[-n] != dfRan$alternative[-1] &
                         dfRan$attribute[-n] != dfRan$attribute[-1])

                   | c(dfRan$trial[-n] == dfRan$trial[-1] &
                         dfRan$alternative[-n] != dfRan$alternative[-1] &
                         dfRan$attribute[-n] == dfRan$attribute[-1] &
                         dfRan$transAtt[-n] != dfRan$transAtt[-1] &
                         dfRan$transAlt[-n] != dfRan$transAlt[-1] &
                         dfRan$transDiff[-n] != dfRan$transDiff[-1] &
                         dfRan$diff1[-n] != dfRan$diff1[-1])

                   | c(dfRan$trial[-n] == dfRan$trial[-1] &
                         dfRan$alternative[-n] != dfRan$alternative[-1] &
                         dfRan$attribute[-n] == dfRan$attribute[-1] &
                         dfRan$transAtt[-n] != dfRan$transAtt[-1] &
                         dfRan$transAlt[-n] != dfRan$transAlt[-1] &
                         dfRan$transDiff[-n] == dfRan$transDiff[-1] &
                         dfRan$diff1[-n] == dfRan$diff1[-1])

                   | c(dfRan$trial[-n] == dfRan$trial[-1] &
                         dfRan$alternative[-n] == dfRan$alternative[-1] &
                         dfRan$attribute[-n] != dfRan$attribute[-1] &
                         dfRan$transAtt[-n] != dfRan$transAtt[-1] &
                         dfRan$transAlt[-n] != dfRan$transAlt[-1] &
                         dfRan$transDiff[-n] != dfRan$transDiff[-1]))
    dfRan$counter = cumsum(dfRan$counter)

    #create alternative- and attribute-wise substrings based on counter variable
    dfRan1 = dfRan[, list(string = paste(attribute, collapse = ""),
                    participant = unique(participant),
                    trial = unique(trial)), by = counter]

    #function that identifies rows with identical elements only
    identElements = function(i){
      length(unique(unlist(strsplit(i, "")))) == 1
    }

    #apply 'identElements' function
    dfRan1$identElements = lapply(dfRan1$string, identElements)

    #subset data with alternative-wise substrings only
    dfRan2 = subset(dfRan1, identElements == "FALSE")

    #function that keeps unique elements and sorts them alphabetically
    relaxedFreqOrder = function(i){
      paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
    }

    #apply 'relaxedFreqOrder' function
    dfRan2$string = lapply(dfRan2$string, relaxedFreqOrder)

    #change class
    dfRan1$string = as.character(dfRan1$string)
    dfRan1$identElements = as.numeric(dfRan1$identElements)
    dfRan2$string = as.character(dfRan2$string)
    dfRan2$identElements = as.numeric(dfRan2$identElements)

    #merge in formatted data set with altwise substrings
    dfRan1 = merge(dfRan1, dfRan2, by = c("counter", "participant", "trial", "identElements"), all.x = T)
    dfRan1$string.y = ifelse(is.na(dfRan1$string.y), dfRan1$string.x, dfRan1$string.y)
    dfRan1$string.x = NULL
    setnames(dfRan1, "string.y", "string")

    #create counter variable for alternative-wise substrings based on string variable within each trial
    dfRan1 = setDT(dfRan1)[, countEqualSubstrings := rleid(string, trial)]

    #if threshold = 4, then compile next two lines, otherwise skip to line 111
    if (threshold == 4) {

      #create variable that assigns 1 to subsequent equal counter variable values
      dfRan1$equalCounter = ifelse(dfRan1$countEqualSubstrings == shift(dfRan1$countEqualSubstrings, 1L) |
                                  dfRan1$countEqualSubstrings == shift(dfRan1$countEqualSubstrings, 1L, type = "lead"), 1, 0)

      #subset substrings
      dfRan1 = dfRan1[dfRan1$equalCounter != 0 | dfRan1$identElements != 0]
    }

    #combine substrings into patterns using counter variable
    #i.e., all substrings with equal count should be collapsed into one pattern
    dfRan2 = dfRan1[,list(pattern = paste(string, collapse = ""),
                    participant = unique(participant),
                    trial = unique(trial)), by = countEqualSubstrings]

    #keep substrings of minimum length four
    dfRan2 = subset(dfRan2, nchar(as.character(pattern)) >= threshold)

    #calculate frequencies for each pattern within each trial and participant
    dfRan3 = setDT(dfRan2)[, list(N = .N), by = c('pattern', 'trial', 'participant')]
    dfRan3
  }

  test1 = patternsSim(dfRan, "participant", "trial", num_alt, num_att, threshold)
  test1

  #replicate 'patternsSim' function n times ####

  patternsSimRep = function(dfRan, participant, trial, num_alt, num_att, threshold, iter) {
    do.call(rbind, lapply(1:iter, function(i) patternsSim(dfRan, participant, trial, num_alt, num_att, threshold)))
  }

  test2 = patternsSimRep(dfRan, "participant", "trial", num_alt, num_att, threshold, iter)
  test2

  #calculate probabilities and probability complements ####

  #function which compares pattern frequencies in original and simulated data sets for each participant, condition and trial
  patternsProb = function(df, df1, iter) {

    probPatterns = function(i) {
      sum(df1$pattern == df$pattern[i]
          & df1$participant == df$participant[i]
          & df1$trial == df$trial[i]
          & df1$N >= df$pattFreq[i])
    }

    #apply 'probPatterns' function
    df$pattFreqSim = sapply(1:nrow(df), probPatterns)

    #calculate probabilities
    df$probability = df$pattFreqSim / iter

    #calculate probability complements
    df$probComplement = 1 - df$probability

    #calculate pattern lenghts
    df$pattLength = nchar(df$pattern)
    df
  }

  test3 = patternsProb(test, test2, iter)
  test3

  #apply SSI equation ####

  applySSIequation = function(df, df1, participant, trial, alternative, attribute) {

    #change class
    df = as.data.table(df)

    #calculate string length for each trial
    df$attributeClean = ifelse(df$trial == shift(df$trial, 1L)
                               & df$alternative == shift(df$alternative, 1L)
                               & df$attribute == shift(df$attribute, 1L), 1, 0)
    df = subset(df, attributeClean != 1 | is.na(attributeClean))
    df[, "attributeClean" := NULL]
    setkey(df, "participant", "trial")
    stringLength = df[, list(N = NROW(attribute)), by = key(df)]

    #calculate numerator for SSI
    df1$numerator = df1$pattLength * df1$pattFreq * df1$probComplement
    setkey(df1, "participant", "trial")
    df2 = df1[, list(patternSum = sum(numerator)), by = key(df1)]

    #format data
    df2 = as.data.table(df2)
    df2$participant = as.numeric(df2$participant)
    df2$trial = as.numeric(df2$trial)
    df2 = df2[order(participant, trial),]

    #merge in trial string leghts
    df2 = merge(df2, stringLength, by = c("participant", "trial"), all = T)
    df2[is.na(df2)] = 0

    #compute SSI
    df2$SSI = df2$patternSum / df2$N
    df2
  }

  test4 = applySSIequation(df, test3, "participant", "trial", "alternative", "attribute")
  test4

  #save table
  write.csv(file = "SSI.csv", x = test4)
}





