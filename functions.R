library(tidyverse)
library(here)
library(fs)


# Create the likelihood function for mutect

model_likelihood <- function(depth, alt_count, qscore=30){
  # Set it to have an average qscore of 30
  ref_count <- depth - alt_count
  f <- alt_count/depth # alt allele frequency
  e <- 10**(-qscore/10)
  ref_prob <- f*e/3. + (1-f)*(1-e)
  alt_prob <- f*(1-e) + (1-f)*e/3.
  l <- ref_prob**ref_count * alt_prob**alt_count
  l
}

null_likelihood <- function(depth, alt_count, qscore=30){
  # Set it to have an average qscore of 30
  ref_count <- depth - alt_count
  f <- 0 # alt allele frequency
  e <- 10**(-qscore/10)
  ref_prob <- f*e/3. + (1-f)*(1-e)
  alt_prob <- f*(1-e) + (1-f)*e/3.
  l <- ref_prob**ref_count * alt_prob**alt_count
  l
}
error_num_expected_wgs <- function(alt,depth){
  log10(sum(dbinom(alt,depth,0.00027))*3e9)
}

error_num_expected_wes <- function(alt,depth){
  p <- sum(dbinom(alt,depth,0.00027))
  log10(p*3e7)
}

truth_expected_wgs <- function(vaf){
  log10(300/vaf)
}

truth_expected_wes <- function(vaf){
  log10(3/vaf)
}

#####
# Create a function for the slope and gradient of the FP/TP ratio.
# TP/FP = #expected False positives at a given depth and vaf/# expected True positives
# Where the expected number of false positivies is 