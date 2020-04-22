#' epiCrop
#'
#' Crop disease epidemiology functions and datasets.
#'
#' @docType package
#' @name epiCrop
"_PACKAGE"

#The way in which you define variables in tidyverse package functions
# can cause confusion for the R CMD check, which sees column names
# and the name of your dataset, and flags them as "undefined global variables".
utils::globalVariables(c("sunset", "doy", "short_date", "sunrise", "."))
