#' List experiment regarding gender and politics
#'
#' A dataset containing responses to a list experiment and direct question
#' regarding the statement ``women are as competent as men in politics.''
#' Data also contain socio-demographics for gender, age, education, region,
#' and political ideology.
#'
#' @format A data frame with 5000 rows and 7 variables:
#' \describe{
#'   \item{y}{response to the list experiment}
#'   \item{treatment}{treatment assignment}
#'   \item{direct}{response to the direct question}
#'   \item{gender}{gender of respondent \{Woman, Man\}}
#'   \item{age}{age of respondent \{18, 19, ..., 94\}}
#'   \item{ageGroup}{age group of respondent \{18-29, 30-39, 40-49, 50-64, 65+\}}
#'   \item{education}{education of respondent \{High school or below, College, University degree\}}
#'   \item{motherTongue}{mother tongue of respondent \{English, French, Other language\}}
#'   \item{region}{region of respondent \{Ontario, Atlantic, Quebec, West\}}
#'   \item{selfPlacement}{political ideology of respondent (0 = right-wing, 10 = left wing) \{0, 1, ..., 10\}}
#'   \item{weight}{survey weight}
#' }
#' @source Eady, Gregory. 2016 "The Statistical Analysis of Misreporting on Sensitive Survey Questions."
"gender"

# setwd("/Users/gregoryeady/Dropbox (Personal)/Software/misreport")
# gender <- read.csv("data-raw/gender.csv", stringsAsFactors = FALSE)
# gender <- subset(gender, attentionCheck == 1 & !is.na(weight))
# names(gender)[names(gender) == "listGender"] <- "y"
# names(gender)[names(gender) == "listGenderTreatment"] <- "treatment"
# names(gender)[names(gender) == "directGender"] <- "direct"
# row.names(gender) <- 1:nrow(gender)
# gender <- gender[, -which(names(gender) %in% c("attentionCheck"))]
# gender$gender <-       factor(gender$gender, levels = c("Man", "Woman"))
# gender$ageGroup <-     factor(gender$ageGroup, levels = c("18-29", "30-39", "40-49", "50-64", "65+"))
# gender$education <-    factor(gender$education, levels = c("High school or below", "College", "University degree"))
# gender$motherTongue <- factor(gender$motherTongue, levels = c("English", "French", "Other language"))
# gender$region <-       factor(gender$region, levels = c("Ontario", "Atlantic", "Quebec", "West"))
# devtools::use_data(gender, pkg = "./data", overwrite = TRUE)