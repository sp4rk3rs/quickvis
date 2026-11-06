#### A real, published dataset for testing the general-factor diagram
#### Source: https://www.personality-project.org/r/psych/help/bifactor.html
#### Thurstone is a nine-test correlation matrix, based on N = 213 people.
#### It is summary data, not 213 rows of individual scores.
library(quickvis)
library(lavaan)
data("Thurstone", package = "psych")

#### Shorten a few test names to keep the diagram readable
thurstone_r <- Thurstone
dimnames(thurstone_r) <- rep(list(c(
  "Sentences", "Vocabulary", "Completion", "FirstLetters", "FourLetter",
  "Suffixes", "Series", "Pedigrees", "Groups"
)), 2)

#### Every test loads on general ability g and one specific ability factor
thurstone_model <- '
  g =~ Sentences + Vocabulary + Completion + FirstLetters + FourLetter +
       Suffixes + Series + Pedigrees + Groups

  Verbal =~ Sentences + Vocabulary + Completion
  Fluency =~ FirstLetters + FourLetter + Suffixes
  Reasoning =~ Series + Pedigrees + Groups
'

#### Fit summary data with sample.cov; these are not raw observations
thurstone_fit <- cfa(
  thurstone_model,
  sample.cov = thurstone_r,
  sample.nobs = 213,
  sample.cov.rescale = FALSE,
  std.lv = TRUE,
  orthogonal = TRUE
)

#### Unit latent variances and zero factor covariances define this bifactor model
stopifnot(
  lavInspect(thurstone_fit, "converged"),
  lavInspect(thurstone_fit, "post.check"),
  all(diag(lavInspect(thurstone_fit, "theta")) > 0)
)
print(fitMeasures(thurstone_fit, c("cfi", "tli", "rmsea", "srmr")))

#### No general or layout argument: g should be detected automatically
thurstone_diagram <- path_diagram(thurstone_fit, std = "all", digits = 2)
stopifnot(thurstone_diagram$layout == "bifactor")

#### Further calls to try interactively after running this example:
# path_diagram(thurstone_fit, variance = "Sentences", tracing = TRUE)
# path_diagram(thurstone_fit, covariance = c("Sentences", "Vocabulary"), tracing = TRUE)
# path_diagram(thurstone_fit, std = "lv", general = "g")
