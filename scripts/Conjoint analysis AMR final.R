#
# R-script belonging to the manuscript
# A summary index for antimicrobial resistance in food animals in the Netherlands
# Arie H. Havelaar, Haitske Graveland, Jan Van De Kassteele, Tizza P. Zomer, Kees Veldman, Martijn Bouwknegt
# Submitted to BMC Vet Res
#
# Author: Jan van de Kassteele
# Last revision date: October 3, 2016
#

#
# Init
#

# Load packages
library(openxlsx)
library(survival)
library(MASS)

#
# Read data
#

# The data are in a folder called 'data'. This folder contains three folders for this script:
# questback: Questback questionnaire data (5 files = 5 versions)
# profielen: link between version and comparison (1 file)
# vergelijkingen: link between comparison and profile (1 file)

# Read Questback questionnaire data
file.list <- list.files(path = "data/questback", full.names = TRUE)
questback.data <- NULL
for (i in 1:5) {
  tmp <- read.xlsx(xlsxFile = file.list[i], sheet = "Ruwe data")
  # Skip if there is nothing there
  if (nrow(tmp) == 0) next
  # Put data into one dataframe named questback.data
  questback.data <- rbind(questback.data, cbind(data.frame(Version = i), tmp))
}

# Read profiles
profiles.data <- read.xlsx(xlsxFile = "data/profielen/Profielen V3 24 juni.xlsx", sheet = "Onder elkaar")

# Read comparisons
comparisons.data <- rbind(
  read.xlsx(xlsxFile = "data/vergelijkingen/Opmaak vergelijkingen V3.xlsx", sheet = "Forward"),
  read.xlsx(xlsxFile = "data/vergelijkingen/Opmaak vergelijkingen V3.xlsx", sheet = "Backward"))

#
# Clean data
#

#
# Questback data
#

# Rename variables
# Remove X in front of column names in columns 4:9
names(questback.data)[4:9] <- substr(names(questback.data)[4:9], start = 5, stop = nchar(names(questback.data)[4:9]))
# Rename columns 10:39 to Q1 to Q30 (question 1 to 30)
names(questback.data)[10:39] <- paste0("Q", 1:30)

# Remove "Profiel " in Q1 t/m Q30
questback.data[, 10:39] <- as.data.frame(
  lapply(X = questback.data[, 10:39], FUN = gsub, pattern = "Profiel ", replacement = ""),
  stringsAsFactors = FALSE)
# Reshape into long format for analysis
questback.data.long <- reshape(
  data = questback.data,
  varying = paste0("Q", 1:30),
  v.names = "Profile",
  timevar = "Question",
  times = 1:30,
  idvar = "Respondent",
  new.row.names = 1:(30*nrow(questback.data)),
  direction = "long")

#
# Profiles data
#

# Give columns informative names
names(profiles.data) <- c("Version", "Comparison", "Direction")
# Add Question 1 to 30
profiles.data <- within(profiles.data, {
  Question <- 1:30
})

#
# Comparisons data
#

# Translate first two columns from Dutch to English
names(comparisons.data)[1:2] <- c("Comparison", "Profile")

# Split variable Comparison into two variables: Comparison and Direction
comparisons.data <- within(comparisons.data, {
  Direction <- substr(Comparison, start = nchar(Comparison), stop = nchar(Comparison))
  Comparison <- as.numeric(substr(Comparison, start = 1, stop = nchar(Comparison)-1))
})

#
# Merge questback data and profiles data
#

# Merge questback data and profiles data, call it amr.data
amr.data <- merge(questback.data.long, profiles.data)
# Re-order records
amr.data <- amr.data[with(amr.data, order(Version, Respondent, Question)), ]
# Check: amr.data should have the same number of rows as questback.data.long
dim(questback.data.long)
dim(amr.data)

#
# Consistency check
#

# Split to Respondent
tmp <- split(amr.data, f = amr.data$Respondent)
# For each respondent, calculate fraction consistent (max = 1)
x <- sapply(tmp, function(x) {
  # Find oud which rows are duplicated (=2*6 = 12 rows)
  index <- with(x, is.element(Comparison, Comparison[duplicated(Comparison)]))
  # Remove those rows
  data <- x[index, c("Comparison", "Direction", "Profile")]
  # Split by Comparison
  data.comparison <- split(data, f = data$Comparison)
  # If chosen profiles within comparision differ -> consistent
  is.consistent <- sapply(data.comparison, function(x) {
    with(x, Profile[1] != Profile[2])
  })
  # Calculate fraction consistent
  mean(is.consistent)
})
# Test consistency (H0: p = 1/2, Ha: p > 1/2)
print(x)
(z.val <- (mean(x)-0.5)/(sd(x)/sqrt(length(x))))
(p.val <- 1-pnorm(z.val))
# Conclusion: there is enough evidence for consistency

#
# Prepare data for conjoint analysis
#

# Conjoint analysis: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1118112/
# Stage 1: Identifying the attributes: Cip, Cef, Tet & Amp
# Stage 2: Assigning levels to the attributes: + / -
# Stage 3: Choice of scenarios: 24 profiles
# Stage 4: Establishing preferences: discrete choices A / B
# Stage 5: Data analysis: Organise data in preferred (1) and not preferred (0). Then do conditional logistic regression
# Also see http://www.xlstat.com/en/products-solutions/feature/conditional-logit-model.html

# Remove the 6 duplicates
# (is allowed to do this using the duplicated function, because of the randomness in the questionnaire data)
# Split to respondent
tmp <- split(amr.data, f = amr.data$Respondent)
# For each respondent, remove duplicates
tmp <- lapply(tmp, function(x) x[with(x, !duplicated(Comparison)), ])
# Put tmp back together in amr.data1
amr.data1 <- do.call(rbind, tmp)
# Check
dim(amr.data1) # 510 records minus 17 participants * 6 questions = 408 records
with(amr.data1, table(Comparison, Direction)) # 24 questions, row totals equal 17

# The data in its current form only shows the chosen option
# For a conditional logistic regression, we also need the not chosen option
# We therefore make our amr.data1 twice as large
amr.data2 <- rbind(
  within(amr.data1, {
    # The first part has been chosen. This is what we have observed
    Chosen <- 1
  }),
  within(amr.data1, {
    # The second part has not been chosen
    Chosen <- 0
    # For the not-chosen part, the profile should of course be swapped A -> B and B -> A
    Profile <- ifelse(Profile == "A", yes = "B", no = "A")
  })
)

# Add resistance profiles to amr.data2
amr.data2 <- merge(amr.data2, comparisons.data)

# Create strata on which conditionling will take place
# In this case each combination of Respondent and Comparison
amr.data2 <- within(amr.data2, strata <- interaction(Respondent, Comparison))

#
# Model
#

# Model: conditional logistic regression on Chosen
amr.mod <- clogit(Chosen ~ Ciprofloxacine + Cefotaxim + Tetracycline + Ampicilline + strata(strata), data = amr.data2)
summary(amr.mod)

# Determine weights
# The coefficients are the contribution to the utility score
# The weights are the normalized coefficients (add up to 1)
# Determine confidence intervals by simulation
coef.sim <- mvrnorm(n = 1000, mu = coef(amr.mod), Sigma = vcov(amr.mod))
w <- coef.sim/rowSums(coef.sim)
w.mean <- colMeans(w)
w.conf <- apply(w, MARGIN = 2, FUN = quantile, prob = c(0.025, 0.975))
# The result (Table 3 in the paper)
round(cbind(w.mean, t(w.conf)), digits = 3)

# Barplot of weights
par(mar = c(2.5, 2.5, 0.5, 0.5))
bar.mid <- barplot(w.mean, width = 1, space = 0.1, ylim = c(0, max(w.conf)))
arrows(bar.mid, w.conf[1, ], bar.mid, w.conf[2, ], code = 3, length = 0.1, angle = 90)

