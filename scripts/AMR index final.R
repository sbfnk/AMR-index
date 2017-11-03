# R-script belonging to the manuscript A summary index for antimicrobial resistance in food animals in the Netherlands Arie Hendrik Havelaar,
# Arie H. Havelaar, Haitske Graveland, Jan Van De Kassteele, Tizza P. Zomer, Kees Veldman, Martijn Bouwknegt Submitted to BMC Vet Res Author: Arie Havelaar
# Last revision date: October 3, 2016
#
# Init
#

# Load packages
library(openxlsx)
library(reshape2)
library(ggplot2)

#
# Read data
#

# The data are in a folder called 'data/simulatie'
# Monitoring data are in an Excel file called Resistance data.xlsx
res.data <- read.xlsx(xlsxFile = "data/simulatie/Resistance data.xlsx", sheet = "Data")

# Weights are in a data frame called w_sim.csv
weights <- read.table("data/simulatie/w_sim.csv", header = TRUE, sep = ",")

#
# Some data modifications
#

# Set antibiotic names
ab.names <- c("Ampicillin", "Cefotaxime", "Ciprofloxacin", "Tetracycline")

# In res.data, rename AmpR, CefR, CipR abd TetR to ab.names
names(res.data) # Check before
names(res.data)[4:7] <- ab.names # Rename
names(res.data) # Check after: OK

# In weights idem dito. Mind the order of the columns!
names(weights) # Check before
weights <- weights[, c(4, 2, 1, 3)] # Reorder columns
names(weights) # Check intermediate result: OK
names(weights) <- ab.names # Rename
names(weights) # Check after: OK

# In res.data, make Animal.species a categorical variable
# This will also set the sequence of the animal species for plotting
res.data <- within(res.data, Animal.species <- factor(Animal.species,
  levels = c("Broiler chickens", "Fattening pigs", "Veal calves", "Dairy cows"),
  labels = c("Broiler chickens", "Fattening pigs", "Veal calves", "Dairy cows")))

#
# Simulate uncertainty in antibiotics resistance using a Beta distribution
#

# Number of simulations is equal to the number of rows of weight simulations
n.sim <- nrow(weights)
# Number of antibiotics is equal to the number of columns of weight simulations
n.ab <- ncol(weights)

# For each combination of year (6) and animal species (4), create n.sim realisations of the n.ab antibiotic resistance prevelances
# Each row of res.data gets expanded n.sim times
res.data.sim <- with(res.data, expand.grid(
  sim = 1:n.sim,
  Year = unique(Year),
  Animal.species = levels(Animal.species)))

# Allocate empty matrix to fill later on
res.data.sim <- cbind(
  res.data.sim,
  matrix(0, nrow = nrow(res.data.sim), ncol = n.ab, dimnames = list(NULL, ab.names)))

# For each row in res.data
for (i in 1:nrow(res.data)) {
  # Create n.sim simulations of the n.ab antibiotic resistance prevelances
  res.data.sim[((i-1)*n.sim + 1):(i*n.sim), ab.names] <- matrix(
    rbeta(
      n = n.sim*n.ab,
      shape1 = as.numeric(res.data[i, ab.names]) + 1,
      shape2 = as.numeric(res.data[i, "Isolates"]) - as.numeric(res.data[i, ab.names]) + 1),
    nrow = n.sim, byrow = TRUE)
}

# Caclulate means of simulation results
res.data.mean <- with(res.data.sim, aggregate(
  cbind(Ampicillin, Cefotaxime, Ciprofloxacin, Tetracycline) ~ Year + Animal.species,
  FUN = mean))

#
# Plot results for all animal species with separate lines for each antibiotic compound
#

# Convert res.data.mean to long format
res.data.long <- melt(res.data.mean,
  id.vars = c("Animal.species", "Year"),
  variable.name = "Antibiotic",
  value.name = "Prevalence")

# Plot
fig2 <- ggplot(data = res.data.long, mapping = aes(x = Year, y = Prevalence, group = Antibiotic, colour = Antibiotic)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Prevalence of resistance") +
  facet_wrap(facets = ~ Animal.species)
fig2

#
# Combine samples from resistance distributions with samples from weights distribution
#

# Split res.data.sim by Year x Animal species combination
res.data.split <- with(res.data.sim, split(x = res.data.sim, f = list(Year, Animal.species)))

# For each dataframe in res.data.split, compute AMR index
amr.data.split <- lapply(X = res.data.split, FUN = function(data) {
  # Calculate weighted prevelance (AMR index)
  amr.sim <- rowSums(weights*data[, ab.names])
  # Calculate statistics
  data.frame(
    amr.m = mean(amr.sim),
    amr.l = quantile(amr.sim, prob = 0.025),
    amr.u = quantile(amr.sim, prob = 0.975))
})

# Combine results into amr.data dataframe
amr.data <- do.call(what = rbind, args = amr.data.split)
amr.data <- cbind(res.data[, c("Year", "Animal.species")], amr.data)

#
# Plot AMR index for each animal species separately
#

# Plot
fig3 <- ggplot(data = amr.data, mapping = aes(x = Year, y = amr.m, ymin = amr.l, ymax = amr.u)) +
  geom_line(size = 1, colour = "green4") +
  geom_ribbon(alpha = 0.1) +
  labs(x = "Year", y = "Prevalence of resistance") +
  facet_wrap(facets = ~ Animal.species)
fig3

#
# Export figures
#

ggsave(file = "figure2.pdf", plot = fig2, device = "pdf", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(file = "figure3.pdf", plot = fig3, device = "pdf", width = 297, height = 210, units = "mm", dpi = 300)
