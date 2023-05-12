library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

num_to_log10scale <- function (x) {
    parts <- strsplit(x, split = " ")[[1]]
    x <- as.numeric(parts[1])
    # assumes that the 2nd factor is 10 to some power
    y <- as.numeric(substr(parts[3], 5, nchar(parts[3]) - 1))
    log10(x) + y
}

interval_to_log10scale <- function(x) {
    # assumes that x is 10 to some power; the function just extracts the exponent
    # exponent must be between { }
    extremes <- strsplit(substr(x, 2, nchar(x) - 1), split = ",")[[1]]
    start <- as.numeric(gsub("^.*\\{(\\d+)\\}$", "\\1", extremes[1]))
    end <- as.numeric(gsub("^.*\\{(\\d+)\\}$", "\\1", extremes[2]))
    paste("[", start, ",", end, "]", sep = "")
}

data <- read.csv("../data/scenarios.csv", sep = ";")

# Convert Scenarios and Intermediates to log10 scale
data["Intermediates"] <- apply(data["Intermediates"], 1, num_to_log10scale)
data[1:5,]["Scenarios"] <- apply(data[1:5,]["Scenarios"], 1, num_to_log10scale)
data[6:nrow(data),]["Scenarios"] <- apply(data[6:nrow(data),]["Scenarios"], 1, interval_to_log10scale)

# Convert the range values in column "y" to separate minimum and maximum columns
data$scen_min <- as.numeric(gsub("\\[(\\d+),.+", "\\1", data$Scenarios))
data$scen_max <- as.numeric(gsub("\\[.+,(\\d+)]", "\\1", data$Scenarios))

# Sort Genomes by rank distance
data$Genome <- factor(data$Genome, levels = data[order(data$d), "Genome"])


p1 <- ggplot(data, aes(x = Intermediates, y = Genome)) +
    geom_point(size = 2) +
    labs(x = "Number of intermediate genomes vs. Human (log scale)",
         y = "Genome")
p2 <- ggplot(data, aes(x = d, y = Genome)) +
    geom_point(size = 2) +
    labs(x = "Rank distance vs. Human",
         y = "Genome")
p3 <- ggplot(data, aes(x = (scen_min + scen_max) / 2, y = Genome)) +
    geom_point(size = 2) +
    geom_errorbar(aes(xmin = scen_min, xmax = scen_max), width = 0.4) +
    labs(x = "Number of scenarios vs. Human (log scale)",
         y = "Genome")
p4 <- plot_grid(p1, p2, p3, nrow = 3, align = "v")

png("../data/img/ints_vs_human.png")
print(p1)
dev.off()

png("../data/img/d_vs_human.png")
print(p2)
dev.off()

png("../data/img/scen_vs_human.png")
print(p3)
dev.off()

png("../data/img/plots_merged.png")
print(p4)
dev.off()