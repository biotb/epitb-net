# Create box plot of MI scores.

library(ggplot2)

input_full_filepath <- "spydrpick.outliers"
suffix="-outliers"

# Read SpydrPick file
input <- read.csv(input_full_filepath, header = FALSE, sep = " ")
# MI: input[input$V5]
# Genome distance: input[input$V3]
summary(input$V5)
summary(input$V3)

maxmi <- max(input$V5)
minmi <- min(input$V5)

hline = 0
extra_y = 0

png(file = paste("boxplot-MI-", suffix, ".png", sep=""), width = 1600, height = 1400, res = 300)

# Plot the chart.
ggplot(input, aes(x="", y=V5)) + 
  geom_boxplot(outlier.size = 0.5) +
  labs(x="SpydrPick outliers", y = "Mutual information") +
  theme_light()
  if (hline > 0) {
    + geom_hline(yintercept=hline, color = "red")
    + scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                         breaks = c(seq(minmi, maxmi + extra_y, by=0.05), hline),
                         expand = c(0, 0), 
                         limits=c(minmi, maxmi + extra_y))
  } else {
    + scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                         breaks = c(seq(minmi, maxmi + extra_y, by=0.05)),
                         expand = c(0, 0), 
                         limits=c(minmi, maxmi + extra_y))
  }

dev.off()

