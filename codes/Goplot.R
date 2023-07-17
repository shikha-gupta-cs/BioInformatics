# to plot go terms
library(ggplot2)

# with columns 'Term', 'Count', '%', and 'PValue'
setwd("/Users/rishav_raj/Desktop/bioinformatics/geo/DEG2")
go_data <- read.csv("3678genesMF.csv",header = TRUE)
dim(go_data)
go_data <- go_data[go_data$Count <= 60, ]
dim(go_data)
# Transform the PValue column to -log10(p value)
go_data$NegLogPValue <- -log10(go_data$PValue)

# Sort the data by the NegLogPValue column in ascending order
go_data <- go_data[order(-go_data$NegLogPValue), ]

# Create the plot using ggplot2
plot <- ggplot(go_data, aes(x = Count, y = Term, fill = NegLogPValue)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "blue", na.value = NA) +
  labs(x = "Gene Count", y = "GO Terms", fill = "-log10(p value)")

#plot <- plot + theme(axis.text.y = element_text(angle = 0))
plot <- plot + facet_wrap(~ Category, ncol = 1,scales="free_y")


# Display the plot
print(plot)

