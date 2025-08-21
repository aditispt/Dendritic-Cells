# Load necessary libraries
library(ggplot2)
library(reshape2)
library(GGally)
library(ggpubr)

install.packages("gridExtra")
library(gridExtra)

##XCR1##
# Data preparation
data <- data.frame(
  Phases = c("Early follicular", "Pre Ovulatory", "Mid Luteal"),
  XCR1_Fold_Change = c(0.3104687225, 0.2349100396, 0.3925018908),
  Progesterone_mg = c(1, 4, 25),
  Estrogen = c(36, 380, 250),
  Testosterone = c(144, 171, 126)
)

# Calculate Pearson correlation coefficients
correlation_matrix <- cor(data[, -1])

# Print correlation matrix
print(correlation_matrix)

##To find significance##
# Perform correlation tests and print results
cor_test_results <- list()
variables <- colnames(data)[-1]

for (var in variables[-1]) {
  test_result <- cor.test(data$XCR1_Fold_Change, data[[var]], method = "pearson")
  cor_test_results[[var]] <- test_result
  cat("Correlation between XCR1_Fold_Change and", var, ":\n")
  print(test_result)
  cat("\n")
}

##A coefficient close to 1 indicates a strong positive correlation, 
#close to -1 indicates a strong negative correlation, 
#and close to 0 indicates no correlation.


##Interpretation
#Strong positive correlation: XCR1 Fold Change with Progesterone (0.8168632)
#Strong negative correlation: XCR1 Fold Change with Testosterone (-0.9903998), Progesterone with Testosterone (-0.7292846)
#Weak to moderate correlations: XCR1 Fold Change with Estrogen (-0.3521110), Progesterone with Estrogen (0.2522636), Estrogen with Testosterone (0.4781103)

# Heatmap
cor_melted <- melt(correlation_matrix)
ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  ggtitle("Correlation Heatmap")


