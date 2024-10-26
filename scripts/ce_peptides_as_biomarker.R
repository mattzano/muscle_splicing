# Install packages if not already installed
install.packages(c("pROC", "caret", "glmnet"))

# Load necessary libraries
library(caret)     # For model training and cross-validation
library(glmnet)    # For elastic net regression (or lasso/ridge)
library(pROC)      # For ROC curve analysis

# Sample Data (Gene expression and binary outcome)
set.seed(123) # For reproducibility
data <- data.frame(
  Gene1 = rnorm(100, 10, 5),  # Simulated gene expression data
  Gene2 = rnorm(100, 12, 6),
  Gene3 = rnorm(100, 8, 4),
  Outcome = sample(c(0, 1), 100, replace = TRUE) # Binary outcome (e.g., disease/control)
)

data <- normed_counts[,c(2:41,43)] %>% 
  #mutate(Outcome = ifelse(grepl("IBM", V4) | grepl("PT", V4), "IBM", 
  #                                 ifelse(grepl("CTR", V4))
  filter(symbol != "") %>% 
  distinct(symbol, .keep_all = T) %>% 
  column_to_rownames("symbol") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("sample_name") %>% 
  mutate(Outcome = ifelse(grepl("IBM", sample_name) | grepl("PT", sample_name), "IBM", 
                          ifelse(grepl("CTR", sample_name) | grepl("Ctrl", sample_name) | grepl("Control", sample_name), "Control",NA))) 
  

data <- splicing_bedparse_invivo_normalised %>% 
  mutate(value = V5/paired.reads,
         name = V7,
         Outcome = condition) %>%
  dplyr::select(c(value,name,Outcome,V4)) %>% 
  pivot_wider() %>% 
  column_to_rownames("V4") %>% 
  relocate("Outcome", .after = last_col())

data[is.na(data)] = 0

data
  

# Split into training and test sets
trainIndex <- createDataPartition(data$Outcome, p = .8, list = FALSE, times = 1)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# Step 1: Model Building using Elastic Net
x_train <- as.matrix(trainData[, 1:3])  # Genes as predictors
y_train <- trainData$Outcome            # Binary outcome

# Fit elastic net (combination of lasso and ridge)
model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5)

# Step 2: Predictions on the test data
x_test <- as.matrix(testData[, 1:3])  # Genes as predictors in test set
y_test <- testData$Outcome            # Binary outcome in test set

# Predict the probabilities for the test set
pred_probs <- predict(model, newx = x_test, type = "response", s = "lambda.min")

# Step 3: ROC Curve for gene combination
roc_curve <- roc(y_test, pred_probs)

# Plot the ROC curve
plot(roc_curve, col = "blue", lwd = 2, main = "ROC Curve for Gene Panel")

# Calculate and display AUC
auc_value <- auc(roc_curve)
cat("AUC for the model:", auc_value, "\n")


# Step 4: Evaluate different gene combinations
# Gene list
genes <- c("SLC24A3", "HDGFL2", "AARS1", "MYO18A", "XPO4", "STMN2", "PHF2", "NECAB2", "ZNF423", 
           "EPB41L4A", "PXDN", "SYNJ2", "PTPRZ1", "ARHGAP22")

# Function to generate all combinations for a given size
get_combinations <- function(genes, size) {
  combn(genes, size, simplify = FALSE)
}

# Generate all possible combinations for sizes from 1 to the total number of genes
gene_combinations <- unlist(lapply(1:length(genes), function(size) get_combinations(genes, size)), recursive = FALSE)

# Print the combinations (optional)
print(gene_combinations)

# Gene list
genes <- c("SLC24A3", "HDGFL2", "AARS1", "MYO18A", "XPO4", "STMN2", "PHF2", "NECAB2", "ZNF423", 
           "EPB41L4A", "PXDN", "SYNJ2", "PTPRZ1", "ARHGAP22")

# Function to generate all combinations for a given size
get_combinations <- function(genes, size) {
  combn(genes, size, simplify = FALSE)
}

# Generate combinations with size >= 2
gene_combinations <- unlist(lapply(2:length(genes), function(size) get_combinations(genes, size)), recursive = FALSE)

# Print the combinations (optional)
print(gene_combinations)


for (genes in gene_combinations) {
  #genes = "HDGFL2"
  # Select genes for the current combination
  x_train <- as.matrix(trainData[, genes, drop = FALSE])
  x_test <- as.matrix(testData[, genes, drop = FALSE])
  
  # Train elastic net model
  model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5)
  
  # Predict probabilities for the test set
  pred_probs <- predict(model, newx = x_test, type = "response", s = "lambda.min")
  
  # Plot ROC curve for this combination
  roc_curve <- roc(y_test, pred_probs)
  plot(roc_curve, col = sample(colors(), 1), lwd = 2, add = TRUE, main = "ROC Curves for Different Gene Combinations")
  
  # Calculate and print AUC
  auc_value <- auc(roc_curve)
  cat("AUC for combination", paste(genes, collapse = ", "), ":", auc_value, "\n")
}




