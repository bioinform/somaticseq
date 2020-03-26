#!/usr/bin/env Rscript

library(xgboost)
library(Matrix)
library(caret)
library(data.table)

#sessionInfo()

args <- commandArgs(TRUE)

training_data_filename = args[1]

##### Main (entry point)

train_filename = paste(training_data_filename)
train_data = read.table(train_filename, na.strings=c("NaN", "nan", "<NA>"), header=TRUE, stringsAsFactors=FALSE)
train_data[is.na(train_data)] <- 0

if (!(1 %in% train_data$TrueVariant_or_False && 0 %in% train_data$TrueVariant_or_False)) {
stop("In training mode, there must be both true positives and false positives in the call set.")
}

# Use substitution identity for training
train_data$GC2CG = 0
train_data$GC2TA = 0
train_data$GC2AT = 0
train_data$TA2AT = 0
train_data$TA2GC = 0
train_data$TA2CG = 0

train_data$GC2CG[ (train_data$REF=='G' & train_data$ALT=='C') | (train_data$REF=='C' & train_data$ALT=='G') ] = 1
train_data$GC2TA[ (train_data$REF=='G' & train_data$ALT=='T') | (train_data$REF=='C' & train_data$ALT=='A') ] = 1
train_data$GC2AT[ (train_data$REF=='G' & train_data$ALT=='A') | (train_data$REF=='C' & train_data$ALT=='T') ] = 1
train_data$TA2AT[ (train_data$REF=='T' & train_data$ALT=='A') | (train_data$REF=='A' & train_data$ALT=='T') ] = 1
train_data$TA2GC[ (train_data$REF=='T' & train_data$ALT=='G') | (train_data$REF=='A' & train_data$ALT=='C') ] = 1
train_data$TA2CG[ (train_data$REF=='T' & train_data$ALT=='C') | (train_data$REF=='A' & train_data$ALT=='G') ] = 1

# Do not use these for training
train_data$CHROM      <- NULL
train_data$POS        <- NULL
train_data$ID         <- NULL
train_data$REF        <- NULL
train_data$ALT        <- NULL
train_data$if_COSMIC  <- NULL
train_data$COSMIC_CNT <- NULL
train_data$T_VAF_REV  <- NULL
train_data$T_VAF_FOR  <- NULL

train_data$Strelka_QSS  <- NULL
train_data$Strelka_TQSS <- NULL

for (var_i in tail(args, -1) ) {
    train_data[, var_i] <- NULL
    cat("Remove", var_i, "\n")
}

# Model building
print("Fitting model...")

train_label <- train_data[, "TrueVariant_or_False"]
train_data  <- train_data[, names(train_data) != "TrueVariant_or_False"]

train_label <- lapply(train_label, as.numeric)
train_data  <- sapply(train_data, as.numeric)

train_data.matrix <- as.matrix(train_data)

# Set deterministic seed for reproducibility on variable removal test
boosting_iters = 500
seed_value     = floor(runif(1, min=100, max=50000))

print( paste("Seed =", seed_value) )


set.seed(seed_value)
bst <- xgboost(data = train_data.matrix, label = train_label, nround=boosting_iters, objective="binary:logistic", verbose=0)
class(bst)

saveFilename = paste(training_data_filename, ".xgboost.Classifier.RData", sep="")
xgb.save(bst, fname=saveFilename)

#xgb_bst <- xgb.load(paste(saveFilename))
#print("Most important features (look at column Gain):")
#imp_matrix<- xgb.importance(model=bst)
#print(head(imp_matrix))
