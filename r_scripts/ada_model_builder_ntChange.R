#!/usr/bin/env Rscript

require("ada")

args <- commandArgs(TRUE)

training_data_filename = args[1]

##### Main (entry point)
train_filename = paste(training_data_filename)
train_data = read.table(train_filename, header=TRUE)

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
train_data$TA2CG[ (train_data$REF=='T' & train_data$ALT=='T') | (train_data$REF=='A' & train_data$ALT=='G') ] = 1

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

model_formula <- as.formula(TrueVariant_or_False ~ .)

print("Fitting model...")
ada.model <- ada(model_formula, data = train_data, iter = 500)

save(ada.model, file = paste(training_data_filename, ".ntChange.Classifier.RData", sep="") )

print(ada.model)
