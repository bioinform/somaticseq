#!/usr/bin/env Rscript

require("ada")

args <- commandArgs(TRUE)

training_data_filename = args[1]

##### Main (entry point)
train_filename = paste(training_data_filename)
train_data = read.table(train_filename, header=TRUE)

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

for (var_i in tail(args, -1) ) {
    train_data[, var_i] <- NULL
    cat("Remove", var_i, "\n")
}

model_formula <- as.formula(TrueVariant_or_False ~ .)

print("Fitting model...")

boosting_iters = 500

seed_value = floor(runif(1, min=100, max=50000))
print( paste("Seed =", seed_value) )
set.seed(seed_value)

ada.model <- ada(model_formula, data = train_data, iter = boosting_iters, control=rpart.control(cp=-1, maxdepth=16, minsplit=0, xval=0))
save(ada.model, file = paste(training_data_filename, ".ada.Classifier.RData", sep="") )

print(ada.model)
