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

for (var_i in tail(args, -1) ) {
    train_data[, var_i] <- NULL
    cat("Remove feature:", var_i, "\n")
}


model_formula <- as.formula(TrueVariant_or_False ~ .)

# Cross validation:

for (ith_try in 1:10)

{
    # split test/train 50-50
    sample <- sample.int(n = nrow(train_data), size = floor(.5*nrow(train_data)), replace = F)
    train  <- train_data[sample, ]
    test   <- train_data[-sample, ]
    
    # do model
    ada.model <- ada(model_formula, data = train, iter = 500)
    # print(ada.model)

    ada.pred <- predict(ada.model, newdata = test, type="both", n.iter=350)
    
    # probability > 0.5
    pass_calls <- ada.pred$prob[,2] > 0.5
    reject_calls <- ada.pred$prob[,2] < 0.1
    
    # Counting
    num_pass_calls <- sum( pass_calls )
    num_reject_calls <- sum( reject_calls )
    num_pass_true_positives <- sum( pass_calls[pass_calls == test$TrueVariant_or_False] )
    num_true_positives <- sum(test$TrueVariant_or_False)
    
    # Calculate results
    precision <- num_pass_true_positives/num_pass_calls
    sensitivity <- num_pass_true_positives/num_true_positives
    F1_score <- 2 * num_pass_true_positives / ( num_true_positives + num_pass_calls )
    
    # Print out
    cat (ith_try, 'th_try', '\n')
    
    cat("PASS_Calls =",          num_pass_calls, "\n")
    cat("REJECT_Calls =",        num_reject_calls, "\n")
    
    cat("PASS_TruePositives =",  num_pass_true_positives, "\n")
    cat("PASS_FalsePositives =", num_pass_calls - num_pass_true_positives, "\n")
    
    cat("Sensitivity =",         sensitivity, "\n")
    cat("Precision =",           precision, "\n")
    cat("F1 =",                  F1_score, "\n")

}
