#!/usr/bin/env R

require("ada")

args <- commandArgs(TRUE)

training_data_filename = args[1]

##### Main (entry point)

train_filename = paste(training_data_filename)

train_data = read.table(train_filename, header=TRUE)

if (FALSE) {
print("Updating missing values...")
train_data <- set_missing_values(train_data)

}else {

train_data <- train_data[,-c(1, 2, 3, 4, 5)]

}


train_data[,'REF'] <- NULL
train_data[,'ALT'] <- NULL


if (FALSE) {
model_formula <- as.formula( TrueVariant_or_False ~
                                ((if_MuTect_ + if_VarScan2_ + if_JointSNVMix2_ + if_SomaticSniper_ +
                                if_MuTect_if_JointSNVMix2 + if_MuTect_if_SomaticSniper + if_JointSNVMix2_if_SomaticSniper +
                                if_MuTect_if_JointSNVMix2_if_SomaticSniper) +
                                (SNVMix2_Score + Sniper_Score +
                                if_dbsnp + BAF + COMMON + G5 + G5A +   # Probably no need for G5/G5A
                                N_VAQ +
                                T_VAQ + T_MQ0  + T_MLEAF +
                                N_StrandBias + N_BaseQBias + N_MapQBias + N_TailDistBias +
                                T_StrandBias + T_BaseQBias + T_MapQBias + T_TailDistBias +
                                N_AMQ_REF + N_AMQ_ALT + N_BQ_REF + N_BQ_ALT + N_MQ +
                                T_AMQ_REF + T_AMQ_ALT + T_BQ_REF + T_BQ_ALT + T_MQ +
                                #N_DP_large +
                                T_DP_small + T_DP_large +
                                N_ALT_FREQ_FOR + N_ALT_FREQ_REV + N_ALT_STRAND_BIAS +
                                T_ALT_FREQ_FOR + T_ALT_FREQ_REV + T_ALT_STRAND_BIAS )))
} else {
model_formula <- as.formula(TrueVariant_or_False ~ .)
}


print("Fitting model...")
ada.model <- ada(model_formula, data = train_data, iter = 1000)

save(ada.model, file = paste(training_data_filename, ".Classifier.RData", sep="") )


print(ada.model)

pdf( paste(training_data_filename, ".varplot.pdf", sep = "") )
varplot(ada.model)
dev.off()

pdf( paste(training_data_filename, ".iterplot.pdf", sep = "") )
plot(ada.model, TRUE, TRUE)
dev.off()

#print("Computing prediction values...")
#ada.pred <- predict(ada.model, newdata = test_data, type="both")
