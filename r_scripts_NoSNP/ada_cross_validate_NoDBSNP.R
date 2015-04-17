### To start, set the input filenames in "Main (entry point)" section.

require("ada")

args <- commandArgs(TRUE)

type = args[1]
numTrueCalls = as.integer( args[2] )


##### Main (entry point)
# Train and test filenames

train_filename = paste("DC", type, "tsv", sep = ".")
test_filename  = paste("DC", type, "tsv", sep = ".")

#train_filename <- ""
test_filename <- ""


# If one filename is set, data is splitted to test/train randomly
#data_filename <- "/net/kodiak/volumes/delta/shared/prj/dream_challenge/stage3.5/take3/SETTING_pileup.q25.Q20_SNVMix2p0.05_sniper.q25.s1e-4_noCOSMIC/Variant_Calls_Merged/DC_Comparison/data_for_training_noindel_norescale_0628_AtLeastOneCaller.tsv"
data_filename <- train_filename
train_filename = ""
test_filename = ""

print(test_filename)
print(train_filename)


if (data_filename != "") {
    data <- read.table(data_filename, header=TRUE)
    
    index <- 1:nrow(data)
    train_index <- sample(index, trunc(length(index)/2))
    
    test_data_ <- data[-train_index,]
    train_data_ <- data[train_index,]
    
}

if (test_filename != "") {
    test_data_ = read.table(test_filename, header=TRUE)
}

if (train_filename != "") {
    train_data_ = read.table(train_filename, header=TRUE)
}

train_data <- train_data_
test_data <- test_data_

if (FALSE) {
print("Updating missing values...")
train_data <- set_missing_values(train_data)
test_data <- set_missing_values(test_data)

}else {

train_data <- train_data[,-c(1, 2, 3, 4, 5)]
test_data <- test_data[,-c(1, 2, 3, 4, 5)]

}
test_data[,'TrueVariant_or_False'] <- NULL

if (FALSE) {
print("Updating features...")
train_data <- update_features(train_data, 0.3)
test_data <- update_features(test_data, 0.3)
}


train_data[,'REF'] <- NULL
train_data[,'ALT'] <- NULL

test_data[,'REF'] <- NULL
test_data[,'ALT'] <- NULL


# Do not use dbsnp information
train_data[,'if_dbsnp'] <- NULL
train_data[,'BAF'] <- NULL
train_data[,'COMMON'] <- NULL
train_data[,'G5'] <- NULL
train_data[,'G5A'] <- NULL

test_data[,'if_dbsnp'] <- NULL
test_data[,'BAF'] <- NULL
test_data[,'COMMON'] <- NULL
test_data[,'G5'] <- NULL
test_data[,'G5A'] <- NULL



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
ada.model <- ada(model_formula, data = train_data, iter = 500)

print(ada.model)

#pdf("varplot.pdf")
#varplot(ada.model)
#dev.off()
#pdf("iterplot.pdf")
#plot(ada.model, TRUE, TRUE)
#dev.off()

print("Computing prediction values...")
ada.pred <- predict(ada.model, newdata = test_data, type="both")



# Print results out:
if (TRUE) {
for (threshold in seq(0,1, .01)) {

    cat("threshold: ", threshold, "\t")
    ada_predicted_call <- ada.pred$prob[,2] > threshold

    # Sensitivity

    # numTrueCalls <- 14194  # stage4 indel
    # numTrueCalls <- 8292   # stage3 indel

    # numTrueCalls <- 16268  # stage4 snv
    # numTrueCalls <- 7903/2+0.5   # stage3 snv
    # numTrueCalls <- 4332   # stage2 snv
    # numTrueCalls <- 3537   # stage1 snv

    num_true_positives_predicted <- sum(ada_predicted_call[ada_predicted_call == test_data_[,'TrueVariant_or_False']])
    num_all_positive_predictied  <- sum(ada_predicted_call)

    Sensitivity <- num_true_positives_predicted / numTrueCalls
    cat("Recall: ", Sensitivity , "\t")

    # Specificity
    Specificity <- num_true_positives_predicted / num_all_positive_predictied
    cat("Precision: ", Specificity, "\t")


    cat("DREAM_Accuracy: ", (Specificity + Sensitivity)/2, "\t")
    cat("F1: ", (  2 * num_true_positives_predicted / ( numTrueCalls + num_all_positive_predictied )  ), "\n")
}
}





# Write predicted values to output file
if(FALSE){
    test_data_output <- cbind(test_data_, SCORE = ada.pred$prob[,2])
    write.table(test_data_output, row.names = FALSE, sep="\t", na = "nan", file = paste( train_stage, ".tsv", sep = ""), quote=FALSE)
}
#source("script/log_regression.R")
