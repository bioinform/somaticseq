#!/usr/bin/env R

require("ada")

args <- commandArgs(TRUE)

trained_model = args[1]
test_filename = args[2]
output_filename = args[3]

test_data_ = read.table(test_filename, header=TRUE)

test_data <- test_data_

# Remove from training identities, i.e., CHROM, POS, etc.
test_data <- test_data[,-c(1, 2, 3, 4, 5)]

test_data[,'TrueVariant_or_False'] <- NULL
test_data[,'REF'] <- NULL
test_data[,'ALT'] <- NULL

# A bug where sometimes SOR gets classified as "types"
test_data$SOR <- as.numeric(test_data$SOR)

# Handle empty input data
if ( nrow(test_data)>=1 ) {
    load( trained_model )
    ada.pred <- predict(ada.model, newdata = test_data, type="both")
    test_data_output <- cbind(test_data_, SCORE = ada.pred$prob[,2])

}   else {

        test_data_output <- test_data_
}

write.table(test_data_output, row.names = FALSE, sep="\t", na = "nan", file = output_filename, quote=FALSE)
