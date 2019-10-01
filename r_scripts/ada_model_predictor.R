#!/usr/bin/env Rscript

require("ada")

args <- commandArgs(TRUE)

trained_model = args[1]
test_filename = args[2]
output_filename = args[3]

# Make a copy of the input data since it will be modified, but don't output the modification into the output file
test_data_ = read.table(test_filename, header=TRUE)
test_data <- test_data_

# Create 6 features based on base substitution types, just in case those are training features. Otherwise, doesn't take much time. 
test_data$GC2CG = 0
test_data$GC2TA = 0
test_data$GC2AT = 0
test_data$TA2AT = 0
test_data$TA2GC = 0
test_data$TA2CG = 0

test_data$GC2CG[ (test_data$REF=='G' & test_data$ALT=='C') | (test_data$REF=='C' & test_data$ALT=='G') ] = 1
test_data$GC2TA[ (test_data$REF=='G' & test_data$ALT=='T') | (test_data$REF=='C' & test_data$ALT=='A') ] = 1
test_data$GC2AT[ (test_data$REF=='G' & test_data$ALT=='A') | (test_data$REF=='C' & test_data$ALT=='T') ] = 1
test_data$TA2AT[ (test_data$REF=='T' & test_data$ALT=='A') | (test_data$REF=='A' & test_data$ALT=='T') ] = 1
test_data$TA2GC[ (test_data$REF=='T' & test_data$ALT=='G') | (test_data$REF=='A' & test_data$ALT=='C') ] = 1
test_data$TA2CG[ (test_data$REF=='T' & test_data$ALT=='C') | (test_data$REF=='A' & test_data$ALT=='G') ] = 1


# Handle empty input data
if ( nrow(test_data)>=1 ) {
    load( trained_model )
    ada.pred <- predict(ada.model, newdata = test_data, type="both", n.iter=300)
    test_data_output <- cbind(test_data_, SCORE = ada.pred$prob[,2])

}   else {

        test_data_output <- test_data_
}

write.table(test_data_output, row.names = FALSE, sep="\t", na = "nan", file = output_filename, quote=FALSE)
