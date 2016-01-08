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
train_data[,'if_COSMIC'] <- NULL
train_data[,'COSMIC_CNT'] <- NULL


model_formula <- as.formula(TrueVariant_or_False ~ .)

print("Fitting model...")
ada.model <- ada(model_formula, data = train_data, iter = 500)

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
