#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

gene = args[1]
v = args[2]

system("source ~/.bashrc")

library(h2o)

train_v <- read.csv(paste(gene,v,'train.fold.csv',sep='_'))
test_v <- read.csv(paste(gene,v,'test.fold.csv',sep='_'))

models = read.csv('sources/base_thresholds', header = FALSE, sep="\t", stringsAsFactors = FALSE)
COLUMNS = models$V1

h2o.init(min_mem_size = "8g", max_mem_size = "16g", nthreads = 1)
train = as.h2o(train_v)
all_h2o_models = h2o.automl(x = COLUMNS, y = "Call",training_frame = train,keep_cross_validation_predictions = FALSE, nfolds = 5,  max_runtime_secs = 3600,sort_metric = 'mean_per_class_error',stopping_metric = 'mean_per_class_error',seed = 100)
test = as.h2o(test_v)
pred = as.data.frame(h2o.predict(all_h2o_models@leader, newdata = test))
write.csv(paste(gene,v,'res.fold.csv',sep='_'))

