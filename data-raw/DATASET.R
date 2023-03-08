## code to prepare `DATASET` dataset goes here

d<-read.csv(file.choose())
perfectionism<-d
usethis::use_data(perfectionism, overwrite = TRUE)
