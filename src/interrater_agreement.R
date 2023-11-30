#install.packages("irr")
library("irr")
#install.packages("readr")
library(readr)

data1 <- read_csv("pubmed.csv")
data2 <- read_csv("scopus.csv")
combined_data <- rbind(data1, data2)

# binarize it
combined_data <- combined_data + 1
combined_data <- sign(combined_data)


kappam.fleiss(combined_data, detail = TRUE)
agree(combined_data)
