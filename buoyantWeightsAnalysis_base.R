# set working directory (where you want files to be saved)
setwd("/Volumes/wrightworkshop/rachel/")

# load files (this can be from your working directory or any other directory, 
# if you know the path)
weights <- read.csv("../shared_files/buoyantWeight.csv")

# look at the file
head(weights)
tail(weights)

# what type of data are contained in this file?
str(weights)

# convert "genet" to a factor, not an integer
weights$genet <- as.factor(weights$genet)
# in english, "overwrite the column of "weights" called "genet" with the former
# contents of the column of weights called genet, but recast as a factor

# did it work?
str(weights)

# summarize your data
# look for missing data or strange outliers
summary(weights)

hist(weights$weight_4_15)
hist(weights$weight_2_16)

# start analyzing! ----
# notice that I just made my first "code section" by adding 
# >=4 symbols (-, =, or #) after the title

# calculate the difference between the first and 
# second time points
weights$weightdiff <- weights$weight_2_16-weights$weight_4_15

hist(weights$weightdiff, breaks = 10)

# calculate the difference as a percentage of starting weight
weights$percDiff <- (weights$weightdiff/weights$weight_4_15)*100
hist(weights$percDiff)

# plot differences in weight gain between variables ----
plot(percDiff ~ treat, data=weights)
plot(percDiff ~ treat, data=weights,
     main = "Skeletal Growth by Treatment",
     ylab = "Skeletal Growth (%)",
     xlab = "Bacterial Treatment (Control or Vibrio-challenged)")

plot(percDiff ~ bank, data=weights,
     main = "Skeletal Growth by Sampling Location",
     ylab = "Skeletal Growth (%)",
     xlab = "Sampling Location (East or West)")

plot(percDiff ~ pheno, data=weights,
     main = "Skeletal Growth by Phenotype",
     ylab = "Skeletal Growth (%)",
     xlab = "Phenotype (Susceptible or Resistant)")

# perform statistics -----
ttest_treat <- t.test(percDiff ~ treat, data=weights)
ttest_treat
ttest_treat$p.value

ttest_pheno <- t.test(percDiff ~ pheno, data=weights)
ttest_pheno 
ttest_pheno$p.value

ttest_bank <- t.test(percDiff ~ bank, data=weights)
ttest_bank 
ttest_bank$p.value

# plot differences in weight gain between variables, add p-values to figures ----
par(mfrow=c(1,3))
plot(percDiff ~ treat, data=weights,
     main = "Skeletal Growth by Treatment",
     ylab = "Skeletal Growth (%)",
     xlab = "Bacterial Treatment (Control or Vibrio-challenged)")
text(x = 2, y = 32, col = "red",
     labels = paste("p-value = ", round(ttest_treat$p.value,2)))

plot(percDiff ~ bank, data=weights,
     main = "Skeletal Growth by Sampling Location",
     ylab = "Skeletal Growth (%)",
     xlab = "Sampling Location (East or West)")
text(x = 2, y = 32, col = "red",
     labels = paste("p-value = ", round(ttest_bank$p.value,2)))

plot(percDiff ~ pheno, data=weights,
     main = "Skeletal Growth by Phenotype",
     ylab = "Skeletal Growth (%)",
     xlab = "Phenotype (Susceptible or Resistant)")
text(x = 1, y = 32, col = "red",
     labels = paste("p-value = ", round(ttest_pheno$p.value,2)))

# save objects as .Rdata
save(weights, file="myWeightData.Rdata")
        
