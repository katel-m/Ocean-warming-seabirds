##%######################################################%##
#                                                          #
####            Breeding success SST models             ####
#                                                          #
##%######################################################%##

#' Script models breeding season SST effect on breeding success with a quasibinomial model to account for overdispersion 
#' Predictions are made for current (2015:2020) and future (2035:2040) average SSTs based on global warming forecast 

##%######################################################%##

library(effects)
library(performance)
library(DHARMa)
library(ggplot2)
library(stringr)

# Load SST data 
sst.data <- read.csv("Data_BS.csv")
sst.data$Unique.id <- paste(sst.data$Species, sst.data$Colony)
names(sst.data)[2] = "year"

# Load breeding success data 
bs.data <- read.csv("Data_SSTbs.csv")
names(bs.data) = c("Species","Year","Colony","BS","N","Successes","SST")

# Control for >1 clutch size for kittiwakes
bs.data$Failures <- ifelse(bs.data$Species == "Black legged kittiwake", (bs.data$N*2)-bs.data$Successes, bs.data$N - bs.data$Successes)

bs.data$Unique.id <- paste(bs.data$Species, bs.data$Colony)


########################
#
# Model running for each population seperately 
#
########################

################################################ # 1 ################################################
# Run model
df1 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[1])

df1$Scaled.SST <- scale(df1$SST)

m1.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST), data = df1, family = quasibinomial())

summary(m1.noN)

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[1])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m1.noN)$coefficients[2,4] < 0.05){
  df1_pred=data.frame(Species=unique(df1$Species),
                      Colony=unique(df1$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df1_pred=data.frame(Species=unique(df1$Species),
                      Colony=unique(df1$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df1$Scaled.SST),
                                   mean(df1$Scaled.SST)))
}

df1_pred$mean=predict(m1.noN,newdata = df1_pred,se.fit = T,type = "response")$fit
df1_pred$se=predict(m1.noN,newdata = df1_pred,se.fit = T,type = "response")$se.fit

################################################ # 2 ################################################
# Model run 
df2 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[2])

df2$Scaled.SST <- scale(df2$SST)

m2.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df2, family = quasibinomial())

summary(m2.noN)

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[2])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m2.noN)$coefficients[2,4] < 0.05){
  df2_pred=data.frame(Species=unique(df2$Species),
                      Colony=unique(df2$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df2_pred=data.frame(Species=unique(df2$Species),
                      Colony=unique(df2$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df2$Scaled.SST),
                                   mean(df2$Scaled.SST)))
}
df2_pred$mean=predict(m2.noN,newdata = df2_pred,se.fit = T,type = "response")$fit
df2_pred$se=predict(m2.noN,newdata = df2_pred,se.fit = T,type = "response")$se.fit

################################################ # 3 ################################################
# Model run 
df3 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[3])

df3$Scaled.SST <- scale(df3$SST)

m3.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST) , data = df3, family = quasibinomial())

summary(m3.noN)

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[3])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m3.noN)$coefficients[2,4] < 0.05){
  df3_pred=data.frame(Species=unique(df3$Species),
                      Colony=unique(df3$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df3_pred=data.frame(Species=unique(df3$Species),
                      Colony=unique(df3$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df3$Scaled.SST),
                                   mean(df3$Scaled.SST)))
}

df3_pred$mean=predict(m3.noN,newdata = df3_pred,se.fit = T,type = "response")$fit
df3_pred$se=predict(m3.noN,newdata = df3_pred,se.fit = T,type = "response")$se.fit
print(df3_pred)

################################################ # 4 ################################################
# Model run 
df4 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[4])

df4$Scaled.SST <- scale(df4$SST)

m4.noN <- glm(cbind(Successes, Failures) ~  as.numeric(Scaled.SST) , data = df4, family = quasibinomial())

summary(m4.noN)

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[4])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m4.noN)$coefficients[2,4] < 0.05){
  df4_pred=data.frame(Species=unique(df4$Species),
                      Colony=unique(df4$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df4_pred=data.frame(Species=unique(df4$Species),
                      Colony=unique(df4$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df4$Scaled.SST),
                                   mean(df4$Scaled.SST)))
}
df4_pred$mean=predict(m4.noN,newdata = df4_pred,se.fit = T,type = "response")$fit
df4_pred$se=predict(m4.noN,newdata = df4_pred,se.fit = T,type = "response")$se.fit
print(df4_pred)

################################################ # 5 ################################################
# Model run 
df5 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[5])

df5$Scaled.SST <- scale(df5$SST)

m5.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df5, family = quasibinomial())

summary(m5.noN)

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[5])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m5.noN)$coefficients[2,4] < 0.05){
  df5_pred=data.frame(Species=unique(df5$Species),
                      Colony=unique(df5$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df5_pred=data.frame(Species=unique(df5$Species),
                      Colony=unique(df5$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df5$Scaled.SST),
                                   mean(df5$Scaled.SST)))
}

df5_pred$mean=predict(m5.noN,newdata = df5_pred,se.fit = T,type = "response")$fit
df5_pred$se=predict(m5.noN,newdata = df5_pred,se.fit = T,type = "response")$se.fit
print(df5_pred)

################################################ # 6 ################################################
# Model run 
df6 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[6])

df6$Scaled.SST <- scale(df6$SST)

m6.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df6, family = quasibinomial())

summary(m6.noN)

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[6])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m6.noN)$coefficients[2,4] < 0.05){
  df6_pred=data.frame(Species=unique(df6$Species),
                      Colony=unique(df6$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df6_pred=data.frame(Species=unique(df6$Species),
                      Colony=unique(df6$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df6$Scaled.SST),
                                   mean(df6$Scaled.SST)))
}

df6_pred$mean=predict(m6.noN,newdata = df6_pred,se.fit = T,type = "response")$fit
df6_pred$se=predict(m6.noN,newdata = df6_pred,se.fit = T,type = "response")$se.fit
print(df6_pred)

################################################ # 7 ################################################
# Model run 
df7 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[7])

df7$Scaled.SST <- scale(df7$SST)

m7.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df7, family = quasibinomial())

summary(m7.noN)

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[7])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m7.noN)$coefficients[2,4] < 0.05){
  df7_pred=data.frame(Species=unique(df7$Species),
                      Colony=unique(df7$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df7_pred=data.frame(Species=unique(df7$Species),
                      Colony=unique(df7$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df7$Scaled.SST),
                                   mean(df7$Scaled.SST)))
}

df7_pred$mean=predict(m7.noN,newdata = df7_pred,se.fit = T,type = "response")$fit
df7_pred$se=predict(m7.noN,newdata = df7_pred,se.fit = T,type = "response")$se.fit
print(df7_pred)

################################################ # 8 ################################################
# Model run 
df8 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[8])

df8$Scaled.SST <- scale(df8$SST)

m8.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df8, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[8])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m8.noN)$coefficients[2,4] < 0.05){
  df8_pred=data.frame(Species=unique(df8$Species),
                      Colony=unique(df8$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df8_pred=data.frame(Species=unique(df8$Species),
                      Colony=unique(df8$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df8$Scaled.SST),
                                   mean(df8$Scaled.SST)))
}

df8_pred$mean=predict(m8.noN,newdata = df8_pred,se.fit = T,type = "response")$fit
df8_pred$se=predict(m8.noN,newdata = df8_pred,se.fit = T,type = "response")$se.fit
print(df8_pred)

################################################ # 9 ################################################
# Model run 
df9 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[9])

df9$Scaled.SST <- scale(df9$SST)

m9.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df9, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[9])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m9.noN)$coefficients[2,4] < 0.05){
  df9_pred=data.frame(Species=unique(df9$Species),
                      Colony=unique(df9$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df9_pred=data.frame(Species=unique(df9$Species),
                      Colony=unique(df9$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df9$Scaled.SST),
                                   mean(df9$Scaled.SST)))
}

df9_pred$mean=predict(m9.noN,newdata = df9_pred,se.fit = T,type = "response")$fit
df9_pred$se=predict(m9.noN,newdata = df9_pred,se.fit = T,type = "response")$se.fit
print(df9_pred)

################################################ # 10 ################################################
# model run 
df10 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[10])

df10$Scaled.SST <- scale(df10$SST)

m10.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df10, family = quasibinomial())

summary(m10.noN)

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[10])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m10.noN)$coefficients[2,4] < 0.05){
  df10_pred=data.frame(Species=unique(df10$Species),
                      Colony=unique(df10$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df10_pred=data.frame(Species=unique(df10$Species),
                      Colony=unique(df10$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df10$Scaled.SST),
                                   mean(df10$Scaled.SST)))
}

df10_pred$mean=predict(m10.noN,newdata = df10_pred,se.fit = T,type = "response")$fit
df10_pred$se=predict(m10.noN,newdata = df10_pred,se.fit = T,type = "response")$se.fit
print(df10_pred)

################################################ # 11 ################################################
# Model run 
df11 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[11])

df11$Scaled.SST <- scale(df11$SST)

m11.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df11, family = quasibinomial())

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[11])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m11.noN)$coefficients[2,4] < 0.05){
  df11_pred=data.frame(Species=unique(df11$Species),
                      Colony=unique(df11$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df11_pred=data.frame(Species=unique(df11$Species),
                      Colony=unique(df11$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df11$Scaled.SST),
                                   mean(df11$Scaled.SST)))
}

df11_pred$mean=predict(m11.noN,newdata = df11_pred,se.fit = T,type = "response")$fit
df11_pred$se=predict(m11.noN,newdata = df11_pred,se.fit = T,type = "response")$se.fit
print(df11_pred)

################################################ # 12 ################################################
# Model run 
df12 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[12])

df12$Scaled.SST <- scale(df12$SST)

m12.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df12, family = quasibinomial())

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[12])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m12.noN)$coefficients[2,4] < 0.05){
  df12_pred=data.frame(Species=unique(df12$Species),
                      Colony=unique(df12$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df12_pred=data.frame(Species=unique(df12$Species),
                      Colony=unique(df12$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df12$Scaled.SST),
                                   mean(df12$Scaled.SST)))
}

df12_pred$mean=predict(m12.noN,newdata = df12_pred,se.fit = T,type = "response")$fit
df12_pred$se=predict(m12.noN,newdata = df12_pred,se.fit = T,type = "response")$se.fit
print(df12_pred)

################################################ # 13 ################################################
# Model run 
df13 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[13])

df13$Scaled.SST <- scale(df13$SST)

m13.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df13, family = quasibinomial())

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[13])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m13.noN)$coefficients[2,4] < 0.05){
  df13_pred=data.frame(Species=unique(df13$Species),
                      Colony=unique(df13$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df13_pred=data.frame(Species=unique(df13$Species),
                      Colony=unique(df13$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df13$Scaled.SST),
                                   mean(df13$Scaled.SST)))
}

df13_pred$mean=predict(m13.noN,newdata = df13_pred,se.fit = T,type = "response")$fit
df13_pred$se=predict(m13.noN,newdata = df13_pred,se.fit = T,type = "response")$se.fit
print(df13_pred)

################################################ # 14 ################################################
# Model run 
df14 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[14])

df14$Scaled.SST <- scale(df14$SST)

m14.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df14, family = quasibinomial())

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[14])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))
if(summary(m14.noN)$coefficients[2,4] < 0.05){
  df14_pred=data.frame(Species=unique(df14$Species),
                      Colony=unique(df14$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df14_pred=data.frame(Species=unique(df14$Species),
                      Colony=unique(df14$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df14$Scaled.SST),
                                   mean(df14$Scaled.SST)))
}

df14_pred$mean=predict(m14.noN,newdata = df14_pred,se.fit = T,type = "response")$fit
df14_pred$se=predict(m14.noN,newdata = df14_pred,se.fit = T,type = "response")$se.fit
print(df14_pred)

################################################ # 15 ################################################
# Model run 
df15 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[15])

df15$Scaled.SST <- scale(df15$SST)

m15.noN<- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df15, family = quasibinomial())

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[15])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m15.noN)$coefficients[2,4] < 0.05){
  df15_pred=data.frame(Species=unique(df15$Species),
                      Colony=unique(df15$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df15_pred=data.frame(Species=unique(df15$Species),
                      Colony=unique(df15$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df15$Scaled.SST),
                                   mean(df15$Scaled.SST)))
}

df15_pred$mean=predict(m15.noN,newdata = df15_pred,se.fit = T,type = "response")$fit
df15_pred$se=predict(m15.noN,newdata = df15_pred,se.fit = T,type = "response")$se.fit
print(df15_pred)

################################################ # 16 ################################################
# Model run 
df16 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[16]) 

df16$Scaled.SST <- scale(df16$SST)

m16.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df16, family = quasibinomial())

# Prediction if significant effect

sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[16])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m16.noN)$coefficients[2,4] < 0.05){
  df16_pred=data.frame(Species=unique(df16$Species),
                      Colony=unique(df16$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df16_pred=data.frame(Species=unique(df16$Species),
                      Colony=unique(df16$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df16$Scaled.SST),
                                   mean(df16$Scaled.SST)))
}

df16_pred$mean=predict(m16.noN,newdata = df16_pred,se.fit = T,type = "response")$fit
df16_pred$se=predict(m16.noN,newdata = df16_pred,se.fit = T,type = "response")$se.fit
print(df16_pred)

################################################ # 17 ################################################
# Model run
df17 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[17]) 

df17$Scaled.SST <- scale(df17$SST)

m17.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df17, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[17])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m17.noN)$coefficients[2,4] < 0.05){
  df17_pred=data.frame(Species=unique(df17$Species),
                      Colony=unique(df17$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df17_pred=data.frame(Species=unique(df17$Species),
                      Colony=unique(df17$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df17$Scaled.SST),
                                   mean(df17$Scaled.SST)))
}

df17_pred$mean=predict(m17.noN,newdata = df17_pred,se.fit = T,type = "response")$fit
df17_pred$se=predict(m17.noN,newdata = df17_pred,se.fit = T,type = "response")$se.fit
print(df17_pred)

################################################ # 18 ################################################
# Model run 
df18 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[18])

df18$Scaled.SST <- scale(df18$SST)

m18.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df18, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[18])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m18.noN)$coefficients[2,4] < 0.05){
  df18_pred=data.frame(Species=unique(df18$Species),
                      Colony=unique(df18$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df18_pred=data.frame(Species=unique(df18$Species),
                      Colony=unique(df18$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df18$Scaled.SST),
                                   mean(df18$Scaled.SST)))
}

df18_pred$mean=predict(m18.noN,newdata = df18_pred,se.fit = T,type = "response")$fit
df18_pred$se=predict(m18.noN,newdata = df18_pred,se.fit = T,type = "response")$se.fit
print(df18_pred)

################################################ # 19 ################################################
# Model run 
df19 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[19])

df19$Scaled.SST <- scale(df19$SST)

m19.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df19, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[19])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m19.noN)$coefficients[2,4] < 0.05){
  df19_pred=data.frame(Species=unique(df19$Species),
                      Colony=unique(df19$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df19_pred=data.frame(Species=unique(df19$Species),
                      Colony=unique(df19$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df19$Scaled.SST),
                                   mean(df19$Scaled.SST)))
}

df19_pred$mean=predict(m19.noN,newdata = df19_pred,se.fit = T,type = "response")$fit
df19_pred$se=predict(m19.noN,newdata = df19_pred,se.fit = T,type = "response")$se.fit
print(df19_pred)

################################################ # 20 ################################################
# Model run 
df20 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[20])

df20$Scaled.SST <- scale(df20$SST)

m20.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df20, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[20])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m20.noN)$coefficients[2,4] < 0.05){
  df20_pred=data.frame(Species=unique(df20$Species),
                      Colony=unique(df20$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df20_pred=data.frame(Species=unique(df20$Species),
                      Colony=unique(df20$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df20$Scaled.SST),
                                   mean(df20$Scaled.SST)))
}

df20_pred$mean=predict(m20.noN,newdata = df20_pred,se.fit = T,type = "response")$fit
df20_pred$se=predict(m20.noN,newdata = df20_pred,se.fit = T,type = "response")$se.fit
print(df20_pred)

################################################ # 21 ################################################
# Model run 
df21 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[21])

df21$Scaled.SST <- scale(df21$SST)

m21.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df21, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[21])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m21.noN)$coefficients[2,4] < 0.05){
  df21_pred=data.frame(Species=unique(df21$Species),
                      Colony=unique(df21$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df21_pred=data.frame(Species=unique(df21$Species),
                      Colony=unique(df21$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df21$Scaled.SST),
                                   mean(df21$Scaled.SST)))
}

df21_pred$mean=predict(m21.noN,newdata = df21_pred,se.fit = T,type = "response")$fit
df21_pred$se=predict(m21.noN,newdata = df21_pred,se.fit = T,type = "response")$se.fit
print(df21_pred)

################################################ # 22 ################################################
# Model run 
df22 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[22])

df22$Scaled.SST <- scale(df22$SST)

m22.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df22, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[22])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m22.noN)$coefficients[2,4] < 0.05){
  df22_pred=data.frame(Species=unique(df22$Species),
                      Colony=unique(df22$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df22_pred=data.frame(Species=unique(df22$Species),
                      Colony=unique(df22$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df22$Scaled.SST),
                                   mean(df22$Scaled.SST)))
}

df22_pred$mean=predict(m22.noN,newdata = df22_pred,se.fit = T,type = "response")$fit
df22_pred$se=predict(m22.noN,newdata = df22_pred,se.fit = T,type = "response")$se.fit
print(df22_pred)

################################################ # 23 ################################################
# Model run 
df23 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[23])

df23$Scaled.SST <- scale(df23$SST)

m23.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df23, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[23])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m23.noN)$coefficients[2,4] < 0.05){
  df23_pred=data.frame(Species=unique(df23$Species),
                      Colony=unique(df23$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df23_pred=data.frame(Species=unique(df23$Species),
                      Colony=unique(df23$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df23$Scaled.SST),
                                   mean(df23$Scaled.SST)))
}

df23_pred$mean=predict(m23.noN,newdata = df23_pred,se.fit = T,type = "response")$fit
df23_pred$se=predict(m23.noN,newdata = df23_pred,se.fit = T,type = "response")$se.fit
print(df23_pred)

################################################ # 24 ################################################
# Model run 
df24 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[24])

df24$Scaled.SST <- scale(df24$SST)

m24.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df24, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[24])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m24.noN)$coefficients[2,4] < 0.05){
  df24_pred=data.frame(Species=unique(df24$Species),
                      Colony=unique(df24$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df24_pred=data.frame(Species=unique(df24$Species),
                      Colony=unique(df24$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df24$Scaled.SST),
                                   mean(df24$Scaled.SST)))
}

df24_pred$mean=predict(m24.noN,newdata = df24_pred,se.fit = T,type = "response")$fit
df24_pred$se=predict(m24.noN,newdata = df24_pred,se.fit = T,type = "response")$se.fit
print(df24_pred)

################################################ # 25 ################################################
# Model run 
df25 <- subset(bs.data, Unique.id == unique(bs.data$Unique.id)[25])

df25$Scaled.SST <- scale(df25$SST)

m25.noN <- glm(cbind(Successes, Failures) ~ as.numeric(Scaled.SST)  , data = df25, family = quasibinomial())

# Prediction if significant effect
sst_sub <- subset(sst.data, Unique.id == unique(bs.data$Unique.id)[25])
sst_sub$SST <- scale(sst_sub$SST, center = attr(df1$Scaled.SST,"scaled:center"), scale = attr(df1$Scaled.SST,"scaled:scale"))

if(summary(m25.noN)$coefficients[2,4] < 0.05){
  df25_pred=data.frame(Species=unique(df25$Species),
                      Colony=unique(df25$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(sst_sub[sst_sub$year %in% (2015:2020),2]),
                                   mean(sst_sub[sst_sub$year %in% (2035:2040),2]))) 
} else {
  df25_pred=data.frame(Species=unique(df25$Species),
                      Colony=unique(df25$Colony),
                      Scenario=c("Current","Future"),
                      Scaled.SST=c(mean(df25$Scaled.SST),
                                   mean(df25$Scaled.SST)))
}

df25_pred$mean=predict(m25.noN,newdata = df25_pred,se.fit = T,type = "response")$fit
df25_pred$se=predict(m25.noN,newdata = df25_pred,se.fit = T,type = "response")$se.fit
print(df25_pred)

################################################ # End ################################################

coeffs <- NULL
pvals <- NULL
lower.ci <- NULL
upper.ci <- NULL
Signif <- NULL

for(i in 1:25){
  m.name <- paste0("m",i,".noN")
  mod <- get(m.name)
  coeffs <- append(coeffs, coef(summary(mod))[2, 'Estimate'])
  pvals <- append(pvals, coef(summary(mod))[2,'Pr(>|t|)'])
  lower.ci <- append(lower.ci, confint(mod)[2,1])
  upper.ci <- append(upper.ci, confint(mod)[2,2])
}

# Combine prediction dataframes 
df.list=list()

for(i in 1:25){
  df.name <- paste0("df",i,"_pred")
  df <- get(df.name)
  df.list[[i]]=df
}
pred_df_all=do.call(rbind,df.list)

effs.df <- data.frame(unique(bs.data$Unique.id), coeffs, pvals, lower.ci, upper.ci,Sig=ifelse(pvals<0.05,1,0))

saveRDS(pred_df_all,file=paste("output/sstPred_formanus_BS.RDS",sep=""))
