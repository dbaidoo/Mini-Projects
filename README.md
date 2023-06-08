# Mini-Projects
ANOVA and Design of Experiments Projects
```{r, echo = F, include=F}
library(SDAResources)
library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)
require(graphics)
library(reshape2)
library(xtable)
library(tidyr)
library(MASS)
library(car)
library(additivityTests)
library(asbio)
```

**Appendix 1.1 Reading the data**

```{r, echo = T, eval = TRUE}
# reading the potato txt file
potato <- read.table(file ="C://Users//baide//Desktop//STAT 545//potato.txt", header = TRUE) 

n <- nrow(potato)


head(potato)

potato$regime <- factor(potato$regime, label =c("R","C"))
potato$variety <- factor(potato$variety, label = c("s1","s2"))
potato$temp <- factor(potato$temp, label = c("-4","-8"))
attach(potato)
```

**Appendix 1.2 Redefining variables of the data and computation of factor level treatment means**
```{r, echo = T, eval = TRUE}
potato <- read.table(file = "C://Users//baide//Desktop//STAT 545//potato.txt", header = TRUE) 
potato$regime <- factor(potato$regime, label =c("R","C"))
potato$variety <- factor(potato$variety, label = c("s1","s2"))
potato$temp <- factor(potato$temp, label = c("-4","-8"))

regime_mean <- tapply(leak, regime, mean,data=potato) 
variety_mean <- tapply(leak, variety, mean, data=potato) 
temp_mean <- tapply(leak, temp, mean)
regime_mean
variety_mean
temp_mean

aggregate(leak~regime+variety+temp, data = potato, mean) # cell means
cell_sizes <- aggregate(leak~regime+variety+temp, data = potato, length) # unequal size
```


**Appendix 1.3 Fitting the ANOVA model**
```{r, echo = T, eval = TRUE}
myfit <- lm(leak~variety*regime*temp, contrasts= c(variety = contr.sum, regime = contr.sum, temp =contr.sum))
table <- Anova(myfit, type = 3)
table
```

**Appendix 1.3 Residual analysis and outlier checking**
```{r,echo=T,results='hold', fig.show='hide',eval=T} 
#### Residual analysis and outliers checking
par(mfrow=c(1,2)) 
qqPlot(myfit$residuals, las=1, main ="QQplot")
plot(myfit$fitted, myfit$resid, xlab ="Fitted", ylab="Residuals")
rstudent <- rstudent(myfit)
bcritical <- qt(1-0.05/(2*n),75-8-1)
outliers <- which(abs(rstudent) > bcritical)
outliers
leveneTest(leak~variety*temp*regime, data = potato)
shapiro.test(myfit$resid)
outlierTest(myfit)
potato[57:58,]

```


** Appendix 1.4 Appropriate Transformation**

```{r, echo = FALSE, eval = TRUE}
min(leak) #-1.11
leak2 <- leak + 2 #make response positive: box-cox needs a positive response
myfit2 <- lm(leak2~variety*regime*temp,contrasts = c(variety=contr.sum, regime=contr.sum,temp=contr.sum))

par(mfrow=c(1,1))
boxcox(myfit2,lambda =seq(-3,3, length = 10))

new_leak <- (leak+2)^(0.2)
```
**Fig 1.4.3.3 Transformation of the data  **

**Appendix 1.5 Fitting new ANOVA after an appropriate transformation using regression approach**

```{r, echo = T, eval = TRUE}
# regression approach: setting the factor

# x1 is the indicator variable for variety factor (need 1 indicator)
x1<-NULL

for(i in 1:n){
  if(potato$variety[i] == "s1"){
    x1[i] <- 1
  }
  else if(potato$variety[i] == "s2"){
          x1[i] <- -1

  }
}

# x2 is the indicator variable for regime factor
x2<-NULL

for(i in 1:n){
  if(potato$regime[i] == "R"){
    x2[i] <- 1
  }
  else if(potato$regime[i] == "C"){
          x2[i] <- -1
}
}


# x3 is the indicator variable for the temperature factor
x3<-NULL

for(i in 1:n){
  if(potato$temp[i] == -4){
    x3[i] <- 1
  }
  else if(potato$temp[i] == -8){
          x3[i] <- -1
}
} 

x1x2 <- x1*x2  
x1x3 <- x1*x3
x2x3 <- x2*x3
x1x2x3 <- x1*x2*x3

```

```{r, echo = T, eval = TRUE}
myfit_full <- lm(new_leak~variety*regime*temp,contrasts = c(variety=contr.sum, regime=contr.sum,temp=contr.sum))
myfit_full
# plot diagnostics
par(mfrow=c(2,2))
plot(myfit_full)
```
**Fig 1.5.1 Residual plots for transformed data**

```{r, echo = FALSE, eval = TRUE}
# see the anova table of full model and check for 3 factor interaction
# from ANOVA table, remove three interaction factor
table_1 <- Anova(myfit_full, type = 3)

# Reduced model with 3-interaction factor removed
myfit_red <- lm(new_leak~variety+regime+temp+variety:temp+regime:temp+variety:regime,contrast = c(variety=contr.sum,regime=contr.sum, temp = contr.sum))

# p-value = 0.464 > 0.05. Null hypothesis holds and accept the reduced model 
anova(myfit_full, myfit_red)

```

**Appendix 1.6 Reduced model with 3-interaction factor removed and interaction plots for the two-factor interaction**
```{r, echo = T, eval = TRUE}
# ANOVA table for the reduced model
# check 2 factor interaction


table_2 <- Anova(myfit_red,type = 3)

full_red <- anova(myfit_full, myfit_red)


## add the leak_new to the data
potato$new_leak <- new_leak



## plot for variety:regime interaction
potato.mean.1 <- ddply(potato, .(variety,regime), summarise, m = mean(new_leak))
p <- ggplot(potato, aes(x = variety , y =new_leak, colour = regime, shape = regime))
p <- p + geom_hline(aes(yintercept = 0),colour ="black", linetype = "solid",
                    size =0.2, alpha = 0.3)
p <- p + geom_boxplot(alpha = 0.25, outlier.size=0.1)
p <- p + geom_point(alpha = 0.5, position=position_dodge(width=0.75))
p <- p + geom_point(data = potato.mean.1, aes(y = m), size = 4)
p <- p + geom_line(data = potato.mean.1, aes(y = m, group = regime), size = 1.5)
p <- p + labs(title = "Potato interaction plot, Regime by Variety")
p

# interaction plot of regime:temperature. NOT significant
potato.mean.2 <- ddply(potato, .(temp,regime), summarise, m = mean(new_leak))
p2 <- ggplot(potato, aes(x = temp , y =new_leak, colour = regime, shape = regime))
p2 <- p2 + geom_hline(aes(yintercept = 0),colour ="black", linetype = "solid",
                    size =0.2, alpha = 0.3)
p2 <- p2 + geom_boxplot(alpha = 0.25, outlier.size=0.1)
p2 <- p2 + geom_point(alpha = 0.5, position=position_dodge(width=0.75))
p2 <- p2 + geom_point(data = potato.mean.2, aes(y = m), size = 4)
p2 <- p2 + geom_line(data = potato.mean.2, aes(y = m, group = regime), size = 1.5)
p2 <- p2 + labs(title = "Potato interaction plot, Regime by Temperature")
p2

## remove the variety:temp interaction.
## new reduced model and final model 
## the regime:temp and variety:regime are both significant interaction. keep them
myfit_new_red <- lm(new_leak~variety+regime+temp+regime:temp+variety:regime,contrast = c(variety=contr.sum,regime=contr.sum, temp = contr.sum))

table_3 <- Anova(myfit_new_red, type = 3)

```


```{r, echo = FALSE, eval = TRUE}
# ANOVA table of the final reduced model
table_3 <- Anova(myfit_new_red, type = 3)

# plot diagnostic plots
# according to the plot, the error of variance appears constant across the fitted values. Normality plot shows linearity in the relation between standardized residuals and the theoretical quantiles.
par(mfrow=c(2,2))
plot(myfit_new_red)
```

**Appendix 1.7 Residual Analysis for the transformed data** 

```{r, echo = FALSE, eval = TRUE} 
par(mfrow=c(1,1))
qqPlot(myfit_new_red$resid, las = 1, main ="QQplot") #linearity in plot appears

# p-value = 0.1292 > 0.05 
# accept the Null Hypothesis: Error terms are normal
normality_test <- shapiro.test(myfit_new_red$residuals) 

# Residuals vs fitted values plot shows constancy in error in variance
par(mfrow=c(1,1))
plot(myfit_new_red$fitted,myfit_new_red$resid, xlab= "Fitted Values", ylab="Residuals", main = "Resiudals vs Fitted Values")+abline(h=0,col = "gray75")


## Check outliers
## there is no outliers
rstudent_new <- rstudent(myfit_new_red)
bcritical2 <- qt(1-0.05/(2*n), 75-6-1)
outliers_new <- which(abs(rstudent_new)>bcritical2)

# OR use the outlier test
# p-value = 0.0024 < 0.05. No outliers 
outlier_test <- outlierTest(myfit_new_red)
```

```{r, echo = FALSE, eval = TRUE} 
# find the affects coefficients
linear_model <- lm(new_leak~x1+x2+x3+x1x2+x2x3)
# overall mean
potato.mean <- ddply(potato, .(variety,temp,regime), summarise, m = mean(new_leak))
overall_mean <- sum(potato.mean$m)/8
```

**Appendix 1.8 Pairwise Comparisons of Factor Level Means**


```{r, echo = FALSE, eval = TRUE} 
library(lsmeans)
regime_temp <- lsmeans(myfit_new_red, pairwise ~regime:temp, adjust = "tukey")
regime_var <- lsmeans(myfit_new_red, pairwise~regime:variety,
                       adjust ="tukey")
plot(regime_temp$contrasts)
```
**Fig 1.5.1.2a Pairwise Comparisons of Factor Level Means* 

```{r, echo = FALSE, eval = TRUE} 

plot(regime_var$contrasts)
```
**Fig 1.5.1.2b Pairwise Comparisons of Factor Level Means* 


**Appendix 1.9 Graphical display for the difference between factor level means.**

```{r, echo = FALSE, eval = TRUE}
potato <- read.table(file = "C://Users//baide//Desktop//STAT 545//potato.txt", header = TRUE) 

n <- nrow(potato)

potato$regime <- factor(potato$regime, label =c("R","C"))
potato$variety <- factor(potato$variety, label = c("s1","s2"))
potato$temp <- factor(potato$temp, label = c("-4","-8"))
attach(potato)
par(mfrow=c(1,3))
boxplot(leak ~ variety, data=potato, xlab= "Variety", 
ylab = "Mean damage score for ion Leakage")

boxplot(leak~regime, data = potato,
        xlab ="Regime",ylab = "Mean damage for ion leakage")
boxplot(leak~temp, data=potato,
        xlab="Temp",ylab = "Mean damage for ion leakage")
```
