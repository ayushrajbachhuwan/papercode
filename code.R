#packages
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(etasFLP)
library(tidyverse)
library(devtools)
library(cluster)
library(neuralnet)
library(sqldf)
library(cowplot)

#data loading cleaning
data <- read.csv("eq1.csv")

#spatial analysis
spatial_data <- filter(data, latitude > 0.00000 & latitude < 43.00000)
spatial_data <- filter(spatial_data, longitude > 56.00000 & longitude < 100.00000)
spatial_data$time <- str_replace(spatial_data$time,"T"," ")
spatial_data$time <- strptime(spatial_data$time, format = "%Y-%m-%d")
spatial_data$time <- as.POSIXct(spatial_data$time)
spacefreq_table <- count(spatial_data, 'mag')

ggplot(spacefreq_table, aes(mag,log10(n))) + geom_line() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
fit <- lm(log10(n)~mag, data = spacefreq_table)
summary(fit)

temp_sp <- filter(spacefreq_table, mag > 4.3)
ggplot(temp_sp, aes(mag,log10(n))) + geom_line() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
fit <- lm(log10(n)~mag, data = temp_sp)
summary(fit)

temp <- scale(spatial_data[,2:3])
set.seed(1)
temp <- kmeans(spatial_data[,2:3], 6)
clusplot(spatial_data[,2:3], temp$cluster, color = T, shade = T)
temp <- cbind(spatial_data[,1:5], cluster = temp$cluster)
temp <- sqldf('SELECT * FROM temp ORDER BY cluster')

#cluster1 spatio-temporal analysis
clus1 <- sqldf('SELECT * FROM temp WHERE cluster = 1')
clus1freq <- count(clus1, mag)
fit <- lm(log10(n)~mag, data = clus1freq)
summary(fit)
c1 <- ggplot(clus1freq, aes(mag,log10(n))) + geom_point() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
tempclus <- sqldf('SELECT AVG(mag) as avg, MIN(mag) as min FROM clus1')
b <- (log10(2.7))/(tempclus$avg-tempclus$min)
a <- log10(262) + 4.2*b

temp <- sqldf('SELECT AVG(depth) as avg, mag FROM clus1 GROUP BY mag')
d1 <- ggplot(temp, aes(mag,avg)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Average Depth")

clus1$y <- format(clus1$time, format = "%Y")
temp_tpc1 <- sqldf('SELECT y, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number, AVG(depth) as depth FROM clus1 GROUP BY y')
ggplot(temp_tpc1, aes(y,number, group = 1)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Frequency")
ggplot(clus1, aes(time,mag)) + geom_line() + theme_gray()
ggplot(clus1, aes(time,depth)) + geom_line() + theme_gray()

b <- (1/temp_tpc1$Diff)*(log10(2.7))
temp <- sqldf('SELECT mag, COUNT(mag) as n, y FROM clus1 group by y, mag')
temp <- sqldf('SELECT mag, MAX(n) as n, y FROM temp group by y')
a <- log10(temp$n) + temp$mag*b
temp_tpc1$b <- b
temp_tpc1$a <- a
temp_tpc1[28,8:9] <- c(1,log10(765) + 4.5*1)

p1 <- ggplot(temp_tpc1) + geom_line(aes(y,maxmag, group = 1), color = "dark green") + theme_gray() + xlab("Time") + ylab("Max Magnitude")
p2 <- ggplot(temp_tpc1) + geom_line(aes(y,b, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("b")
p3 <- ggplot(temp_tpc1) + geom_line(aes(y,a, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("a")

plot_grid(p1, p2, p3, ncol = 1)

finalc1 <- temp_tpc1[,c(1,2,3,4,6,7,8,9)]
finalc1[,2:8] <- scale(finalc1[,2:8])

nnc1 <- neuralnet(formula = a ~ avgmag+minmag+maxmag+number+depth+b, data = finalc1[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
predict <- compute(nnc1, finalc1[13,2:7], rep = 5)
plot(nnc1)
error <- predict$net.result - finalc1[13,8]
error <- abs(error)

for (i in 1:100) {
  
  nnc1 <- neuralnet(formula = a ~ avgmag++minmag+maxmag+number+depth+b, data = finalc1[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
  predict <- compute(nnc1, finalc1[13,2:7], rep = 5)
  y <- c(y,predict$net.result)
  error <- predict$net.result - finalc1[13,8]
  error <- abs(error)
  error <- (error*100)/finalc1[13,8]
  x <- c(x,error)
  
}

y <- mean(y)
x <- mean(x)

#cluster2 spatio-temporal analysis
clus2 <- sqldf('SELECT * FROM temp WHERE cluster = 2')
clus2freq <- count(clus2, mag)
fit <- lm(log10(freq)~mag, data = clus2freq)
summary(fit)
c2 <- ggplot(clus2freq, aes(mag,log10(n))) + geom_point() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
tempclus <- sqldf('SELECT AVG(mag) as avg, MIN(mag) as min FROM clus2')
b <- (log10(2.7))/(tempclus$avg-tempclus$min)
a <- log10(562) + 4.2*b

temp <- sqldf('SELECT AVG(depth) as avg, mag FROM clus2 GROUP BY mag')
d2 <- ggplot(temp, aes(mag,avg)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Average Depth")

clus2$y <- format(clus2$time, format = "%Y")
temp_tpc2 <- sqldf('SELECT y, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number, AVG(depth) as depth FROM clus2 GROUP BY y')
ggplot(temp_tpc2, aes(y,number, group = 1)) + geom_line() + theme_gray()
ggplot(clus2, aes(time,mag)) + geom_line() + theme_gray()
ggplot(clus2, aes(time,depth)) + geom_line() + theme_gray()

b <- (1/temp_tpc2$Diff)*(log10(2.7))
temp <- sqldf('SELECT mag, COUNT(mag) as n, y FROM clus2 group by y, mag')
temp <- sqldf('SELECT mag, MAX(n) as n, y FROM temp group by y')
a <- log10(temp$n) + temp$mag*b
temp_tpc2$a <- a
temp_tpc2$b <- b

p1 <- ggplot(temp_tpc2) + geom_line(aes(y,maxmag, group = 1), color = "dark green") + theme_gray() + xlab("Time") + ylab("Max Magnitude")
p2 <- ggplot(temp_tpc2) + geom_line(aes(y,b, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("b")
p3 <- ggplot(temp_tpc2) + geom_line(aes(y,a, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("a")

plot_grid(p1, p2, p3, ncol = 1)

finalc2 <- temp_tpc2[,c(1,2,3,4,6,7,8,9)]
finalc2[,2:8] <- scale(finalc2[,2:8])

nnc2 <- neuralnet(formula = a ~ avgmag+minmag+maxmag+number+depth+b, data = finalc2[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
predict <- compute(nnc2, finalc2[13,2:7])
plot(nnc2)
error <- predict$net.result - finalc2[13,8]
error <- abs(error)

for (i in 1:100) {
  
  nnc2 <- neuralnet(formula = a ~ avgmag++minmag+maxmag+number+depth+b, data = finalc2[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
  predict <- compute(nnc2, finalc2[13,2:7])
  y <- c(y,predict$net.result)
  error <- predict$net.result - finalc2[13,8]
  error <- abs(error)
  error <- (error*100)/finalc2[13,8]
  x <- c(x,error)
  
}

y <- mean(y)
x <- mean(x)

#cluster3 spatio-temporal analysis
clus3 <- sqldf('SELECT * FROM temp WHERE cluster = 3')
clus3freq <- count(clus3, mag)
fit <- lm(log10(freq)~mag, data = clus3freq)
summary(fit)
c3 <- ggplot(clus3freq, aes(mag,log10(n))) + geom_point() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
tempclus <- sqldf('SELECT AVG(mag) as avg, MIN(mag) as min FROM clus3')
b <- (log10(2.7))/(tempclus$avg-tempclus$min)
a <- log10(176) + 4.4*b

temp <- sqldf('SELECT AVG(depth) as avg, mag FROM clus3 GROUP BY mag')
d3 <- ggplot(temp, aes(mag,avg)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Average Depth")

clus3$y <- format(clus3$time, format = "%Y")
temp_tpc3 <- sqldf('SELECT y, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number, AVG(depth) as depth FROM clus3 GROUP BY y')
ggplot(temp_tpc3, aes(y,number, group = 1)) + geom_line() + theme_gray()
ggplot(clus3, aes(time,mag)) + geom_line() + theme_gray()
ggplot(clus3, aes(time,depth)) + geom_line() + theme_gray()

b <- (1/temp_tpc3$Diff)*(log10(2.7))
temp <- sqldf('SELECT mag, COUNT(mag) as n, y FROM clus3 group by y, mag')
temp <- sqldf('SELECT mag, MAX(n) as n, y FROM temp group by y')
a <- log10(temp$n) + temp$mag*b
temp_tpc3$b <- b
temp_tpc3$a <- a

p1 <- ggplot(temp_tpc3) + geom_line(aes(y,maxmag, group = 1), color = "dark green") + theme_gray() + xlab("Time") + ylab("Max Magnitude")
p2 <- ggplot(temp_tpc3) + geom_line(aes(y,b, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("b")
p3 <- ggplot(temp_tpc3) + geom_line(aes(y,a, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("a")

plot_grid(p1, p2, p3, ncol = 1)

finalc3 <- temp_tpc3[,c(1,2,3,4,6,7,8,9)]
finalc3[,2:8] <- scale(finalc3[,2:8])

nnc3 <- neuralnet(formula = a ~ avgmag+minmag+maxmag+number+depth+b, data = finalc3[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
predict <- compute(nnc3, finalc3[13,2:7])
plot(nnc3)
error <- predict$net.result - finalc3[13,8]
error <- abs(error)

for (i in 1:100) {
  
  nnc3 <- neuralnet(formula = a ~ avgmag++minmag+maxmag+number+depth+b, data = finalc3[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
  predict <- compute(nnc3, finalc3[13,2:7])
  y <- c(y,predict$net.result)
  error <- predict$net.result - finalc3[13,8]
  error <- abs(error)
  error <- (error*100)/finalc3[13,8]
  x <- c(x,error)
  
}

y <- mean(y)
x <- mean(x)

#cluster4 spatio-temporal analysis
clus4 <- sqldf('SELECT * FROM temp WHERE cluster = 4')
clus4freq <- count(clus4, mag)
fit <- lm(log10(freq)~mag, data = clus4freq)
summary(fit)
c4 <- ggplot(clus4freq, aes(mag,log10(n))) + geom_point() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
tempclus <- sqldf('SELECT AVG(mag) as avg, MIN(mag) as min FROM clus4')
b <- (log10(2.7))/(tempclus$avg-tempclus$min)
a <- log10(174) + 4.3*b

temp <- sqldf('SELECT AVG(depth) as avg, mag FROM clus4 GROUP BY mag')
d4 <- ggplot(temp, aes(mag,avg)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Average Depth")

clus4$y <- format(clus4$time, format = "%Y")
temp_tpc4 <- sqldf('SELECT y, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number, AVG(depth) as depth FROM clus4 GROUP BY y')
ggplot(temp_tpc4, aes(y,number, group = 1)) + geom_line() + theme_gray()
ggplot(clus4, aes(time,mag)) + geom_line() + theme_gray()
ggplot(clus4, aes(time,depth)) + geom_line() + theme_gray()

b <- (1/temp_tpc4$Diff)*(log10(2.7))
temp <- sqldf('SELECT mag, COUNT(mag) as n, y FROM clus4 group by y, mag')
temp <- sqldf('SELECT mag, MAX(n) as n, y FROM temp group by y')
a <- log10(temp$n) + temp$mag*b
temp_tpc4$b <- b
temp_tpc4$a <- a

p1 <- ggplot(temp_tpc4) + geom_line(aes(y,maxmag, group = 1), color = "dark green") + theme_gray() + xlab("Time") + ylab("Max Magnitude")
p2 <- ggplot(temp_tpc4) + geom_line(aes(y,b, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("b")
p3 <- ggplot(temp_tpc4) + geom_line(aes(y,a, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("a")

plot_grid(p1, p2, p3, ncol = 1)

finalc4 <- temp_tpc4[,c(1,2,3,4,6,7,8,9)]
finalc4[,2:8] <- scale(finalc4[,2:8])

nnc3 <- neuralnet(formula = a ~ avgmag+minmag+maxmag+number+depth+b, data = finalc3[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
predict <- compute(nnc3, finalc3[13,2:7])
plot(nnc3)
error <- predict$net.result - finalc3[13,8]
error <- abs(error)

for (i in 1:100) {
  
  nnc4 <- neuralnet(formula = a ~ avgmag++minmag+maxmag+number+depth+b, data = finalc4[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
  predict <- compute(nnc4, finalc4[13,2:7])
  y <- c(y,predict$net.result)
  error <- predict$net.result - finalc4[13,8]
  error <- abs(error)
  error <- (error*100)/finalc4[13,8]
  x <- c(x,error)
  
}

y <- mean(y)
x <- mean(x)

#cluster5 spatio-temporal analysis
clus5 <- sqldf('SELECT * FROM temp WHERE cluster = 5')
clus5freq <- count(clus5, mag)
fit <- lm(log10(freq)~mag, data = clus5freq)
summary(fit)
c5 <- ggplot(clus5freq, aes(mag,log10(n))) + geom_point() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
tempclus <- sqldf('SELECT AVG(mag) as avg, MIN(mag) as min FROM clus5')
b <- (log10(2.7))/(tempclus$avg-tempclus$min)
a <- log10(314) + 4.3*b

temp <- sqldf('SELECT AVG(depth) as avg, mag FROM clus5 GROUP BY mag')
d5 <- ggplot(temp, aes(mag,avg)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Average Depth")

clus5$y <- format(clus5$time, format = "%Y")
temp_tpc5 <- sqldf('SELECT y, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number, AVG(depth) as depth FROM clus5 GROUP BY y')
ggplot(temp_tpc5, aes(y,number, group = 1)) + geom_line() + theme_gray()
ggplot(clus5, aes(time,mag)) + geom_line() + theme_gray()
ggplot(clus5, aes(time,depth)) + geom_line() + theme_gray()

b <- (1/temp_tpc5$Diff)*(log10(2.7))
temp <- sqldf('SELECT mag, COUNT(mag) as n, y FROM clus4 group by y, mag')
temp <- sqldf('SELECT mag, MAX(n) as n, y FROM temp group by y')
a <- log10(temp$n) + temp$mag*b
temp_tpc5$b <- b
temp_tpc5$a <- a

p1 <- ggplot(temp_tpc5) + geom_line(aes(y,maxmag, group = 1), color = "dark green") + theme_gray() + xlab("Time") + ylab("Max Magnitude")
p2 <- ggplot(temp_tpc5) + geom_line(aes(y,b, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("b")
p3 <- ggplot(temp_tpc5) + geom_line(aes(y,a, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("a")

plot_grid(p1, p2, p3, ncol = 1)

finalc5 <- temp_tpc5[,c(1,2,3,4,6,7,8,9)]
finalc5[,2:8] <- scale(finalc5[,2:8])

nnc5 <- neuralnet(formula = a ~ avgmag+minmag+maxmag+number+depth+b, data = finalc5[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
predict <- compute(nnc5, finalc5[13,2:7])
plot(nnc3)
error <- predict$net.result - finalc3[13,8]
error <- abs(error)

for (i in 1:100) {
  
  nnc5 <- neuralnet(formula = a ~ avgmag++minmag+maxmag+number+depth+b, data = finalc5[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
  predict <- compute(nnc5, finalc5[13,2:7])
  y <- c(y,predict$net.result)
  error <- predict$net.result - finalc5[13,8]
  error <- abs(error)
  error <- (error*100)/finalc5[13,8]
  x <- c(x,error)
  
}

y <- mean(y)
x <- mean(x)

#cluster6 spatio-temporal analysis
clus6 <- sqldf('SELECT * FROM temp WHERE cluster = 6')
clus6freq <- count(clus6, mag)
fit <- lm(log10(freq)~mag, data = clus6freq)
summary(fit)
c6 <- ggplot(clus6freq, aes(mag,log10(n))) + geom_point() + geom_smooth(se=FALSE) + theme_gray() + xlab("Magnitude") + ylab("Frequency")
tempclus <- sqldf('SELECT AVG(mag) as avg, MIN(mag) as min FROM clus6')
b <- (log10(2.7))/(tempclus$avg-tempclus$min)
a <- log10(138) + 4.5*b

temp <- sqldf('SELECT AVG(depth) as avg, mag FROM clus6 GROUP BY mag')
d6 <- ggplot(temp, aes(mag,avg)) + geom_line() + theme_gray() + xlab("Magnitude") + ylab("Average Depth")

clus6$y <- format(clus6$time, format = "%Y")
temp_tpc6 <- sqldf('SELECT y, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number, AVG(depth) as depth FROM clus6 GROUP BY y')
ggplot(temp_tpc6, aes(y,number, group = 1)) + geom_line() + theme_gray()
ggplot(clus6, aes(time,mag)) + geom_line() + theme_gray()
ggplot(clus6, aes(time,depth)) + geom_line() + theme_gray()

b <- (1/temp_tpc6$Diff)*(log10(2.7))
temp <- sqldf('SELECT mag, COUNT(mag) as n, y FROM clus4 group by y, mag')
temp <- sqldf('SELECT mag, MAX(n) as n, y FROM temp group by y')
a <- log10(temp$n) + temp$mag*b
temp_tpc6$b <- b
temp_tpc6$a <- a

p1 <- ggplot(temp_tpc6) + geom_line(aes(y,maxmag, group = 1), color = "dark green") + theme_gray() + xlab("Time") + ylab("Max Magnitude")
p2 <- ggplot(temp_tpc6) + geom_line(aes(y,b, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("b")
p3 <- ggplot(temp_tpc6) + geom_line(aes(y,a, group = 2) , color = "steel blue") + theme_gray() + xlab("Time") + ylab("a")

plot_grid(p1, p2, p3, ncol = 1)

finalc6 <- temp_tpc6[,c(1,2,3,4,6,7,8,9)]
finalc6[,2:8] <- scale(finalc6[,2:8])

nnc6 <- neuralnet(formula = a ~ avgmag+minmag+maxmag+number+depth+b, data = finalc6[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
predict <- compute(nnc5, finalc5[13,2:7])
plot(nnc3)
error <- predict$net.result - finalc3[13,8]
error <- abs(error)

for (i in 1:100) {
  
  nnc6 <- neuralnet(formula = a ~ avgmag++minmag+maxmag+number+depth+b, data = finalc6[1:12,2:8], hidden = 10, rep = 10, linear.output = T, err.fct = "sse")
  predict <- compute(nnc6, finalc6[13,2:7])
  y <- c(y,predict$net.result)
  error <- predict$net.result - finalc6[13,8]
  error <- abs(error)
  error <- (error*100)/finalc6[13,8]
  x <- c(x,error)
  
}

y <- mean(y)
x <- mean(x)

plot_grid(c1, c2, c3, c4, c5, c6, ncol = 3, labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"), label_size = 8, hjust = -2.25, vjust = 2.25)
plot_grid(d1, d2, d3, d4, d5, d6, ncol = 3, labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"), label_size = 8, hjust = -2.25, vjust = 2.25)

#temporal analysis
temporal_data <- filter(data, latitude > 0.00000 & latitude < 43.00000)
temporal_data <- filter(temporal_data, longitude > 56.00000 & longitude < 100.00000)
temporal_data$time <- str_replace(temporal_data$time,"T"," ")
temporal_data$time <- strptime(temporal_data$time, format = "%Y-%m-%d")
temporal_data$time <- as.POSIXct(temporal_data$time)
timefreq_table <- count(temporal_data, time)
ggplot(timefreq_table, aes(time,n)) + geom_line() + theme_gray() + xlab("Time") + ylab("Frequency")
ggplot(temporal_data, aes(time,mag)) + geom_line() + theme_gray()  + xlab("Time") + ylab("Magnitude")
ggplot(temporal_data, aes(time,depth)) + geom_line() + theme_gray() + xlab("Time") + ylab("Depth")

temporal_data$year <- format(temporal_data$time, format = "%Y")
temp_tp <- sqldf('SELECT YEAR, AVG(mag) as avgmag, MIN(mag) as minmag, MAX(mag) as maxmag, AVG(mag)-MIN(mag) as Diff, COUNT(*) as number FROM temporal_data GROUP BY year')
b <- (1/temp_tp$Diff)*(log10(2.7))
temp_tp$b <- b

p1 <- ggplot(temp_tp) + geom_line(aes(year,maxmag, group = 1), color = "dark green") + theme_gray()
p2 <- ggplot(temp_tp) + geom_line(aes(year,b, group = 2) , color = "steel blue") + theme_gray()

plot_grid(p1, p2, labels = "AUTO", ncol = 1)