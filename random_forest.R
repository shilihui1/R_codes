library(randomForest)
library(ggplot2)
library(caret)
library(e1071)
library(reshape)
library(party)

data("iris")

dim(iris)
str(iris)

head(iris)
tail(iris)

par(mfrow=c(1,3))
boxplot(iris$sepal.length[which(iris$species == 'setosa')])
boxplot(iris$sepal.length[which(iris$species == 'versicolor')])
boxplot(iris$sepal.length[which(iris$species == 'virginica')])
par(mfrow=c(1,1)
    
    
DF <- melt(iris, vars = c("iris$species"), vars = c("iris$sepal.Length", "iris$sepal.Width"))
head(DF)
boxplot(value ~ species + variable, DF)
    
    
inTrain
    
dim(training)
head(training)

dim(training)
head(training)
    
modFit <- train(species ~ ., data=training, method="rf", prox= TRUE)
print(modFit)
    
getTree(modFit$finalModel, k=2)
    
    
irisP <- classCenter(training[,c(3,4)],training$species, modFit$finalModel$prox)
    
    
    
    
    
    
    
    
