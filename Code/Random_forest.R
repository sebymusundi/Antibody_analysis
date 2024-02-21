
# clear working environment 
rm(list = ls())

#  load packages 
library(glmnet)
library(randomForest)
library(tidyverse)
library(GGally)
library(ROCR)

# Load dataset 

iris_data <- datasets::iris

# Check data 
 head(iris_data)
 
 # check associations 
 ggpairs(iris_data)

 # Random Forest classification 
 set.seed(123)

 # split training and test data 
 
 index_rows <- sample(1:2, 
                      nrow(iris), 
                      replace = T, 
                      prob = c(0.7, 0.3))
 
 train_data <- iris_data[index_rows==1, ]
 test_data <- iris_data[index_rows==2, ]
 
 
 # create model using the test data 
 iris_clasisfier <- randomForest(Species ~ ., data = train_data, 
                                 importance=T)
 
 

iris_clasisfier

# plot the training model 
plot(iris_clasisfier)

# importance of every feature
#This tells the predicting power of each predictor. 
#In consistent with the ggpair plot we’ve seen, Petal variables # 
#have greater predicting power. This is visualised in the plot above

importance(iris_clasisfier)

#  confirm importance by the plots below 
qplot(Petal.Width, Petal.Length, data=iris, color = Species)

qplot(Sepal.Width, Sepal.Length, data=iris, color = Species)

# Making predictions 
predicted_table <- predict(iris_clasisfier, test_data[,-5] )

table(observed=test_data[,5], predicted=predicted_table)

# Making predictions
predicted_probs <- predict(iris_clasisfier, test_data[, -5], type = "prob")

# Create a prediction object for ROC curve calculation
prediction_obj <- prediction(predicted_probs[, "versicolor"], test_data$Species)
