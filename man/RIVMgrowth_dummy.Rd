\name{RIVMgrowth_dummy}
\alias{RIVMgrowth_dummy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latent Growth Curve models with dummy variables
}
\description{
Function for non-linear LGC model with dummy variable predictors
}
\usage{
RIVMgrowth_dummy(data, type = c("non-linear", "linear"), independent, time, ordered = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame
}
  \item{type}{
\code{"non-linear"} for a non-linear model and \code{"linear"} for a linear model
}
  \item{independent}{
A vector of strings indicating variables used as independent predictor variables (categorical)
}
  \item{time}{
A vector encoding relative time of the measurements
}
  \item{ordered}{
Logial, are independent variables ordered-categorical or continuous? Can also be a vector with a value per independent variable.
}
  \item{\dots}{
Encoding of the dependent variables. Currenlty only one dependent variable is supported. USse \code{name_of_dependent = c("name_of_wave_1","name_of_wave_2",...)}
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
# Load needed packages:

# Let's first simulate some data:
library("lavaan")
simmodel <- '
int_eenzaam =~ 1*eenzaam_1 + 1*eenzaam_2 + 1*eenzaam_3 + 1*eenzaam_4 + 1*eenzaam_5 + 1*eenzaam_6 + 1*eenzaam_7
slope_eenzaam =~ 0*eenzaam_1 + 1*eenzaam_2 + 2*eenzaam_3 + 3.4*eenzaam_4 + 4.4*eenzaam_5 + 4.4*eenzaam_6 + 3.7*eenzaam_7
slope_eenzaam ~ -0.079*1
int_eenzaam ~ 2.7*1

int_mental =~ 1*mental_1 + 1*mental_2 + 1*mental_3 + 1*mental_4 + 1*mental_5 + 1*mental_6 + 1*mental_7
slope_mental =~ 0*mental_1 + 1*mental_2 + 28*mental_3 + 53*mental_4 + 72*mental_5 + 70*mental_6 + 52*mental_7
slope_mental ~ 0.015 * 1
int_mental ~ 84 * 1
'

# Simulate:
data <- simulateData(simmodel, model.type = "growth", sample.nobs = 10000)

# Add some random grouping variables:
data$leeftijd_1 <- sample(1:5,nrow(data),TRUE)
data$geslacht_1 <- sample(0:1,nrow(data),TRUE)


# Better to make nice factors with labels for predictors:
data$leeftijd <- factor(
  data$leeftijd_1,
  levels = 1:5,
  labels = c("< 25","25 - 39","40 - 54","55 - 69","70+")
)

data$geslacht <- factor(
  data$geslacht_1,
  levels = c(0,1),
  labels = c("man", "vrouw")
)

# Call the function (dependent = eenzaam, independent = geslacht + leeftijd):
result <- RIVMgrowth_dummy(
  data = data, # Dataset
  type = "non-linear", # Type of analysis to do
  independent = c("geslacht","leeftijd"), # Variable name of predictor (only one supported at the moment)
  time =  c(1, 2, 3, 4, 5, 7, 9), # time vector encoding distances in time
  
  # To indicate the dependent variable, we can use this notation:
  eenzaam = c("eenzaam_1","eenzaam_2","eenzaam_3","eenzaam_4","eenzaam_5","eenzaam_6","eenzaam_7")
)

# Now we can print the results:
print(result)

# The expected values are stored in:
result$expected_values

# And for the marginals:
result$marginal_expected_values

# We can now easily plot the marginals (this returns a ggplot2 object you can edit):
plot(result,independent = "leeftijd", wave_label = "Ronde")
plot(result,independent = "geslacht", wave_label = "Ronde")

# We can also look at correlations:
corplot(result)

### Ordered categorical variables:
# Let's make the dependent variable ordered categorical:
data$eenzaam_1 <- as.numeric(as.ordered(cut(data$eenzaam_1, c(-Inf,0,5,Inf))))
data$eenzaam_2 <- as.numeric(as.ordered(cut(data$eenzaam_2, c(-Inf,0,5,Inf))))
data$eenzaam_3 <- as.numeric(as.ordered(cut(data$eenzaam_3, c(-Inf,0,5,Inf))))
data$eenzaam_4 <- as.numeric(as.ordered(cut(data$eenzaam_4, c(-Inf,0,5,Inf))))
data$eenzaam_5 <- as.numeric(as.ordered(cut(data$eenzaam_5, c(-Inf,0,5,Inf))))
data$eenzaam_6 <- as.numeric(as.ordered(cut(data$eenzaam_6, c(-Inf,0,5,Inf))))
data$eenzaam_7 <- as.numeric(as.ordered(cut(data$eenzaam_7, c(-Inf,0,5,Inf))))

# Run the analysis:
result_ordered <- RIVMgrowth_dummy(
  data = data, # Dataset
  type = "non-linear", # Type of analysis to do
  independent = c("geslacht","leeftijd"), # Variable name of predictor (only one supported at the moment)
  time =  c(1, 2, 3, 4, 5, 7, 9), # time vector encoding distances in time
  ordered = TRUE,
  
  # To indicate the dependent variable, we can use this notation:
  eenzaam = c("eenzaam_1","eenzaam_2","eenzaam_3","eenzaam_4","eenzaam_5","eenzaam_6","eenzaam_7")
)

# Now we can print the results:
print(result_ordered)

# The expected values are stored in:
result_ordered$expected_values

# And for the marginals:
result_ordered$marginal_expected_values

# We can now easily plot these (this returns a ggplot2 object you can edit):
plot(result_ordered, independent = "leeftijd", wave_label = "Ronde")
plot(result_ordered, independent = "geslacht", wave_label = "Ronde")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
