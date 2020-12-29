set.seed(6)

# Data values
n_1 <- c(5,6,4)
n_2 <- c(5,4,6)
y <- c(7,5,6)

theta_1_values <- seq(0,1,0.01) # Generate values for theta
theta_2_values <- seq(0,1,0.01) # Generate values for theta

multiNomial <- function(theta_1_values,n_1,theta_2_values,n_2,y){ # Likelihood
  
  first_likelihood <- 0 # i = 1 iteration
  for (j in max(0,y[1]-n_2[1]):min(n_1[1],y[1])){ # summing together the product of dbinom
    first_likelihood  <- first_likelihood  + 
      (dbinom(j, n_1[1],theta_1_values)*dbinom(y[1]-j, n_2[1],theta_2_values))
  }
  
  second_likelihood <- 0 # i = 2 iteration
  for (j in max(0,y[2]-n_2[2]):min(n_1[2],y[2])){ # summing together the product of dbinom
    second_likelihood  <- second_likelihood  + 
      (dbinom(j, n_1[2],theta_1_values)*dbinom(y[2]-j, n_2[2],theta_2_values))
  }
  
  third_likelihood <- 0 # i = 3 iteration
  for (j in max(0,y[3]-n_2[3]):min(n_1[3],y[3])){ # summing together the product of dbinom
    third_likelihood  <- third_likelihood  + 
      (dbinom(j, n_1[3],theta_1_values)*dbinom(y[3]-j, n_2[3],theta_2_values))
  }
  
  # Take the product of the likelihoods. I use the negative log likelihood.
  
  likelihood <- -1*(first_likelihood*second_likelihood*third_likelihood) 
}

# Create an empty matrix to store the likelihoods
likelihood.matrix<-matrix(nrow=length(theta_1_values),ncol=length(theta_2_values))


# Compute the likelihood and add the values to the matrix
for (m in 1:length(theta_1_values))
{
  for (p in 1:length(theta_2_values))
  {
    likelihood.matrix[m,p]<-multiNomial(theta_1_values[m],n_1,theta_2_values[p],n_2,y)
  }
}

mle_element <- which(likelihood.matrix==min(likelihood.matrix),arr.ind=T)
theta_1_values[mle_element[1]]
theta_2_values[mle_element[2]]

# Use a contour plot to plot the 2D likelihood. 
contour(theta_1_values,theta_2_values,likelihood.matrix,nlevels=30,xlab='Theta One',
        ylab='Theta Two',
        main='Likelihood Contour Plot')