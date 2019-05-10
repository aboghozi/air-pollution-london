import pandas as pd  
import numpy as np
import matplotlib.pyplot as pl
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C


np.random.seed(1)

##===================================================================================
## read in data
## collapsed over time
dta = pd.read_csv("../data/kcl_london_model_data_winter_collapsed.csv", sep=',') 
print(dta.head())

#dta_nowinter_col <- readRDS("data/kcl_london_model_data_nowinter_collapsed.rds")

## aggregated over time
#dta_winter_agg <- readRDS("data/kcl_london_model_data_winter_agg_time.rds")
#dta_nowinter_agg <- readRDS("data/kcl_london_model_data_nowinter_agg_time.rds")

#dta_monthly_agg <- readRDS("data/kcl_london_model_data_monthly.rds")
##===================================================================================

## divide into features and variable
X = dta.iloc[:,0:5].values  
y = dta.iloc[:,5].values  

print(y[0:10])
print(X[1:10,:])



def f(x):
    """The function to predict."""
    return x * np.sin(x)

# now the noisy case
X = np.linspace(0.1, 9.9, 20)
X = np.atleast_2d(X).T

# Observations and noise
y = f(X).ravel()
dy = 0.5 + 1.0 * np.random.random(y.shape)
noise = np.random.normal(0, dy)
y += noise

# Instantiate a Gaussian Process model
kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
gp = GaussianProcessRegressor(kernel=kernel, alpha=dy ** 2,
                              n_restarts_optimizer=10)

# Fit to data using Maximum Likelihood Estimation of the parameters
gp.fit(X, y)

# Make the prediction on the meshed x-axis (ask for MSE as well)
#y_pred, sigma = gp.predict(x, return_std=True)




# Test data
#n = 50
#Xtest = np.linspace(-5, 5, n).reshape(-1,1)

# Define the kernel function
#def kernel(a, b, param):
#    sqdist = np.sum(a**2,1).reshape(-1,1) + np.sum(b**2,1) - 2*np.dot(a, b.T)
#    return np.exp(-.5 * (1/param) * sqdist)

#param = 0.1
#K_ss = kernel(Xtest, Xtest, param)
# Get cholesky decomposition (square root) of the
# covariance matrix
#L = np.linalg.cholesky(K_ss + 1e-15*np.eye(n))
# Sample 3 sets of standard normals for our test points,
# multiply them by the square root of the covariance matrix
#f_prior = np.dot(L, np.random.normal(size=(n,3)))

# Now let's plot the 3 sampled functions.
#pl.plot(Xtest, f_prior)
#pl.axis([-5, 5, -3, 3])
#pl.title('Three samples from the GP prior')
#pl.show()
