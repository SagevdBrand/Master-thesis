# Master-thesis
*Title: The Battle of Internal Validation within Medical Prediction Models: Bootstrap vs. Cross-Validation*

**Supervisors: Maarten Van Smeden & Ben Van Calster**

The initial research proposal and a report based on the work done during the first half of the thesis can be found in their respective folders.
In the R folder the code to run simulations can be found.


### To do list:

**Within scenario definitions**
Concerning the `Scenarios.R` script:
[ ] Come up with solution, maybe a function or something, to get sample size and beta coefficients into scenario matrix
[ ] Save scenario matrices once they have been finalized


**Within Data generation**
Concerning the `Data generation scenario 1.R`:
[ ] Save the generated data in a folder per sub-scenario

Concerning the `Data generation functions.R`:
[ ] Create a test script for beta optimization
[ ] Add options for different predictor effects in beta optimization
[ ] Add options for different number of predictors within beta optimization
[ ] Add options to sample from different distributions, other than multivariate normal
[ ] Check convergence of optimization, especially when changing effect and number of predictors
[ ] Add a check for how many events have been sampled in `generate_data` function
