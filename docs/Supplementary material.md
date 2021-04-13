# Supplementary material
This document contains the following material.
1. Numerical estimation procedure for the regression coefficients
used within the data generating mechanism
2. Error handling
3. Scenario settings
4. Details on the software and code used

## 1. Numerical estimation procedure for regression coefficients used within the data generating mechanism
This approach is based on the supplemental material of Van Smeden et al., 2019.
Mathematical calculations are represented in Rcode.
1.	Determine preferred AUC and prevalence
    -  Result	= pref_auc \& pref_prev
2.	Approximate R2_CS  given pref_auc and pref_prev, using the script provided Riley et al., 2020.
    -  Result = pref_R2
3.	Simulate candidate predictors
    -	 Result	= dm (data matrix)
4.	Provide an initial guess for the regression coefficients,
while constraining the coefficients such that 30% is stronger, 50% weak and
20% noise (except if otherwise specified within a scenario). Given this
constraint, only 2 coefficients need optimization: one for the intercept and one
 for all candidate predictors. A strong effect is defined as 3 * &beta;  , weak = &beta;  ,
noise = 0.
    - Result 	= dgm_par (data generating mechanism parameters)
5.	Use logistic regression to determine probability of an event.
This probability will then be used as input for a uniform distribution to
determine the outcomes: <br>
`p <- 1/(1+exp(-dm %*% dgm_par))` <br>
`y <- as.numeric(p > runif(length(p)))`
6.	 Calculate observed AUC and prevalence.
    -  Results = obs_R2 \& obs_prev
7.	Define a loss function: which minimizes the sum of 1) the squared distance
    between preferred R2and observed R2 and 2) the squared distance between
    the preferred prevalence of events and average predicted probability for an
    event in generated data `(intercept = mean(y))`
8.	Use `optim()` to arrive at the best coefficients:
     - Execute the optimizaton 20 times for each combination of prevalence and AUC
   as described in the scenarios. Where R2_CS is approximated for each
   combination as well.
     - At each iteration; the optimized betas provide new predicted probabilities p
   (step 5) and hence a new outcome vector y, on which the AUC and prevalence
   can be calculated as an extra check.
     - Take the median coefficients
     - Validate results on independent validation set N = 100,000
   (generated under same dgm).
9.	Using the just defined regression coefficients we can now simulate the
datasets for each scenario and each iteration. The size of the dataset is
defined by using the ‘pmsampsize’ package, which uses pref_prev and pref_R2 as
defined above as input (Riley et al., 2020).

## 2. Error handling
Below I describe how anticipated errors in estimation and validation on the
simulation data are detected and handled.

### No events sampled:
No or few events may be present in generated simulation data, i.e. degenerate
outcome distribution. It is anticipate that there will be none or a negligible
number of instances where zero events will be sampled in the development data.
For the tuning of lasso and ridge regression, 8 events/non-events are used as a
criterion, following van Smeden et al. (2019).
1. How to detect: Build in checks throughout the script:
  - While sampling the development dataset:
      - Check for `var(Y) = 0`.
  - While resampling the data (bootstrap methods):
      - Check for `var(Y) = 0` in each bootstrap sample. Within cross-validation
      stratified folds are used, which makes sure that there are events within
      each fold. However, an extra check is implemented as well.
  -  While tuning ridge and lasso regression:
      - Check for `sum_i(y_i) > 7` and `sum_i(1-y_i) > 7`

2. How to handle:
  - Development dataset when `var(Y) = 0`:
     - No models will be developed on these data. Performance measures over all
     simulation iterations will be based on complete case analysis. The number of
     instances is documented and reported per simulation scenario.
  - CV fold or bootstrap sample with `var(Y) = 0`:
     - validation will not be performed on this resampling/fold, the remaining
    folds/resamples will be used to calculate performance. The number of
    instances will be documented and reported per study and simulation scenario.
  - When tuning of hyperparameters with not: `sum_i(y_i) > 7` and
      `sum_i(1-y_i) > 7`:
     - When using Lasso or Ridge regression use LOOCV to tune the shrinkage parameter.
     - When using CART or RF, use LOOCV to tune hyperparameters

### No predictors selected:
For models in which predictors are selected (maximum likelihood with backwards
  elimination, Lasso) it is possible that zero predictors are chosen in the
  final model. This leaves the model with only the intercept. This will result
  in problems for the estimation of the calibration slope. With CART it is
  possible that the final model contains no splits (i.e. it bases the predicted
    probabilities on the prevalence of the majority class). New observations
    will have the same predicted probability, yielding an AUC of 0.5 and again
    estimation problems for the calibration slope.
1.	How to detect:
  - Check for each model with either backwards elimination or Lasso whether
    variance of the linear predictor = 0.
  - For CART check the number of splits in the final model.
2.	How to handle:
  - Following Van Smeden et al. (2019), replace the calibration slope by the
    highest estimated calibration slope within the scenario (replacement is
    done after all simulation iterations are completed). The number of instances
    will be documented and reported per model and simulation scenario.

### Separation:
Data separation might occur (probability increases when decreasing sample size
  and increasing dimensionality). We anticipate separation will not occur often
  due to the continuous predictor data and moderate discrimination, AUC of 0.75,
  for all simulation scenarios.
1.	How to detect:
  - Supplemental material of (van Smeden et al., 2019) assumes separation when
    any estimated standard error of maximum likelihood model was > 70 (study 3).
  - Warnings or errors might be returned that indicate separation.
2.	How to handle:
  - Report an error when separation occurs within the results. The number of
    separated datasets will be registered per scenario. No action is needed, as
    models developed on separated data can still be used to make predictions.
    With the low prevalence of separation that is anticipated, we assume it will
    have a negligible influence on performance.

### Probabilities of exactly 0 or 1:
It might occur that predicted probability are exactly 1 or 0. In those cases the
 calibration slope cannot be calculated.
1.	How to detect:
  - Check for predicted probabilities of exactly 1 or 0.
2. How to handle:
  - Report an error when probabilities of 0 or 1 occurred. Change the values
    of these probabilities to 0.000001 or 0.999999, respectively.

## 3. Scenario settings
For a clear overview of all scenarios considered in the thesis, see the table below.
- dim = dimensionality, or number of candidate predictors
- n_setting = the setting for the sample size, either half, double or the minimum sample size required
- noise = the setting for the percentage of noise variables, were the default is 20%, none is 0%, and half 50%
- prev = prevalence or event-fraction
- pred_selection = predictor selection, <0.157 represents Backwards selection using the AIC criterion.
- R2 = R2_CS as obtained through the code by Riley et al., 2020. This is used for calculating the sample size.
- par1 = The coefficient for the intercept that was obtained through the numerical estimation procedure explained above.
- par2 = The value of &beta;, which represent a weak candidate predictor.

|  AUC| dim|n_setting |noise   | prev|model |pred_selection |        R2|    n|       par1|       par2|scenario    |study   |
|----:|---:|:---------|:-------|----:|:-----|:--------------|---------:|----:|----------:|----------:|:-----------|:-------|
| 0.75|  10|n/2       |default | 0.05|ML    |none           | 0.0413198| 1064| -3.3705418|  0.1222722|Scenario_1  |Study_1 |
| 0.75|  10|n         |default | 0.05|ML    |none           | 0.0413198| 2128| -3.3705418|  0.1222722|Scenario_2  |Study_1 |
| 0.75|  10|n*2       |default | 0.05|ML    |none           | 0.0413198| 4256| -3.3705418|  0.1222722|Scenario_3  |Study_1 |
| 0.75|  10|n/2       |default | 0.20|ML    |none           | 0.1265814|  330| -1.6740826|  0.1280781|Scenario_4  |Study_1 |
| 0.75|  10|n         |default | 0.20|ML    |none           | 0.1265814|  660| -1.6740826|  0.1280781|Scenario_5  |Study_1 |
| 0.75|  10|n*2       |default | 0.20|ML    |none           | 0.1265814| 1320| -1.6740826|  0.1280781|Scenario_6  |Study_1 |
| 0.75|  10|n/2       |default | 0.50|ML    |none           | 0.1840704|  219|  0.0140326|  0.1314477|Scenario_7  |Study_1 |
| 0.75|  10|n         |default | 0.50|ML    |none           | 0.1840704|  438|  0.0140326|  0.1314477|Scenario_8  |Study_1 |
| 0.75|  10|n*2       |default | 0.50|ML    |none           | 0.1840704|  876|  0.0140326|  0.1314477|Scenario_9  |Study_1 |
| 0.75|  10|n/2       |default | 0.05|ML    |<0.157         | 0.0413198| 1064| -3.3705418|  0.1222722|Scenario_10 |Study_1 |
| 0.75|  10|n         |default | 0.05|ML    |<0.157         | 0.0413198| 2128| -3.3705418|  0.1222722|Scenario_11 |Study_1 |
| 0.75|  10|n*2       |default | 0.05|ML    |<0.157         | 0.0413198| 4256| -3.3705418|  0.1222722|Scenario_12 |Study_1 |
| 0.75|  10|n/2       |default | 0.20|ML    |<0.157         | 0.1265814|  330| -1.6740826|  0.1280781|Scenario_13 |Study_1 |
| 0.75|  10|n         |default | 0.20|ML    |<0.157         | 0.1265814|  660| -1.6740826|  0.1280781|Scenario_14 |Study_1 |
| 0.75|  10|n*2       |default | 0.20|ML    |<0.157         | 0.1265814| 1320| -1.6740826|  0.1280781|Scenario_15 |Study_1 |
| 0.75|  10|n/2       |default | 0.50|ML    |<0.157         | 0.1840704|  219|  0.0140326|  0.1314477|Scenario_16 |Study_1 |
| 0.75|  10|n         |default | 0.50|ML    |<0.157         | 0.1840704|  438|  0.0140326|  0.1314477|Scenario_17 |Study_1 |
| 0.75|  10|n*2       |default | 0.50|ML    |<0.157         | 0.1840704|  876|  0.0140326|  0.1314477|Scenario_18 |Study_1 |
| 0.75|   6|n/2       |none    | 0.20|ML    |none           | 0.1265814|  198| -1.6705880|  0.2089430|Scenario_1  |Study_2 |
| 0.75|  30|n/2       |none    | 0.20|ML    |none           | 0.1265814|  990| -1.6670190| -0.0534409|Scenario_2  |Study_2 |
| 0.75|   6|n         |none    | 0.20|ML    |none           | 0.1265814|  396| -1.6705880|  0.2089430|Scenario_3  |Study_2 |
| 0.75|  30|n         |none    | 0.20|ML    |none           | 0.1265814| 1980| -1.6670190| -0.0534409|Scenario_4  |Study_2 |
| 0.75|   6|n*2       |none    | 0.20|ML    |none           | 0.1265814|  792| -1.6705880|  0.2089430|Scenario_5  |Study_2 |
| 0.75|  30|n*2       |none    | 0.20|ML    |none           | 0.1265814| 3960| -1.6670190| -0.0534409|Scenario_6  |Study_2 |
| 0.75|   6|n/2       |half    | 0.20|ML    |none           | 0.1265814|  198| -1.6716030|  0.2760249|Scenario_7  |Study_2 |
| 0.75|  30|n/2       |half    | 0.20|ML    |none           | 0.1265814|  990| -1.6234870|  0.0778107|Scenario_8  |Study_2 |
| 0.75|   6|n         |half    | 0.20|ML    |none           | 0.1265814|  396| -1.6716030|  0.2760249|Scenario_9  |Study_2 |
| 0.75|  30|n         |half    | 0.20|ML    |none           | 0.1265814| 1980| -1.6234870|  0.0778107|Scenario_10 |Study_2 |
| 0.75|   6|n*2       |half    | 0.20|ML    |none           | 0.1265814|  792| -1.6716030|  0.2760249|Scenario_11 |Study_2 |
| 0.75|  30|n*2       |half    | 0.20|ML    |none           | 0.1265814| 3960| -1.6234870|  0.0778107|Scenario_12 |Study_2 |
| 0.75|   6|n/2       |none    | 0.20|ML    |<0.157         | 0.1265814|  198| -1.6705880|  0.2089430|Scenario_13 |Study_2 |
| 0.75|  30|n/2       |none    | 0.20|ML    |<0.157         | 0.1265814|  990| -1.6670190| -0.0534409|Scenario_14 |Study_2 |
| 0.75|   6|n         |none    | 0.20|ML    |<0.157         | 0.1265814|  396| -1.6705880|  0.2089430|Scenario_15 |Study_2 |
| 0.75|  30|n         |none    | 0.20|ML    |<0.157         | 0.1265814| 1980| -1.6670190| -0.0534409|Scenario_16 |Study_2 |
| 0.75|   6|n*2       |none    | 0.20|ML    |<0.157         | 0.1265814|  792| -1.6705880|  0.2089430|Scenario_17 |Study_2 |
| 0.75|  30|n*2       |none    | 0.20|ML    |<0.157         | 0.1265814| 3960| -1.6670190| -0.0534409|Scenario_18 |Study_2 |
| 0.75|   6|n/2       |half    | 0.20|ML    |<0.157         | 0.1265814|  198| -1.6716030|  0.2760249|Scenario_19 |Study_2 |
| 0.75|  30|n/2       |half    | 0.20|ML    |<0.157         | 0.1265814|  990| -1.6234870|  0.0778107|Scenario_20 |Study_2 |
| 0.75|   6|n         |half    | 0.20|ML    |<0.157         | 0.1265814|  396| -1.6716030|  0.2760249|Scenario_21 |Study_2 |
| 0.75|  30|n         |half    | 0.20|ML    |<0.157         | 0.1265814| 1980| -1.6234870|  0.0778107|Scenario_22 |Study_2 |
| 0.75|   6|n * 2     |half    | 0.20|ML    |<0.157         | 0.1265814|  792| -1.6716030|  0.2760249|Scenario_23 |Study_2 |
| 0.75|  30|n*2       |half    | 0.20|ML    |<0.157         | 0.1265814| 3960| -1.6234870|  0.0778107|Scenario_24 |Study_2 |
| 0.75|  20|n/2       |default | 0.20|ML    |none           | 0.1265814|  660| -1.6497583| -0.0708246|Scenario_1  |Study_3 |
| 0.75|  20|n         |default | 0.20|ML    |none           | 0.1265814| 1320| -1.6497583| -0.0708246|Scenario_2  |Study_3 |
| 0.75|  20|n*2       |default | 0.20|ML    |none           | 0.1265814| 2640| -1.6497583| -0.0708246|Scenario_3  |Study_3 |
| 0.75|  20|n/2       |default | 0.20|Firth |none           | 0.1265814|  660| -1.6497583| -0.0708246|Scenario_4  |Study_3 |
| 0.75|  20|n         |default | 0.20|Firth |none           | 0.1265814| 1320| -1.6497583| -0.0708246|Scenario_5  |Study_3 |
| 0.75|  20|n*2       |default | 0.20|Firth |none           | 0.1265814| 2640| -1.6497583| -0.0708246|Scenario_6  |Study_3 |
| 0.75|  20|n/2       |default | 0.20|Ridge |none           | 0.1265814|  660| -1.6497583| -0.0708246|Scenario_7  |Study_3 |
| 0.75|  20|n         |default | 0.20|Ridge |none           | 0.1265814| 1320| -1.6497583| -0.0708246|Scenario_8  |Study_3 |
| 0.75|  20|n*2       |default | 0.20|Ridge |none           | 0.1265814| 2640| -1.6497583| -0.0708246|Scenario_9  |Study_3 |
| 0.75|  20|n/2       |default | 0.20|Lasso |none           | 0.1265814|  660| -1.6497583| -0.0708246|Scenario_10 |Study_3 |
| 0.75|  20|n         |default | 0.20|Lasso |none           | 0.1265814| 1320| -1.6497583| -0.0708246|Scenario_11 |Study_3 |
| 0.75|  20|n*2       |default | 0.20|Lasso |none           | 0.1265814| 2640| -1.6497583| -0.0708246|Scenario_12 |Study_3 |
| 0.75|  20|n/2       |default | 0.20|CART  |none           | 0.1265814|  660| -1.6497583| -0.0708246|Scenario_13 |Study_3 |
| 0.75|  20|n         |default | 0.20|CART  |none           | 0.1265814| 1320| -1.6497583| -0.0708246|Scenario_14 |Study_3 |
| 0.75|  20|n*2       |default | 0.20|CART  |none           | 0.1265814| 2640| -1.6497583| -0.0708246|Scenario_15 |Study_3 |
| 0.75|  20|n/2       |default | 0.20|RF    |none           | 0.1265814|  660| -1.6497583| -0.0708246|Scenario_16 |Study_3 |
| 0.75|  20|n         |default | 0.20|RF    |none           | 0.1265814| 1320| -1.6497583| -0.0708246|Scenario_17 |Study_3 |
| 0.75|  20|n*2       |default | 0.20|RF    |none           | 0.1265814| 2640| -1.6497583| -0.0708246|Scenario_18 |Study_3 |

## 4. Details on the software and code used
The following code was used and/or adapted:
-	For generating the data: Approximate R2_CS function Riley et al., 2020
-	For ridge and Lasso: Lambda tuning approach by Van Calster et al., 2020
-	For creating [stratified cross-validation folds and training and testing the model](https://github.com/ledell/cvAUC).
-	For obtaining the AUC: [FastAUC](https://gist.github.com/traversc/1446ebe1dcc2d84bccdca781d4f1fa2a)
-	Multiple functions from Van Smeden's [Beyond EPV simulation study](https://github.com/MvanSmeden/Beyond-EPV)
-	Parts from [here](https://github.com/easystats/performance/blob/master/R/r2_tjur.R) to obtain Tjur’s R2

The following describes all information about the HPC used,
and packages and dependencies that were used to obtain the results:

**R version 4.0.2 (2020-06-22) <br>
Platform: x86_64-pc-linux-gnu (64-bit) <br>
Running under: CentOS Linux 7 (Core)** <br>

**attached base packages:**
- stats
- graphics
- grDevices utils
- datasets
- methods
- base     

**other attached packages:**
- kableExtra_1.3.1 
- ranger_0.12.1
- rpart_4.1-15
- caret_6.0-86
- glmnetUtils_1.1.8
- glmnet_4.1-1
- Matrix_1.2-18
- logistf_1.24
- glue_1.4.2
- rms_6.2-0
- SparseM_1.81
- Hmisc_4.4-1
- Formula_1.2-3
- survival_3.1-12
- lattice_0.20-41   
- MASS_7.3-53.1
- pmsampsize_1.0.3
- forcats_0.5.0
- stringr_1.4.0
- dplyr_1.0.2
- purrr_0.3.4
- readr_1.3.1
- tidyr_1.1.2
- tibble_3.0.3
- ggplot2_3.3.2    
- tidyverse_1.3.0  