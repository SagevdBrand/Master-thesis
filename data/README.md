# Data folder

## Folder structure
Within this data there are two folders, `simulation settings` and
`validation data`.
The former is created using the `study scenarios.R` script as found in the
[`src` folder](https://github.com/SagevdBrand/Master-thesis/tree/master/src).
The latter is created by the `create validation data.R` again present with all
other scripts in the
[`src` folder](https://github.com/SagevdBrand/Master-thesis/tree/master/src).

In case you would like to reproduce the results, you would first have to
run the `study scenarios.R` script and then the `create validation data.R`
script. The simulation scenario settings will be saved as `studies.RDS`
within the simulation settings folder.
The creation of the validation data requires this file to create the validation
data.

## Your own simulation study
If you would like to investigate some other scenarios than what has been studied
in the thesis, you can adapt the study scenarios to your liking. However,
please note that scenarios with a different AUC setting or event-fraction
require new regression coefficients for the data generating mechanism. Please
take a look at the `beta optimization script.R` for an example of how I arrived
at the estimates for my regression coefficients.
