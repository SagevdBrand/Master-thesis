# Output

### Errors
The erros folder contains .RDS files that were created within the
`error handling.R`. They contain the information on which errors occured during
the situation, both per scenarios specifically, `all_errors_per_scenario.Rds`,
and overall, `all_errors_together.Rds`.

### Estimands

Within this folder one can find the results from the simulation studies.
The zip file contains all the loose results, which are also present in the `all`
folder. Using the `error handling.R` script, this latter folder is turned into
the single file: `all_estimands_batch_7.RDS`.
This is the RDS file used for data analyses and visualization
(`check_results.R`, `visualization.R`, and `performance visualization.R`).

In the `old` folder, one can find the results that were downloaded while the
simulations were still running. Meaning that not all runs for all scenarios are
present within the files.

### Performance
A folder containing the performance measures as obtined from `check_results.R`.
Again, in the old folder, older versions of the results are present. The
`internal_external_by_columns.RDS` is created for nested loop plots in
`performance visualization.R`.
