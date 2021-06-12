# Master Thesis Sofie van den Brand MSBBSS research master

Version 1.2.
Date 12-06-2021

*Update*: Working on a newer version of the paper. The figure and error counts are no longer in line with the numbers and figures of the manuscript.

*Title*: Clinical Prediction Models: A Comparison Between Cross-Validation and Bootstrap Approaches for Internal Validation

*Abstract*: Clinical prediction models are important tools to support medical decision making. Statistical evaluation through internal validation is an essential step before such a model can be implemented in practice. However, with many internal validation strategies available, in particular various bootstrap and cross-validation variants, it is currently unclear which approach achieves the best out-of-sample performance under varying circumstances. In this thesis three extensive simulation studies are performed, evaluating the impact of outcome prevalence, strength of predictors, and type of modeling algorithm (both regression and tree-based). The ability of the approaches to estimate out-of-sample performance is evaluated in terms of Area under ROC curve, calibration slope and intercept, Mean Absolute Prediction Error, Root Mean Squared Prediction Error, R2_Cox-Snell, R2_Tjur and Estimated Calibration Index. The results illustrate that all internal validation approaches show more variability, compared to large sample external validation, for the AUC and R2 measures, but less for the calibration intercept. Moreover, no single internal validation approach outperforms all others across all investigated scenarios and  performance measures. This indicates that the optimal strategy is dependent on the specific setting and evaluation metric used. The results show that Harrell's bootstrap works generally well when using regression-based models. However, when using tree-based models, which showed a high amount of overfitting, Harrell's, .632 and .632+ bootstrap display large overestimation, while 5, 10 and 10x10 fold cross-validation approach the external validation estimates.

## Project organization
* RO = Read only
* PG = Project generated
* HW = Human work

```
.
├── data                    <- Data about the simulation settings (PG)
│   └── simulation settings <- The simulation studies as defined in the thesis (RO)
├── docs                    <- Documentation notebook for users + supplementary material (HW)
│   ├── manuscript          <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
|   ├── suppl. material     <- A document containing the supplementary information (RO)
│   └── reports             
|       ├── Proposal        <- Research proposal (HW)
|       └── Research report <- Research report (HW)            
├── results                 
│   ├── figures             <- Figures for the manuscript (PG)
│   └── output              <- Other output for the manuscript or reports (PG)
|       ├── estimands       <- Estimands for each simulation run (PG)
|       └── performance     <- Performance measures for each simulation run (PG)           
├── src                     <- Source code for this project (HW)
├── .gitignore              <- Files that are not taken into account when commiting (RO)
├── CITATION.md             <- How to cite my work (RO)
├── LICENSE.md              <- License for my project (RO)
├── README.md               <- The document you are reading right now (RO)
└── requirements.txt        <- Which packages are needed to run the scripts? (RO)

```


## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

You are free to:

    Share — copy and redistribute the material in any medium or format
    Adapt — remix, transform, and build upon the material
    for any purpose, even commercially.

Under the following terms:

   Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

   No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

## Citation

Please [cite this project as described here](/CITATION.md).

## Maintenance
The creation and maintenance of this archive are the responsibility of the author: Sofie van den Brand.
All files are created and saved by the author.
In line with the license, this archive is completely Open Access, so anyone can access the archive for an unspecified amount of time.
