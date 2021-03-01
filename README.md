# Master thesis internal validation

Version 0.1.0

The battle of internal validation The Battle of Internal Validation within Medical Prediction Models: Bootstrap vs. Cross-Validation


## Project organization
* RO = Read only
* PG = Project generated
* HW = Human work

```
.
├── .gitignore              <- Files that are not taken into account when commiting (RO)
├── CITATION.md             <- How to cite my work (RO)
├── LICENSE.md              <- License for my project (RO)
├── README.md               <- The document you are reading right now (RO)
├── requirements.txt        <- Which packages are needed to run the scripts? (RO)
├── bin                     <- Compiled and external code, ignored by git (PG)
│   └── external            <- Any external source code, ignored by git (RO)
├── config                  <- Configuration files (HW)
├── data                    <- All project data, ignored by git
│   ├── simulation data     <- Folder where the simulated data is temporarily 
|   |                          stored for each run per study (PG) 
│   ├── simulation settings <- The simulation studies as defined in the thesis methods (RO)
│   └── validation data     <- Folder where the simulated validation data is temporarily
|                              stored for each run (PG) 
├── docs                    <- Documentation notebook for users (HW)
│   ├── manuscript          <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports             
|       ├── Proposal        <- Research proposal (HW)
|       └── Research report <- Research report (HW)            
├── results                 
│   ├── figures             <- Figures for the manuscript (PG)
│   └── output              <- Other output for the manuscript or reports (PG)
|       ├── estimands       <- Estimands for each simulation run (PG)
|       └── performance     <- Performance measures for each simulation run (PG)           
└── src                     <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)

## Citation

Please [cite this project as described here](/CITATION.md).
