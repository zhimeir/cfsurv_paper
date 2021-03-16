# conformalized_survival_analysis_paper
This repository contains the code to reproduce the numerical results in the [Conformalized Survival Analysis]() paper. 

## Overview
We provide the code to <em>exactly</em> reproduce the numerical results in Section 4 of the paper. Based on the code,
we further develop an R package `cfSurvival` that implements the procedure. The package will be constantly improved and udpated. 
We recommend the users to check out the R package if they want to want to apply our procedure.

## Folders
- `R/`: contains the main functions that implement the variable selection procedures.
- `utils/`: contains the utility functions that supports the experiments. 
- `simulation/`: contains the scripts to carry out the simulations.
- `results/`: stores the result.
- `bash/`: bash files to run the simulations.

## Usage
### Single run
Each script in the `simulation/` folder implements one run of the simulation. The users can specify the random seed when running the script. 
For example, to implement one run of the low-dimensional-homoscedastic-noise experiment in Section 4 with random seed 1, run the following command in your terminal:
```{r}
Rscript ./simulation/ld_homosc.R 1
```

### Multiple runs
The results presented in the paper are averaged over multiple runs. The user can use the bash file in `bash/` to automatically
implement multiple runs of the simulation. But note that it may take a long time if it is run on a laptop. To use the bash file,
run the following code in  your terminal:
```{r}
bash ./bash/run_all.sh
```
