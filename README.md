# ComBatCorr
This repository contains code for the paper *Overcoming the impacts of two-step batch effect correction on gene expression estimation and inference*.
## Example 
### 1. TB example
- `tbbatchsim.rds`: Simulated data for the TB example.
- `tbdata.R`: R processing code.
- `tbdata.rds`: Gene expression.
- `tbsample.rds`: Sample information.
### 2. Bladderbatch example
- `balanced.rds`: Output for simulation with balanced design.
- `bladderbatch.R`: R processing code.
- `unbalanced.rds`: Output for simulation with unbalanced design. 
### 3. Johnson example
- `dataexample2.txt`: Gene expression. 
- `sampleInfoExample2.txt`: Sample information.
- `johnson.R`: R processing code.
### 4. Towfic example
- `sampleannotation.txt`: Sample information.
- `towfic.R`: R processing code.
### 5. Helper functions
### 6. `simple_example.R`: R code for generating the introduction example. 

## R
- `design_adj.R`: R function for doing ComBatCorr adjustment. 
