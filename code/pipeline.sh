# Run from project root directory, typing "bash code/pipeline.sh" at the command line
Rscript code/generate_webs.R 8354201 # Random seed for generating food webs
Rscript code/run_simulations.R 100000 # Number of Monte Carlo iterations
Rscript code/fit_parameter.R
