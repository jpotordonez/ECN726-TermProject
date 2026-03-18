The file data.csv contains the following variables:

1. Y: log of real per capita GDP
2. G: log of real per capita federal government spending
3. T: log of real per capital federal tax revenue
4. GSHK_R: instrument for government spending shocks, from Ramey and Zubairy (2018)
5. T90_MMO: instrument for tax revenue shocks, from Mertens and Montiel Olea (2018)
6. TFPSHK_FRBNY_P: instrument for non-fiscal shocks, from the NY Fed DSGE model
7. TBILL3: interest rate on three-month Treasury bills

The file main.R replicates the empirical results in the paper. This draws on auxiliary files stored in the folder Aux_files as well as the data.csv file.

The file monteCarlo.R replicates the Monte Carlo experiments from the paper.
