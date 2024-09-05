Performance of model selection methods for the SCAD penalty on Least Squares Regression

Applied various tuning parameter selectors for the smoothly clipped absolute deviation (SCAD) penalty on least squares estimation. 
Empirical data was used to test the model selection methods of generalized cross validation and BIC, with multiple simulations considered in order to assess the performance of these selectors through metrics such as model error and number of correct predictors selected. 
Simulations were conducted in R, and the process of obtaining the tuning parameter using BIC through minimization was also coded manually as no such package providing BIC with this method exists (at this date).
Simulations are based of off population data provided by the following article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2663963/

Article Citation:
Wang H, Li R, Tsai CL. Tuning parameter selectors for the smoothly clipped absolute deviation method. Biometrika. 2007 Aug 1;94(3):553-568. doi: 10.1093/biomet/asm053. PMID: 19343105; PMCID: PMC2663963.
