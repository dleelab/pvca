PVCA (Principal Variance Component Analysis) for sequencing read counts
-----------------------------------------------------------------------

Author: Donghyung Lee

Note: The function is written based on the 'pvcaBatchAssess' function of the PVCA R package and slightly changed to make it more efficient and flexible. (<http://watson.nci.nih.gov/bioc_mirror/packages/release/bioc/manuals/pvca/man/pvca.pdf>)

Description: Function for Principal Variance Component Analysis

Input:

      1. counts: normalized(e.g. TMM) or log-transformed reads count matrix from sequencing data (row:gene/feature, col:sample) 
               
      2. meta: meta data matrix containing predictor variables (row:sample, col:predictor)
      
      3. threshold: proportion of the variation in read counts explained by k top PCs. This value determines the number of PCs to be used in pvca. 
      
      4. inter: TRUE/FALSE - include/do not include pairwise interactions of predictors

Output:

       a vector of proportions of the variation in read counts data explained by each predictor.
