##PVCA (Principal Variance Component Analysis) for sequencing read counts


Author: Donghyung Lee


Note: The function is written based on the 'pvcaBatchAssess' function of the PVCA R package 
       and slightly changed to make it more efficient and flexible. 
       (http://watson.nci.nih.gov/bioc_mirror/packages/release/bioc/manuals/pvca/man/pvca.pdf)


Description: Function for Principal Variance Component Analysis


Release Date: 04/18/2016
 
 
Input:

      1. counts: Normalized(e.g. TMM)/raw reads count matrix from sequencing data (row:gene/feature, col:sample) 
               
      2. meta: Meta data matrix containing predictor variables (row:sample, col:predictor)
      
      3. threshold: Minimum amount of variation needs to be explained by top leading PCs
      
      4. inter: TRUE/FALSE - include/do not include pairwise interactions of predictors


Output: 

       a vector of proportions of variation explained by each predictor.
