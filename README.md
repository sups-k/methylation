# methylation

Pipeline:

1. Run the code "detP.R" for all the samples.
  a) Make sure the following 3 libraries with dependencies are installed:
    1) minfi
    2) data.table
    3) ewastools
  
  b) Change the base directory to the directory where the samples and the sample sheet are stored. This is on line 10 of the code. The variable name is "baseDir".
  
  c) On line 48, specify the path (file name is given, don't change it, only change the path) of the PDF to be generated. The line looks like "pdf(...)".
  
2. After you run "detP.R", you will obtain a PDF file that contains a plot of how many Y chromosome probes are excluded for both males and females (Y- axis) for a given detection p-value (X-axis). Based on the graph, choose an appropriate p-value cutoff and enter it in step 5.

3. Run the code "sex_pred.R" for all the samples. Make sure the samples contain both sexes otherwise the code will not work.
  a) Make sure the following 3 libraries with dependencies are installed:
    1) minfi
    2) IlluminaHumanMethylation450kmanifest
    3) IlluminaHumanMethylation450kanno.ilmn12.hg19
    
  b) Change the base directory to the directory where the samples and the sample sheet are stored. This is on line 18 of the code. The variable name is "baseDir".
  
  c) Change the path of the CSV file on line 31 (last line). Enter the full path.
  
4. Make sure that the predicted sex matches the given sex. In case of a mismatch, change the value of the sex in the sample sheet with the value predicted.

5. Before running the main code "meth_pipeline_v1.R", make the following changes:
  a) Make sure the following libraries with dependencies are installed:
    1) missMethyl
    2) limma
    3) minfi
    4) IlluminaHumanMethylation450kmanifest
    5) IlluminaHumanMethylation450kanno.ilmn12.hg19
    
  b) Change the base directory to the directory where the samples and the sample sheet are stored. This is on line 62 of the code. The variable name is "baseDir".
  
  c) Set the paths (the file names are given approriately, just change the path) in:
    1) line 100 that starts with "qcReport(..., pdf="file name with full path")"
    2) line 107 that calls function "performQC"
    3) line 122 that calls function "performQC"
    4) line 126 calls function "performQC"
    5) line 129 that starts with "pdf(....)"
    6) line 148 that starts with "pdf(....)"
    6) line 174 that starts with "pdf(....)"
    7) line 246 that starts with "pdf(....)"
    8) line 282 that starts with "pdf(....)"
    9) line 292 that starts with "write.csv(....)"
    
  d) IN LINE 116, ENTER THE P-VALUE AS CHOSEN FROM THE PDF IN STEP 2.
  
  For entire pipeline, visit https://docs.google.com/document/d/1SQIegq0TmIWoD1OtN66pI_apQJ3qLzUdrfwHtjXe4mQ/edit?usp=sharing
