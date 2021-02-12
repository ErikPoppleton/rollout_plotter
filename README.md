# rollout_plotter
plotting scripts for outputs from the rollout/RNAfold algorithm

These scripts provide analysis and visualization capabilities for [ExpertRNA](https://github.com/MenghanLiu212/RL-RNA), a dynamic programming approach to RNA secondary structure prediction that integrates multiple folding and evaluation algorithms.  The ExpertRNA implementation linked produces as output a csv with the following columns:


| column number | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 |
|--------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
|**Data**|Number|Name|Sequence|Actual Structure|RNAfold Structure|ExpertRNA Structure|RNAfold Hamming Distance|ExpertRNA Hamming Distance|ent_3|GC Perentage|Ensemble Diversity|Expected Accuracy|fe_per|Running Time(sec)|Actual Foldability MFE|Actual Foldability NFE|Actual abs Foldability difference|RNAfold Foldability MFE|RNAfold Foldability NFE|RNAfold abs Foldability difference|ExpertRNA Foldability MFE|ExpertRNA Foldability NFE|ExpertRNA abs Foldability difference|Actual FE|RNAfold FE|ExpertRNA FE

## Data sanitization
These scripts assume that ExpertRNA is run with 4 branches and that all 4 branches successfully predicted a structure.  However, this does not always happen.  For MCC.py to not get offset (since it is reading every 4 lines of the output files), the .csv files must be padded with "NONE,NONE" in any cases where insufficient structures were produced.  The following shell command will highlight the lines where the NONE,NONE must be added:

`cat name.csv | cut -d ',' -f2 | uniq -c | tr -s ' ' | awk '{if ($1 != 4) {print}}'`

In the Mathews dataset, there is one sequence, srp_Crit.fasc.\_AY781797, which does not have a secondary structure associated.  This crashes MCC.py and the lines corresponding to this structure must be deleted.

## Script descriptions

**MCC.py**: Reads all .csv files in the current directory and calculates the [Matthews Correlation Coefficent](https://en.wikipedia.org/wiki/Matthews_correlation_coefficient) for every prediction compared with the actual structure.  This script also produces a large number of figures and statistics on the aggregate MCC data.  The user is encouraged to look through the plots and comment out those that they do not need.  The following files are produced:<br>
Three .fasta files for each line in the input .csv files: These files contain the structure name, the sequence, and each of the actual structure, the RNAfold prediction and the ExpertRNA prediction.  
Two .col files for each line in the input .csv files: These files contain a [Forna](http://rna.tbi.univie.ac.at/forna/) colormap marking which nucleotides were correctly predicted and which were incorrect.  
Five score files: One for each branch and one with the aggregate data (labeled 0) These have the name, and the MCC scores for the RNAfold prediction, and the ExpertRNA prediction and finally the difference between the two predictions.  This file is sorted from worst to best difference.

**forna_generator.py**: Reads all the .fasta and .col files produced by MCC.py in the present directory and spits out [Forna](http://rna.tbi.univie.ac.at/forna/) links for every structure.  These can be pasted in a web browser to view the actual, and predicted structures.  The predicted structures are annotated with red for incorrectly-predicted nuclotides and green for correct predictions.

**extract_RNAfold_failures.py**, **fix_naming.py**, and **ranking.py** are old scripts that were not used in the final paper.

