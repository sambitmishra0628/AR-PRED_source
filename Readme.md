## Active and Regulatory site PREDiction (AR-PRED)
AR-PRED is a prediction method to detect putative regulatory (allosteric) and active site residues in a given protein. The AR-PRED package contains separate models for predicting these two residue types: active and allosteric. Based on their requirement, the user can use either or both models to make predictions for a given protein.

## Workflow
![ar-pred-approach](https://user-images.githubusercontent.com/6353495/63666447-9d0a4680-c79d-11e9-9012-b51d811d32f7.png)

## Instructions for Running AR-PRED

The code for AR-PRED is written in Perl and MatLab and tested on a 64-bit Linux machine (RHEL, release 6.8). It is designed to be run on a Linux machine of similar or higher configuration. A set of specific instructions is provided below to execute AR-PRED. 

### A. Software
Download and install the following software
1.	Perl version 5 (https://www.activestate.com/products/activeperl/downloads/)
2.	MatLab 2017 or later (https://www.mathworks.com/downloads/)
3.	Naccess (http://wolf.bms.umist.ac.uk/naccess/)
4.	Fpocket version 2 (http://fpocket.sourceforge.net/)
5.	DSSP (https://swift.cmbi.umcn.nl/gv/dssp/)
6.	CD-HIT (http://weizhongli-lab.org/cd-hit/download.php)
7.	BLAST 2.7 and nr protein sequence database (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/)
8.	Rate4Site (https://m.tau.ac.il/~itaymay/cp/rate4site.html)
9.	Clustal Omega (http://www.clustal.org/omega/). Rename the executable to clustalomega.
10.	Download and unzip the MAVEN source code (MAVEN_source.v121.zip)  from https://sourceforge.net/projects/maven/files/Release%201.21/. Store all the source code files in a directory named ‘MAVEN’ (rename the default directory after unzipping from MAVEN_source.v121 to MAVEN).

### B. Environment variables
Download the AR-PRED source code and unzip to a local directory on your Linux work station. Make sure all the AR-PRED code is under 	the directory named ‘AR-PRED_source’ (by default when you unzip, it should be unzipped into a folder with that name). Copy the MAVEN 	directory to the AR-PRED_source directory. Then set the required environment variables and paths using the following steps.

Open the .bashrc file in your home directory using any text editor and add the following lines to the file after excluding the quotes and any special instructions given in parenthesis.

`export ARPRED_HOME='path to AR-PRED_source directory'`

`export MATLAB_PATH='path to matlab executable' (the binary/executable matlab file must be included in the path)`

`export PATH=$PATH:'path to BLAST 2.7 bin directory`

`export BLASTDB='path to BLAST nr database'`

`export NACESS_PATH='path to NACCESS executable/binary'`

`export PATH=$PATH:'path to naccess binary' (the binary file should be named as naccess)`

`export PATH=$PATH:'path to rate4site binary' (the binary executable should be named as rate4site)`

`export PATH=$PATH:'path to cd-hit binary' (the default cd-hit binary is named as cd-hit and it should be left unchanged)`

`export PATH=$PATH:'path to fpocket binary' (the default binary name is fpocket and it should be left unchanged)`

`export PATH=$PATH:'path to dssp' (the binary should be named as dssp)`

### C. Compiling MAVEN \*.c codes with mex
We will need to compile the MAVEN \*.c codes for generating Hessian and adjacency matrices before running AR-PRED. Here are the instructions.

	a. Start Matlab
	b. Within Matlab, navigate to the MAVEN directory inside the ARPRED_HOME path (as set above)
	c. To compile the \*.c codes, we will need to set up mex (if it is not already installed). Type “doc mex” on the matlab command line for instructions on how to set up mex
	d. Compile the c codes with the “mex adjacency.c” and “mex ANMHess.c”

### D. Running AR-PRED
AR-PRED requires an all-atom PDB file as its input. The user should make sure that non-standard or modified amino acids (other than the 20 standard amino acids) are be replaced by their standard names in the PDB file before executing AR-PRED. Otherwise, AR-PRED will report error.

AR-PRED has separate prediction models for active site and allosteric site residues in a given protein chain. Following is the generic command.

	`$ARPRED_HOME/master_script.pl pdbfile  chain outputdir  pred_type activesite_file`
	
	`pdbfile`: Name of the all-atom PDB file. Must be present inside the output directory.

	`chain`: Name of the chain in the PDB file to make predictions for

	`outputdir`: Name of the output directory. Must be created before running the master script.

	`pred_type`: 1 for active site predictions and 2 for allosteric site predictions

	`activesite_file`: Comma separated values (csv) file containing residue IDs in the PDB forming the active site. 
	The residue IDs must be separated by commas and must correspond to the residue IDs in the PDB file. The file name must include
	the path to the file.

For predicting active site residues change directory into the `outputdir` (`cd outputdir`) and run the following command.

	`$ARPRED_HOME/master_script.pl pdbfile  chain outputdir  1`

For predicting regulatory or allosteric site residues however, AR-PRED requires prior knowledge of active site residues. If you already know the residues which form the active site (or a part of the active site) in the chain of interest, then write the residue ids as comma separated values into a .csv (e.g., activesite_res.csv) file and pass the file as an argument to the master_script. In case of no prior knowledge of active site residues, the AR-PRED active site prediction tool or any other active site prediction software may be used. Use the following command to predict allosteric site residues.

`$ARPRED_HOME/master_script.pl pdbfile  chain outputdir  2 activesite_res.csv`

Tip: When running the master script from the output directory, make sure to include the path alongwith the dir name in the `outputdir` argument.

### E. AR-PRED output
AR-PRED outputs its predictions into a .csv file: allostericsite_predictions.csv and activesite_predictions.csv for the allosteric and active site residues, respectively. Each file has two columns: column 1 lists the residue IDs in the PDB file and column 2 lists the weighted probability scores corresponding to each residue ID. Residues having higher scores are more likely to be active or allosteric site residues, depending upon the type of prediction made. 

Note: Depending on the database size and the number of sequences used for calculation of residue conservation scores and due to the random sampling of perturbation points for calculations of dynamic flexibility and perturbation response, there may be minor inconsistency in the residue ranking and scoring by AR-PRED in independent runs of the same protein. However, we have usually observed that the ranking of the top 50 highly probably active or allosteric site residues usually remain the same. To obtain more consistent results, it is advisable to get the average predictions from independent runs of AR-PRED.

In case of any questions or suggestions, please feel free to contact Sambit Mishra (skmishra@iastate.edu). Please include the following citation (Mishra, S. K., Kandoi, G., & Jernigan, R. L. (2019). Coupling Dynamics and Evolutionary Information with Structure to Identify Protein Regulatory and Functional Binding Sites. Proteins: Structure, Function, and Bioinformatics.) if you use AR-PRED for your work.
