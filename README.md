# Protein function and subcellular localization prediction

Proteins directed toward the **cell secretory pathway** (ER-Golgi-Membrane-Extracellular) are endowed with a **signal sequence** in the N-terminal region.
The **signal sequence is cleaved** after the protein reaches its final destination
*In silico* recognition of the presence of the signal peptide is a key step for the characterization of protein function and subcellular localization.

## Signal Peptide (SP) prediction ##
* Recognized to direct the protein toward the secretory pathway
* Modular “structure” comprising three separate regions with different properties
* Hydrophobic core similar to a transmembrane helix
* Weakly conserved sequence motifs can be observed at cleavage sites

 ![image](https://github.com/ibojovic/LB2_2022/assets/62520977/0a24407f-ad81-4c2c-a869-9962f4897eba)

Two distinct prediction tasks:
1. Discrimination/detection: Detect the presence of the SP sequence
2. Labelling: Identification of the precise position of the cleavage site along the sequence

# Project setting #

* Collection and organization of training/benchmark data
* Visualization of statistics on both datasets
* Feature extraction
* Implementation of the von Heijne algorithm
* Implementation of the SVM-based approach (using sklearn)
* In cross-validation and blind test

## Training and benchmark datasets ##

* Derived from **SignalP-5.0 datasets** (Almagro Armenteros et al, 2019)
-Extracted from UniProt Knowledgebase release 2018_04
-Only reviewed entries (SwissProt), hypothetical proteins were not included
-Protein sequences shorter than 30 AAs were discarded
-Only signal peptides with experimental evidence (ECO: 0000269) for the cleavage site included

* Training dataset: **1723 eukaryotic sequences**
-258 positive examples i.e., sequences endowed with N-terminal secretory signal peptides
-1465 negative examples i.e, proteins with a subcellular location annotated as cytosolic, nuclear, mitochondrial, plastid, and/or peroxisomal in Eukarya and not belonging to the secretory pathway with experimental evidence
-A randomly selected subset of the original eukaryotic SignalP-5.0 training dataset

* Benchmark set: **7456 eukaryotic sequences**
- 209 and 7247 positive and negative examples, respectively
- Exactly the same dataset used for benchmarking SignalP-5.0 and other approaches
* **All sequences were shortened to the first 50 N-terminal residues**

***These datasets are provides as files named: training_set.tsv and benchmark_set.tsv***

### The training set file (tab-separated) ###
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/009e5e89-843e-4017-a019-1a6116a1dbe4)

Fields are:
- UniProtKB accession
- Taxa: the species to which the protein belongs to
- Kingdom: the kingdom of the species, can be Metazoa, Plants, Fungi or Other
- Class: whether the protein has the SP, can be SP or NO_SP
- Cross-validation fold: the cross-validation subset to which the protein belongs to (better clarified later)
- The amino acid sequence of the first 50 N-terminal positions
- The annotation of the SP along the sequence. S represents position belonging to the signal sequence, N other positions. For sequences in class NO_SP, all positions are labelled as N

### The benchmark set file (tab-separated) ###
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/fd3ea7c0-c021-461f-b844-c65427af9c64)

Fields are the same as for the training_set, except for cross-validation fold which is missing here

## Cross validation split ##
* Data were randomly split into 5 equally-sized different subsets
* Training data are already non-redundant
-No pair of sequences above similarity threshold is present in the dataset
* Split has been performed randomly on non-redundant data
* This avoids biases during execution of CV, ensuring no similarity is present between training and testing sequences

# Statistical analysis of the datasets #
* Basic statistics on training/benchmark datasets
* Assessing the adequacy of the dataset 

Following plots are generated :
* The distribution of SP lengths
* Comparative amino-acid composition of SPs against some background distribution (amino acid composition of SwissProt available at https://web.expasy.org/docs/relnotes/relstat.html) 
* Taxonomic classification (at kingdom and species levels)
* Sequence logos of SP cleavage sites created with [WebLogo](https://weblogo.berkeley.edu/logo.cgi)

# The vonHeijne method for SP detection #


![image](https://github.com/ibojovic/LB2_2022/assets/62520977/1edb70e1-a48e-45e6-b81e-a78a61a35bb1)

## Position-Specific Weight Matrices (PSWM) ##
* A way of representing patterns or motifs in biological sequences (DNA/RNA or proteins)
* The number of rows is equal to the number of different characters in the alphabet (20 for proteins, 4 for nucleotide sequences)
* The number of columns is equal to the length of the motif
* The starting point is a set of aligned/stacked sequence fragments of length *L*

Given a set *S* of *N* aligned sequences of length *L*, the PSPM *M* is computed as follows:
M<sub>k,j</sub> = (1/N) ∑<sub>i=1</sub><sup>N</sup> I(s<sub>i,j</sub> = k)

Where:
* s<sub>i,j</sub> is the observed residue of aligned sequence *i* at position *j*
* *k* is the residue corresponding to the *k-th* row in the matrix
* I(s<sub>i,j</sub> = k) is an indicator function (1 if the condition is met, 0 otherwise)

From the PSPM *M*, the PSWM *W* is computed as follows: 
W<sub>k,j</sub> = log(M<sub>k,j</sub> / b<sub>k</sub>)

Where:
* b<sub>k</sub> is the frequency of residue type k in the background model

  ## PSWM: meaning of the values ##
Values W<sub>k,j</sub> can be:
* Positive, when M<sub>k,j</sub>> b<sub>k</sub>: the probability of residue *k* at position *j* in the motif differs from the background and it is higher (more likely to be an important/functional site than random)
* Zero or negative, when b<sub>k</sub>>= M<sub>k,j</sub>: the probability of residuekat positionjin themotif differs from the background and it is lower (more likely to be a random site than a functional one)

  ## PSWM: background model ##
* The simplest background model is the uniform distribution i.e. assuming all the residues as equally frequent with probability 1/20=0.05
* In general, any background distribution can be adopted e.g., the overall AA composition computed on the SwissProt database (https://web.expasy.org/docs/relnotes/relstat.html)

  ## PSWM: pseudocounts ##
* When dealing with datasets of finite size some of the counts may be zero i.e., not all residues are observed in some positions along the motifs
* In order to avoid zero probabilities in the PSPM and hence the impossibility of computing the log-odds, pseudocountsare added during computation of PSPM
* In the most simple setting, the count matrix is initialized assuming each residue is observed at least once in all positions
* Formally, the formula for computing the PSPM *M* becomes:
M<sub>k,j</sub> = (1 / (N + 20)) * (1 + ∑<sub>i=1</sub><sup>N</sup> I(s<sub>i,j</sub> = k))

  ## PSWM: scoring motifs on new sequences ##
* Given any piece of sequence X=(x<sub>1</sub>,..,x<sub>L</sub>) of length *L*,  compute the (log-likelihood) score of *X* given the PSWM *W* as:
score(X|W) = ∑<sub>i=1</sub><sup>L</sup> W<sub>x<sub>i</sub>, i</sub>
* The values of the scores can be used to scan a protein sequence in order to find regions where the likelihood of occurrence of the motif represented by *W* is higher

## Scoring the presence of SPs in new sequences ##
* Extract the first 50 N-terminal positions from the complete protein sequence
* Scan positions from 1 to 36 (=50-15+1) with the PSWM (extracting all subsequence of 15 residues), obtaining a score for each position
* Compute the global score for SP detection as the maximum positional score along the sequence: the higher the global score the higher is the likelihood that the protein has a SP

### Classification ###
Optimizing the value of the threshold in cross-validation for each training run:
* At each training run, the value of the threshold is estimated using the training set and fixed for predicting the testing set
* The optimal value for each training set is obtained as the one maximizing some performance metric on training data e.g. using precision-recall values at varying prediction thresholds

###Implement the code of the vonHeijne method for SP detection###
* Cleavage-site context: (-13,+2)
* Pseudocounts:1
* Background model: SwissProt
* Threshold selection: cross-validated (using PR curve)

# Support Vector Machine for SP detection #

## SVMs with scikit-learn ##
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/3f2d4ac2-9b79-4865-ac4b-a897c00ff00d)
SVC constructor parameters
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/d84833a8-079c-4d14-b684-a7fec59d512e)
SVC class methods
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/7a8898b9-5f37-47c5-b3e1-a9cbf9c89746)
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/aa9d0f2b-e7b2-4fbc-9c5d-3c2c58d14a74)
SVC attributes
![image](https://github.com/ibojovic/LB2_2022/assets/62520977/6dadfeb7-abeb-4e2d-ade2-0aae4b8f1e2a)
**WARNING:** attributes have defined values only after training

## Input data ##
Training data consist of:
* A matrix X containing all example vectors (in rows)
* A vectory containing target classes (one for each vector in X)
* The number of rows in X = the number of elements in y

Trained models are saved/loaded to/from  compressed files using the pickle and gzip Python modules.

## SVM for SP detection ##
* Each sequence is  represented as a point in some D-dimensional space
* By combining the *SP length distribution* and *SP sequence composition analysis*  we can encode each protein by the composition of its N-terminus

  ## SP length distribution ##
   * Median length of SPs in the training set is 21
   * Mean length of SPs in the training set is about 22±5
 
  **We expect most of the SPs having a length between 19 to 22**

  ## SP sequence composition ##
   * SP compositions differ from the background SwissProt distribution
   * Non-polar residues are abundant in SP
   * Polar and charged residues are less abundant
 
  **Composition of the N-terminal region can be used to distinguish SPs from non-SPs**


* Each protein is encoded with a 20-dimensional vector corresponding to the composition of the first K residues in the protein
* K is an hyperparameter that need to be optimized
- Using the expected SP length we test different values around the mean length (e.g. all the values from 20 to 24)

* Use non-linear SVMs with RBF kernel
- C parameter is optimized
- RBF kernel parameter gamma is optimized

* For the three hyperparameters (K, C and gamma)  a grid-search with cross-validation is run
* Considered parameters:
   - K= {20, 22, 24}
  - C = {1, 2, 4}
  - gamma = {0.5, 1, “scale”}  ( Gamma set to “scale” in SKlearn means setting the value to: 1 / (n_features * var(X)) )
* Each value in the grid is obtained after a complete cross-validation with SVM models trained using the corresponding combination of K, C and gamma.
* After CV the average MCC is computed
* The best combination is selected as the one providing the highest MCC:
  (<i>K'</i>, <i>C'</i>, <i>γ'</i>) = argmax<sub><i>K</i>,<i>C</i>,<i>γ</i></sub> MCC<sub><i>K</i>,<i>C</i>,<i>γ</i></sub>

## Evaluating binary classifiers (quantitatively) ##
* Binary confusion matrix
* Main binary scoring indexes: Sensitivity, Precision, Accurarcy, PPV, MCC

## Cross-validation ##
* After each CV  performance scores are comouted  on the current testing set
* Once the entire cross-validation procedure is completed the final value for each score is computed as the average value over the five CV runs
* Average values are accompanied by standard errors, defined as:
  SE = σ / √n
where:
    - σ is the standard deviation
    - n is the number of samples (5 in this case)

## Qualitative evaluation of results ##
Using the benchmarking data analyses:
* False positive (FP) predictions: non-SP sequences predicted as SPs
* False negative (FN) predictions: SPs predicted as non-SPs

  ### False positive analysis ###
  To understand the impact of this situation we isolated FP predictions in the benchmark dataset
   * Many errors are due to the presence of an hydrophobic transmembrane helix at the protein N-terminus
  ### False negative analysis ###
  * For the von-Heijne method, FNs are due to a different composition of the cleavage-site context than expected
   - For all FN in the benchmark, we computed the corresponding LOGO and compared it with the expected one
 
  * For the SVM method
    - FN may be due to a different composition of the N-terminal region (with the selected K) than expected
    - FN may be due to the presence of a longer/shorter SP than expected



