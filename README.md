# Predicting antimicrobial resistance phenotypes using a rules-based algorithm and logistic regression 

Antimicrobial resistance phenotypes of clinical *E. coli* and *P. aeruginosa* isolates were predicted using two different algorithms, 1) rules based and 2) logistic regression.

## Metadata of clinical isolates

Raw FASTQs of clinical isolates can be found in [NCBI BioProject PRJNA532924](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA532924)

[Antibiotic susceptibility test results](https://github.com/karatsang/rulesbased_logisticregression/tree/master/AST) coded as Resistant (R), Resistant interpretation by the VITEK 2 Advanced Expert System, even though the MIC should have been intermediate or susceptible (RA), Intermediate (I), Susceptible (S). Please note that ertapenem was not tested on *P. aeruginosa* isolates. 

## Clonality and resistomes of clinical isolates

### Multi-locus sequence typing (MLST)
[Multi-locus sequence typing (MLST) results](https://github.com/karatsang/rulesbased_logisticregression/tree/master/MLST) to test whether either dataset conform to the clonal model of population structure. MLSTs were performed via [comparison to the reference sequences using pubMLST](https://github.com/agmcarthur/pubMLST). Novel MLSTs are new combinations of known alleles and unresolved MLSTs indicate that there is one or more allele of a known clonal complex that was not identified. Clonal complexes are existing combinations of known alleles. Missing MLSTs indicate that one or more alleles could not be identified. 

### Resistance Gene Identifier (RGI)
[Resistance Gene Identifier v4.1.0](https://card.mcmaster.ca/analyze/rgi) and [CARD v2.0.2](https://card.mcmaster.ca/home) were used to identify the collection of all the antibiotic resistance genes in the [*E. coli*](https://github.com/karatsang/rulesbased_logisticregression/tree/master/ResistanceGeneIdentifier/Ecoli_RGI_4.1.0) and [*P. aeruginosa*](https://github.com/karatsang/rulesbased_logisticregression/tree/master/ResistanceGeneIdentifier/Paeruginosa_RGI_4.1.0) datasets. 

Visual overviews of the prevalence of resistance genes and uniqueness of resistomes can be viewed as [heatmaps](https://github.com/karatsang/rulesbased_logisticregression/tree/master/ResistanceGeneIdentifier/Heatmaps_RGI_5.1.0). Heatmaps feature AMR genes categorized by gene family and only unique resistome profiles are displayed with their frequency (histogram at top). Yellow represents a Perfect RGI hit, teal represents a Strict RGI hit, and purple represents no RGI hit. Genes with asterisks (*) appear multiple times because they belong to more than one Drug Class category in the Antibiotic Resistance Ontology (ARO).

### Efflux Pump Identifier (EPI)
[Efflux Pump Identifier v1.0.0](https://github.com/karatsang/rulesbased_logisticregression/tree/master/rulesbased/EffluxPumpIdentifier) was used to identify overexpressed, multi-component efflux pumps in the clinical isolates. 

`python effluxpumpidentifier.py` should be run in a folder of RGI results in .json format. A folder called "EffluxPumpIdentifierResults" will be created with resulting outputs in .json format. 

`python ConvertEPIJsonToTSV.py` should be run in a folder of Efflux Pump Identifier results in .json format to be converted into tab-delimited (TSV) format.

## Predicting antibiotic resistance phenotypes

### Rules-based method
The rules-based method uses RGI and EPI results and traverses CARD's Antibiotic Resistance Ontology to predict antibiotic resistance phenotypes. 

`python 1_predict_phen.py phenotype.json` should be run in a folder with RGI and EPI results in tab-delimited format and a file called "phenotype.json" which has all of the 'confers_resistance_to_drug' relationships from CARD. This will create a *.gp.txt file for every isolate that describes the RGI and EPI results in relation to which antibiotic they confer resistance towards.

`python 2_predict_phen.py ast.tsv` should be run in a folder with all *.gp.txt files created from the previous script (1_predict_phen.py) and an (antibiotic susceptibility testing result file)[https://github.com/karatsang/rulesbased_logisticregression/tree/master/AST]. This will create a file called "phenotype_comparison.txt" that will have a description of AST result and predicted phenotype for every clinical isolate and every tested antibiotic. 

`python 3_phen_summary.py phenotype_comparison.txt` should be run in the same folder as the phenotype_comparison.txt file that was created in the previous step. This will print the true positive, true negative, false positive, and false negative predictions for every antibiotic tested onto the (standard out)[https://github.com/karatsang/rulesbased_logisticregression/blob/master/rulesbased/rulesbased_algorithm/Paeruginosa_rulesbased/rules_based_results.txt]. 

### Logistic regression

`python LR_data_formatter.py rgi_out [folder of RGI results in TSV format] ast_tsv [path to ast.tsv file]` will format the data required for logistic regression.

`python LR_genotypephenotype.py` should be run to perform logistic regression statistical modelling on the clinical isolates. File paths to the ast.tsv and rgi_encoded.pkl (created in previous step) are required to be changed in the script. Area under the curve, receiver operating characteristic, and average precision curves can be generated using the script. Coefficients / weights of importance and their respective p-values generated by logistic regression to predict antibiotic resistance phenotypes can be viewed [here](https://github.com/karatsang/rulesbased_logisticregression/tree/master/logisticregression/weights_pvalues).

## Conda environments

[Conda](https://docs.conda.io/en/latest/) environments were used to run all analyses and they can be found [here](https://github.com/karatsang/rulesbased_logisticregression/tree/master/condaenv). 

`rgi4.yaml`: to run RGI v4.1.0 analyses, but RGI v4.1.0 mut be downloaded from https://card.mcmaster.ca/download and installed as this version is not available on Conda.

`rgi5.yaml`: to generate the heatmaps for visualization of resistomes. 

`rulesbased_logreg.yaml`: to run all remaining analyses, including rules-based and logistic regression statistical modelling.
