# Welcome to Enrichment analysis project

In this GitLab project, you will find all files and documentation related to E.coli enrichment analysis project.

:warning: You are on the branch "Vic", you will find here the documentation regarding the work done by Victoria FATHI :warning: 

## Organization

- **benchmark_data/** : all data related to the benchmark
- **benchmark_results/**: results from bencharmking analysis
- **raw_data**/ : all data related to E.coli retrieved form UniProt and EcoCyc 
- **import/**: all files needed to create the data base
- **query_sets/**: query sets for enrichment analysis
- **results_set20/**: results from set20 analysis
- **data_integration_procedure.hmtl** : procedure used to integrate E.coli data in Neo4j database.
- **get_data.Rmd** : procedure used to prepare data for integration
- **get_benchmark_data.r** : procedure used to prepare benchmark data for analysis
- **benchmark_analysis.Rmd** : procedure used for benchmark analysis
- **set.M2.20.txt**:: unknown set to analyze
- **idh_env.yml** : to recreate conda environment

## How to navigate between branches

:warning: By switching branches, the documentation (README.md) will change according to work done by each student. :warning:  

### Directly from GitLab

You can switch between branches directly from GitLab. This allows you to preview selected branch files and access the documentation (README.md). 
![](https://i.imgur.com/1C1CQD0.png)

### From terminal

Once the project has been cloned, you can between switch branches directly from a terminal.
```
git checkout <branch_name>
```

You can also look in witch branch you are.
```
git branch -a
```

This is an example of all the steps to access the files of 'vic' branch for Part 2.
```
git clone git@gitlab.com:best-partnership/m2idh.git
cd m2idh/Part_2
git checkout vic
```
## Environment 
Note that we used conda 4.10.1 during this project and **idh_env.yml** is available on the gitlab project to recreate the environment.

## Documentation

# Database creation - get_data.Rmd documentation

Required packages:
- Tidyverse 1.3.1

To reproduce data parsing you will need to change the working directory: setwd('/<path to Part2>/Part_2/')


# Enrichement analysis - Enrichment.py documentation

Required packages:
- numpy 1.19.5
- scipy 1.5.3 
- py2neo 2021.2.3 
- pandas 1.1.5


## How to use
```tiddlywiki
enrichment.py -p [PASSWORD] -q [QUERY] -t [NODE_TYPE] -m [METRICS] [OPTION]
                     
-p, --password [required] password neo4j database
-q, --query [required] path to query set
-t --target_type [required] Target sets node type [GOTerm, Interpro, Keyword,Pathway, PubMed, TU"]
-m --measure [required] Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage
-s --species taxon id (default=511145)
-a --alpha [optional] Significance threshold (default = 0.05)
-c, --adjust [optional] Adjust for multiple testing (FDR)
-l --limit [optional] Maximum number of results to report.
-w, --write [optional] path and name for results (tsv file)
-v, --verbose [optional] print intermediary results/queries
```

## Example of use

```typescript
./enrichment.py -p password -q query_sets/set01.txt -t GOTerm -s 511145 alpha 0.2 -c -m chi2 -l 10
```
This will print the first 10 enriched GOTerm in set01.txt, using Chi-Square Test of Independence , with a significance threshold of 0.2 and using adjust for mutliple testing

```typescript
./enrichissement.py -p password -q query_sets/set01.txt -t Keyword -s 511145 -e -l 5
```
This will print the first 5 enriched Keyword in set01.txt with the 4 methods.

```
./enrichment.py -p pwd -q benchmark_data/big_set.txt -t GOTerm -m binomial -c -w benchmark_results/big_GOTerm_binomial.tsv
```
