Welcome to Enrichment analysis project (branch : lou)
===

**Authors : Lou Duron**

# Descritpion
In this GitLab project, you will find all files and documentation related to E.coli enrichment analysis project.

:warning: You are on the branch "lou", you will find here the documentation regarding the work done by Lou Duron :warning:

# Organisation

- **benchmark_data** : all data related to the benchmark
- **import** : all files needed to create the data base
- **query_sets** : query sets for enrichment analysis
- **raw_data** : all data related to E.coli retrieved form UniProt and EcoCyc 
- **CR_projet_IDH_Lou_Duron.pdf** : project report
- **Enrichment.py** : modified script for enrichment analysis
- **Intégration_des_données_Neo4j.hmtl** : procedure used to integrate E.coli data in Neo4j database.
- **get_data.Rmd** : procedure used to prepare data for integration
- **get_benchmark_data.r** : procedure used to prepare benchmark data for analysis
- **benchmark_analysis.Rmd** : procedure used for benchmark analysis

# Enrichment.py documentation

Required packages :
- scipy 1.5.3 
- py2neo 2021.2.3

## How to use
```tiddlywiki
Enrichment.py -p [PASSWORD] -q [QUERY] -t [NODE_TYPE] -s [TAXON_ID] [OPTION]
-p    --password    [required] Password to connect to the DBMS.
-q    --query    [required] Query set.
-t    --type    [required] Target sets node type.
-s    --spieces    [required] Taxon id (int).
-a    --alpha    [optional] Significance threshold (float).
-c    --adjust    [optional] Adjust for multiple testing (FDR).
-m    --metric    [optional] Dissimilarity method: binomial (default), hypergeometric, chi2, coverage.
-l    --limit    [optional] Maximum number of results to report (int).
-e    --eval    [optional] Metrics evaluation.
-v    --verbose    [optional] Talk a lot.
```

## Example of use
```typescript
./Enrichment.py -p password -q query_sets/set01.txt -t GOTerm -s 511145 alpha 0.2 -c -m chi2 -l 10
```
This will print the first 10 enriched GOTerm in set01.txt, using chi 2 of independance, with a significance threshold of 0.2 and using adjust for mutliple testing

```typescript
./recherche_enrichissement.py -p password -q query_sets/set01.txt -t Keyword -s 511145 -e -l 5
```
This will print the first 5 enriched Keyword in set01.txt with the 4 methods.