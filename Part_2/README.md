# Welcome to Enrichment analysis project

This Gitlab summarizes the work done during IDH (Heterogenous Data Integration) class as part of the Master Sc in Bioinformatics of Toulouse.

**Authors :**
- **Victoria Fathi**
- **Lou Duron**

## Organisation

- **benchmark_data** : all data related to the benchmark
- **raw_data** : all data related to E.coli retrieved form UniProt and EcoCyc 
- **import** : all files needed to create the data base
- **query_sets** : query sets for enrichment analysis
- **Neo4j_Roland.py** : initial script (this has been modified in each branch)
- **data_integration_procedure.hmtl** : procedure used to integrate E.coli data in Neo4j database.
- **get_data.Rmd** : procedure used to prepare data for integration
- **get_benchmark_data.r** : procedure used to prepare benchmark data for analysis
- **benchmark_analysis.Rmd** : procedure used for benchmark analysis


## Description

This project is divided in two parts.
- Data intregration
- Enrichment analysis

The first part has been done in collaboration as part of IDH class. The second part had to be done individualy, therefore this GitLab project has been divided in two branches (lou and vic).

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

