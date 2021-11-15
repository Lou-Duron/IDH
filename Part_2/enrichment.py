#!/usr/bin/env python

import argparse
from os.path import isfile
from scipy.stats import binom, hypergeom, chi2_contingency
from py2neo import *
from py2neo.matching import *
import numpy as np
import pandas as pd
from termcolor import colored


# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-p', '--password', required=True, help='password neo4j database')
parser.add_argument('-q', '--query', required=True, help='path to query set')
parser.add_argument('-t', '--target_type', required=True, choices=["GOTerm", "Interpro", "Keyword", "Pathway", "PubMed", "TU"], help='Target sets node type [GOTerm, Interpro, Keyword, Pathway, PubMed, TU"]')
parser.add_argument('-s', '--species', required=False, type=int, help='Taxon id', default=511145)
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR)')
parser.add_argument('-m', '--measure', required=False, choices = ["binomial", "hypergeometric", "coverage", "chi2"], default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-r', '--revigo', required=False, action="store_true", help='only if --write, write id and pvalue for revigo')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='print intermediary results/queries')
parser.add_argument('-w', '--write', required=False, help='path and name for results (tsv file)')
param = parser.parse_args()


# LOAD QUERY
text = param.query
query = set()
if isfile(text):
	with open(text) as f:
		content = ' '.join(f.read().split('\n')).split()
		query |= set(content)
else: # parse string
	query |= set(text.split())

if param.verbose:
  print(f'query set: {query}')

# CONNECT TO Neo4J
neo = Graph("bolt://localhost:7687", auth=("neo4j", param.password))
nodes = NodeMatcher(neo)

# COMPUTE POPULATION SIZE
population_size = nodes.match('Gene', taxon_id=param.species).count()
if param.verbose:
	print("population_size: ", population_size)

# RETRIEVE TARGET SETS FOR QUERY
path = '[:is_a|part_of|annotates*]' if param.target_type =='GOTerm' else ''
cypher = f"MATCH (t:{param.target_type})-{path}->(n:Gene {{taxon_id:{param.species} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN DISTINCT t"
if param.verbose:
	print(cypher)

target_sets = neo.run(cypher).to_table()
target_sets_size = len(target_sets)
if target_sets_size == 0:
	print(f"{cypher} returned an empty results.\nCheck if genes are in the database or in the correct format (bnumber)")
	exit(1)
	
# EVALUATE SETS
results = []
query_size = len(query)

for i, target in enumerate(target_sets): #_ids:
	# PROGRESSION STATUS
	print(f"In progress : {i}/{target_sets_size} ({round(i/target_sets_size*100,1)}%)", end = "\r")

	# RETRIEVE TARGET SET ELEMENTS
	target_id = target[0]['id']
	path = '[:is_a|part_of|annotates*1..20]' if param.target_type=='GOTerm' else ''
	cypher = "MATCH (t:{})-{}->(n:Gene {{taxon_id:{} }}) WHERE t.id='{}' RETURN n.id".format(param.target_type, path, param.species, target_id)
	
	if param.verbose:
		print(cypher)
	
	table = neo.run(cypher).to_table()
	target_elements = set(map(lambda x: x[0], table))
	common_elements = target_elements.intersection(query)

	common_elements_size = len(common_elements)
	target_elements_size = len(target_elements)

	# ENRICHMENT MEASURES
	if len(common_elements) < 2:
		next
	if param.measure == 'coverage':
		# based of shared elements
		measure = (common_elements_size/query_size * (common_elements_size/target_elements_size))

	if param.measure =='binomial': # binom.cdf(>=success, attempts, proba)
		measure = binom.cdf(query_size - len(common_elements), query_size, 1 - float(target_elements_size)/population_size)

	if param.measure =='hypergeometric': #hypergeometric.cdf(>=sucess, total number of object ,proba, attempts)
		inverse_freq = int((1 - float(target_elements_size)/population_size)*population_size)
		measure = hypergeom.cdf(query_size - common_elements_size, population_size, inverse_freq, query_size)

	if param.measure =='chi2': 
		unique_query = query_size - common_elements_size
		unique_target = target_elements_size - common_elements_size
		residue =  population_size - query_size - target_elements_size + common_elements_size
		contengency_table = np.array([
			[common_elements_size, unique_query],
			[unique_target, residue]
			])
		measure = chi2_contingency(contengency_table)[1]
	
	# RESULTS
	r = { 'id': target_id, 'desc': target[0]["desc"], 'common_size':len(common_elements), 'target_size': len(target_elements), 'measure': measure}
	results.append(r)

if param.verbose:
  print(results)

# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['measure'])
adjusted_results = []
for r in results:
	# FDR
	if param.adjust and r['measure'] > param.alpha * i / len(results) and param.measure != 'coverage': break
	
	# limited output
	if param.limit > 0 and i+1>param.limit: break
	# alpha threshold
	
	elif r['measure'] > param.alpha: break
	# OUTPUT
	pval = "{:.4f}".format(r['measure']) if r['measure']>0.01 else "{:.2e}".format(r['measure'])
	adjusted_results.append({"id": r['id'], "desc": r['desc'], "pval":pval, "common_size":r['common_size'], "target_size":r['target_size']})

print("")
df_adjusted_results = pd.DataFrame(adjusted_results)	
print(df_adjusted_results)

if param.write:
	if param.revigo:
			df_adjusted_results[['id', 'pval']].to_csv(param.write, sep = "\t", index= False)
	else:
		df_adjusted_results.to_csv(param.write, sep = "\t", index= False)
	print(colored(f"file written: {param.write}", 'green')) 
