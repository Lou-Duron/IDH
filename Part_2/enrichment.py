#!/usr/bin/env python

import argparse
from os.path import isfile
from scipy.stats import binom, hypergeom, chi2_contingency
from pprint import pprint
from py2neo import *
from py2neo.matching import *
import numpy as np

# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--target_type', required=True, help='Target sets node type [GOTerm, Interpro, Keyword, Pathway, PubMed, TU"].')
parser.add_argument('-s', '--species', required=True, type=int, help='Taxon id.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
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
neo = Graph("bolt://localhost:7687", auth=("neo4j", "a"))
nodes = NodeMatcher(neo)

# COMPUTE POPULATION SIZE
population_size = nodes.match('Gene', taxon_id=param.species).count()
if param.verbose:
	print("population_size: ", population_size)

# RETRIEVE TARGET SETS FOR QUERY
path = '[:is_a|part_of|annotates*]' if param.target_type =='GOTerm' else ''
cypher = f"MATCH (t:{param.target_type})-{path}->(n:Gene {{taxon_id:{param.species} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN DISTINCT t.id"
if param.verbose:
	print(cypher)
target_sets = set(map(lambda x: x[0], neo.run(cypher).to_table()))
target_sets_size = len(target_sets)

# EVALUATE SETS
results = []
query_size = len(query)

for i, target in enumerate(target_sets): #_ids:
	# PROGRESSION STATUS
	print(f"In progress : {i}/{target_sets_size} ({round(i/target_sets_size*100,1)}%)", end = "\r")

	# RETRIEVE TARGET SET ELEMENTS
	path = '[:is_a|part_of|annotates*1..20]' if param.target_type=='GOTerm' else ''
	cypher = "MATCH (t:{})-{}->(n:Gene {{taxon_id:{} }}) WHERE t.id='{}' RETURN n.id".format(param.target_type, path, param.species, target)
	
	if param.verbose:
		print(cypher)
	
	table = neo.run(cypher).to_table()
	target_elements = set(map(lambda x: x[0], table))
	common_elements = target_elements.intersection(query)

	common_elements_size = len(common_elements)
	target_elements_size = len(target_elements)

	if len(common_elements) < 2:
		next
	if param.measure == 'coverage':
		measure = (common_elements_size/query_size * (common_elements_size/target_elements_size))
	if param.measure =='binomial': # binom.cdf(>=success, attempts, proba)
		measure = binom.cdf(query_size - len(common_elements), query_size, 1 - float(target_elements_size)/population_size)
	if param.measure =='hypergeometric':
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
	r = { 'id': target, 'common_size':len(common_elements), 'target_size': len(target_elements), 'measure': measure}
	results.append(r)

if param.verbose:
  print(results)

# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['measure'])

if param.measure != "coverage":
	for r in results:
		# FDR
		if param.adjust and r['measure'] > param.alpha * i / len(results): break
		
		# limited output
		if param.limit > 0 and i+1>param.limit: break
		# alpha threshold
		
		elif r['measure'] > param.alpha: break
		# OUTPUT
		pval = "{:.4f}".format(r['measure']) if r['measure']>0.01 else "{:.2e}".format(r['measure'])
		print(f"{r['id']}\t{pval}\t{r['common_size']}\t{r['target_size']}")

