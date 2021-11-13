#!/usr/bin/env python

import argparse
from os.path import isfile
from scipy.stats import binom, hypergeom, chi2_contingency
from py2neo import *
from py2neo.matching import *
import numpy as np
import pprint
import time

start_time = time.time() # Execution time

# Parameters
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-p', '--password', required=True, help='Password to connect to the DBMS')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--type', required=True, help='Target sets node type.')
parser.add_argument('-s', '--species', required=True, type=int, help='Taxon id.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--metric', required=False, default='binomial', help='Dissimilarity method: binomial (default), hypergeometric, chi2, coverage.')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-e', '--eval', required=False, action="store_true", help='Metrics evaluation')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
arg = parser.parse_args()

# DBMS connection
try:
	neo = Graph("bolt://localhost:7687", auth=("neo4j", arg.password))
	nodes = NodeMatcher(neo)
	print('Connection to the DBMS successful')
except Exception:
	print('Error : An error occurred while trying to connect to the DBMS please check your password')
	exit(1)

# Parameters checks
# Type validity check
if arg.type not in ['Alias', 'GOTerm', 'Interpro', 'Keyword', 'Pathway', 'PubMed', 'TU'] :
	print("Error : Please choose a valid node type : Alias, GOTerm, Interpro, Keyword, Pathway, PubMed, TU")
	exit(1)
# Taxon id check
if neo.run(f"MATCH (n:Gene {{taxon_id:{arg.species}}}) RETURN COUNT(n)").to_table()[0][0] == 0:
	print("Error : Taxon id not in the Database")
	exit(1)

# Query loading
text = arg.query
query = set()
if isfile(text):
	with open(text) as f:
		content = ' '.join(f.read().split('\n')).split()
		query |= set(content)
else: # parse string
	query |= set(text.split())

if arg.verbose:
  print(f'query set: {query}')

# Population size recovery
population_size = nodes.match('Gene', taxon_id=arg.species).count()
if arg.verbose:
	print("Population_size: ", population_size)

# Retrieve target sets for query

path = '[:is_a|part_of|annotates*]' if arg.type=='GOTerm' else ''
cypher = f"MATCH (t:{arg.type})-{path}->(n:Gene {{taxon_id:{arg.species} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN COUNT(DISTINCT t)"
sample_size = neo.run(cypher).to_table()[0][0] # Get sample size for feedback purpose
if sample_size == 0:
	print("Error : It seems that the query is not valid, please check the documentation for more informations")
	exit(1)
cypher = f"MATCH (t:{arg.type})-{path}->(n:Gene {{taxon_id:{arg.species} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN DISTINCT t"
sets = neo.run(cypher)
if arg.verbose:
	print(cypher)

# Results recovery
if arg.eval: # Stores the result from all methods for metrics evaluation
	results = {'Binomial':[],'Coverage':[],'Hypergeometric':[],'Chi2':[]}
else:
	results = []

# Metrics
def binomial(c, q, t, g):
	return binom.cdf(q - c, q, 1 - float(t) / g)

def coverage(c, q, t):
	return (c/q) * (c/t)

def hypergeometric(q, t, g):
	return hypergeom.cdf(q - t, g, int((1 - float(t)/g)*g), q)

def chi2(c, q, t, g):
	return chi2_contingency(np.array([[c, q-c],[t-c, g-q-t+c]]))[1]

# Enrichment evaluation
for count, s in enumerate(sets):
	print(f"In progress : {count}/{sample_size} ({round(count/sample_size*100,1)}%)", end = "\r") # Feedback
	sid = s['t']['id']
	path = '[:is_a|part_of|annotates*1..20]' if arg.type=='GOTerm' else '' # Get correct path
	cypher = f"MATCH (t:{arg.type})-{path}->(n:Gene {{taxon_id:{arg.species} }}) WHERE t.id='{sid}' RETURN n.id"
	if arg.verbose:
		print(cypher)
	table = neo.run(cypher).to_table()
	elements = set( map(lambda x: x[0], table))
	common_elements = elements.intersection(query)
	c = len(common_elements)
	q = len(query)
	t = len(elements)
	g = population_size
	if c < 2:
		next
	if not arg.eval: # If we are not evaluating the metrics
		if arg.metric=='binomial':
			measure = binomial(c, q, t, g)
		elif arg.metric == 'coverage':
			measure = coverage(c, q, t)
		elif arg.metric =='hypergeometric':
			measure = hypergeometric(q, t, g)
		elif arg.metric =='chi2':
			measure = chi2(c, q, t, g)
		else:
			print(f'Sorry, {arg.metric} not implemented')
			print("Error : Please choose an implemented method : binomial, coverage, hypergeometric, chi2")
			exit(1)
		r = { 'id': sid, 'desc': s['t']['desc'], 'common.n':c, 'target.n': t, 'measure': measure, 'elements.target': elements, 'elements.common': common_elements }
		results.append(r)

	else: # Metrics evaluation
		temp_results = [] # Stores results from all methods
		temp_results.append(binomial(c, q, t, g)) # Binomial
		temp_results.append(coverage(c, q, t)) # Coverage
		temp_results.append(hypergeometric(q, t, g)) # Hypergeometric
		temp_results.append(chi2(c, q, t, g)) # Chi2
		for count, key in enumerate(results.keys()):
			r = { 'id': sid, 'desc': s['t']['desc'], 'common.n':c, 'target.n': t, 'measure': temp_results[count], 'elements.target': elements, 'elements.common': common_elements }
			results[key].append(r)
	
if arg.verbose:
	pprint.pprint(results, width=1)

# Print significant results
print(f"Finished in {round(time.time() - start_time, 3)} seconds !   ")
print("Significant results :")
if not arg.eval:
	results.sort(key=lambda an_item: an_item['measure'])
	print("ID   p_value   ratio   Description   Common_elements") #################
	for count, r in enumerate(results):
		if arg.adjust and r['measure'] > arg.alpha * count+1 / len(results): break # FDR
		if arg.limit > 0 and count+1>arg.limit: break # Limited output
		elif r['measure'] > arg.alpha: break # Alpha threshold
		# OUTPUT
		pval = "{:.4f}".format(r['measure']) if r['measure']>0.01 else "{:.2e}".format(r['measure'])
		print(f"{r['id']}\t{pval}\t{r['common.n']}/{r['target.n']}\t{r['desc']}", end=' ')
		if arg.verbose : 
			print(r['elements.common'])
		else :
			print(len(r['elements.common']))
else:
	for key in results.keys():
		results[key].sort(key=lambda an_item: an_item['measure'])
		print(f"\n{key} :")
		print("ID   p_value   ratio   Description   Common_elements") #################
		for count, r in enumerate(results[key]):
			if arg.adjust and r['measure'] > arg.alpha * count+1 / len(results): break # FDR
			if arg.limit > 0 and count+1>arg.limit: break # Limited output
			elif r['measure'] > arg.alpha: break # Alpha threshold
			# OUTPUT
			pval = "{:.4f}".format(r['measure']) if r['measure']>0.01 else "{:.2e}".format(r['measure'])
			print(f"{r['id']}\t{pval}\t{r['common.n']}/{r['target.n']}\t{r['desc']}", end=' ')
			if arg.verbose : 
				print(r['elements.common'])
			else :
				print(len(r['elements.common']))
