#!/usr/bin/env python

import argparse
from os.path import isfile
from scipy.stats import binom, hypergeom, chi2_contingency
from py2neo import *
from py2neo.matching import *
import numpy as np

# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets node type.')
parser.add_argument('-s', '--species', required=True, type=int, help='Taxon id.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2, coverage or all.')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
arg = parser.parse_args()

# LOAD QUERY
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

# CONNECT TO Neo4J
neo = Graph("bolt://localhost:7687", auth=("neo4j", "a"))
nodes = NodeMatcher(neo)

# COMPUTE POPULATION SIZE
population_size = nodes.match('Gene', taxon_id=arg.species).count()
if arg.verbose:
	print("population_size: ", population_size)

# RETRIEVE TARGET SETS FOR QUERY
path = '[:is_a|part_of|annotates*]' if arg.sets=='GOTerm' else ''
cypher = f"MATCH (t:{arg.sets})-{path}->(n:Gene {{taxon_id:{arg.species} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN COUNT(DISTINCT t)"
sample_size = neo.run(cypher).to_table()[0][0] # Get sample size for feedback purpose
cypher = f"MATCH (t:{arg.sets})-{path}->(n:Gene {{taxon_id:{arg.species} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN DISTINCT t"
sets = neo.run(cypher)
if arg.verbose:
	print(cypher)
# EVALUATE SETS
cpt = 0
results = []
for s in sets: #_ids:
	cpt += 1
	print(f"In progress : {cpt}/{sample_size}", end = "\r") # Feedback
	sid = s['t']['id']
	# RETRIEVE TARGET SET ELEMENTS
	path = '[:is_a|part_of|annotates*1..20]' if arg.sets=='GOTerm' else ''
	cypher = f"MATCH (t:{arg.sets})-{path}->(n:Gene {{taxon_id:{arg.species} }}) WHERE t.id='{sid}' RETURN n.id"
	if arg.verbose:
		print(cypher)
	table = neo.run(cypher).to_table()
	elements = set( map(lambda x: x[0], table))
	common_elements = elements.intersection( query ) # C
	c = len(common_elements)
	q = len(query)
	t = len(elements)
	g = population_size
	if c < 2:
		next
	if arg.measure=='binomial': # binom.cdf(>=success, attempts, proba)
		measure = binom.cdf(q - c, q, 1 - float(t) / g)
	elif arg.measure == 'coverage':
		measure = (c/q) * (c/t)
	elif arg.measure =='hypergeometric':
		measure = hypergeom.cdf(q - t, g, int((1 - float(t)/g)*g), q)
	elif arg.measure =='chi2':
		contengency_table = np.array([[c, q-c],[t-c, g-q-t+c]])
		g, measure, dof, expctd = chi2_contingency(contengency_table)
	elif arg.measure =='all':
		break
	else:
		print(f'sorry, {arg.measure} not (yet) implemented')
		exit(1)
	r = { 'id': sid, 'desc': s['t']['desc'], 'common.n':c, 'target.n': t, 'measure': measure, 'elements.target': elements, 'elements.common': common_elements }
	results.append( r )

if arg.verbose:
  print(results)

# PRINT SIGNIFICANT RESULTS
i = 0
results.sort(key=lambda an_item: an_item['measure'])
print("Significant results :")
for r in results:
	i += 1
	# FDR
	if arg.adjust and r['measure'] > arg.alpha * i / len(results): break
	# limited output
	if arg.limit > 0 and i>arg.limit: break
	# alpha threshold
	elif r['measure'] > arg.alpha: break
	# OUTPUT
	pval = "{:.4f}".format(r['measure']) if r['measure']>0.01 else "{:.2e}".format(r['measure'])
	print(f"{r['id']}\t{pval}\t{r['common.n']}/{r['target.n']}\t{r['desc']}")
	if arg.verbose : print(f"Common elements : {r['elements.common']}")
		

