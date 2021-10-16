#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2neo import Graph
import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom

################################################################################################
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets node type.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-s', '--taxon_id', required=True, type=int, help='Taxon id.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
param = parser.parse_args()
################################################################################################

# Connection to Neo4J
try:
    graph = Graph("bolt://localhost:7687", auth=("neo4j", "a"))
except:
    print("Connection to the database failed!\nPlease check your password and try again.\n--help for more information")
    exit()

# LOAD QUERY
text = param.query
query = set()
if isfile(text):
    with open(text) as f:
        content = ' '.join(f.read().split('\n')).split()
        query |= set(content)
else: # parse string
    query |= set(text.split())


# COMPUTE POPULATION SIZE
q =(f"MATCH (n:Gene) return count(n)")
population_size = graph.run(q).to_table()[0][0]
print(f"population size is {population_size}")

# Get keyword associated to the genes query 
path = '[:is_a|part_of|annotates*]' if param.sets=='GOTerm' else ''
cypher = f"MATCH (t:{param.sets})-{path}->(n:Gene {{taxon_id:{param.taxon_id} }}) WHERE n.id IN ['"+"', '".join(query)+"'] RETURN distinct t"
sets = graph.run(cypher)

if param.verbose:
    print(cypher)
    print(sets)

#EVALUATE SETS
results = []
query_size = len(query)

for s in sets:
    sid = s['t']['id']
    print(sid)
	# RETRIEVE TARGET SET ELEMENTS
    path = '[:is_a|part_of|annotates*]' if param.sets=='GOTerm' else ''
    cypher = "MATCH (t:{})-{}->(n:Gene {{taxon_id:{} }}) WHERE t.id='{}' RETURN distinct n.id".format(param.sets, path, param.taxon_id, sid)
    table = graph.run(cypher).to_table()
    elements = set( map(lambda x: x[0], table))
	#print("elements:", elements)
    common_elements = elements.intersection(query)
    if param.measure=='binomial': # binom.cdf(>=success, attempts, proba)
    	pvalue = binom.cdf( query_size - len(common_elements), query_size, 1 - float(len(elements))/population_size)
    else:
    	print(f'sorry, {param.measure} not (yet) implemented')
    	exit(1)
    r = { 'id': sid, 'desc': s['t']['desc'], 'common.n':len(common_elements), 'target.n': len(elements), 'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements }
    results.append( r )


# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['p-value'])

for i, r in enumerate(results):
    if param.adjust:
        if r['p-value'] > ((i+1)/len(results))*param.alpha:
            break
    else:
        if r['p-value'] > param.alpha: 
            break
    # OUTPUT
    print("{}\t{}\t{}/{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], ', '.join(r['elements.common'])))