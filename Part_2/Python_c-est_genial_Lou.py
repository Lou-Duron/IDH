#!/usr/bin/env python

import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom

# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets filename.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
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

# LOAD REFERENCE SETS
sets = json.loads(open(param.sets).read())
if param.verbose:
    print('first target sets: ', sets[0:2])

# COMPUTE POPULATION SIZE
population = set()
pop = []
for dict in sets:
  for id in dict['elements']:
    population.add(id)
popSize = len(population)
if param.verbose:
  print("Population size : ", popSize)


# EVALUATE SETS
results = []
query_size = len(query)
for s in sets:
    elements = set(s['elements' ])
    common_elements = elements.intersection(query)
    if param.measure=='binomial': # binom.cdf(>=success, attempts, proba)
        pvalue = binom.cdf( query_size - len(common_elements), query_size, 1 - float(len(elements))/popSize)
    else:       
        print(f'sorry, {param.measure} not (yet) implemented')
        exit(1)
    r = { 'id': s['id'], 'desc': s['desc'], 'common.n':len(common_elements), 'target.n': len(elements), 'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements }
    results.append( r )
if param.verbose:
  print(results)


# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['p-value'])

if param.adjust:
    for i, r in enumerate(results):
        if r['p-value'] > ((i+1) / len(results)) * param.alpha:
            break
        print("{}\t{}\t{}/{}\t{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))
else:
    for r in results:
        if r['p-value'] > param.alpha: 
            break
            # OUTPUT
        print("{}\t{}\t{}/{}\t{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))



