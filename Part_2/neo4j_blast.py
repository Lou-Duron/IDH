#!/usr/bin/env python

import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom
from py2neo import Graph

# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--target', required=True, help='PubMed, Keyword or GO')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial',
                    help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
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
else:  # parse string
    query |= set(text.split())

if param.verbose:
    print(f'query set: {query}')

# LOAD REFERENCE SETS
if param.verbose:
    print("Conncetion to graph bolt://localhost:7687...")
graph = Graph("bolt://localhost:7687", auth=("neo4j", "a"))

if param.verbose:
    print("Connected")
q = "match (p:Gene) return count(p)"

population_size = graph.run(q).to_table()[0][0]

if param.verbose:
    print("Population size :", population_size)


# EVALUATE SETS


def get_Keyword_sets():
    if param.verbose:
        print("Target Type KeyWord Processing...")
    q2 = "match (k:Keyword) return k.keyword"

    keywords = graph.run(q2).to_table()
    if param.verbose:
        print("Nombre d'eléments : ", len(keywords))

    sets = {}
    for k in keywords:
        sets[k[0]] = []
        # print(k[0])
        q3 = "match (g:Gene)-[]-(:Keyword {keyword : \"" + k[0] + "\"}) return g.gene_id"
        for g in graph.run(q3).to_table():
            sets[k[0]].append(g[0])
    return sets


def get_Keyword_sets_v2(query):
    if param.verbose:
        print("Target Type KeyWord Processing...")
    # q2 = "match (k:Keyword) return count(k.keyword)"

    # keywords = graph.run(q2).to_table()
    q = "match (g:Gene)-[]-(k:Keyword) where g.gene_id IN [" + "','".join(query) + "] return distinct(k.keyword)"
    common = graph.run(q).to_table()
    if param.verbose:
        print("Nombre d'eléments : ", len(common))
    for k in common:
        sets[k[0]] = []
        q2 = "match (g:Gene)-[]-(:Keyword {Keyword: \"" + k[0] + "\"}) return g.gene_id"
        for g in graph.run(q2).to_table():
            sets[k[0]].append(g[0])
    return sets


def get_GO_sets():
    if param.verbose:
        print("Target Type GOTerm Processing...")
    q = "match (g:Gene)<-[:is_a|part_of|annotates*]-(t:GOTerm) where g.gene_id IN ['" + "','".join(
        query) + "'] return distinct(t.id)"
    common = graph.run(q).to_table()
    if param.verbose:
        print("Nombre d'eléments : ", len(common))
    sets = {}
    for t in common:
        sets[t[0]] = []
        q2 = "match (g:Gene)<-[:is_a|part_of|annotates*]-(t:GOTerm {id: \"" + t[0] + "\"}) return g.gene_id"
        for g in graph.run(q2).to_table():
            sets[t[0]].append(g[0])
    return sets


def get_Pubmed_sets():
    if param.verbose:
        print("Target Type PubMed Processing...")
    req = "match (g:Gene)-[]-(p:PubMed) where g.gene_id IN ['" + "','".join(query) + "'] return distinct(p.PMID)"
    articles = graph.run(req).to_table()
    if param.verbose:
        print("Nombre d'eléments : ", len(articles))
    sets = {}
    for ar in articles:
        sets[ar[0]] = []
        q3 = "match (g:Gene)-[]-(:PubMed {PMID: \"" + ar[0] + "\"}) return g.gene_id"
        for g in graph.run(q3).to_table():
            sets[ar[0]].append(g[0])
    return sets


def compute_pval(sets, query):
    results = []
    query_size = len(query)
    if param.verbose:
        print(query_size)
        print(population_size)
    for s in sets:
        elements = set(sets[s])
        common_elements = elements.intersection(query)
        if param.measure == 'binomial':  # binom.cdf(>=success, attempts, proba)
            count_common = len(common_elements)
            pvalue = binom.cdf(query_size - count_common, query_size, 1 - len(elements) / population_size)
        else:
            print(f'sorry, {param.measure} not (yet) implemented')
            exit(1)
        r = {'id': s, 'common.n': len(common_elements), 'target.n': len(elements),
             'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements}
        results.append(r)
    if param.verbose:
        print(results)

    # PRINT SIGNIFICANT RESULTS
    results.sort(key=lambda an_item: an_item['p-value'])
    k = 1
    m = len(results)
    for r in results:
        if r['p-value'] > (k / m) * param.alpha:
            break
        # OUTPUT
        k += 1
        print("{}\t{}\t{}/{}\t{}".format(r['id'], r['p-value'], r['common.n'], r['target.n'], ', '.join(r['elements'
                                                                                                          '.common'])))


sets = None

if param.target == "Keyword":
    sets = get_Keyword_sets()
if param.target == "PubMed":
    sets = get_Pubmed_sets()
if param.target == "GOTerm":
    sets = get_GO_sets()

if param.verbose:
    print("\n\n\n\n************************ RESULT ******************\n\n\n\n\n")
compute_pval(sets, query)
