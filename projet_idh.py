#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2neo import Graph
import argparse

################################################################################################
parser = argparse.ArgumentParser(description='truc')
parser.add_argument('--username', '-u', type=str, required=True, help="Name of the user")
parser.add_argument('--track', '-t', type=str, required=True, help="Name of the track")
parser.add_argument('--alpha', '-a', type=float, required=False, default=0.25, help="Alpha value (0.25 by default")
parser.add_argument('--beta', '-b',  type=float, required=False, default=0.75, help="Beta value (0.75 by default")
args = parser.parse_args()
################################################################################################

#Connection to Neo4J
graph = Graph("bolt://localhost:7687", auth=("neo4j", "a"))

# Adéquation au genre musical
def fb(u,t):
    listGenre = []
    q = f"match (:Person{{name:'{u}'}})-[]->(:Track )-[]->(:Album)-[]->(:Artist)-[]->(g:Genre) return distinct g.name"
    res = graph.run(q).to_table()
    for el in res : 
        listGenre.append(el[0])
    q = f"match (:Track{{name:'{t}'}})-[]->(:Album)-[]->(:Artist)-[]->(g:Genre) return distinct g.name"
    res = graph.run(q).to_table()
    for el in res:
        if el[0] in listGenre:
            return 2
    return 1

# Facteur social
def fs(u,t):
    q = f"match path=allShortestPaths((u:Person{{name:'{u}'}})-[*]->(t:Track{{name:'{t}'}})) RETURN length(path), count(*)"
    res = graph.run(q).to_table()
    print(res)
    return res[0][0] * res[0][1]
    """
    Nombre de PPC * la taille des PPC
    """

def ft(u,t):
    """
    Cohérence thématique (Jaccard mes couilles)
    """

def score(u, t,a=0.25, b=0.75):
    return fb(u,t) * (a * fs(u,t) + b * ft(u,t))



print(fs(test, test2))

