#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2neo import Graph
import argparse

################################################################################################
parser = argparse.ArgumentParser(description='truc')
parser.add_argument('--password', '-p', type=str, required=True, help="Password to connect to the database")
parser.add_argument('--username', '-u', type=str, required=True, help="Name of the user")
parser.add_argument('--alpha', '-a', type=float, required=False, default=0.25, help="Alpha value (0.25 by default")
parser.add_argument('--beta', '-b',  type=float, required=False, default=0.75, help="Beta value (0.75 by default")
args = parser.parse_args()
################################################################################################

#Connection to Neo4J
graph = Graph("bolt://localhost:7687", auth=("neo4j", args.password))

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
    return res[0][0] * res[0][1]

def ft(u,t):
    """
    Cohérence thématique
    """
    #Get tracks liked by user u
    tracksliked = []
    q=(f"MATCH (p:Person{{name:'{u}'}})-[:LIKE]->(t:Track)"
       f"RETURN t.name")
    res = graph.run(q).to_table()
    for el in res : 
        tracksliked.append(el[0])
    #computes jaccard coefficient for each track pair
    jaccard = []
    for track in tracksliked:
        q=(f"MATCH (t1:Track{{name:'{track}'}})<-[:WORD_OF]-(w:Word)-[:WORD_OF]->(t2:Track{{name:'{t}'}})"
           f"RETURN COUNT (distinct w.name)")
        res = graph.run(q).to_table()
        intercept = res[0][0] 
        union = []
        q=(f"MATCH (p:Person{{name:'{u}'}})-[:LIKE]->(t1:Track{{name:'{track}'}})<-[:WORD_OF]-(w1:Word)"
           f"RETURN distinct w1.name AS word "
           f"UNION "
           f"MATCH (w2:Word)-[:WORD_OF]->(t2:Track{{name:'{t}'}})"
           f"RETURN distinct w2.name AS word"
        )
        res = graph.run(q).to_table()
        union = len(res)
        coeff = intercept/union
        jaccard.append(coeff)
    #Get max jaccard coefficient
    return max(jaccard)

def score(u, t, a, b):
    return fb(u,t) * (a * fs(u,t) + b * ft(u,t))

if __name__ == "__main__":
    q = f"match (:Person{{name:'{args.username}'}})-[]->(t:Track) with collect(distinct t) as test match (t2:Track) where  NOT t2 IN test  return distinct t2.name"
    res = graph.run(q).to_table()
    tracks = []
    cpt = 0
    for el in res:
        cpt += 1
        print(f"In progress : {cpt}/{len(res)} tracks", end = "\r") 
        tracks.append([el[0],score(args.username, el[0], args.alpha, args.beta)])
        
    tracks = sorted(tracks, key=lambda x:x[1], reverse = True)
    print(tracks[0:20])
   



