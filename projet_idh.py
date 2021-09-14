#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2neo import Graph
import argparse

################################################################################################
parser = argparse.ArgumentParser(description='Tracks recommendation')
parser.add_argument('--password', '-p', type=str, required=True, help="Password to connect to the database")
parser.add_argument('--username', '-u', type=str, required=True, help="Username, (Case sensitive)")
parser.add_argument('--alpha', '-a', type=float, required=False, default=0.25, help="Alpha value (0.25 by default)")
parser.add_argument('--beta', '-b',  type=float, required=False, default=0.75, help="Beta value (0.75 by default)")
args = parser.parse_args()
################################################################################################

# Connection to Neo4J
try:
    graph = Graph("bolt://localhost:7687", auth=("neo4j", args.password))
except:
    print("Connection to the database failed!\nPlease check your password and try again.\n--help for more information")
    exit()

# Musical genre match
def fb(u,t):
    # Get genres liked by user u
    listGenre = []
    q =(f"MATCH (:Person{{name:'{u}'}})-[]->(:Track)-[]->(:Album)-[]->(:Artist)-[]->(g:Genre) "
         "RETURN distinct g.name")
    res = graph.run(q).to_table()
    for el in res : 
        listGenre.append(el[0])
    # Get track t genre
    q =(f"MATCH (t:Track)-[]->(:Album)-[]->(:Artist)-[]->(g:Genre) "
        f"WHERE id(t) = {t} "
        f"RETURN distinct g.name")
    res = graph.run(q).to_table()
    # Check musical genre match
    for el in res:
        if el[0] in listGenre:
            return 2
    return 1

# Social factor
def fs(u,t):
    # Get shortest path length and number 
    q =(f"MATCH path=allShortestPaths((u:Person{{name:'{u}'}})-[*]->(t:Track)) "
        f"WHERE id(t) = {t} "
        f"RETURN length(path), count(*)")
    res = graph.run(q).to_table()
    return res[0][0] * res[0][1]

# Thematics coherence
def ft(u,t):
    #Get tracks liked by user u
    tracksliked = []
    q=(f"MATCH (p:Person{{name:'{u}'}})-[:LIKE]->(t:Track)"
       f"RETURN id(t)")
    res = graph.run(q).to_table()
    for el in res : 
        tracksliked.append(el[0])
    #computes jaccard coefficient for each track pair
    jaccard = []
    for track in tracksliked:
        q=(f"MATCH (t1:Track)<-[:WORD_OF]-(w:Word)-[:WORD_OF]->(t2:Track) "
           f"WHERE id(t1) = {track} AND id(t2) = {t} "
           f"RETURN COUNT (distinct w.name)")
        res = graph.run(q).to_table()
        intercept = res[0][0] 
        q=(f"MATCH (w1:Word)-[:WORD_OF]->(t1:Track)"
           f"WHERE id(t1) = {track} RETURN distinct lower(w1.name) AS word " # Case insensitive
           f"UNION "
           f"MATCH (w2:Word)-[:WORD_OF]->(t2:Track)"
           f"WHERE id(t2) = {t} RETURN distinct lower(w2.name) AS word"
        )
        res = graph.run(q).to_table()
        union = len(res)
        coeff = intercept/union
        jaccard.append(coeff)
    #Get max jaccard coefficient
    return max(jaccard)

# Track score 
def score(u, t, a, b):
    return fb(u,t) * (a*fs(u,t) + b*ft(u,t))

if __name__ == "__main__":
    # Checks if user in database
    q = f"MATCH (p:Person{{name:'{args.username}'}}) RETURN count(*)"
    res = graph.run(q).to_table()
    if res[0][0] == 0:
        print("Username not in the database!\nPlease check the username and try again.\n--help for more information")
        exit()
    
    # Tracks retrieval (already liked tracks are excluded)
    q =(f"MATCH (:Person{{name:'{args.username}'}})-[]->(t:Track) "
        f"WITH collect(distinct id(t)) AS likedTracks "
        f"MATCH (t2:Track) "
        f"WHERE  NOT id(t2) IN likedTracks "
        f"RETURN distinct id(t2)")
    res = graph.run(q).to_table()
    tracks = []
    cpt = 0
    for el in res:
        cpt += 1
        print(f"In progress : {cpt}/{len(res)} track's scores computed", end = "\r") 
        tracks.append([el[0],score(args.username, el[0], args.alpha, args.beta)])
    tracks = sorted(tracks, key=lambda x:x[1], reverse = True)
    print("\nTracks recommendation : ")
    for n in range(20):
        q = f"MATCH (t:Track)-[]->(:Album)-[]->(a:Artist) WHERE id(t) = {tracks[n][0]} RETURN t.name, a.name"
        res = graph.run(q).to_table()
        print(f"#{n+1} : " + res[0][0]+ " by "+ res[0][1] )
        print("Score = " + str(round(tracks[n][1], 2)))


