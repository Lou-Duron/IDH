from py2neo import Graph

graph = Graph("bolt://localhost:7687", auth=("neo4j", "a"))

test = "Newark Chesterton"
test2 = "Going Home"

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
    """
    Si t est dans la liste des genres de u alors return 2
    sinon return 1

    !!!
    Regler les  f string de merde
    Meilleur retour de querry ? to_data ??
    """

def fs(u,t):
    """
    Nombre de PPC * la taille des PPC
    """

def ft(u,t):
    """
    Cohérence thématique (Jaccard mes couilles)
    """
    #Get tracks liked by user u
    tracksliked = []
    q=(f"MATCH (p:Person{{name:{u}}})-[:LIKE]->(t:Track)"
       f"RETURN t.name")
    res = graph.run(q).to_table()
    for el in res : 
        trackliked.append(el[0])

    #computes jaccard coefficient for each track pair
    jaccard = []
    for track in tracksliked:
        q=(f"MATCH (p:Person{{name:'{u}'}})-[:LIKE]->(t1:Trac{{name:'{track}'}})<-[:WORD_OF]-(w:Word)-[:WORD_OF]->(t2:Track{{name:'{t}'}})"
           f"RETURN COUNT (distinct w.name)"
        )

        res = graph.run(q).to_table()
        intercept = res[O] 
        print(intercept)
        
        union = []
        q=(f"MATCH (p:Person{{name:'{u}'}})-[:LIKE]->(t1:Track{{name:'{track}'}})<-[:WORD_OF]-(w1:Word)"
           f"RETURN distinct w1.name AS word"
           f"UNION"
           f"MATCH (w2:Word)-[:WORD_OF]->(t2:Track{{name:'{t}'}})"
           f"RETURN distinct w2.name AS word"
        )
        res = graph.run(q).to_table()
        union = length(res)
        print(union)
        coeff = intercept/union
        jaccard.append(coeff)

    #Get max jaccard coefficient
    return max(jaccard)
    

def score(u, t,a=0.25, b=0.75):
    return fb(u,t) * (a * fs(u,t) + b * ft(u,t))

print(fb(test, test2))

aze = ft('Bogota Worcestershire', 'Going Home')
print(aze)