from py2neo import Graph

graph = Graph("bolt://localhost:7687", auth=("neo4j", "a"))

test = "{name:'Newark Chesterton'}"
test2 = "{name:'Going Home'}"

def fb(u,t):
    listGenre = []
    q = f"match (:Person{u})-[]->(:Track )-[]->(:Album)-[]->(:Artist)-[]->(g:Genre) return distinct g.name"
    res = graph.run(q).to_table()
    for el in res : 
        listGenre.append(el[0])
    q = f"match (:Track{t})-[]->(:Album)-[]->(:Artist)-[]->(g:Genre) return distinct g.name"
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

def score(u, t,a=0.25, b=0.75):
    return fb(u,t) * (a * fs(u,t) + b * ft(u,t))

print(fb(test, test2))

