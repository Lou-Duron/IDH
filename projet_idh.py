from py2neo import Graph
graph = Graph("bolt://localhost:7687", auth=("neo4j", "votremotdepassepourlabase"))
q = "match (t:Track) return t.name as track"
res = graph.run(q).to_table()
for record in res :
    print(record[0])