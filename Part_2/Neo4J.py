#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from py2neo import Graph

# Connection to Neo4J
try:
    graph = Graph("bolt://localhost:7687", auth=("neo4j", "a"))
except:
    print("Connection to the database failed!\nPlease check your password and try again.\n--help for more information")
    exit()