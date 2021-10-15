---
title: E. coli - Gestion de données non structurées - Applications post-génomiques
  - Mise en oeuvre
author: "Roland Barriot"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_notebook:
    toc: TRUE
    toc_float: TRUE
    toc_depth: 4
    theme: paper
    highlight: tango
  html_document:
    toc: TRUE
    toc_float: TRUE
    toc_depth: 4
    df_print: paged
---


```{r setup, include=TRUE}
#knitr::opts_chunk$set(echo = TRUE)
options(width = 120)
# knitr::opts_chunk$set(cache=F)
# knitr::opts_chunk$set(fig.width=6, fig.height=6) # for rstudio
# knitr::opts_chunk$set(fig.width=14, fig.height=14) # for html
``` 


```{r}
library(reticulate)
```


Illustration sur un organisme. Constitution d'une base dédiée à un organisme et son exploitation. 

# Approche par ensembles : Enrichissement 

Il a été identifié un ou quelques ensembles de gènes d'intérêt chez *E. coli*. A quels processus biologiques, fonctions moléculaires ou localisations sub-cellulaires peut-on les relier ?

Ensembles : 

  * [set 1](http://silico.biotoul.fr/enseignement/m2bbs/idh/2021/set.01.txt)
  * [set 2](http://silico.biotoul.fr/enseignement/m2bbs/idh/2021/set.02.txt)
  * [set 3](http://silico.biotoul.fr/enseignement/m2bbs/idh/2021/set.03.txt)


Utiliser l'interface Web http://silico.biotoul.fr/enrichment/ afin d'analyser ces ensembles. 

## Ensembles cibles et script de recherche d'enrichissement

Un premier exercice, relativement simple, va être de reproduire ce type d'analyse. Nous allons utiliser les mots-clés associer aux protéines dans UniProt. 

### Génération des ensembles cibles (UniProt keywords, Interpro domains, EcoCyc pathways, EcoCyc transcription units)

Aller sur [UniProt](https://www.uniprot.org) et télécharger le protéome correspondant à *E. coli K-12 MG1655* ; étapes :  
→ sélectionner `Proteomes` puis rechercher *E. coli* K-12 MG1655. Identifier le bon résultat et cliquer sur son `Proteome ID`  
→ Cliquer ensuite sur `View all proteins`  
→ Dans l'onglet `Columns`, sélectionner celles que nous allons utiliser par la suite : `Entry name`, `Gene names (ordered locus)`, `Gene Ontology IDs`, `Interpro`  
→ Puis télécharger au format `Tab-separated`  

Pour générer un fichier texte au format simple qui donne pour chaque mot-clé, la liste des protéines associées, nous allons utiliser R et les librairies `tidyverse` et `jsonlite`. 

#### Choix des identifiants de référence

Un des problèmes souvent rencontrés dans ce type d'analyse est d'arriver à identifier : les génomes et les gènes, protéines, *etc.* En effet, chaque source de données peut utiliser des méthodes de référencement et d'identification qui lui est propre. Pour le protéome d'*E. coli*, nous allons utiliser un type d'identifiant très souvent utilisé : les `bnumbers`.

A partir du fichier télécharger, vous allez donc le charger sous R et le reformater pour avoir un *mapping* des identifiants UniProt vers les *bnumbers*.

```{r}
library(tidyverse)
```


```{r}
uniprot = read_tsv("uniprot/uniprot.proteome.UP000000625.tab")
uniprot
```


A partir de ce *tibble*, il s'agit d'extraire les *bnumbers* de la colonne `Gene names (ordered locus)`. La fonction `str_extract` permet d'extraire un *pattern* à partir d'une expression régulière (ici ce sera `'b\\d+'`) d'une chaîne de caractères. 

```{r}
mapping = uniprot %>% 
  select(Entry, names = `Gene names  (ordered locus )`) %>%
  mutate(bnumber=str_extract(names, 'b\\d+')) %>%
  select(bnumber, uniprotID=Entry) %>%
  filter(!is.na(bnumber)) %>% # 2021 → P0DQD7 and P0A6D5 are lost (no bnumber)
  arrange(bnumber)
mapping
```

#### Reformatage des données

A présent, il s'agit de générer la liste des protéines associées à chaque mot-clé (colonne `Keywords`). Pour cela, on ne gardera que les colonnes `Entry` et `Keywords`, une jointure avec *mapping* va permettre de faire le mapping avec les *bnumbers*. La fonction `separate_rows` va nous servir à découper la colonne `Keywords` et répartir les mots-clés sur des lignes du tibble :
```{r}
keywords = uniprot %>% 
  select(uniprotID=Entry, keyword=Keywords) %>%
  right_join(mapping) %>% # right join to remove those without bnumber
  separate_rows(keyword, sep=';') %>%
  select(bnumber, keyword) %>%
  arrange(bnumber)
keywords
```

Maintenant, il ne s'agit plus que de regrouper les *bnumbers* par mot-clé :

```{r}
ref_sets = keywords %>% 
  group_by(keyword) %>%
  summarise(count=n(), elements = list(bnumber)) %>%
  ungroup %>%
  filter(count>1)  %>%
  select(id=keyword, desc=count, elements)
ref_sets
```


Pour utiliser "facilement" ces données avec python, nous allons sauvegarder le tibble au format JSON :
```{r}
library(jsonlite)
```

```{r}
ref_sets %>% toJSON %>% write("reference.sets/uniprot.keywords.sets.json")
```

### Script de recherche d'enrichissement

La recherche d'enrichissement se résume à comparer chacun des ensembles cibles (pour ce premier exemple, constitués des protéines associées à un même mot-clé) au moyen d'un test statistique. Nous utiliserons pour commencer une loi binomiale.

Ecrire un script python qui prendra en paramètre :

  * la liste ou le nom de fichier contenant les identifiants composant l'ensemble requête `query`
  * le nom de fichier contenant les ensembles cibles `target`

On vous fournit le début du script :
```{python eval=F}
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
```

Il s'agit ensuite de charger les identifiants de l'ensemble requête :
```{python eval =F}
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
```

Les ensembles cibles peuvent être chargés avec la librairie `json` (importée au début):
```{python eval=F}
# LOAD REFERENCE SETS
sets = json.loads(open(param.sets).read())
if param.verbose:
    print('first target sets: ', sets[0:2])
```

Afin d'appliquer un test avec la loi binomiale, nous avons besoin de connaître la taille de la population pour calculer la probabilité de succès ou d'échec à chaque tentative (nombre d'éléments pris au hasard dans la population, avec remise).


Fonction disponible dans `scipi.stats` → https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html

Rappel sur la loi binomiale : $k$ succès obtenus sur $n$ tentatives, chaque tentative ayant une probabilité $p$ de succès. Pour chaque ensemble cible $T$, ayant $C$ éléments en commun avec l'ensemble requête $Q$, il s'agit donc de considérer la probabilité d'avoir au moins $C$ succès en effectuant $Q$ tirages aléatoires avec remise, ayant chacun une probabilité de $T$/$G$ de succès, $G$ étant le nombre de gènes/protéines considérés pour le tirage.

Pour appliquer le test, il nous faut déterminer la probabilité de succès et donc le nombre total d'identifiants possibles. Ajouter une partie dans le script qui calcule le nombre d'identifiants "disponibles" parmi les ensembles cibles :
```{python eval=F}
# COMPUTE POPULATION SIZE
population = set()
for s in sets:
  # TO DO
```

Ensuite, il faut compléter la partie *TO DO* ci-dessous pour stocker la p-valeur obtenue pour chaque test.

```{python eval=F}
# EVALUATE SETS
results = []
query_size = len(query)
for s in sets:
	elements = set(s['elements' ])
	common_elements = elements.intersection( query )
	if param.measure=='binomial': # binom.cdf(>=success, attempts, proba)
		pvalue = 100000 # TO DO 
	else:
		print(f'sorry, {param.measure} not (yet) implemented')
		exit(1)
	r = { 'id': s['id'], 'desc': s['desc'], 'common.n':len(common_elements), 'target.n': len(elements), 'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements }
	results.append( r )
if param.verbose:
  print(results)
```

Une fois l'ensemble des ensembles cibles comparés et évalués, il s'agit d'afficher uniquement les résultats statistiquement significatifs :
```{python eval=F}
# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['p-value'])
for r in results:
	if r['p-value'] > param.alpha: 
	  break
	# OUTPUT
	print("{}\t{}\t{}/{}\t{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))
```

Exemple d'utilisation :
```
./enrichment/blastsets.json.py -q ecoli/query.sets/set.01.txt  -t ecoli/reference.sets/uniprot.keywords.sets.json 
Amino-acid biosynthesis	1.524947628271214e-115	102/102	102	b2020, b2026, b2021, b0077, b3672, b1264, b0388, b3958, b3960, b3828, b1262, b0073, b3237, b4388, b1693, b2024, b0273, b0261, b0243, b3281, b3940, b0075, b0242, b0071, b4254, b3671, b4488, b0166, b0908, b1622, b2413, b3213, b2019, b3766, b0529, b3773, b3008, b2472, b3809, b2599, b3744, b0386, b2329, b2601, b2764, b3390, b3212, b4019, b4054, b2421, b3829, b0754, b3772, b0072, b1275, b2913, b3959, b2818, b3607, b2763, b2018, b3941, b3939, b0032, b4024, b0031, b3433, b3771, b1704, b3774, b4243, b3769, b1263, b3359, b1261, b0004, b2022, b2414, b0078, b3938, b0002, b3670, b2551, b3172, b1265, b4013, b0003, b2838, b0907, b2478, b2025, b0001, b2600, b0674, b2023, b1260, b0074, b0033, b3770, b0159, b3389, b3957
Amino-acid transport	9.940384887710163e-87	79/79	79	b2131, b0810, b2308, b0813, b1798, b2310, b0652, b1453, b1917, b2306, b1918, b0576, b4208, b3269, b3110, b0861, b0863, b1533, b2679, b2309, b3460, b0199, b2678, b0112, b4077, b1729, b2670, b0197, b1296, b2923, b3116, b3271, b2663, b2307, b3268, b0874, b2578, b2677, b2365, b0654, b3270, b0486, b0653, b3161, b3370, b0692, b4115, b3456, b0811, b3089, b1920, b3795, b0862, b1492, b2130, b4132, b2128, b4141, b3458, b0864, b0809, b2156, b0655, b3457, b2014, b1015, b3454, b0401, b0198, b2129, b3709, b1907, b3653, b1336, b3455, b0402, b1605, b1473, b0860
Aromatic amino acid biosynthesis	1.5397106646021513e-20	19/19	19	b4054, b0754, b0908, b1264, b0388, b1265, b1262, b2600, b2599, b2329, b1704, b1260, b2601, b1693, b3390, b1263, b3281, b1261, b3389
Branched-chain amino acid biosynthesis	1.5397106646021513e-20	19/19	19	b0071, b3671, b4488, b0078, b3772, b0077, b0072, b3670, b3672, b3766, b0073, b3773, b3771, b3774, b4243, b0074, b3769, b3770, b0075
Transport	3.59923218981218e-18	81/747	747	b2131, b0810, b2308, b0813, b1798, b2310, b0652, b1453, b1917, b2306, b1918, b0576, b4208, b3269, b3110, b0861, b0863, b1533, b2679, b2309, b3460, b0199, b2678, b0112, b4077, b1729, b0197, b2413, b1296, b2923, b3116, b3271, b2663, b2307, b3268, b0874, b2578, b2677, b2365, b0654, b2764, b3270, b0486, b0653, b3161, b3370, b0860, b0692, b4115, b3456, b0811, b3089, b1920, b3795, b0862, b1492, b2130, b4132, b2128, b4141, b3458, b0864, b0809, b2156, b0655, b3457, b2014, b1015, b3454, b0401, b0198, b2129, b3709, b1907, b3653, b1336, b3455, b0402, b1605, b1473, b2670
Methionine biosynthesis	1.904483923982984e-15	14/14	14	b3941, b1622, b3939, b3008, b0261, b0159, b3938, b3829, b3433, b4019, b4013, b3828, b3940, b0529
Arginine biosynthesis	2.0255230660176893e-13	12/12	12	b4254, b0033, b0273, b3958, b3959, b3172, b0032, b3237, b2818, b3359, b3960, b3957
Histidine biosynthesis	2.1364884961775092e-11	10/10	10	b2018, b2021, b2023, b2022, b2025, b2020, b2024, b2026, b2019, b0529
Lysine biosynthesis	2.1913479560176622e-10	9/9	9	b2478, b2472, b4024, b3359, b3809, b0031, b2838, b3433, b0166
Direct protein sequencing	5.269908912363197e-09	73/934	934	b0242, b0071, b2020, b4254, b0655, b2414, b3938, b1015, b3829, b0754, b2310, b0002, b0166, b0077, b0072, b0860, b3213, b1264, b0388, b2913, b2551, b3959, b3172, b0811, b3089, b2306, b1920, b4013, b1262, b2763, b0529, b0073, b0907, b2478, b2025, b3939, b3008, b0032, b4024, b0031, b3458, b3744, b0386, b3433, b2329, b3771, b3774, b4243, b0674, b2601, b2764, b2023, b1693, b3769, b1260, b0074, b3390, b0273, b1263, b0033, b0863, b2156, b0243, b3281, b2679, b3212, b4019, b3359, b1261, b2309, b3940, b0004, b3460
Cysteine biosynthesis	2.3127031349291697e-08	7/7	7	b2413, b2764, b1275, b2414, b2421, b3607, b2763
Pyridoxal phosphate	3.212141355952543e-08	15/59	59	b1622, b2021, b0907, b3770, b2551, b2414, b3939, b3008, b4054, b2421, b3359, b1261, b2838, b3772, b0004
Tryptophan biosynthesis	2.3887776556004005e-07	6/6	6	b1260, b1264, b1263, b1265, b1261, b1262
Lyase	6.071889359377002e-07	23/169	169	b2022, b0071, b2414, b3772, b1622, b0072, b1264, b3960, b2838, b1262, b2478, b2025, b3008, b2599, b2329, b3771, b2023, b1260, b1693, b1263, b3389, b1261, b0004
NADP	1.4770687101029291e-06	16/91	91	b3941, b3774, b2329, b0386, b2764, b3213, b3958, b0243, b3281, b3212, b3433, b0031, b0002, b3940, b2763, b0529
Diaminopimelate biosynthesis	2.4866486372087227e-06	5/5	5	b2478, b2472, b0031, b3433, b0166
Leucine biosynthesis	2.4866486372087227e-06	5/5	5	b0073, b0072, b0071, b0074, b0075
Threonine biosynthesis	2.4866486372087227e-06	5/5	5	b0001, b0003, b0002, b3433, b0004
Allosteric enzyme	3.37483004700035e-05	9/38	38	b2478, b4254, b0199, b1264, b2414, b1263, b4024, b0002, b3772
Leader peptide	6.115758805811655e-05	6/16	16	b2018, b3672, b0001, b1265, b3766, b0075
Cell inner membrane	7.902115006894048e-05	63/977	977	b0810, b2308, b0199, b0486, b0653, b3457, b3161, b2014, b2678, b0813, b0112, b4077, b1015, b1729, b2670, b3370, b0652, b2413, b1453, b0401, b0692, b1296, b4115, b1917, b3456, b2923, b3116, b3089, b3271, b2663, b2306, b1918, b2307, b0874, b3795, b0862, b1492, b2130, b2578, b0198, b4132, b0576, b2677, b2128, b4141, b3709, b1907, b3653, b1336, b2365, b3269, b3110, b0861, b0654, b0864, b1473, b0809, b0402, b1605, b2156, b1533, b3270, b4208
Glutamine amidotransferase	0.00015695166263730433	5/12	12	b2023, b1263, b0032, b3212, b0674
Multifunctional enzyme	0.00016687777995440557	9/47	47	b2022, b2026, b1263, b2600, b2599, b0002, b3940, b1262, b0529
Proline biosynthesis	0.00028340060169450157	3/3	3	b0242, b0243, b0386
Serine biosynthesis	0.00028340060169450157	3/3	3	b0907, b2913, b4388
Cytoplasm	0.0006556873426611716	45/677	677	b0242, b2022, b4254, b2026, b4054, b3938, b0166, b0908, b1275, b0388, b2551, b3958, b3959, b3172, b2019, b3960, b2818, b4013, b3607, b3828, b0003, b0073, b0907, b2478, b2025, b3773, b3939, b3008, b2600, b3237, b3809, b0031, b2599, b3744, b0386, b4243, b3774, b2023, b2024, b3390, b0273, b0243, b3359, b3389, b3957
Symport	0.001283103138374326	7/39	39	b0692, b3653, b4077, b1015, b3089, b3116, b1729
Cell membrane	0.001283619289086577	65/1121	1121	b0810, b2308, b0199, b0486, b0653, b3457, b3161, b2014, b2678, b0813, b1798, b0112, b4077, b1015, b1729, b2670, b3370, b0652, b0197, b2413, b1453, b0401, b0692, b1296, b4115, b1917, b3456, b2923, b3116, b3089, b3271, b2663, b2306, b1918, b2307, b0874, b3795, b0862, b1492, b2130, b2578, b0198, b4132, b0576, b2677, b2128, b4141, b3709, b1907, b3653, b1336, b2365, b3269, b3110, b0861, b0654, b0864, b1473, b0809, b0402, b1605, b2156, b1533, b3270, b4208
Aminotransferase	0.0015503613834804647	5/20	20	b0907, b2021, b3770, b4054, b3359
Antiport	0.0015503613834804647	5/20	20	b4132, b4115, b0692, b1605, b1492
Transmembrane helix	0.0017279367747848216	56/939	939	b0810, b2308, b3457, b0486, b0653, b2014, b3161, b2678, b0813, b1798, b0112, b4077, b1015, b1729, b3370, b2413, b1453, b0401, b0692, b1296, b4115, b3456, b2923, b3116, b3089, b2663, b2307, b1918, b3795, b0874, b0862, b1492, b2130, b0198, b2578, b4132, b0576, b2128, b4208, b4141, b3709, b1907, b3653, b1336, b2365, b3269, b3110, b0861, b0654, b1473, b0402, b1605, b2156, b1533, b3270, b2670
Transferase	0.0026930137944873025	38/583	583	b0242, b4254, b2414, b3671, b0078, b4054, b2421, b3829, b0754, b0002, b0166, b2021, b0077, b0908, b3670, b2551, b0388, b3959, b2019, b2818, b4013, b3607, b0003, b0907, b3939, b4024, b1704, b0074, b2601, b3769, b3770, b3390, b0273, b0261, b1263, b3359, b4019, b3940
Asparagine biosynthesis	0.003204130883915414	2/2	2	b3744, b0674
Glutamate biosynthesis	0.003204130883915414	2/2	2	b3213, b3212
Isoleucine biosynthesis	0.003204130883915414	2/2	2	b4243, b3772
Transmembrane	0.005826055531125065	56/991	991	b0810, b2308, b3457, b0486, b0653, b2014, b3161, b2678, b0813, b1798, b0112, b4077, b1015, b1729, b3370, b2413, b1453, b0401, b0692, b1296, b4115, b3456, b2923, b3116, b3089, b2663, b2307, b1918, b3795, b0874, b0862, b1492, b2130, b0198, b2578, b4132, b0576, b2128, b4208, b4141, b3709, b1907, b3653, b1336, b2365, b3269, b3110, b0861, b0654, b1473, b0402, b1605, b2156, b1533, b3270, b2670
Thiamine pyrophosphate	0.013848266281607343	3/12	12	b3769, b0077, b3671
Membrane	0.01769991657552778	65/1249	1249	b0810, b2308, b0199, b0486, b0653, b3457, b3161, b2014, b2678, b0813, b1798, b0112, b4077, b1015, b1729, b2670, b3370, b0652, b0197, b2413, b1453, b0401, b0692, b1296, b4115, b1917, b3456, b2923, b3116, b3089, b3271, b2663, b2306, b1918, b2307, b0874, b3795, b0862, b1492, b2130, b2578, b0198, b4132, b0576, b2677, b2128, b4141, b3709, b1907, b3653, b1336, b2365, b3269, b3110, b0861, b0654, b0864, b1473, b0809, b0402, b1605, b2156, b1533, b3270, b4208
ATP-binding	0.025342292724889085	26/422	422	b0242, b0199, b2026, b3454, b0002, b0652, b1917, b0388, b3959, b3172, b2019, b3271, b2306, b0003, b2129, b2677, b0032, b4024, b3455, b3744, b0674, b0864, b0809, b0033, b3390, b3940
FAD	0.027004750944179177	7/70	70	b3941, b0077, b2764, b3671, b4488, b3212, b2329
One-carbon metabolism	0.04366622422295997	2/8	8	b2551, b0529
Cobalt	0.04874798966838455	4/33	33	b3957, b3389, b2472, b4019
```


### Autres ensembles d'ensembles cibles

Même chose à faire pour

  * les domaines Interpro à partir du protéome téléchargé UniProt Proteome
  * les GOTerms référencés dans le même fichier d'UniProt
  * les voies métaboliques (pathways) et les unités de transcription (TUs) que vous trouverez sur EcoCyc : 
    * gènes : https://www.ecocyc.org/group?id=biocyc13-55140-3842501533
    * pathways : https://www.ecocyc.org/group?id=biocyc17-55140-3842483872
    * unités de transcription : https://www.ecocyc.org/group?id=biocyc17-55140-3842483291

Une autre source qui peut s'avérer intéressante est la littérature biomédicale relative à chacun·e des gènes/protéines. Pour cela *Entrez* du NCBI va nous permettre de constituer, pour chaque article référençant au moins 2 gènes d'*E. coli*, un ensemble cible dans la section suivante.

### Ficher d'annotation du génome

A partir du site du NCBI, avec Entrez, rechercher *E. coli* K-12 MG1655 : https://www.ncbi.nlm.nih.gov/genome/?term=escherichia+coli+k-12+MG1655

Dans les *NCBI Datasets* obtenus, le problème, encore une fois est d'identifier sur la bonne version des données. Pour cet exemeple, nous utiliserons le génome de référence et son annotation RefSeq ASM584v2 → RefSeq:GCF_000005845.2

Télécharger le ficher au format GBFF. Nous aurons besoin ensuite du module biopython pour analyser le contenu et extraire les gènes codants (pour des protéines).

Nous allons ici ne conserver que les GeneID (pour chercher les références bibliographiques) et la localisation des gènes sur le chromosome :

**Attention :** pour ne pas surcharger le serveur du NCBI (et se faire bloquer), il faut passer par leur API et limiter le nombre de requêtes effectuées par seconde.

#### Mapping des identifiants

Encore une fois, il faut choisir les "bons" identifiants parmi ceux utiliser. Ici, nous utiliserons les *bnumbers*.

Pour extraire le mapping à partir du fichier GenBank : 
```{python eval=F}
#!/bin/env python
# conda install python3 biopython requests
# ./gbff.to.bnumber.GeneID.py > mapping.bnumber_ncbi.tsv

from Bio import SeqIO, Entrez

gbff = 'ncbi/data/GCF_000005845.2/genomic.gbff'
record = SeqIO.read(gbff, "genbank")
print('bnumber\tdbname\tdbid')
for f in record.features:
	if f.type=='CDS':
		f_xref = f.qualifiers['db_xref']
		f_id = f.qualifiers['locus_tag']
		for i in f_id:
			for j in f_xref:
				(dbsource, dbid) = j.split(':')
				print(f'{i}\t{dbsource}\t{dbid}')
```

On obtient à la sortie le *mapping* entre identifiants :
```{bash}
head mapping.bnumber_ncbi.tsv
```

#### Génération des ensembles cibles (PubMed ids)

Pour chaque gène/identifiant/GeneID,, l'utilitaire elink d'Entrez a été utilisé pour récupérer les publications relatives à chacun des gènes :

Le script suivant a été utilisé pour récupérer, à partir des GeneID, les identifiants PubMed relatifs à une séquence :     

Ceci n'est **PAS à réaliser** → récupérer plutôt l'ensemble des résultats sur http://silico.biotoul.fr/enseignement/m2bbs/idh/2021/bnumber.PMID.tsv
```{python eval=F}
#!/bin/env python

import sys
import datetime
import time
import requests
import os.path
import json

sleep_for = 1
n=0
last_sleep = time.time()
GeneIDs = set()
mapping = {}
with open('mapping.bnumber_ncbi.tsv', 'r') as f:
	for line in f.readlines():
		(bnumber, dbname, dbid) = line.strip().split('\t')
		if dbname == 'GeneID':
			GeneIDs.add(dbid)
			mapping[dbid] = bnumber
			resfile = 'ncbi/pubmed/GeneID.'+dbid+'.PMIDlinks.json'
			if os.path.isfile(resfile):
				sys.stderr.write('%s file %s exists for %s skipped\n' % (datetime.datetime.now(), resfile, dbid))
			else:
				url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=json&dbfrom=gene&db=pubmed&id='+dbid
				sys.stderr.write(f'fetching pmid links of {dbid} from {url}\n')
				r = requests.get(url)
				fh = open(resfile, mode='w')
				fh.write(r.text)
				fh.close()
				# timing
				n += 1
				elapsed = (time.time() - last_sleep)
				rps = n / elapsed
				sys.stderr.write("  %3d requests in %2.2f seconds → %2.1f requests per second" % (n, elapsed, rps))
				if rps > 3 : # max 3 requests/second at NCBI = 9 req in 3 sec
					sys.stderr.write("sleeping for %s seconds\n" % (sleep_for))
					time.sleep(sleep_for)
					last_sleep = time.time()
					n=0
				elif n>=10: # reset counters after 10 requests
					last_sleep = time.time()
					n=0

print('bnumber\tPMID')
for dbid in GeneIDs:
	sys.stderr.write("extracing PMIDs for %s\n" % (dbid))
	# extract PMIDs from XML file
	resfile = 'ncbi/pubmed/GeneID.'+dbid+'.PMIDlinks.json'
	f = open(resfile, 'r')
	d = json.load(f)
	for ls in d['linksets']:	
		for lsd in ls['linksetdbs']:
			if lsd['linkname'] == "gene_pubmed":
				for i in lsd['links']:
					print(f'{mapping[dbid]}\t{i}')
```

Nous obtenons donc le fichier associant les identifiants de gène aux identifiants PubMed :
```{bash}
head ncbi/bnumber.PMID.tsv
```

Utiliser la même technique que pour les mots-clés UniProt afin de générer le fichier au format JSON correspondant associant à un PMID la liste des gènes (*bnumbers*) associée.

```{r echo=F, warning=F, message=F}
pmids = read_tsv("ncbi/bnumber.PMID.tsv", col_types = "cc")
pmids %>%
  group_by(PMID) %>%
  summarise(count=n(), elements = str_c(sort(bnumber), collapse=' ')) %>%
  ungroup %>%
  filter(count>1)  %>%
  mutate(desc=paste0("https://pubmed.ncbi.nlm.nih.gov/", PMID)) %>%
  select(id=PMID, desc, elements)
```

Effecctuer un travail similaire à partir du fichier d'annotation GenBank afin d'extraire la localisation des gènes sur le chromosome et obtenir le tibble suivant :
```{r echo=F, message=F, warning=F}
read_tsv("neo.data/ncbi.bnumber.location.tsv")
```

### Intégration dans Neo4J

![Schéma de la base de données](graph_schema.png)



#### Neo4j installation

Sites et documentations :

  * https://neo4j.com/
  * https://neo4j.com/developer/get-started/ 

Récupérer la dernière version : Téléchargement linux https://neo4j.com/download-center/#releases (onglet *Community server*)

Création des répertoires et désarchivage (shell):
```{bash eval=F}
version=4.3.5
tar xf neo4j-community-$version-unix.tar.gz
ln -s neo4j-community-$version neo4j
cd neo4j/
```

Démarrage et arrêt du serveur (shell):
```{bash eval=F}
./bin/neo4j console
```

Le processus est au premier plan donc pour arrêter le serveur il faut faire `Ctrl + C` dans le terminal.

Utilisation depuis le navigateur (vérifier le port renseigné lors de la précédente commande) : http://localhost:7474/

A la première connexion, le mot de passe est `neo4j`, le système demande ensuite de changer le mot de passe. Explorez l'interface Web de Neo4j browser, notamment le côté gauche avec les paramètres, et les informations sur la base de données. 

Utilisation depuis le shell :
```{bash eval=F}
./bin/cypher-shell
```

Aller sur http://localhost:7474 une première fois afin de définir le mot de passe : au début login: neo4j passwd: neo4j

Librairie neo4r https://github.com/neo4j-rstats/neo4r et https://neo4j-rstats.github.io/user-guide/
```{r}
if (!require('neo4r')) { # client neo4j (pas dispo avec conda)
  install.packages('neo4r')
}
library(neo4r)
neodb = neo4j_api$new(
  url = "http://localhost:7474", 
  user = "neo4j", 
  password = "bioinfo"
)
cypher = function(query, neo4j=neodb, ...) call_neo4j(query=query, con=neo4j, ...)
# 'MATCH (p:Protein {name:"dnaJ"}) RETURN p.protein_id, p.annotation' %>% cypher
```

Créer un lien symbolique pour neo4j/import 

##### Création du graphe

**Genes**
```{r}
genes = read_tsv("neo.data/ncbi.bnumber.location.tsv")
genes$organism=511145
genes %>% rename(gene_id=bnumber) %>% write_csv("neo4j.import/ncbi.bnumber.location.csv")
```


```{r}
# IMPORT
'MATCH (n:Gene) DELETE n' %>% cypher
'LOAD CSV WITH HEADERS FROM "file:///ncbi.bnumber.location.csv" AS row 
CREATE (n:Gene)
SET n = row,
 n.id = row.gene_id,
 n.taxon_id = toInteger(row.organism),
 n.rank = toInteger(row.rank),
 n.strand = row.strand,
 n.begin = toInteger(row.begin),
 n.end = toInteger(row.end)
'  %>% cypher
```

Index
```{r}
'CREATE INDEX ON :Gene(id)' %>% cypher
```
**Keywords**
```{r}
keywords %>% 
  select(keyword) %>%
  unique %>%
  write_csv("neo4j.import/uniprot.keywords.csv")
'LOAD CSV WITH HEADERS FROM "file:///uniprot.keywords.csv" AS row 
CREATE (n:Keyword)
SET n = row,
 n.id = row.keyword
'  %>% cypher
```

```{r}
'CREATE INDEX ON :Keyword(id)' %>% cypher
```

Links Keyword → Gene
```{r}
keywords %>% write_csv("neo4j.import/uniprot.keywords.genes.csv")
```

```{r}
"MATCH (:Keyword)-[r:describes]->(:Gene) DELETE r" %>% cypher
"
LOAD CSV WITH HEADERS FROM 'file:///uniprot.keywords.genes.csv' AS line 
MATCH (k:Keyword),(g:Gene) 
WHERE k.id=line.keyword AND g.id=line.bnumber
WITH k,g 
MERGE (k)-[:describes]->(g)
" %>% cypher
```
**TUs**


**Pathways**

**PubMed**
```{r}
pmids %>% select(PMID) %>% unique %>% write_csv("neo4j.import/ncbi.pmid.csv")
pmids %>% rename(gene_id = bnumber) %>% write_csv("neo4j.import/ncbi.pmid.genes.csv")
```

```{r}
'MATCH (n:PubMed) DELETE n' %>% cypher
'LOAD CSV WITH HEADERS FROM "file:///ncbi.pmid.csv" AS row 
CREATE (n:PubMed)
SET n = row,
 n.id = row.PMID
'  %>% cypher

'CREATE INDEX ON :PubMed(id)' %>% cypher

"MATCH (:PubMed)-[r:cites]->(:Gene) DELETE r" %>% cypher
"
LOAD CSV WITH HEADERS FROM 'file:///ncbi.pmid.genes.csv' AS line 
MATCH (p:PubMed),(g:Gene) 
WHERE p.id=line.PMID AND g.id=line.gene_id
WITH p,g 
MERGE (p)-[:cites]->(g)
" %>% cypher

```

### Adaptation du script pour l'utilisation de neo4j

Adapter ou créer un nouveau qui utilisera, non plus les fichiers au format JSON, mais une connexion au serveur Neo4J pour effectuer les recherches d’enrichissement.

