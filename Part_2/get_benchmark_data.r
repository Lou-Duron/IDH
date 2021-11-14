# Benchmark dataset
 
 
####### LIBRARY #####
setwd('/home/victoriafathi/Bureau/M2/m2idh/Part_2/') #to modify
library("tidyverse")
####### LIBRARY #####


# LOAD DATA
benchmark_set = read_tsv('sets/uniprot-helicaseorganism__EscherichiacolistrainK12ECOLI--.tab')
uniprot_mapping = read_tsv('import/Alias/uniprot_mapping.tsv')
ecocyc_mapping = read_tsv('import/Alias/ecocyc_mapping.tsv')


head(benchmark_set)


# Mapping 
# The first step is to keep genes annotated in Uniprot genome annotation of 
# E.coli and annotated id Ecocyc

mapping = benchmark_set %>% 
  select(uniprotID=Entry) %>%
  inner_join(uniprot_mapping) %>%
  select(bnumber, uniprotID) %>%
  inner_join(ecocyc_mapping)

dim(mapping) # check

# Conversion to sets format: space separated
bnumber_bm = paste(mapping$bnumber, sep=" ", collapse=" ")
#write_file(bnumber_bm, "benchmark_data/bnumber_bm.txt")


# GoTerm benchmark
GO_genes_bm = benchmark_set %>% 
  select(uniprotID=Entry, GO="Gene ontology IDs") %>%
  right_join(mapping) %>%
  select(bnumber, GO) %>% 
  separate_rows(GO, sep='; ')%>%
  arrange(bnumber)
length(unique(GO_genes_bm$bnumber)) #check
#write_delim(GO_genes_bm, "benchmark_data/GO_bm.tsv", delim='\t')


# Interpro Domains 
domains_genes_bm = benchmark_set %>% 
  select(uniprotID=Entry, interpro='Cross-reference (InterPro)') %>%
  right_join(mapping) %>% # right join to remove those without bnumber
  separate_rows(interpro, sep=';') %>%
  select(bnumber, interpro) %>%
  filter(interpro != '') %>%
  arrange(bnumber)

length(unique(domains_genes_bm$bnumber)) #check
#write_delim(domains_genes_bm, "benchmark_data/interpro_bm.tsv", delim='\t')


# Keywords
keywords_genes_bm = benchmark_set %>% 
  select(uniprotID=Entry, keyword=Keywords) %>%
  right_join(mapping) %>% # right join to remove those without bnumber
  separate_rows(keyword, sep=';') %>%
  select(bnumber, keyword) %>%
  arrange(bnumber)
length(unique(GO_genes_bm$bnumber)) #check
#write_delim(keywords_genes_bm, "benchmark_data/keyword_bm.tsv", delim='\t')


# TU
TU_genes = read_tsv("import/TU/TU_link.tsv")
TU_genes_bm = mapping %>%
  select(bnumber) %>%
  left_join(TU_genes)
length(unique(TU_genes_bm$bnumber)) #Check
#write_delim(TU_genes_bm, "benchmark_data/TU_bm.tsv", delim='\t')


# Pathway
Pathway_genes = read_tsv('import/Pathway/pathway_link.tsv')
Pathway_genes_bm = mapping %>%
  select(bnumber) %>%
  left_join(Pathway_genes)
length(unique(Pathway_genes_bm$bnumber))
#write_delim(Pathway_genes_bm, "benchmark_data/Pathway_bm.tsv", delim='\t')


# PubMed
PubMed = read_tsv('import/PubMed/pubmed_link.tsv')
Pubmed_bm = mapping %>%
  select(bnumber) %>%
  left_join(PubMed)
length(unique(Pubmed_bm$bnumber))
#write_delim(Pubmed_bm, "benchmark_data/Pubmed_bm.tsv", delim='\t')


