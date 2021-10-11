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