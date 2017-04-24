from Bio import SeqIO

import sys, gzip, re
from itertools import imap


#Creates File to test each sectiona and to view
def file_handle():

    file = gzip.open("e_coli_genome.gbff.gz")

    gb_record = SeqIO.read(file,"genbank")

    return gb_record

#Creates file which we use BioPython Methods on
def file_parse():
    file = gzip.open("e_coli_genome.gbff.gz")

    record = SeqIO.parse(file, "genbank")

    return record


def extract_source(gb_record):
    my_source = gb_record.features[0]
    return my_source

def extract_gene(gb_record):
    my_gene = gb_record.features[1]
    return my_gene

def extract_cds(gb_record):
    my_CDS = gb_record.features[2]
    return my_CDS



def gen_header():
    Header = 'Accession' + '\t' + "Coordinates" + '\t' "Strand" + '\t' +  "Gene Name" + '\t' + "Locus Tag" + '\t' + "Synonyms"+ "\t"+ "Protein Name" + '\t' + "EC - Numbers" +'\t' + "External Refs"
    return Header

#Extract source info -  taxon
def source_info(record):

    rec = next(record)
    for f in rec.features:
        if f.type == 'source':
            taxon = (f.qualifiers['db_xref'])
    for i in taxon:
        return i

#Extract CDs Information -

def cds_info(record):

    #go through each entry
    rec = next(record)
    for f in rec.features:
        if f.type == 'CDS':

            #check for key error if no proein ID exists: return "Pseudo"
            if 'protein_id' in f.qualifiers:
                protein_ids = '\t'.join(f.qualifiers['protein_id'])
            else: protein_ids = "Pseudo"

            #Location and Strand Direction taken from Location NOT qualifiers
            locations_strands = (str(f.location))
            #regular expression to seperate the start, stop, direction
            reg = "\[([0-9]+):([0-9]+)\]\((.)\)"
            loc_str = re.findall(reg, locations_strands)

            #Genes, Locus Tags, Synonyms taken from Qualifiers
            genes = "\t".join(f.qualifiers['gene'])
            locus_tags = '\t'.join(f.qualifiers['locus_tag'])
            gene_synonyms = '\t'.join(f.qualifiers['gene_synonym'])

            #check for key error if no product exists: return "pseduo"
            if 'product' in f.qualifiers:
                products = '\t'.join(f.qualifiers['product'])
            else: products = "Pseudo"

            #Ext Refs taken from Qualifiers
            ext_refs = '\t'.join(f.qualifiers['db_xref'])

            #Join the Start, Stop, Direction into single string
            string = '\t'.join(map(str, loc_str))

            print '\t'.join([protein_ids, string, genes,locus_tags,gene_synonyms,products, ext_refs])


#Output
print "Tax ID: ", source_info(file_parse())
print gen_header()
cds_info(file_parse())


#Tests to view each section in the larger file and export each to new text file

#Gene Section
gene = extract_gene(file_handle())
with open('gene.txt', 'w') as file:
    print >> file, gene
#CDS Section
cds = extract_cds(file_handle())
with open('cds.txt', 'w') as file:
    print >> file, cds
#Source Section
source = extract_source(file_handle())
with open('source.txt', 'w') as file:
    print >> file, source