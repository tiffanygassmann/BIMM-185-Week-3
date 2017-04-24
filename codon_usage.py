#Tiffany Gassmann
#BIMM 185
#Exercise 3: Codon Freqeunce


from collections import defaultdict

import matplotlib.pyplot as plt

#numnber of total codons for each gene, used for CUI
codon_master = {}

def main ():


    out = gen_table()
    with open('count.txt', 'w') as file:
        for row in out:
            print >> file, row

    return


#opens file and gets the gene, sequence pairs for gene file
def get_sequence():
    file = open("genes")
    genes_condon = file.readlines()

    for i in xrange(len(genes_condon)):

        # making table into list of lists without tab/newline
        # [["boo1", "ATG.......TGA],["b002', "ATG........]
        genes_condon[i] = genes_condon[i].strip().split("\t")

    return genes_condon


#function that handles counting  codons
def codon_counter():

    # DNA codon table from http://en.wikipedia.org/wiki/Genetic_code
    degenerated =  ['ATG','GCT', 'GCC', 'GCA', 'GCG',
     'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
     'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
     'AAA', 'AAG', 'AAT', 'AAC', 'GAT', 'GAC',
     'TTT', 'TTC', 'TGT', 'TGC', 'CCT', 'CCC', 'CCA', 'CCG',
     'CAA', 'CAG', 'TCT', 'TCC', 'TGG' , 'TCA', 'TCG', 'AGT', 'AGC',
     'GAA', 'GAG', 'ACT', 'ACC', 'ACA', 'ACG',
     'GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC', 'TAT', 'TAC',
     'ATT', 'ATC', 'ATA', 'GTT', 'GTC', 'GTA', 'GTG',
     'TAA', 'TGA', 'TAG']

    #CUI init
    for codon in degenerated:
        codon_master[codon] = 0

    #get sequences for each gene (gene,sequence)
    gene_sequences = get_sequence()

    #dictionary of codon counts
    gene_codon_counts = {}


    #calcululating count
    for pair in gene_sequences:
        for seq in pair:
            seq = pair[1]
            name = pair[0]
            max_seq = len(pair[1])

            #seperate sequence into codons
            query_codons = [seq[i:i + 3] for i in range(0, max_seq, 3)]

            #prepare dictionary of counts:
            counts = {}

            #initilialization
            for codon in degenerated:
                counts[codon] = 0

            for codon in query_codons:
                #calculate all counts for every codon
                codon_master[codon] += 1
                #counts for gene
                counts[codon] +=1
            gene_codon_counts[name] = counts
    return gene_codon_counts




#table formating function
def gen_table():
    #codon header order
    codons = ['ATG', 'GCT', 'GCC', 'GCA', 'GCG',
              'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
              'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
              'AAA', 'AAG', 'AAT', 'AAC', 'GAT', 'GAC',
              'TTT', 'TTC', 'TGT', 'TGC', 'CCT', 'CCC', 'CCA', 'CCG',
              'CAA', 'CAG', 'TCT', 'TCC', 'TGG', 'TCA', 'TCG', 'AGT', 'AGC',
              'GAA', 'GAG', 'ACT', 'ACC', 'ACA', 'ACG',
              'GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC', 'TAT', 'TAC',
              'ATT', 'ATC', 'ATA', 'GTT', 'GTC', 'GTA', 'GTG',
              'TAA', 'TGA', 'TAG']

    header = "Locus" + '\t'

    for cod in codons:
        header += cod + "\t"

    final_h = header + "Length" + "\t" + "CUI"

    #new dict of codons so all are present in table
    dict = codon_counter()
    codon_list = []

    #append all codons
    for key in dict:
        codon_list.append((key, dict[key]))

    #sort by first in tuple (gene name)
    codon_list = sorted(codon_list, key = lambda x:x[0])

    rows = []
    rows.append(final_h)

    #summing all counts of all codons
    T_g = float(sum(codon_master[codon] for codon in codon_master))
    #table of P_G - probability of codon c in genome g
    P_g = {}
    for key in codon_master:
        P_g[key] = float(codon_master[key]/T_g)


    #actual table
    for pair in codon_list:
        length = 0
        row = pair[0] + '\t'
        for codon in codons:
            number_codons = pair[1][codon]
            row += str(number_codons) + '\t'
            length += number_codons
        row += str(length) + '\t'

        # calculating CUI
        CUI_g = 0.0

        for codon in codons:
            #prob of codon x in gene i/ total number of codons in gene
            q_c = float(pair[1][codon]/float(length))
            CUI_g += q_c * P_g[codon]

        row += str(CUI_g)

        rows.append(row)

    return rows


def make_plots():
    #CUI vs genes in chromosomal order

    rows = gen_table()
    table = [row.split('\t') for row in rows[1:]]

    #genes in order
    x_vals = [int(row[0][1:]) for row in table]
    #CUI- values
    y_vals = [float(row[-1]) for row in table]

    plt.plot(x_vals, y_vals, "ro")
    plt.axis([0,5000,0.0, 0.03])

    plt.xlabel('Genes')
    plt.ylabel('CUI')
    plt.show()

    #CUI vs genes sorted by CUI

    gene_cui = zip(x_vals,y_vals)
    gene_cui = sorted(gene_cui, key = lambda x:x[1])


    plt.plot([i for i in range(len(gene_cui))], [pair[1] for pair in gene_cui] , "ro")
    plt.axis([0, 2700, 0.0, 0.03])

    plt.xlabel('Genes')
    plt.ylabel('CUI')
    plt.show()


    #Historam of CUI numbers

    plt.hist(y_vals)
    plt.xlabel("CUI Value")
    plt.ylabel("Frequency")
    plt.show()

    return

make_plots()


if __name__ == '__main__':
    main()
