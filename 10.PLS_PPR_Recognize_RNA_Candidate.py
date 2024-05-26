# -*- encoding: utf-8 -*-
'''
@File    :   PLS_PPR_Recognize_RNA_Candidate.py
@Time    :   2023/05/12 16:52:14
@Author  :   Ming Chen
@Contact :   chenm@big.ac.cn
'''

# Put the import lib here.
import sys
from Bio.Seq import Seq

def collect_all_sequence(organelle_fasta):
    '''
    Get organelle genome fasta and change multi-line sequence into one-line sequence.
    '''
    mit_seq = []

    for line in open(organelle_fasta, 'r'):
        new_line = line.strip()

        if line[0] == ">":
            continue
        else:
            mit_seq.append(new_line)
    all_mit_seq = ''.join(mit_seq)

    return all_mit_seq


def seize_editing_sites_strand(editing_out):
    '''
    Get the RNA editing site and the corresponding strand information.
    '''
    site_strand = {}
    site_gene_name = {}

    for line in open(editing_out, 'r'):
        new_line = line.strip()
        cols = new_line.split('\t')

        if cols[0] == "editSite":
            continue
        else:
            units = cols[0].split("_")
            site_strand[units[1]] = units[2]
            site_gene_name[units[1]] = cols[1]

    return site_strand, site_gene_name


def ppr_motif_score():
    '''
    Based on the data in (Gutmann et al., 2021, Cells), get each base (A, C, G, U) score acoording to each PPR code.
    '''
    motif_score = {}

    for line1 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\P1.tsv", 'r'):
        new_line1 = line1.strip()
        col1 = new_line1.split('\t')
        if new_line1[0:3] == "For" or new_line1[0:8] == "5th/last":
            continue
        else:
            motif_score["P1-"+col1[0]+"-A"] = col1[1]
            motif_score["P1-"+col1[0]+"-C"] = col1[2]
            motif_score["P1-"+col1[0]+"-G"] = col1[3]
            motif_score["P1-"+col1[0]+"-U"] = col1[4]

    for line2 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\P2.tsv", 'r'):
        new_line2 = line2.strip()
        col2 = new_line2.split('\t')
        if new_line2[0:3] == "For" or new_line2[0:8] == "5th/last":
            continue
        else:
            motif_score["P2-"+col2[0]+"-A"] = col2[1]
            motif_score["P2-"+col2[0]+"-C"] = col2[2]
            motif_score["P2-"+col2[0]+"-G"] = col2[3]
            motif_score["P2-"+col2[0]+"-U"] = col2[4]

    for line3 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\L1.tsv", 'r'):
        new_line3 = line3.strip()
        col3 = new_line3.split('\t')
        if new_line3[0:3] == "For" or new_line3[0:8] == "5th/last":
            continue
        else:
            motif_score["L1-"+col3[0]+"-A"] = col3[1]
            motif_score["L1-"+col3[0]+"-C"] = col3[2]
            motif_score["L1-"+col3[0]+"-G"] = col3[3]
            motif_score["L1-"+col3[0]+"-U"] = col3[4]

    for line4 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\L2.tsv", 'r'):
        new_line4 = line4.strip()
        col4 = new_line4.split('\t')
        if new_line4[0:3] == "For" or new_line4[0:8] == "5th/last":
            continue
        else:
            motif_score["L2-"+col4[0]+"-A"] = col4[1]
            motif_score["L2-"+col4[0]+"-C"] = col4[2]
            motif_score["L2-"+col4[0]+"-G"] = col4[3]
            motif_score["L2-"+col4[0]+"-U"] = col4[4]

    for line5 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\S1.tsv", 'r'):
        new_line5 = line5.strip()
        col5 = new_line5.split('\t')
        if new_line5[0:3] == "For" or new_line5[0:8] == "5th/last":
            continue
        else:
            motif_score["S1-"+col5[0]+"-A"] = col5[1]
            motif_score["S1-"+col5[0]+"-C"] = col5[2]
            motif_score["S1-"+col5[0]+"-G"] = col5[3]
            motif_score["S1-"+col5[0]+"-U"] = col5[4]

    for line6 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\S2.tsv", 'r'):
        new_line6 = line6.strip()
        col6 = new_line6.split('\t')
        if new_line6[0:3] == "For" or new_line6[0:8] == "5th/last":
            continue
        else:
            motif_score["S2-"+col6[0]+"-A"] = col6[1]
            motif_score["S2-"+col6[0]+"-C"] = col6[2]
            motif_score["S2-"+col6[0]+"-G"] = col6[3]
            motif_score["S2-"+col6[0]+"-U"] = col6[4]

    for line7 in open("E:\\Rice_Endosperm_RNA_Editing\\15.PPR_RNA_Editing\\03.PPR_RNA_interaction\\02.PPR_Score\\PPRmatcher-main\\scoring_tables\\Millman\\SS.tsv", 'r'):
        new_line7 = line7.strip()
        col7 = new_line7.split('\t')
        if new_line7[0:3] == "For" or new_line7[0:8] == "5th/last":
            continue
        else:
            motif_score["SS-"+col7[0]+"-A"] = col7[1]
            motif_score["SS-"+col7[0]+"-C"] = col7[2]
            motif_score["SS-"+col7[0]+"-G"] = col7[3]
            motif_score["SS-"+col7[0]+"-U"] = col7[4]

    return motif_score


def seize_editing_site_adjancent_sequence(organelle_fasta, editing_out, ppr_code_file, output_path):
    '''
    Score the upstream adjancent sequence of RNA editing site.
    '''

    motif_score_dict = ppr_motif_score()

    all_organelle_seq = collect_all_sequence(organelle_fasta)
    site_strand = seize_editing_sites_strand(editing_out)[0]
    site_gene_name = seize_editing_sites_strand(editing_out)[1]

    for line in open(ppr_code_file, 'r'):
        new_line = line.strip()
        cols = new_line.split('\t')
        ppr_motif = cols[1].split('-')
        ppr_code = cols[2].split('-')
        rna_base_num = int(cols[3]) + 4
        
        seq_out = open(output_path + "\\" + cols[0]+"_RNA_Binding_Score.txt", 'w')
        
        for key, value in site_strand.items():
            if value == "0":
                # RNA editing site on the minus-strand
                StartPosition = int(key) - 1
                EndPosition = int(key) + int(rna_base_num) - 1
                editing_site_adjacent = all_organelle_seq[StartPosition:EndPosition]
                adjacent_sequence = Seq(editing_site_adjacent)
                real_editing_site_adjacent = str(adjacent_sequence.reverse_complement())
                replace_t_u = real_editing_site_adjacent.replace('T', 'U')

                rna_binding_score = 0
                # cols[3]: The number of PLS.
                for i in range(int(cols[3])):
                    rna_binding_score += float(motif_score_dict[ppr_motif[i]+"-"+ppr_code[i]+"-"+replace_t_u[i]])
                seq_out.write(key + "\t" + site_gene_name[key] + "\t" + replace_t_u + "\t" + str(rna_binding_score) + "\n")

            elif value == "1":
                # RNA editing site on the forward-strand
                StartPosition = int(key) - int(rna_base_num)
                EndPosition = int(key)
                editing_site_adjacent = all_organelle_seq[StartPosition:EndPosition]
                replace_t_u = editing_site_adjacent.replace('T', 'U')

                rna_binding_score = 0
                # cols[3]: The number of PLS.
                for i in range(int(cols[3])):
                    rna_binding_score += float(motif_score_dict[ppr_motif[i]+"-"+ppr_code[i]+"-"+replace_t_u[i]])
                seq_out.write(key + "\t" + site_gene_name[key] + "\t" + replace_t_u + "\t" + str(rna_binding_score) + "\n")

            elif value == "2":
                '''RNA editing site on the minus-strand or/and forward-strand'''

                # minus-strand
                StartPosition = int(key) - 1
                EndPosition = int(key) + int(rna_base_num) - 1
                editing_site_adjacent = all_organelle_seq[StartPosition:EndPosition]
                adjacent_sequence = Seq(editing_site_adjacent)
                real_editing_site_adjacent = str(adjacent_sequence.reverse_complement())
                replace_t_u = real_editing_site_adjacent.replace('T', 'U')

                # forward-strand
                StartPosition2 = int(key) - int(rna_base_num)
                EndPosition2 = int(key)
                editing_site_adjacent2 = all_organelle_seq[StartPosition2:EndPosition2]
                replace_t_u2 = editing_site_adjacent2.replace('T', 'U')

                if replace_t_u[-1] == 'C':
                # if replace_t_u2[-1]:
                    rna_binding_score = 0
                    # cols[3]: The number of PLS.
                    for i in range(int(cols[3])):
                        rna_binding_score += float(motif_score_dict[ppr_motif[i]+"-"+ppr_code[i]+"-"+replace_t_u[i]])

                    if ";" in site_gene_name[key]:
                        gene_names = site_gene_name[key].split('; ')
                        same_strand_site_gene_name = []
                        for name in gene_names:
                            if "-" in name:
                                same_strand_site_gene_name.append(name.split("-")[0])
                        seq_out.write(key + "\t" + ','.join(same_strand_site_gene_name) + "\t" + replace_t_u + "\t" + str(rna_binding_score) + "\n")

                    elif "-" in site_gene_name[key]:
                        seq_out.write(key + "\t" + (site_gene_name[key]).split("-")[0] + "\t" + replace_t_u + "\t" + str(rna_binding_score) + "\n")

                    elif site_gene_name[key] == "Intergenic":
                        seq_out.write(key + "\t" + site_gene_name[key] + "\t" + replace_t_u + "\t" + str(rna_binding_score) + "\n")

                if replace_t_u2[-1] == 'C':
                # if replace_t_u2[-1]:

                    rna_binding_score = 0

                    # cols[3]: The number of PLS.
                    for i in range(int(cols[3])):
                        rna_binding_score += float(motif_score_dict[ppr_motif[i]+"-"+ppr_code[i]+"-"+replace_t_u2[i]])
                    
                    if ";" in site_gene_name[key]:
                        gene_names = site_gene_name[key].split('; ')
                        same_strand_site_gene_name = []
                        for name in gene_names:
                            if "+" in name:
                                same_strand_site_gene_name.append(name.split("+")[0])
                        seq_out.write(key + "\t" + ','.join(same_strand_site_gene_name) + "\t" + replace_t_u2 + "\t" + str(rna_binding_score) + "\n")

                    elif "+" in site_gene_name[key]:
                        seq_out.write(key + "\t" + (site_gene_name[key]).split("+")[0] + "\t" + replace_t_u2 + "\t" + str(rna_binding_score) + "\n")

                    elif site_gene_name[key] == "Intergenic":
                        seq_out.write(key + "\t" + site_gene_name[key] + "\t" + replace_t_u2 + "\t" + str(rna_binding_score) + "\n")

            else:
                print("Check the whole process!")

def main():
    '''
    sys.argv[1]: organelle genome fasta
    sys.argv[2]: editing site indentification output file
    sys.argv[3]: all_ppr.txt, including PPR gene name, PPR motif, PPR code and PLS motif number
    sys.argv[4]: ouput path
    '''
    seize_editing_site_adjancent_sequence(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

if __name__ == '__main__':
    main()