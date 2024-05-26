# -*- encoding: utf-8 -*-
'''
@File    :   12.mitochondrial_editing_cds_synonymous.py
@Time    :   2023/08/28 21:28:02
@Author  :   Ming Chen
@Contact :   chenm@big.ac.cn
'''

# Put the import lib here.
import re, os, argparse, collections
from Bio.Seq import Seq


def CDS_Editing_Info(annotation_gff3_dir):
    """
    Get the start and end position of the target gene.
    """

    # Dict of species name and CDS info.
    species_name_cds_info = collections.defaultdict(list)

    # Dict of species name and RNA editing sites.
    species_editing_info = collections.defaultdict(list)

    # Dict of species name and Reported RNA editing, as we only consider .gff3 files with reported RNA editing events.
    species_name_reported_editing = {}

    # All annotation files (.gff3) of all species.
    all_gff3 = os.listdir(annotation_gff3_dir)

    # Iterate over all annotation files (.gff3)
    for each_gff3 in all_gff3:
        if not os.path.isdir(each_gff3):

            # Seize species name from annotation file (.gff3)
            species_name = each_gff3.split('.')[0] + '.' + each_gff3.split('.')[1]

            # RNA editing sites in a single species.
            editing_sites = 0

            # Keep species fasta and .gff3 files with reported RNA editing events.
            for line in open(annotation_gff3_dir + "\\" + each_gff3, 'r'):

                if re.findall('.*sequence_feature.*', line) and re.findall('.*editing.*', line):

                    editing_sites += 1

            if editing_sites > 0:
                
                species_name_reported_editing[species_name] = "Reported RNA editing events"

                # The array of each coding gene and start position, end position and strand info of CDS.
                cds_start_end_strand = []

                # All RNA editing sites in .gff3, including sites not in CDS.
                candidate_rna_editing_sites = collections.defaultdict(list)

                for line in open(annotation_gff3_dir + "\\" + each_gff3, 'r'):
                    new_line = line.strip()
                    if line[0] == "#" or line[0] == '\n':
                        continue
                    else:
                        cols = new_line.split('\t')
                        gene_annotation = cols[8]

                        # Seize CDS name, start position, end position and strand.
                        if cols[2] == "CDS" and re.findall(".*;gene=(.*?);.*", gene_annotation):

                            gene_name = ''.join(re.findall(".*;gene=(.*?);.*", gene_annotation))

                            # Substitution of synonymous genes.
                            if re.findall('mttB', gene_name, re.IGNORECASE) or re.findall('tatC', gene_name, re.IGNORECASE):
                                reference_gene_name = "orfX"
                                new_gene_name = gene_name

                            elif re.findall('ccmFN1', gene_name, re.IGNORECASE):
                                reference_gene_name = "ccmFn"
                                new_gene_name = gene_name

                            elif re.findall('ccmFN2', gene_name, re.IGNORECASE):
                                reference_gene_name = "ccmFn"
                                new_gene_name = gene_name

                            elif gene_name == "mat-R" or gene_name == "matR":
                                reference_gene_name = "mat-r"
                                new_gene_name = gene_name

                            # Make gene names consistent with those in rice_coding_gene_list.txt
                            elif re.findall('atp6', gene_name, re.IGNORECASE):
                                reference_gene_name = "atp6"
                                new_gene_name = gene_name

                            elif re.findall('atp9', gene_name, re.IGNORECASE):
                                reference_gene_name = "atp9"
                                new_gene_name = gene_name

                            elif re.findall('ccmB', gene_name, re.IGNORECASE):
                                reference_gene_name = "ccmB"
                                new_gene_name = gene_name

                            elif re.findall('ccmC', gene_name, re.IGNORECASE):
                                reference_gene_name = "ccmC"
                                new_gene_name = gene_name

                            elif re.findall('ccmFc', gene_name, re.IGNORECASE):
                                reference_gene_name = "ccmFc"
                                new_gene_name = gene_name

                            elif re.findall('ccmFn', gene_name, re.IGNORECASE):
                                reference_gene_name = "ccmFn"
                                new_gene_name = gene_name

                            elif re.findall('cob', gene_name, re.IGNORECASE):
                                reference_gene_name = "cob"
                                new_gene_name = gene_name

                            elif re.findall('cox2', gene_name, re.IGNORECASE):
                                reference_gene_name = "cox2"
                                new_gene_name = gene_name

                            elif re.findall('nad1', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad1"
                                new_gene_name = gene_name

                            elif re.findall('nad2', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad2"
                                new_gene_name = gene_name

                            elif re.findall('nad4L', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad4L"
                                new_gene_name = gene_name

                            elif re.findall('nad4', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad4"
                                new_gene_name = gene_name

                            elif re.findall('nad5', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad5"
                                new_gene_name = gene_name

                            elif re.findall('nad6', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad6"
                                new_gene_name = gene_name

                            elif re.findall('nad7', gene_name, re.IGNORECASE):
                                reference_gene_name = "nad7"
                                new_gene_name = gene_name

                            elif re.findall('orf183', gene_name, re.IGNORECASE):
                                reference_gene_name = "orf183"
                                new_gene_name = gene_name

                            elif re.findall('orf288', gene_name, re.IGNORECASE):
                                reference_gene_name = "orf288"
                                new_gene_name = gene_name

                            elif re.findall('orfX', gene_name, re.IGNORECASE):
                                reference_gene_name = "orfX"
                                new_gene_name = gene_name

                            elif re.findall('rpl2', gene_name, re.IGNORECASE):
                                reference_gene_name = "rpl2"
                                new_gene_name = gene_name

                            elif re.findall('rps13', gene_name, re.IGNORECASE):
                                reference_gene_name = "rps13"
                                new_gene_name = gene_name

                            elif re.findall('rps19', gene_name, re.IGNORECASE):
                                reference_gene_name = "rps19"
                                new_gene_name = gene_name

                            elif re.findall('rps4', gene_name, re.IGNORECASE):
                                reference_gene_name = "rps4"
                                new_gene_name = gene_name

                            elif re.findall('rps7', gene_name, re.IGNORECASE):
                                reference_gene_name = "rps7"
                                new_gene_name = gene_name

                            elif re.findall('mat-r', gene_name, re.IGNORECASE):
                                reference_gene_name = "mat-r"
                                new_gene_name = gene_name

                            # This part is for rps1.
                            else:
                                reference_gene_name = gene_name
                                new_gene_name = gene_name

                            cds_start_end_strand.append(reference_gene_name + "~" + new_gene_name + "~" + cols[3] + "~" + cols[4] + "~" + cols[6])

                            # Dict of species name and CDS info.
                            species_name_cds_info[species_name].append(reference_gene_name + "~" + new_gene_name + "~" + cols[3] + "~" + cols[4] + "~" + cols[6])

                        # Seize gene name, position of RNA editing sites, strand, and RNA editing type (C-to-U or U-to-C).
                        if cols[2] == "sequence_feature":
                            if re.findall(".*C to U.*", cols[8], re.I) or re.findall(".*c->u.*", cols[8], re.I) or re.findall(".*C-to-U.*", cols[8], re.I) \
                                or re.findall(".*C to T.*", cols[8], re.I) or re.findall(".*c->t.*", cols[8], re.I) or re.findall(".*C-to-T.*", cols[8], re.I):

                                if re.findall(".*gene=.*", gene_annotation.split(";")[-1]):
                                    gene_name = gene_annotation.split(";")[-1].split("=")[-1]
                                    candidate_rna_editing_sites[gene_name].append(cols[3] + "_" + cols[6] + "_C-to-U")

                                else:
                                    gene_name = ''.join(re.findall(".*;gene=(.*?);.*", gene_annotation))
                                    candidate_rna_editing_sites[gene_name].append(cols[3] + "_" + cols[6] + "_C-to-U")

                            if re.findall(".*U to C.*", cols[8], re.I) or re.findall(".*u->c.*", cols[8], re.I) or re.findall(".*U-to-C.*", cols[8], re.I) \
                                or re.findall(".*T to C.*", cols[8], re.I) or re.findall(".*t->c.*", cols[8], re.I) or re.findall(".*T-to-C.*", cols[8], re.I):

                                if re.findall(".*gene=.*", gene_annotation.split(";")[-1]):
                                    gene_name = gene_annotation.split(";")[-1].split("=")[-1]
                                    candidate_rna_editing_sites[gene_name].append(cols[3] + "_" + cols[6] + "_U-to-C")

                                else:
                                    gene_name = ''.join(re.findall(".*;gene=(.*?);.*", gene_annotation))
                                    candidate_rna_editing_sites[gene_name].append(cols[3] + "_" + cols[6] + "_U-to-C")

                # Seize RNA editing sites in CDS.
                for cds_info in cds_start_end_strand:
                    for key, value in candidate_rna_editing_sites.items():
                        for editing_site in value:
                            if key == cds_info.split('~')[1] and int(editing_site.split('_')[0]) >= int(cds_info.split('~')[2]) and \
                                int(editing_site.split('_')[0]) <= int(cds_info.split('~')[3]) and \
                                    editing_site.split('_')[1] == cds_info.split('~')[4]:
                                
                                # Dict of species and editing info.
                                species_editing_info[species_name].append(editing_site)
            else:
                species_name_reported_editing[species_name] = "No reported RNA editing events"

    return species_name_cds_info, species_editing_info, species_name_reported_editing


def MultiOneGenome(genome_fasta_dir, annotation_gff3_dir):
    """
    Convert multi-line genome sequence into one-line.
    """

    # Dict of species name and genome sequence
    species_genome = {}

    # Dict of species name and reported RNA editing events, we only keep genome sequence with reported RNA editing events
    species_reported_rna_editing = CDS_Editing_Info(annotation_gff3_dir)[2]

    # All genome of all species.
    all_genome = os.listdir(genome_fasta_dir)

    # Change each multi-genome into one-line-genome.
    for each_genome in all_genome:
        if not os.path.isdir(each_genome):

            # Seize species name from name of .fasta
            species_name = each_genome.split('.')[0] + '.' + each_genome.split('.')[1]

            # We only keep genome sequence with reported RNA editing events
            if species_reported_rna_editing[species_name] == "Reported RNA editing events":

                # Multi-line in each genome.
                multi_seq_genome = []

                # Merge all genome sequence into one line.
                for line in open(genome_fasta_dir + "\\" + each_genome, 'r'):
                    new_line = line.strip()
                    if line[0] == ">":
                        continue
                    else:
                        multi_seq_genome.append(new_line)

                # Merge multi lines (genome) into one line (genome).
                all_seq = ''.join(multi_seq_genome)

                # Dict of species name and one-line genome.
                species_genome[species_name] = all_seq

    return species_genome


def ProteinBeforeEditing(all_genome_fasta, all_gff3_file, rice_coding_list_file, output_file_dir):
    """
    Translating before RNA editing.
    """

    # Whole genome in one line.
    one_line_genome = MultiOneGenome(all_genome_fasta, all_gff3_file)

    # All CDS info in the whole genome.
    All_CDS_INFO = CDS_Editing_Info(all_gff3_file)[0]

    # Build output files name of genes in rice coding list
    for line in open(rice_coding_list_file, 'r'):
        new_line = line.strip()

        # Output file: one gene including multi-species protein before RNA editing.
        protein_before_editing = open(output_file_dir + "\\" + new_line + ".fasta", 'w')

        # Dict of species and all genes
        species_gene = collections.defaultdict(list)

        for key, value in All_CDS_INFO.items():

            for cds in value:
                species_gene[key].append(cds.split('~')[0])

            if new_line in species_gene[key]:

                # The array of all CDS info (position and strand) per gene.
                cds_info_per_gene = collections.defaultdict(list)

                # All CDS in each gene.
                for each_cds_info in value:

                    if each_cds_info.split('~')[0] == new_line:
                        cds_info_per_gene[each_cds_info.split('~')[0] + '~' + each_cds_info.split('~')[1]].append(each_cds_info.split('~')[2] + "~" + each_cds_info.split('~')[3] + "~" + each_cds_info.split('~')[4])

                # All CDS in one gene.
                all_cds_per_gene = []

                # Translating corresponding gene.
                for gene_name, all_cds_info in cds_info_per_gene.items():

                    for each_cds_info in all_cds_info:

                        cds_start_pos = int(each_cds_info.split('~')[0])
                        cds_end_pos = int(each_cds_info.split('~')[1])

                        if each_cds_info.split('~')[2] == '+':
                            cds_seq = one_line_genome[key][cds_start_pos-1:cds_end_pos]
                            all_cds_per_gene.append(cds_seq)

                        if each_cds_info.split('~')[2] == '-':
                            cds_seq = one_line_genome[key][cds_start_pos-1:cds_end_pos]
                            cds_seq_class = Seq(cds_seq)
                            cds_real = cds_seq_class.reverse_complement()
                            all_cds_per_gene.append(str(cds_real))

                    # Merge all cds
                    merged_cds = ''.join(all_cds_per_gene)

                    # The length of CDS must be an integral multiple of 3.
                    if len(merged_cds) % 3 == 0:

                        # Write gene name and protein sequence after RNA editing.
                        protein_before_editing.write(">" + key + "~" + gene_name + "\n")
                        protein_before_editing.write(merged_cds + "\n")


def ProteinAfterEditing(all_genome_fasta, all_gff3_file, rice_coding_list_file, output_file_dir, paper_editing_list):
    """
    Translating after RNA editing.
    """

    # Except RNA editing sites in annotation file (.gff3), we add RNA editing sites from published papers.
    # Format: e.g. Oryza_sativa (species)   atp6 (edited gene)    5 (position of amino acid)   H (amino acid before editing)   Y (amino acid after editing)

    # Dict of species-gene and amino acid position-edited amino acid.
    species_gene_edited_aa = collections.defaultdict(list)

    for line in open(paper_editing_list, 'r'):

        new_line = line.strip()
        cols = new_line.split('\t')

        # Dict of accession-id-species-gene and position of edited base.
        species_gene_edited_aa[cols[0] + "~" + cols[1]].append(cols[3])

    # Whole genome in one line.
    one_line_genome = MultiOneGenome(all_genome_fasta, all_gff3_file)

    # All CDS info in the whole genome.
    All_CDS_INFO = CDS_Editing_Info(all_gff3_file)[0]
    
    # Build output files name of genes in rice coding list
    for line in open(rice_coding_list_file, 'r'):
        
        new_line = line.strip()

        # Output file: one gene including multi-species protein after RNA editing.
        protein_after_editing = open(output_file_dir + "\\" + new_line + ".fasta", 'w')

        # Dict of species and all genes
        species_gene = collections.defaultdict(list)

        for key, value in All_CDS_INFO.items():

            # All RNA editing sites and strand in CDS in the whole genome per species.
            editing_per_species = CDS_Editing_Info(all_gff3_file)[1][key]

            # Get the whole genome sequence after RNA editing.
            for j in editing_per_species:

                # RNA editing (C-to-U) on plus strand.
                if j.split('_')[1] == '+' and j.split('_')[2] == "C-to-U" and one_line_genome[key][int(j.split('_')[0])-1] == 'C':
                    one_line_genome[key] = one_line_genome[key][:int(j.split('_')[0])-1] + 'T' + one_line_genome[key][int(j.split('_')[0]):]

                # RNA editing (U-to-C) on plus strand.
                if j.split('_')[1] == '+' and j.split('_')[2] == "U-to-C" and one_line_genome[key][int(j.split('_')[0])-1] == 'T':
                    one_line_genome[key] = one_line_genome[key][:int(j.split('_')[0])-1] + 'C' + one_line_genome[key][int(j.split('_')[0]):]

                # RNA editing (C-to-U) on minus strand.
                if j.split('_')[1] == '-' and j.split('_')[2] == "C-to-U" and one_line_genome[key][int(j.split('_')[0])-1] == 'G':
                    one_line_genome[key] = one_line_genome[key][:int(j.split('_')[0])-1] + 'A' + one_line_genome[key][int(j.split('_')[0]):]

                # RNA editing (C-to-U) on minus strand.
                if j.split('_')[1] == '-' and j.split('_')[2] == "U-to-C" and one_line_genome[key][int(j.split('_')[0])-1] == 'A':
                    one_line_genome[key] = one_line_genome[key][:int(j.split('_')[0])-1] + 'G' + one_line_genome[key][int(j.split('_')[0]):]

            for cds in value:
                species_gene[key].append(cds.split('~')[0])

            if new_line in species_gene[key]:

                # The array of all CDS info (position and strand) per gene.
                cds_info_per_gene = collections.defaultdict(list)

                # All CDS in each gene.
                for each_cds_info in value:

                    if each_cds_info.split('~')[0] == new_line:
                        cds_info_per_gene[each_cds_info.split('~')[0] + '~' + each_cds_info.split('~')[1]].append(each_cds_info.split('~')[2] + "~" + each_cds_info.split('~')[3] + "~" + each_cds_info.split('~')[4])

                # All CDS in one gene.
                all_cds_per_gene = []

                # Translating corresponding gene.
                for gene_name, all_cds_info in cds_info_per_gene.items():
                    for each_cds_info in all_cds_info:

                        cds_start_pos = int(each_cds_info.split('~')[0])
                        cds_end_pos = int(each_cds_info.split('~')[1])

                        if each_cds_info.split('~')[2] == '+':
                            cds_seq = one_line_genome[key][cds_start_pos-1:cds_end_pos]
                            all_cds_per_gene.append(cds_seq)

                        if each_cds_info.split('~')[2] == '-':
                            cds_seq = one_line_genome[key][cds_start_pos-1:cds_end_pos]
                            cds_seq_class = Seq(cds_seq)
                            cds_real = cds_seq_class.reverse_complement()
                            all_cds_per_gene.append(str(cds_real))

                    # Merge all cds
                    merged_cds = ''.join(all_cds_per_gene)

                    # The length of CDS must be an integral multiple of 3.
                    if len(merged_cds) % 3 == 0:

                        # Consider RNA editing sites from published papers.
                        for position_aa in species_gene_edited_aa[key + "~" + new_line]:
                            merged_cds = merged_cds[:int(position_aa)-1] + 'T' + merged_cds[int(position_aa):]
                            
                        # Write gene name and protein sequence after RNA editing.
                        protein_after_editing.write(">" + key + "~" + gene_name + "\n")
                        protein_after_editing.write(merged_cds + "\n")


def main():
    """ 
    Pass parameters and enter the matrix generating process.
    """

    # Set parameters
    parser = argparse.ArgumentParser(description = "Please enter the proper dirs and files.")

    # The dirs and files.
    parser.add_argument('--MitFasta', '-mf', type = str, help = "The dir of .fasta files of mitochondrial genome.")
    parser.add_argument('--MitGff', '-mg', type = str, help = "The dir of .gff files of mitochondrial genome.")
    parser.add_argument('--MitCodingGene', '-mcg', type = str, help = "The file containing rice mitochondrial coding genes harbouring RNA editing sites.")
    parser.add_argument('--OutBeforeEditing', '-obe', type = str, help = "The dir of output files of protein sequence before RNA editing.")
    parser.add_argument('--OutAfterEditing', '-oae', type = str, help = "The dir of output files of protein sequence after RNA editing.")
    parser.add_argument('--RiceEditingEvents', '-ree', type = str, help = "The file containing rice editing events we identified, or editing events in other species curated from literatures.")

    args = parser.parse_args()

    mitFasta = args.MitFasta
    mitGff = args.MitGff
    mitCodingGene = args.MitCodingGene
    outBeforeEditing = args.OutBeforeEditing
    outAfterEditing = args.OutAfterEditing
    riceEditingEvents = args.RiceEditingEvents

    ProteinBeforeEditing(mitFasta, mitGff, mitCodingGene, outBeforeEditing)
    ProteinAfterEditing(mitFasta, mitGff, mitCodingGene, outAfterEditing, riceEditingEvents)

if __name__ == '__main__':
    main()