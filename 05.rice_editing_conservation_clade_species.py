# -*- encoding: utf-8 -*-
'''
@File    :   03.rice_editing_conservation_percent.py
@Time    :   2023/07/13 09:10:03
@Author  :   Ming Chen
@Contact :   chenm@big.ac.cn
'''

# Put the import lib here.
import sys, os, re
import collections

"""
Usage
sys.argv[1]: List of rice RNA editing sites.
sys.argv[2]: Dir of multi-alignment results.
sys.argv[3]: Output ratio amino acids produced by rice RNA editing among land plants.
sys.argv[4]: Dir or plant clades
"""

# Dict of gene and position-edited amino acid (list) from rice editing sites.
gene_position_aa = collections.defaultdict(list)

# List of rice RNA editing sites
for line in open(sys.argv[1], 'r'):
    new_line = line.strip()
    cols = new_line.split('\t')

    # Gene name, position of amino acid and amino acid produced by RNA editing.
    gene_position_aa[cols[1]].append(cols[2] + "~" + cols[4])

# All mitochondrial proteins multi_alignment output.
all_multi_alignment = os.listdir(sys.argv[2])

# List of amino acids and Stop codon
aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

# Amino acid produced by rice RNA editing conservation percent output.
aa_conservation = open(sys.argv[3], 'w')

# Output of conservation of RNA editing and the corresponding amino acid.
aa_conservation.write("GeneName\tAAPosition\tAminoAcid\tEditingAAPercentage\n")

# All species clades of plants.
all_clades = os.listdir(sys.argv[4])

# Dict of clades and species name.
clades_species = collections.defaultdict(list)

for clade in all_clades:

    # clades: seed_plants, lycophytes, ferns, mosses, hornworts and liverworts.
    clade_name = clade.split('.')[-1]

    for accession_species in os.listdir(sys.argv[4] + "\\" + clade):
        
        species_name = accession_species.split('-')[-1].split('.')[0]

        # Dict of clade and species name.
        clades_species[clade_name].append(species_name)

# Clades of species in the output of multi-alignment.
clade_species_multi_alignment = collections.defaultdict(list)

# Iterate over all multi alignment output.
for each_multi_alignment in all_multi_alignment:

    if not os.path.isdir(each_multi_alignment):

        # Seize gene name from multi alignment output.
        gene_name = each_multi_alignment.split('.')[0]

        # Dict of species and multi-alignment output.
        species_sequence = {}

        # Species names in multi-alignment output.
        species = []

        # Multi-alignment output.
        for line in open(sys.argv[2] + "\\" + each_multi_alignment, 'r'):

            new_line = line.strip()
            
            if new_line[0] == '\n':
                continue

            elif new_line[0] == ">" and re.findall(".*(Sequence).*", new_line):
                accession_id_species_name = "Sequence_1"
                species_sequence[accession_id_species_name] = ""

            elif new_line[0] == ">":
                accession_id_species_name = '_'.join(new_line.split('>')[1].split(' ')[:-2])
                species_sequence[accession_id_species_name] = ""

                # Species name.
                species_name = '_'.join(new_line.split(' ')[:-2]).split('-')[-1]
                species.append(species_name)

                if species_name in clades_species['seed_plants']:
                    clade_species_multi_alignment['seed_plants'].append(species_name)

                if species_name in clades_species['lycophytes']:
                    clade_species_multi_alignment['lycophytes'].append(species_name)

                if species_name in clades_species['ferns']:
                    clade_species_multi_alignment['ferns'].append(species_name)

                if species_name in clades_species['mosses']:
                    clade_species_multi_alignment['mosses'].append(species_name)

                if species_name in clades_species['hornworts']:
                    clade_species_multi_alignment['hornworts'].append(species_name)

                if species_name in clades_species['liverworts']:
                    clade_species_multi_alignment['liverworts'].append(species_name)

            else:
                species_sequence[accession_id_species_name] = new_line

        species_sequence.pop("Sequence_1", None)

        # Number of multi-alignment species
        species_number = len(list(set(species)))

        # Dict of clade and number of species.
        clade_num_species = {}
        
        clade_num_species['seed_plants'] = len(list(set(clade_species_multi_alignment['seed_plants'])))
        clade_num_species['lycophytes'] = len(list(set(clade_species_multi_alignment['lycophytes'])))
        clade_num_species['ferns'] = len(list(set(clade_species_multi_alignment['ferns'])))
        clade_num_species['mosses'] = len(list(set(clade_species_multi_alignment['mosses'])))
        clade_num_species['hornworts'] = len(list(set(clade_species_multi_alignment['hornworts'])))
        clade_num_species['liverworts'] = len(list(set(clade_species_multi_alignment['liverworts'])))

        print("seed_plants: " + str(len(list(set(clade_species_multi_alignment['seed_plants'])))))
        print("lycophytes: " + str(len(list(set(clade_species_multi_alignment['lycophytes'])))))
        print("ferns: " + str(len(list(set(clade_species_multi_alignment['ferns'])))))
        print("mosses: " + str(len(list(set(clade_species_multi_alignment['mosses'])))))
        print("hornworts: " + str(len(list(set(clade_species_multi_alignment['hornworts'])))))
        print("liverworts: " + str(len(list(set(clade_species_multi_alignment['liverworts'])))))

        print("\n")

        # amino acid position in rice RNA editing sites list.
        amino_acid_position = 0

        # amino acid index in rice RNA editing sites list.
        amino_acid_index = 0

        # Calculate conservation percent of rice RNA editing site.
        for index, base in enumerate(species_sequence["NC_011033.1-Oryza_sativa"]):

            # amino acid position 
            if base in aa_list:

                amino_acid_position += 1

                for gene, position_aa in gene_position_aa.items():

                    # Because of multi accessions of a gene in a species, it is necessary to construct a dict of species name and amino acids at the same RNA editing site
                    species_amino_acid = collections.defaultdict(list)

                    if gene == gene_name and str(amino_acid_position)+'~'+base in position_aa:

                        amino_acid_index = index

                        for accession_id_species_name, protein_sequence in species_sequence.items():

                            if accession_id_species_name == "NC_011033.1-Oryza_sativa":
                                continue

                            else:
                                # Seize amino acid in other species corresponding to RNA editing events in rice.
                                # species_amino_acid[accession_id_species_name.split('-')[-1].split('_')[0]].append(protein_sequence[amino_acid_index])

                                species_amino_acid[accession_id_species_name.split('-')[-1]].append(protein_sequence[amino_acid_index])

                        # The RNA editing site and the amino acid caused by RNA editing.
                        editing_amino_acid_num = 0

                        # The amino acid on the position of RNA editing
                        amino_acid_num = 0

                        # Dict of plant clade and number of species containing amino acid corresponding to RNA editing in rice
                        plant_clade_num_species_editing = collections.defaultdict(int)

                        plant_clade_num_species_editing['seed_plants'] = 0
                        plant_clade_num_species_editing['lycophytes'] = 0
                        plant_clade_num_species_editing['ferns'] = 0
                        plant_clade_num_species_editing['mosses'] = 0
                        plant_clade_num_species_editing['hornworts'] = 0
                        plant_clade_num_species_editing['liverworts'] = 0

                        for species_name, amino_acid_list in species_amino_acid.items():

                            if species_name in clades_species['seed_plants'] and base in amino_acid_list:
                                plant_clade_num_species_editing['seed_plants'] += 1

                            if species_name in clades_species['lycophytes'] and base in amino_acid_list:
                                plant_clade_num_species_editing['lycophytes'] += 1

                            if species_name in clades_species['ferns'] and base in amino_acid_list:
                                plant_clade_num_species_editing['ferns'] += 1

                            if species_name in clades_species['mosses'] and base in amino_acid_list:
                                plant_clade_num_species_editing['mosses'] += 1

                            if species_name in clades_species['hornworts'] and base in amino_acid_list:
                                plant_clade_num_species_editing['hornworts'] += 1

                            if species_name in clades_species['liverworts'] and base in amino_acid_list:
                                plant_clade_num_species_editing['liverworts'] += 1

                            # convert set of amino acids into character string.
                            all_merged_amino_acid = ''.join(set(amino_acid_list))

                            amino_acid_stop_codon = re.compile(r'[A-Za-z*]', re.S)

                            if re.findall(amino_acid_stop_codon, all_merged_amino_acid):
                                amino_acid_num += 1

                        # print(gene_name + "_" + str(amino_acid_position) +  ": ")
                        # print("seed_plants: " + str(plant_clade_num_species_editing['seed_plants']))
                        # print("lycophytes: " + str(plant_clade_num_species_editing['lycophytes']))
                        # print("ferns: " + str(plant_clade_num_species_editing['ferns']))
                        # print("mosses: " + str(plant_clade_num_species_editing['mosses']))
                        # print("hornworts: " + str(plant_clade_num_species_editing['hornworts']))
                        # print("liverworts: " + str(plant_clade_num_species_editing['liverworts']) + "\n")

                        num_clades_containing_species = 0

                        editing_amino_acid_percentage = 0

                        # Percentage of amino acid caused by RNA editing at per editing site.
                        for i in ['seed_plants', 'lycophytes', 'ferns', 'mosses', 'hornworts', 'liverworts']:

                            if plant_clade_num_species_editing[i] > 0:

                                editing_amino_acid_percentage += plant_clade_num_species_editing[i] / clade_num_species[i]

                                num_clades_containing_species += 1

                        if editing_amino_acid_percentage > 0:

                            averaged_editing_amino_acid_percentage =  editing_amino_acid_percentage / num_clades_containing_species

                        # Percentage of amino acid at per editing site.
                        amino_acid_percentage = amino_acid_num / species_number

                        # Threshold of sequency conservation.
                        if amino_acid_percentage >= 0.8:

                            # Output of conservation of RNA editing and the corresponding amino acid.
                            aa_conservation.write(gene_name + '\t' + str(amino_acid_position) + '\t' + base + '\t' + str(averaged_editing_amino_acid_percentage) + '\n')

aa_conservation.close()