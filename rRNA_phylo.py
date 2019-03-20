from Bio import SearchIO, SeqIO, Entrez, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Tkinter import *
from tkFileDialog import askopenfilenames
from Bio.Align.Applications import MuscleCommandline
import operator
import subprocess

def phylo_tree(muscle_align, output_file):
    result = open(output_file, "w") 
    result.write('\nPhylogenetic Tree\n***********************\n')
    align = AlignIO.read(muscle_align, 'clustal')
    
    # Generate distance matrix
    calc = DistanceCalculator('identity')
    dist_matrix = calc.get_distance(align)
    
    # Generate tree
    makeTree = DistanceTreeConstructor()
    tree = makeTree.upgma(dist_matrix)
    Phylo.draw_ascii(tree, file=result)

    # Format output
    result.write('Distance Matrix\n***********************\n' + str(dist_matrix) + '\n')
    result.close() 
    return str(result)

def muscle_align(cat_file, align_file):
    muscle_exe = r"muscle.exe"
    muscle_cline = MuscleCommandline(muscle_exe, input=cat_file, out=align_file)
    muscle_cline()
    return str(align_file)

# Ask for Fasta files to align
root = Tk().withdraw()
Input = askopenfilenames(filetypes = [("FASTA files", "*.fasta")])
file_names = ''

for val in Input:
    file_names += str(val) + ' '

# Concatonate Input fasta files in bash
subprocess.call('cat {} > cat_seqs.fasta'.format(file_names), shell=True)

# Write muscle alignment to output file
clust_file = muscle_align('cat_seqs.fasta', 'muscle_out.clw')

# Write phylogeny tree to output file    
phylo_tree('muscle_out.clw', 'out.txt')