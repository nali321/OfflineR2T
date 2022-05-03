import os
import gtdbtk_extractor
import summary
import gtotree_text
import argparse
import switch_maker

parser = argparse.ArgumentParser()

parser.add_argument("-r1", "--read1", type=str,
                    help="Filepath of first half of read pair")

parser.add_argument("-r2", "--read2", type=str,
                    help="Filepath of second half of read pair")

parser.add_argument("-s", "--summary", type=str,
                    help="Filepath to assembly_summary.txt")

parser.add_argument("-r", "--references", type=str,
                    help="Filepath to references folder")

parser.add_argument("-l", "--leaves", type=int,
                    help="Maximum number of leaves for each tree")

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go")

#make assembly summary filepath and references filepath flags
#make cpus flag

args = parser.parse_args()

#the input the user has to provide are the paired reads
read1 = args.read1
read2 = args.read2

#get filepath to assembly_summary
summary_path = args.summary

#get filepath to references folder
ref_path = args.references

#sets max number of leaves
if args.leaves is None:
    max = 10

else:
    max = args.leaves

#if output directory already exists, send error message
path = args.outdir
try:
    os.mkdir(path)
except OSError as error:
    print(error)

#makes switches directory
switches = os.path.join(path, "switches").replace("\\", "/")
os.mkdir(switches)

#PIPELINE: trimmomatic > spades.py > quast.py > prokka > gtdbtk > gtotree

#trimmomatic
trim = os.path.join(path, "trim").replace("\\", "/")
os.mkdir(trim)
os.system(f"trimmomatic PE -phred33 {read1} {read2} {trim}/forward_paired.fq.gz {trim}/forward_unpaired.fq.gz {trim}/reverse_paired.fq.gz {trim}/reverse_unpaired.fq.gz ILLUMINACLIP://mmfs1//groups//HPC-Marshall//miniconda3//pkgs//trimmomatic-0.39-hdfd78af_2//share//trimmomatic-0.39-2//adapters//NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70")

#spades
spades = os.path.join(path, "spades").replace("\\", "/")
os.mkdir(spades)
os.system(f"spades.py -o {spades} -1 {trim}/forward_paired.fq.gz -2 {trim}/reverse_paired.fq.gz -s {trim}/forward_unpaired.fq.gz -m 128 -k 33,55,77,99 --careful -t 16")

#quast
quast = os.path.join(path, "quast").replace("\\", "/")
os.mkdir(quast)
os.system(f"/mmfs1/groups/HPC-Marshall/build/quast-5.1.0rc1/quast.py {spades}/contigs.fasta -o {quast}")

#prokka
prokka = os.path.join(path, "prokka").replace("\\", "/")
os.mkdir(prokka)
switch_maker.prokka_switch(switches, prokka, spades)
os.system(f"bash {switches}/prokka.sh")

#gtdbtk
#create a folder with contigs.fna in it
genome_dir = os.path.join(path, "genome_dir").replace("\\", "/")
os.mkdir(genome_dir)
os.system(f"cp {spades}/contigs.fasta {genome_dir}/contigs.fna")

gtdbtk = os.path.join(path, "gtdbtk").replace("\\", "/")
os.mkdir(gtdbtk)
switch_maker.gtdbtk_switch(switches, gtdbtk, genome_dir)
os.system(f"bash {switches}/gtdbtk.sh")

#gtotree
#call gtdbtk_extractor
filepath = f"{gtdbtk}/classify/gtdbtk.bac120.summary.tsv"

#get the species name and h-flag name
species_name, h_flag, species_acc = gtdbtk_extractor.flag_extractor(filepath)

#create the dictionary structure for species name lookup
name_structure, acc_structure = summary.structure(summary_path)

#obtain leaves by looking up species name
leaves = summary.leaves(name_structure, acc_structure, max, species_name)

#take accessions from leaves and feed it into other_related to avoid outgroup
#accession matching one of the accessions of the leaves
leaf_accessions = []
for leaf in leaves:
    leaf_accessions.append(leaf[0])

#get the other related species name
other_related = gtdbtk_extractor.other_related(filepath, species_name, species_acc, leaf_accessions)

#obtain outgroup by searching its accession
outgroup = summary.outgroup(name_structure, acc_structure, other_related[1], other_related[0])

#combine outgroup tuple with leaves to get final list
leaves.append(outgroup)

#create the text files needed for gtotree
gtotree_text_files = os.path.join(path, "gtotree_text_files").replace("\\", "/")
os.mkdir(gtotree_text_files)
gtotree_text.fasta_files(gtotree_text_files, leaves)
gtotree_text.map_id(gtotree_text_files, leaves)

#EXTRACT FILEPATH NAME FROM COLUMN 20: "FTP_PATH", AND ADD "_genomic.fna.gz" AT THE END TO FIND
#IT IN REFERENCES: /mmfs1/groups/HPC-Marshall/database/genbank_3-2022/references
#need to rename .fa.gz file as accession.fa.gz file only

#move accessions and rename them in the main directory
for x in leaves:
    acc = x[0] + ".fa.gz"
    ftp = x[2] + "_genomic.fna.gz"
    ftp_path = os.path.join(ref_path, ftp).replace("\\", "/")
    acc_path = os.path.join(path, acc).replace("\\", "/")
    os.system(f"cp {ftp_path} {acc_path}")

#move contigs.fasta from ./spades to main dir, renamed contigs.fa
os.system(f"cp {spades}/contigs.fasta {path}/contigs.fa")

#create the gtotree switch and gtotree output folder
gtotree = os.path.join(path, "gtotree").replace("\\", "/")

#use gtotree switch
switch_maker.gtotree_switch(path, gtotree_text_files, switches, gtotree, h_flag)
os.system(f"bash {switches}/gtotree.sh")

#TREE MADE