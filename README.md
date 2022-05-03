# OfflineR2T
OfflineR2T (Offline Read to Tree) is a tool for the Marshall Lab to automate the read pair to phylogenetic tree pipeline through slurm.

Inputs:

-r1: read pair one

-r2: read pair two

-s: filepath to assembly_summary.txt

-r: filepath to references folder

-l: max number of leaves on the tree (Optional flag - default is 10)

-o: output directory 

Example (-s and -r flags are actual filepaths): 

    python offliner2t.py -r1 [read pair 1] -r2 [read pair 2] -s /mmfs1/groups/HPC-Marshall/database/genbank_3-2022/assembly_summary.txt -r /mmfs1/groups/HPC-Marshall/database/genbank_3-2022/references -l 10 -o ./output

Pipeline Details (in order): 
1. Trimmomatic - hardcoded: ILLUMINACLIP://mmfs1//groups//HPC-Marshall//miniconda3//pkgs//trimmomatic-0.39-hdfd78af_2//share//trimmomatic-0.39-2//adapters//NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
2. Spades - hardcoded: -m 128 | -k 33,55,77,99 | --careful | -t 16
3. Quast - harcoded: /mmfs1/groups/HPC-Marshall/build/quast-5.1.0rc1/quast.py
4. Prokka - hardcoded: --centre X | --force | --prefix genome | --gffver 3 | --cpus 8
5. Gtdbtk - hardcoded: classify_wf | --cpus 32
6. assembly_summary.txt & references folder usage
7. GToTree - harcoded: -t -L Species,Strain | -T IQ-TREE | -j 16

Script Details:

offliner2t.py - Main python script that runs the pipeline; manages folder creation and data movement. Calls upon other helper scripts below:

switch_maker.py - When automating the tool usage of things with conda environments that need to be activated beforehand, there is a workaround needed to "conda activate" a tool properly on the server. You cannot just os.system("conda activate prokka") with Python because the current shell has no conda environment access. The solution is that switch_maker.py creates a bash script with requisite flags for Prokka, Gtdbtk, and GToTree. The conda profile filepath is also hardcoded in each of these bash scripts: "source /mmfs1/groups/HPC-Marshall/miniconda3/etc/profile.d/conda.sh". This is needed to allow the current shell, in either tmux or slurm, the access to activate any conda environment. Then, offliner2t.py runs the bash script through os.system("bash [switch file]") to "switch" on the environment and run the tool properly.

gtdbtk_extractor.py - Extracts the 2nd and 16th columns of the gtdbtk.bac120.summary.tsv file output from gtdbtk. The 2nd column contains the species name of the organism, which is needed to search for leaves of the tree. It also contains the h-flag needed for GToTree. The 16th column finds the outgroup for the tree, by taking the very first name in the column that does not match the species name from the 2nd column. Hardcoded:
    - If the phylum is proteobacteria, take the class instead for the h-flag
    - If the h-flag is "Actinobacteriota", rename it to "Actinobacteria"

summary.py -

summary.structure - Reads in each line of assembly_summary.txt, and creates a dictionary file structure where the name of a species can be used to search for every single accession under that name stored as a list of lists. Ex: ["Aeromonas popoffii"] = [[information about accession 1], [information about accession 2], etc...]. Create a second dictionary structure where dictionary entries are accession numbers, Ex: [GCA_001006765.1] = [information about accession] in order to search for outgroup information.
    
summary.leaves - Uses the created dictionary structure to take the species name from gtdbtk_extractor.py and searches it. Takes no more than the max amount of leaves defind by user. 
        NOTE: Opportunity to create code to prioritize accessions with a complete genome assembly type instead of just take the top x accessions
        
summary.outgroup - Searches the second dictionary structure with the outgroup's accession, and takes the first outgroup that matches the accession taken from the 16th column

Other:

Check tree success -

    for i in {1..x}; do
        file="run_1/read_pair_"$i"/gtotree/gtotree.tre"
        if [[ -f "$file" ]]; then
            echo "Tree exists for read pair "$i"" >> ~/tree_report.txt
        else
            echo "NO TREE FOR "$i"" >> ~/tree_report.txt
        fi
    done
