#make prokka switch
def prokka_switch(switches, outdir, spades):
    f = open(f"{switches}/prokka.sh", "w+")

    #shebang statement
    f.write("#!/bin/bash")
    f.write("\n")

    #conda profile
    f.write("source /mmfs1/groups/HPC-Marshall/miniconda3/etc/profile.d/conda.sh")
    f.write("\n")

    #activate tool
    f.write("conda activate prokka")
    f.write("\n")

    #prokka code
    f.write(f"prokka --centre X --force --locustag MB_test121 --outdir {outdir} --prefix genome --gffver 3 --cpus 8 {spades}/contigs.fasta")
    f.write("\n")

    f.write("conda deactivate")

def gtdbtk_switch(switches, outdir, genome_dir):
    f = open(f"{switches}/gtdbtk.sh", "w+")

    #shebang statement
    f.write("#!/bin/bash")
    f.write("\n")

    #conda profile
    f.write("source /mmfs1/groups/HPC-Marshall/miniconda3/etc/profile.d/conda.sh")
    f.write("\n")

    #activate tool
    f.write("conda activate gtdbtk")
    f.write("\n")

    #gtdbtk code
    f.write(f"gtdbtk classify_wf --genome_dir {genome_dir} --out_dir {outdir} --cpus 32")
    f.write("\n")

    f.write("conda deactivate")

#make gtotree switch
#path: to text files
def gtotree_switch(path, dir1, switches, outdir, h_flag):
    f = open(f"{switches}/gtotree.sh", "w+")

    #shebang statement
    f.write("#!/bin/bash")
    f.write("\n")

    #conda profile
    f.write("source /mmfs1/groups/HPC-Marshall/miniconda3/etc/profile.d/conda.sh")
    f.write("\n")

    #activate tool
    f.write("conda activate gtotree")
    f.write("\n")

    #move to directory
    f.write(f"cd {path}")
    f.write("\n")

    #gtotree code
    f.write(f'GToTree -f {dir1}/fasta_files.txt -H {h_flag} -t -L Species,Strain -m {dir1}/map_id.txt -T IQ-TREE -j 16 -o {outdir}')
    f.write("\n")

    f.write("conda deactivate")