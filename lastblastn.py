import os
import subprocess
import sys

blast_db = "/localdisk01/backup/home/common/blast/blast_db_teste/core_nt"

def run_blastn(genome_file, blast_db):
    """Executa BLASTN no genoma especificado."""

    outputs_folder = "/data/nfs/home/viruses/brenno/a"
    
    # Pasta de saída para blastn
    blastn_folder = os.path.join(outputs_folder, "blastn")
    os.makedirs(blastn_folder, exist_ok=True)

    genome_name = os.path.basename(genome_file).replace(".fa", "").replace(".fasta", "")
    output_file = os.path.join(blastn_folder, f"{genome_name}_blastnResults.tsv")

    # Comando para executar BLASTN
    blastn_command = [
        "blastn",
        "-query", genome_file,
        "-db", blast_db,
        "-out", output_file,
        "-evalue", "0.001",
        "-outfmt", "6 qseqid sseqid qlen slen pident evalue bitscore stitle",
        "-max_target_seqs", "10",
        "-num_threads", "4"
    ]

    subprocess.run(blastn_command, check=True)
    print(f"BLASTN concluído para {genome_name}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python run_blastn.py <genome_file.fa>")
        sys.exit(1)

    genome_file = sys.argv[1]
    run_blastn(genome_file, blast_db)
