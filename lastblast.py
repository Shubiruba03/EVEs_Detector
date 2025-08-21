import os
import subprocess
import sys

db_path = "/localdisk01/backup/home/common/nr/nr.dmnd"

def run_diamond_blastx(fasta_file, db_path):
    """Executa DIAMOND BLASTX e filtra os resultados."""

    outputs_folder = "/data/nfs/home/viruses/brenno/a"
    
    # Definir a pasta blastx dentro de Outputs
    blastx_folder = os.path.join(outputs_folder, "blastx")
    os.makedirs(blastx_folder, exist_ok=True)

    genome_name = os.path.basename(fasta_file).replace(".fa", "")  
    diamond_output = os.path.join(blastx_folder, f"{genome_name}_blastResults.tsv")  

    # Comando para executar DIAMOND BLASTX
    cmd = [
        "diamond", "blastx",
        "-d", db_path,
        "-q", fasta_file,
        "--outfmt", "6", "qseqid", "sseqid", "qlen", "slen", "pident", "length", "mismatch", 
        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle", "qtitle", 
        "full_qseq", "--max-target-seqs", "5",
        "-o", diamond_output
    ]
    subprocess.run(cmd, check=True)
    print(f"DIAMOND blastx concluído para {genome_name}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python run_diamond.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    run_diamond_blastx(fasta_file, db_path)
