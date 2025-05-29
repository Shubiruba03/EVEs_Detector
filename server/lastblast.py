import os
import subprocess

fasta_file = "/home/brenno/a/GCA_910591515.1_iyCerRyby1.1_genomic_Cicadasta.fa"
db_path = "/home/common/nr/nr.dmnd"

def run_diamond_blastx(fasta_file, db_path):
    """Executa DIAMOND BLASTX e filtra os resultados."""

    outputs_folder = "/home/brenno/a"
    
    # Definir a pasta Blastn dentro de Outputs
    blastn_folder = os.path.join(outputs_folder, "BlastN")
    os.makedirs(blastn_folder, exist_ok=True)
    genome_name = os.path.basename(fasta_file).replace(".fa", "")  
    diamond_output = os.path.join(blastn_folder, f"{genome_name}_blastResults.tsv")  

    # Comando para executar DIAMOND BLASTX
    cmd = [
        "diamond", "blastx",
        "-d", db_path,
        "-q", fasta_file,
        "--outfmt", "6", "qseqid", "sseqid", "qlen", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle", "qtitle", "full_qseq",
        "--max-target-seqs", "5",
        "-o", diamond_output
    ]
    subprocess.run(cmd, check=True)
    print(f"DIAMOND blastx concluído para {genome_name}")


if __name__ == "__main__":
    run_diamond_blastx(fasta_file, db_path)