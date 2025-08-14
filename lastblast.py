import os
import subprocess

fasta_file = "/data/nfs/home/viruses/brenno/a/GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_genomic_blastxResultssta.fa"
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
    run_diamond_blastx(fasta_file, db_path)



import subprocess
import sys
import os

# Verifica se o nome do genoma foi fornecido
if len(sys.argv) != 2:
    print("Uso: python run_blastn.py <arquivo_genoma.fasta>")
    sys.exit(1)

# Parâmetros
genome_file = sys.argv[1]
genome_name = os.path.splitext(os.path.basename(genome_file))[0]

# Caminhos fixos
blast_db = "/home/common/blast/blast_db_teste/core_nt"
output_file = f"{genome_name}_blastn.tabular"

# Comando BLASTN
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

# Executa o comando
print(f"Rodando BLASTn para: {genome_file}")
try:
    subprocess.run(blastn_command, check=True)
    print(f"Arquivo de saída gerado: {output_file}")
except subprocess.CalledProcessError as e:
    print("Erro ao executar o BLASTn:", e)
