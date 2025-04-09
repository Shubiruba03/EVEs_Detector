from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import csv
import time
import random
from concurrent.futures import ThreadPoolExecutor, as_completed

# Função para calcular a cobertura
def calculate_coverage(query_length, alignment_length):
    return (alignment_length / query_length) * 100

# Função para executar o BLAST e processar os resultados
def run_blast(seq_record, program, database):
    query_id = seq_record.id
    query_sequence = str(seq_record.seq)
    query_length = len(seq_record.seq)

    print(f"Running {program} for: {query_id}")

    attempts = 3  # Número máximo de tentativas
    for attempt in range(attempts):
        try:
            # Rodando o BLAST online
            result_handle = NCBIWWW.qblast(program, database, query_sequence)
            blast_records = NCBIXML.parse(result_handle)

            results = []
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        per_identity = (hsp.identities / hsp.align_length) * 100
                        coverage = calculate_coverage(query_length, hsp.align_length)
                        results.append([
                            query_id,
                            alignment.hit_id,
                            alignment.title,
                            hsp.expect,
                            per_identity,
                            coverage,
                            alignment.accession,
                            query_sequence
                        ])
            
            result_handle.close()
            return results  # Retorna os resultados ao invés de escrever diretamente no arquivo

        except Exception as e:
            print(f"Erro no BLAST ({query_id}, tentativa {attempt + 1}): {e}")
            time.sleep(random.uniform(2, 5))  # Espera antes de tentar novamente

    print(f"Falha ao obter BLAST para {query_id} após {attempts} tentativas.")
    return []  # Retorna lista vazia em caso de falha

# Função para salvar os resultados de forma eficiente
def save_results_to_csv(results, output_csv):
    with open(output_csv, "a", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerows(results)

# Arquivo de entrada e saída
input_fasta = "/home/brenno/IC/Projeto2/scripts/Pipeline/GCF_000503995.2_genoma/Outputs/CDHIT/GCF_000503995.2_CerSol_1.0_genomic_Cicadasta.fa"
output_blastn_csv = "/home/brenno/IC/Projeto2/scripts/Pipeline/GCF_000503995.2_genoma/Outputs/BlastN_BlastX/blastn_results.csv"
output_blastx_csv = "/home/brenno/IC/Projeto2/scripts/Pipeline/GCF_000503995.2_genoma/Outputs/BlastN_BlastX/blastx_results.csv"

# Parâmetros do BLAST
blast_configs = [
    {"program": "blastn", "database": "nt", "output_csv": output_blastn_csv},
    {"program": "blastx", "database": "nr", "output_csv": output_blastx_csv},
]

# Carregar as sequências do arquivo FASTA
with open(input_fasta, "r") as fasta_file:
    fasta_sequences = list(SeqIO.parse(fasta_file, "fasta"))

# Processar cada configuração de BLAST
for config in blast_configs:
    program = config["program"]
    database = config["database"]
    output_csv = config["output_csv"]

    print(f"Running {program} for all sequences...")

    # Criar o arquivo de saída e escrever o cabeçalho
    with open(output_csv, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["query_id", "subject_id", "subject_title", "evalue", "per_identity", "coverage", "accession", "query_sequence"])

    # Usando ThreadPoolExecutor para paralelizar as requisições BLAST
    with ThreadPoolExecutor(max_workers=2) as executor:  # Reduzi para 2 para evitar sobrecarga no servidor
        future_to_seq = {executor.submit(run_blast, seq_record, program, database): seq_record.id for seq_record in fasta_sequences}
        
        all_results = []
        for future in as_completed(future_to_seq):
            results = future.result()
            if results:
                all_results.extend(results)

    # Salvar os resultados no arquivo CSV
    save_results_to_csv(all_results, output_csv)

    print(f"{program} completed. Results saved in {output_csv}.")
