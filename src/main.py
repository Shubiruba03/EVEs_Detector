import sys
import os
from dotenv import load_dotenv
load_dotenv()
from PipelineFunctions import(
    baixar_genoma, 
    run_getorf, 
    run_diamond_blastx, 
    filtrar_evalue_phage, 
    ncbi_taxon_filter, 
    converter_xlsx_para_fasta,
    executar_cd_hit_est)

# Solicita o Assembly ID ao usuário
if len(sys.argv) != 2:
    print("Uso: python main.py <Assembly_ID>")
    sys.exit(1)

assembly_id = sys.argv[1]

# Caminho do arquivo de banco de dados de proteínas virais
db_path = os.getenv('DB_FILE')

#Caminho pra o arquivo ICTV
ictv_file = os.getenv('ICTV_FILE')

# Executa as todas as funções em sequência
genome_file = baixar_genoma(assembly_id) #Baixa genoma
orf_file = run_getorf(genome_file) #Executa o getorf
blast_file = run_diamond_blastx(orf_file, db_path) #Executa o Diamond blastx
filtered_file= filtrar_evalue_phage(blast_file) #Filtra com base no evalue e retira fagos
ncbi_file = ncbi_taxon_filter(filtered_file, ictv_file) #Busca taxonomia e retira baseado no material genético
fasta_file = converter_xlsx_para_fasta(ncbi_file) #Converte para fasta
CDHIT_file = executar_cd_hit_est(fasta_file)
#ultimo_blast_file
print("Pipeline concluído com sucesso!")
 