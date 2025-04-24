import sys
from PipelineFunctions import(
    download_genoma, 
    run_getorf, 
    run_diamond_blastx, 
    filtrar_evalue_phage, 
    ncbi_taxon_filter, 
    converter_xlsx_para_fasta,
    executar_cd_hit_est,
    run_blastx_blastn)


# Solicita o Assembly ID ao usuário
if len(sys.argv) != 2:
    print("Uso: python main.py <Assembly_ID>")
    sys.exit(1)

assembly_id = sys.argv[1]


# Caminho do arquivo de banco de dados de proteínas virais
db_path = "/mnt/c/Users/deyvi/OneDrive/Documentos/IC/eves/EVEs_Detector/Input/ICTV_Master_Species_List_2022_MSL38.v2 (1).xlsx"

#Caminho pra o arquivo ICTV
ictv_file = "/mnt/c/Users/deyvi/OneDrive/Documentos/IC/eves/EVEs_Detector/Input/bancodedados.dmnd"

# Executa as todas as funções em sequência
genome_file = download_genoma(assembly_id) #Baixa genoma
orf_file = run_getorf(genome_file) #Executa o getorf
blast_file = run_diamond_blastx(orf_file, db_path) #Executa o Diamond blastx
filtered_file= filtrar_evalue_phage(blast_file) #Filtra com base no evalue e retira fagos
ncbi_file = ncbi_taxon_filter(filtered_file, ictv_file) #Busca taxonomia e retira baseado no material genético
fasta_file = converter_xlsx_para_fasta(ncbi_file) #Converte para fasta
CDHIT_file = executar_cd_hit_est(fasta_file) # Executa cdhit para retirar redundâcias
blastX_N = run_blastx_blastn(CDHIT_file) #Executa blastx e blastn
print("Pipeline concluído com sucesso!")




