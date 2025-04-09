#importação das bibliotecas 
import os
import sys
import zipfile
import subprocess
import re
import openpyxl as op
import csv
import time
import random
from Bio import Entrez
import pandas as pd
from time import sleep
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from dotenv import load_dotenv
load_dotenv()


def download_genoma(assembly_id):
    """Baixa e extrai o genoma a partir do Assembly ID usando NCBI datasets."""
    
    # Nome do arquivo de saída
    zip_filename = f"{assembly_id}.zip"

    # Comando para baixar o genoma
    print(f"Baixando genoma para {assembly_id}...")
    cmd_download = f"datasets download genome accession {assembly_id} --include genome --filename {zip_filename}"
    os.system(cmd_download)

    # Verifica se o download foi bem-sucedido
    if not os.path.exists(zip_filename):
        print("Erro no download. Verifique o Assembly ID.")
        return

    # Extrair o arquivo ZIP
    print("Extraindo os arquivos...")
    with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
        zip_ref.extractall(f"{assembly_id}_genoma")

    # Encontrar o arquivo .fna extraído
    extracted_path = f"{assembly_id}_genoma/ncbi_dataset/data/{assembly_id}"
    fna_file = None

    for root, dirs, files in os.walk(extracted_path):
        for file in files:
            if file.endswith("genomic.fna"):
                fna_file = os.path.join(root, file)
                break

    if fna_file:
        print(f" Download concluído! Arquivo salvo em: {fna_file}")
    else:
        print("Arquivo .fna não encontrado após extração.")

    # Opcional: Remover o ZIP para economizar espaço
    os.remove(zip_filename)

    return fna_file

#*******************************************************************************************
def run_getorf(input_file):
    """
    Executa o programa GETORF para encontrar regiões de ORFs no genoma de entrada.
    - input_file: Caminho do arquivo de genoma (.fna) a ser analisado.
    - Retorna o caminho do arquivo de saída gerado pelo GETORF.
    """
    # Caminho base para a pasta Outputs
    base_folder = os.path.abspath(os.path.join(input_file, "../../../../Outputs"))
    os.makedirs(base_folder, exist_ok=True)
    print(f"Pasta de saídas criada em {base_folder}")

    # Criar a pasta ORFs dentro de Outputs
    orf_folder = os.path.join(base_folder, "ORFs")
    os.makedirs(orf_folder, exist_ok=True)
    
    genome_name = os.path.basename(input_file).replace(".fna", "")  # Obtém o nome do arquivo sem extensão
    orf_output = os.path.join(orf_folder, f"{genome_name}_ORF.fasta")  # Define o nome do arquivo de saída
    
    # Comando para executar o GETORF
    cmd = [
        "getorf",
        "-sequence", input_file,  # Arquivo de entrada
        "-outseq", orf_output,  # Arquivo de saída com as ORFs
        "-minsize", "100",  # Tamanho mínimo das ORFs
        "-maxsize", "6000",  # Tamanho máximo das ORFs
        "-find", "3"  # Busca apenas ORFs na fita direta
    ]
    
    # Executa o comando e verifica se ocorreu algum erro
    subprocess.run(cmd, check=True)
    print(f"GETORF concluído para {genome_name}")
    
    return orf_output  # Retorna o caminho do arquivo gerado pelo GETORF

#*******************************************************************************************
#*******************************************************************************************

def run_diamond_blastx(orf_file, db_path):
    """Executa DIAMOND BLASTX e filtra os resultados."""
    
    # Voltar um nível a partir do arquivo ORF para encontrar a pasta "Outputs"
    outputs_folder = os.path.abspath(os.path.join(os.path.dirname(orf_file), ".."))
    
    # Definir a pasta Diamond dentro de Outputs
    blastx_folder = os.path.join(outputs_folder, "Diamond")
    os.makedirs(blastx_folder, exist_ok=True)

    genome_name = os.path.basename(orf_file).replace("_ORF.fasta", "")  
    diamond_output = os.path.join(blastx_folder, f"{genome_name}_blastResults.tsv")  

    # Comando para executar DIAMOND BLASTX
    cmd = [
        "diamond", "blastx",
        "-d", db_path,
        "-q", orf_file,
        "--outfmt", "6", "qseqid", "sseqid", "qlen", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle", "qtitle", "full_qseq",
        "--max-target-seqs", "5",
        "--out", diamond_output
    ]
    subprocess.run(cmd, check=True)
    print(f"DIAMOND BLASTX concluído para {genome_name}")

    return diamond_output  # Retorna o caminho do arquivo blastx
#*******************************************************************************************
#*******************************************************************************************

def filtrar_evalue_phage(blastx_file):
    """
    Filtra um arquivo BLASTX removendo 'phage' e aplicando filtros de E-value e Bitscore.
    
    Args:
        blastx_file (str): Caminho do arquivo de saída do BLASTX (.tsv).
    
    Returns:
        str: Caminho do arquivo filtrado salvo.
    """
    
    # Voltar um nível para acessar "Outputs"
    outputs_folder = os.path.abspath(os.path.join(os.path.dirname(blastx_file), ".."))
    
    # Criar a pasta "E-value" dentro de Outputs
    output_folder = os.path.join(outputs_folder, "E-value")
    os.makedirs(output_folder, exist_ok=True)

    # Verificar se o arquivo de entrada existe
    if not os.path.exists(blastx_file):
        print(f" Erro: Arquivo {blastx_file} não encontrado!")
        return None

    # Definir colunas do arquivo BLASTX
    column_names = ["qseqid", "sseqid", "qlen", "slen", "pident", "length", "mismatch", 
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                    "stitle", "qtitle", "full_qseq"]

    # Ler arquivo BLASTX
    df = pd.read_csv(blastx_file, sep="\t", names=column_names, dtype=str)

    # Filtro 1: Remover linhas que contêm 'phage' na coluna 'stitle'
    df = df[~df["stitle"].str.contains("phage", case=False, na=False)]

    # Filtro 2: Selecionar apenas a melhor sequência para cada 'qseqid' com base em E-value e Bitscore
    df["evalue"] = df["evalue"].astype(float)
    df["bitscore"] = df["bitscore"].astype(float)
    df = df.sort_values(by=["qseqid", "evalue", "bitscore"], ascending=[True, True, False])
    df = df.drop_duplicates(subset=["qseqid"], keep="first")

    # Criar nome do arquivo de saída substituindo "blastResults" por "filtered"
    file_base = os.path.basename(blastx_file).replace("blastResults", "filtered").replace(".tsv", ".xlsx")
    output_path = os.path.join(output_folder, file_base)

    # Salvar em Excel
    df.to_excel(output_path, index=False, header=True)

    print(f" Fagos removidos e melhores seq. selecionadas para: {file_base}")
    return output_path

#*******************************************************************************************

def ncbi_taxon_filter(filtered_file, ictv_file):
    """
    Filtra um arquivo do NCBI com base na taxonomia e genoma, removendo certas categorias.
    
    Args:
        filtered_file (str): Caminho do arquivo filtrado (.xlsx).
        ictv_file (str): Caminho do arquivo ICTV.
    
    Returns:
        str: Caminho do arquivo processado salvo.
    """
    
    # Configuração do e-mail e chave da API do NCBI
    #Entrez.email = 'jp.uesc17@gmail.com'
    #Entrez.api_key = 'ee7ccfdfa22559163c2bd8f3c822157ae108'
    
    # Configuração do e-mail e chave da API do NCBI
    Entrez.email = os.getenv('NCBI_EMAIL')
    Entrez.api_key = os.getenv('NCBI_API_KEY')

    # Voltar um nível para acessar "Outputs"
    output_folder = os.path.abspath(os.path.join(os.path.dirname(filtered_file), ".."))
    
    # Criar a pasta "NCBI" dentro de Outputs
    output_folder = os.path.join(output_folder, "NCBI")
    os.makedirs(output_folder, exist_ok=True)

    # Ler o arquivo ICTV com dados de vírus
    print("\nCarregando arquivo ICTV...")
    lt_ICTV = pd.read_excel(ictv_file, sheet_name='MSL')
    coluna_familia = lt_ICTV['Family']
    coluna_genoma = lt_ICTV['Genome Composition']
    print("Arquivo ICTV carregado com sucesso!")
    
    # Dicionários para armazenar resultados e evitar buscas repetitivas
    data = {}
    data_genome = {}
    
    def get_ncbi_tax(taxon, max_retries=5):
        """Obtém a taxonomia do NCBI para um determinado taxon"""
        retries = 0
        while retries < max_retries:
            try:
                print(f"Buscando taxonomia para: {taxon}")
                
                if not re.match(r'\d+', taxon):  # Se não for um ID numérico
                    taxon2 = '"' + taxon + '"'
                    handle = Entrez.esearch(db='taxonomy', term=taxon2, rettype='gb', retmode='text')
                    record = Entrez.read(handle, validate=False)
                    handle.close()
                    
                    if not record['IdList']:
                        print(f"Erro: Taxon inválido - {taxon}")
                        return "Desconhecido"
                    tax_id = record['IdList'][0]
                else:
                    tax_id = taxon
    
                handle2 = Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml')
                record2 = Entrez.read(handle2, validate=False)
                handle2.close()
    
                tax_list = record2[0]['LineageEx']
                
                apresenta = []
                valor, atualiza = 0, 0
                for tax_element in tax_list:
                    if tax_element['Rank'] == 'family':
                        apresenta.append(tax_element['ScientificName'])
                        valor += 1
                    elif (tax_element['Rank'] == 'no rank' or 'unclassified' in tax_element['ScientificName']) and valor == atualiza:
                        apresenta.append(tax_element['ScientificName'])
                        valor += 1
                        atualiza += 1
                    elif valor == 0 and tax_element == tax_list[-1]:
                        apresenta.append(tax_element['ScientificName'])
    
                print(f"Taxonomia encontrada: {apresenta[-1]}")
                return apresenta[-1]
    
            except:
                retries += 1
                print("Erro na conexão com NCBI. Aguardando 5 segundos e tentando novamente...")
                sleep(5)
    
        print("Falha após múltiplas tentativas.")
        return "Erro na busca"
    
    print(f"\nProcessando arquivo: {filtered_file}")
    print("Carregando dados...")
    Tabela = pd.read_excel(filtered_file, header=None)
    print("Arquivo carregado!")
    
    lista_esp, lista_analisada, lista_adict, lista_adict2 = [], [], [], []
    
    print("Extraindo espécies da coluna 14...")
    for especie in Tabela[14]:
        especie = re.findall(r'\[.*\]', str(especie))
        lista_esp.append(str(especie))
    
    for i in lista_esp:
        i = i.replace('[', '').replace(']', '').replace("'", '')
        lista_analisada.append(i)
    
    print("Buscando famílias para cada espécie...")
    for j in lista_analisada:
        if j in data:
            lista_adict.append(data[j])
        else:
            sleep(0.5)
            variavel = get_ncbi_tax(j)
            lista_adict.append(variavel)
            data[j] = variavel
    
    print("Associando famílias ao genoma correspondente...")
    for j in lista_adict:
        if str(j) in data_genome:
            lista_adict2.append(data_genome[j])
        else:
            a = 0
            for i in range(len(coluna_familia)):
                if str(coluna_familia[i]) == j:
                    data_genome[j] = coluna_genoma[i]
                    a -= 1
                else:
                    a += 1
                    if a == len(coluna_familia):
                        data_genome[j] = j
            lista_adict2.append(data_genome[j])
    
    print("Dados processados com sucesso!")
    
    # Adicionar novas colunas na tabela
    print("Adicionando novas colunas ao arquivo...")
    Tabela['adicao'] = lista_adict
    Tabela['material'] = lista_adict2
    
    # Remover linhas onde a coluna 'material' tenha valores indesejados
    Tabela = Tabela[~Tabela['material'].str.contains('dsDNA|RNA-RT|DNA-RT', na=False)]
    
    # Gerar nome do arquivo de saída
    nome_arquivo = os.path.basename(filtered_file)
    if not nome_arquivo.endswith(".xlsx"):
        nome_arquivo += ".xlsx"
    caminho_saida = os.path.join(output_folder, nome_arquivo)
    
    # Salvar o arquivo processado
    print(f"Salvando arquivo processado como: {caminho_saida}")
    Tabela.to_excel(caminho_saida, sheet_name='famílias')
    print("Arquivo salvo com sucesso!")
    
    return caminho_saida

#*******************************************************************************************

def converter_xlsx_para_fasta(ncbi_file):
    """
    Converte um arquivo Excel do NCBI para o formato FASTA.
    
    Args:
        ncbi_file (str): Caminho do arquivo NCBI de entrada (.xlsx).
    
    Returns:
        str: Caminho do arquivo FASTA gerado ou None em caso de erro.
    """
    # Voltar um nível para acessar "Outputs"
    outputs_folder = os.path.abspath(os.path.join(os.path.dirname(ncbi_file), ".."))
    
    # Criar a pasta "fastas" dentro de Outputs
    fasta_folder = os.path.join(outputs_folder, "fastas")
    os.makedirs(fasta_folder, exist_ok=True)
    
    print(f"\nProcessando arquivo: {ncbi_file}")
    
    # Verificar se o arquivo de entrada existe
    if not os.path.isfile(ncbi_file):
        print(f"Erro: O arquivo {ncbi_file} não foi encontrado.")
        return None
    
    # Nome do arquivo de saída
    output_filename = os.path.splitext(os.path.basename(ncbi_file))[0] + ".fasta"
    output_path = os.path.join(fasta_folder, output_filename)
    
    try:
        # Ler o arquivo Excel
        df = pd.read_excel(ncbi_file, header=None)
    except Exception as e:
        print(f"Erro ao ler o arquivo {ncbi_file}: {e}")
        return None
    
    # Verificar se há dados suficientes
    if df.shape[0] < 2:
        print(f"Erro: O arquivo {ncbi_file} não possui dados suficientes.")
        return None
    
    # Extrair IDs (coluna 3, índice 2) e sequências (coluna 18, índice 17)
    query = df.iloc[1:, 2].astype(str).str.replace(' ', '', regex=True).tolist()
    seq = df.iloc[1:, 17].astype(str).tolist()
    
    print(f"Convertendo {len(query)} sequências para FASTA...")
    
    # Criar o arquivo FASTA
    with open(output_path, "w") as fasta_file:
        for j, sequence in zip(query, seq):
            if sequence.strip():  # Ignorar sequências vazias
                fasta_file.write(f">{j}\n{sequence}\n")
    
    print(f"Arquivo FASTA salvo para: {output_filename}")
    return output_path

    #******************************************************************

def executar_cd_hit_est(fasta_file):
    """
    Executa o CD-HIT-EST para agrupar sequências semelhantes.

    Args: input_file (str): Caminho do arquivo de entrada (.fa).
        output_folder (str): Pasta onde o arquivo processado será salvo.

    Returns: str: Caminho do arquivo de saída do CD-HIT-EST.
    """
    # Voltar um nível para acessar "Outputs"
    output_folder = os.path.abspath(os.path.join(os.path.dirname(fasta_file), ".."))
    
    # Criar a pasta "CDHIT" dentro de Outputs
    output_folder = os.path.join(output_folder, "CDHIT")
    os.makedirs(output_folder, exist_ok=True)

    # Nome base do arquivo sem extensão
    nome_base = os.path.basename(fasta_file).replace(".fa", "")
    nome_cicada = nome_base.replace("filtered", "Cicada")  # Ajuste no nome

    # Caminho do arquivo de saída do cd-hit
    cdhit_output = os.path.join(output_folder, f"{nome_cicada}.fa")

    # Comando do CD-HIT-EST
    comando = [
        "cd-hit-est",
        "-i", fasta_file,
        "-o", cdhit_output,
        "-c", "0.9",
        "-aS", "0.9",
        "-n", "5",
        "-M", "4000",
        "-d", "0",
        "-T", "2"
    ]

    print(f"Executando CD-HIT-EST para {nome_base}...")
    subprocess.run(comando, check=True)
    print(f"CD-HIT-EST concluído. Resultado salvo em {cdhit_output}")

    return cdhit_output


def executar_diamond_blastx(CDHIT_file, db_path):
    """
    Executa o DIAMOND BLASTX para buscar sequências correspondentes.

    Args:
        input_file (str): Caminho do arquivo de entrada (.fa) já processado pelo CD-HIT-EST.
        output_folder (str): Pasta onde o arquivo processado será salvo.
        db_path (str): Caminho do banco de dados DIAMOND.

    Returns:
        str: Caminho do arquivo de saída do DIAMOND BLASTX.
    """
    if not os.path.exists(CDHIT_file):
        raise FileNotFoundError(f"Arquivo de entrada não encontrado: {CDHIT_file}")

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Banco de dados DIAMOND não encontrado: {db_path}")

     # Voltar um nível para acessar "Outputs"
    output_folder = os.path.abspath(os.path.join(os.path.dirname(CDHIT_file), ".."))
    
    # Criar a pasta "LastDiamond" dentro de Outputs
    output_folder = os.path.join(output_folder, "LastDiamond")
    os.makedirs(output_folder, exist_ok=True)

    # Nome base do arquivo sem extensão
    nome_base = os.path.basename(CDHIT_file).replace(".fa", "")
    diamond_output = os.path.join(output_folder, f"{nome_base}_Diamond.tsv")

    # Comando do DIAMOND BLASTX
    comando = [
        "diamond", "blastx",
        "-d", db_path,
        "-q", CDHIT_file,
        "--outfmt", "6",
        "qseqid", "sseqid", "qlen", "slen", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle", "qtitle", "full_qseq",
        "--max-target-seqs", "1",
        "--out", diamond_output
    ]

    print(f"Executando DIAMOND BLASTX para {CDHIT_file}...")
    subprocess.run(comando, check=True)
    print(f"DIAMOND BLASTX concluído. Resultado salvo em {diamond_output}")

    return diamond_output

#*******************************************************************************************

def run_blastx_blastn(input_fasta):
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

    # Criar as pastas BlastN e BlastX dentro de Outputs
    outputs_folder = os.path.abspath(os.path.join(os.path.dirname(input_fasta), ".."))
    blastn_folder = os.path.join(outputs_folder, "BlastN")
    blastx_folder = os.path.join(outputs_folder, "BlastX")
    os.makedirs(blastn_folder, exist_ok=True)
    os.makedirs(blastx_folder, exist_ok=True)

    # Definir os caminhos de saída
    output_blastn_csv = os.path.join(blastn_folder, "blastn_results.csv")
    output_blastx_csv = os.path.join(blastx_folder, "blastx_results.csv")

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
        with ThreadPoolExecutor(max_workers=2) as executor:
            future_to_seq = {executor.submit(run_blast, seq_record, program, database): seq_record.id for seq_record in fasta_sequences}
            all_results = []
            for future in as_completed(future_to_seq):
                results = future.result()
                if results:
                    all_results.extend(results)

        # Salvar os resultados no arquivo CSV
        save_results_to_csv(all_results, output_csv)
        print(f"{program} completed. Results saved in {output_csv}.")

    return output_blastn_csv, output_blastx_csv
