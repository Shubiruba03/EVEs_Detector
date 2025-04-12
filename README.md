# EVE Detector

Ferramenta em Python para detecção e caracterização de Elementos Virais Endógenos (EVEs).

## 📌 Objetivo

Este pipeline realiza uma série de etapas automáticas para identificar possíveis EVEs em genomas, utilizando:

- Predição de ORFs (GETORF)
- Alinhamento com proteínas virais (DIAMOND)
- Filtragem com base no material genético (NCBI + ICTV)
- Remoção de redundância (CD-HIT)
- BLASTx/BLASTn contra bancos nt/nr

---

## ⚙️ Requisitos

> ⚠️ Esta ferramenta é compatível apenas com **Linux** ou **Windows com WSL**.

### 📦 Dependências Python

- Python 3.9+
- Biopython
- pandas
- openpyxl
- python-dotenv

Instale todos os pacotes Python necessários com o seguinte comando:

```bash
pip install biopython pandas openpyxl python-dotenv
```

---

### 🛠️ Instalação das ferramentas externas (Linux / WSL / Ubuntu)

- [DIAMOND](https://github.com/bbuchfink/diamond)
- [CD-HIT](https://github.com/weizhongli/cdhit)
- [GETORF (EMBOSS)](http://emboss.sourceforge.net/)

Execute os seguintes comandos para instalar as ferramentas necessárias:

```bash
sudo apt update && sudo apt install -y diamond-aligner cd-hit emboss
```

Verifique se tudo está funcionando com:

```bash
diamond --version
cd-hit-est -h
getorf -h
```

---

## 🔐 Configuração da API do NCBI

Para consultar a taxonomia via NCBI, é necessário fornecer um e-mail e uma chave de API:

1. Crie uma conta gratuita em: https://www.ncbi.nlm.nih.gov/account/
2. Gere sua chave em: https://www.ncbi.nlm.nih.gov/account/settings/

Depois crie um arquivo chamado `.env`

Edite o arquivo `.env`:

```
NCBI_EMAIL = 'seu_email@exemplo.com'
NCBI_API_KEY = 'sua_chave_api'
```

---

## 📥 Entrada obrigatória

Coloque os seguintes arquivos na pasta `Input/`.

- `bancodedados.dmnd` → banco de dados de proteínas virais no formato DIAMOND
- `ICTV_Master_Species_List_2022_MSL38.v2.xlsx` → planilha oficial do ICTV

⚠️ No arquivo `src/main.py`, você **deve editar manualmente** os caminhos desses arquivos, ajustando as variáveis:

```python
db_path = "/caminho/para/Input/bancodedados_proteinas_virais.dmnd"
ictv_file = "/caminho/para/Input/ICTV_Master_Species_List_2022_MSL38.v2.xlsx"
```

Substitua pelos caminhos corretos no seu sistema local.

## ▶️ Uso

Execute a ferramenta com:

```bash
python src/main.py <Assembly_ID>
```

> Os resultados serão salvos automaticamente na pasta `Outputs/`, criada dinamicamente dentro da pasta de download do genoma.

---
