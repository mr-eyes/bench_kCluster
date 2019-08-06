# Data Preparation

## 1- [Human EST database](https://www.ncbi.nlm.nih.gov/genbank/dbest)

Running the following script will download all the [dbEST reports](https://ftp.ncbi.nih.gov/repository/dbEST/), parse and covert them to FASTA files with applied filter to extract only "Homo Sapiens".

**Output:** `debEST/dbEST_fasta/`

```bash
cd dbEST && bash prepare_dbEST.sh
```
