# scripts_mass_spectrometry
Used for MS Data processing and analysis.

## Notes
### yeast_protein_table_generation.py
Generate maps linking UniProtID to various protein annotations, output is a json file including following information:<br>
  - UniProt_ID
  - Description
  - Gene_Name
  - OLN
  - SGD_ID
  - Structure
  - GO
  - Cluster_Name (from KEGG & Complex Portal)

#### Sources of annotation<br>
UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz<br>
-  UniProt_ID
-  Description<br> 
UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/yeast.txt<br>
- Gene_Name(if multiple names exist, record the 1st one)
- OLN(Ordered locus name)
- SGD_ID
- Structure(if related structures solved)
- Residu_Number<br> 
UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping_selected.tab.gz<br>
- GO<br> 
KEGG: see class KEGG_DataBase()<br>
- Cluster_Name -> KEGG<br>
Complex Portal: https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv<br>
- Cluster_Name -> Complex Portal
#### Notice<br>
- write_name_oln_sgdid_structure_resinum(): protein P10081，P02994，P32324，P02309 have 2 gene locus, so there are 2 OLN records.
- class KEGG_DataBase() is used to generate local KEGG database using KEGG API, but it will not be updated unless manually operated.
