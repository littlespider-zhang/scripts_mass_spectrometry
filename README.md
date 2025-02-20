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
  - Cluster_KEGG
  - Cluster_Complex_Portal
  - Cluster_Custom

#### Sources of annotation<br>
##### UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz<br>
-  UniProt_ID
-  Description
##### UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/yeast.txt<br>
- Gene_Name(if multiple names exist, record the 1st one)
- OLN(Ordered locus name)
- SGD_ID
- Structure(if related structures solved)
- Residu_Number
##### UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping_selected.tab.gz<br>
- GO<br>
##### KEGG: see class KEGG_DataBase()<br>
- Cluster_KEGG
##### Complex Portal: https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv<br>
- Cluster_Complex_Portal
##### Custom defined clusters: db_custom_cluster.txt
- Cluster_Custom
#### Notice<br>
- class KEGG_DataBase() is used to generate local KEGG database using KEGG API, but it will not be updated unless manually operated.
- write_name_oln_sgdid_structure_resinum(): protein P10081，P02994，P32324，P02309 have 2 gene locus, so there are 2 OLN records.
- write_Complex_Portal(): some proteins can be assigned to several complexes, so they have several Cluster_Complex_Portal records.
