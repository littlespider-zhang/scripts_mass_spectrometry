# scripts_mass_spectrometry
Used for MS Data processing and analysis.

## Notes
### yeast_protein_table_generation.py
Generate maps linking UniProtID to various protein annotations, output is a json file:<br>
Key<br>
  - UniProt_ID<br>
Values<br>
  - Description
  - Gene_Name
  - OLN
  - SGD_ID
  - Structure
  - GO
  - Cluster_Name

#### Sources of annotation<br>
UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz<br>
-  UniProt_ID
-  Description
UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/yeast.txt<br>
- Gene_Name(if multiple names exist, record the 1st one)
- OLN(Ordered locus name)
- SGD_ID
- Structure
- Residu_Number
UniProtKB: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping_selected.tab.gz<br>
- GO
Custom files:<br>
- Cluster_Name<br>
#### Notice<br>
- write_name_oln_sgdid_structure_resinum(): protein P10081，P02994，P32324，P02309 have 2 gene locus, so there are 2 OLN records.

