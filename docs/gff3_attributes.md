# GFF3 attribute reference

This is the full attribute list emitted by `dante_ltr` in the output GFF3.
Most users only need the core fields (`Rank`, `ID`, `Final_Classification`,
`LTR_Identity`, `TSD`); the remainder are diagnostic, carried over from
DANTE, or used by downstream tools.

## Element-level attributes

- **Rank**: Rank of the element (`D`, `DL`, `DLT`, `DLP`, `DLTP`).
- **ID**: A unique identifier for the feature.
- **Parent**: Indicates the parent feature of the current feature, here a
  transposable element ID.
- **Ndomains**: The number of protein domains found within a transposable
  element.
- **LTR_Identity**: The percentage of sequence identity between 5' and 3'
  LTR sequences.
- **LTR5_length** / **LTR3_length**: The lengths of the 5' and 3' LTRs,
  respectively.
- **TSD** (Target Site Duplication): The sequence of the target site
  duplication.
- **Final_Classification**: A hierarchical classification of the
  transposable element based on the REXdb classification system.

## Protein-domain attributes (carried from DANTE)

- **Name**: The attribute is part of DANTE output and corresponds to the
  name of the protein domain (RT, RH, PROT, ...).
- **Best_Hit**: Information about the best match of the protein domain to
  a known database entry.
- **Best_Hit_DB_Pos**: Position of the best hit within the database
  sequence.
- **DB_Seq**: The database sequence that corresponds to the best hit.
- **Region_Seq**: The sequence of the region where the query sequence was
  aligned to the database sequence.
- **Query_Seq**: The sequence of the query used to find the best hit.
- **Identity**: The percentage identity of the best hit match.
- **Similarity**: The similarity score of the best hit match.
- **Relat_Length**: The relative length of the match compared to the
  database sequence.
- **Relat_Interruptions**: Indicates the relative number of interruptions
  in the domain sequence. Interruption could be either a stop codon or a
  frameshift.
- **Hit_to_DB_Length**: The length of the hit compared to the database
  sequence length.

## Primer binding site attributes

- **trna_id**: The identifier for the tRNA related to the primer binding
  site.
- **PBS_evalue**: The E-value associated with the primer binding site.
