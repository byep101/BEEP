# TamI literature note v0.1

## System
- System ID: TamI
- BGC/product: tirandamycin biosynthetic gene cluster / tirandamycin
- Organism: Streptomyces sp. 307-9
- MIBiG: BGC0001052.5
- Protein IDs: NCBI GenPept ADC79647.1; UniProt D3Y1J3
- Structure: PDB 6XA2 (tirandamycin C-bound), 6XA3 (ligand-free)

## Why TamI is the first benchmark
TamI has BGC information, P450 function, crystal structures, mechanistic papers, and experimentally tested key residues. It is suitable for checking whether residue-level coevolution analysis can recover known functional residues.

## Known benchmark residues
- Leu101
- Leu244
- Leu295

These residues line the bicyclic ketal binding site and provide residue-level benchmark positions for future coevolution or structural mapping analyses. Evidence tiers differ by residue: Leu101 and Leu295 were included in targeted double mutants, while Leu244 is currently a structure-supported pocket residue without a direct point-mutant record in the checked TamI paper.

## Key residue mutation table

Primary evidence: Newmister et al., 2020, ACS Catalysis, "Molecular Basis of Iterative C-H Oxidation by TamI, a Multifunctional P450 monooxygenase from the Tirandamycin Biosynthetic Pathway" (PMID: 33569241; PMCID: PMC7869845; DOI: 10.1021/acscatal.0c03248).

| Benchmark residue | Structural context | Tested TamI variant | Experimental observation | Benchmark interpretation |
| --- | --- | --- | --- | --- |
| Leu101 | Bicyclic ketal binding-site residue near the C14 methyl/ketal region of tirandamycin C. | L101V/H102S double mutant | Product profile did not show relaxed site selectivity, but TirC to TirE conversion decreased substantially. | Activity-sensitive active-site contact; useful positive benchmark residue. |
| Leu244 | Hydrophobic residue lining the bicyclic ketal binding site. | No direct Leu244 point mutant found in the checked 2020 TamI text. | Structure-supported pocket contact; locatable in both 6XA2 and 6XA3. | Structural benchmark residue; mutation evidence should be treated as indirect unless another source is added. |
| Leu295 | Bicyclic ketal binding-site residue on the opposite side of the active site from Leu101/His102. | L295A/T299A double mutant | Double mutant displayed wild-type-like TirC to TirE conversion. | Nearby neutral/tolerant pocket residue; useful contrast against activity-sensitive contacts. |

## Known reaction logic
TamI is an iterative P450 involved in a defined oxidative sequence, including C10 hydroxylation, C11/C12 epoxidation, and C18 hydroxylation. MITE further records pathway steps covering tirandamycin C/E/D/A intermediates and final tirandamycin B formation.

## Data to compare against future predictions
- Do EVcouplings/MIp top-ranked residues include Leu101/Leu244/Leu295 or their local residue network?
- Are predicted residues mapped to the substrate binding pocket or hydrophobic channel?
- Does the signal appear around SRS/loop/channel residues rather than only the heme-proximal conserved motif?

## Local raw data package
- Collection date: 2026-06-09
- MIBiG JSON: `data/raw/mibig/BGC0001052.5_annotations.json`
- NCBI FASTA: `data/raw/ncbi/ADC79647.1_TamI_Streptomyces_307-9.fasta`
- UniProt FASTA: `data/raw/uniprot/D3Y1J3_TamI_UniProt.fasta`
- PDB/mmCIF structures: `data/raw/pdb/6XA2_TamI_tirandamycinC_bound.pdb`, `data/raw/pdb/6XA3_TamI_ligand_free.pdb`, `data/raw/pdb/6XA2.cif`, `data/raw/pdb/6XA3.cif`
- Data collection log: `metadata/data_collection_log_v0.1_TamI.csv`

The official NCBI and UniProt FASTA downloads both contain a 413 aa TamI sequence with matching N- and C-terminal sequence checks. The MIBiG JSON parsed successfully and records BGC0001052, tirandamycin, Streptomyces sp. 307-9, GU385216.1, quality `questionable`, status `active`, and completeness `complete`. The RCSB PDB headers distinguish 6XA2 as tirandamycin C-bound and 6XA3 as ligand-free.
