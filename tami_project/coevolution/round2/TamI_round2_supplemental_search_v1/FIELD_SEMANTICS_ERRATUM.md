# ClusteredNR hit-table field correction

The downloaded 13-column tables were previously interpreted as ending in `qcovs`. Direct reconstruction with local BLAST shows that column 13 is `ppos` (percent positives): for the TamI branch, 250/250 rows match local `ppos` within 0.01, whereas 0/250 match local `qcovs` within 0.01. The median absolute difference from true qcovs is 29.505 percentage points.

No raw file or QC threshold was changed. True qcovs was recomputed by aligning every complete protein sequence to TamI reference ADC79647.1. This recovers 152 unique sequences that pass the frozen TamI-reference rules and were previously excluded only because percent positives had been mislabeled as coverage. Previously rejected post-MSA motif failures remain blocked from re-entry.
