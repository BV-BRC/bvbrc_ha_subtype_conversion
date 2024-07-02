# HA Subtype Numbering Conversion

## Overview

The HA Subtype Numbering Conversion tool takes influenza HA protein sequence(s) and converts their existing HA position numbering to a different HA numbering scheme using David Burke and Derek Smith’s method that uses both sequence and structure information to propose positions of functional equivalence across different HA subtypes. The analysis starts with the user inputting protein sequence(s). The sequence(s) are BLASTed against the Burke Reference sequences, which will determine the best reference subtype to use in the HA numbering pipeline. The HA numbering pipeline will generate a pairwise multiple sequence alignment using the reference protein sequence selected and the user inputted protein sequence(s). This alignment will generate a mapping between the user input sequence(s) and the BLAST reference sequence. Then this mapping is used to align the input sequence(s) to the selected HA subtype positions.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

There is one application service specification defined here:
1. [HA Subtype Numbering Conversion](app_specs/HASubtypeNumberingConversion.md): The HA Subtype Numbering Conversion service allows user to renumber Influenza HA sequences according to a cross-subtype numbering scheme proposed by Burke and Smith in [Burke DF, Smith DJ.2014. A recommended numbering scheme for influenza A HA subtypes. PLoS One 9:e112302.](https://pubmed.ncbi.nlm.nih.gov/25391151/)

The code in this module provides the BV-BRC application service wrapper scripts for the genome annotation service as well
as some backend utilities:

| Script name | Purpose |
| ----------- | ------- |
| [App-HASubtypeNumberingConversion.pl](service-scripts/App-HASubtypeNumberingConversion.pl) | App script for the [HA Subtype Numbering Conversion Service](https://www.bv-brc.org/docs/quick_references/services/ha_numbering_service.html) |

## See also

* [HA Subtype Numbering Conversion Service](https://www.bv-brc.org/app/HASubtypeNumberingConversion)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/ha_numbering_service.html)
* [HA Subtype Numbering Conversion Service Tutorial](https://www.bv-brc.org/docs/tutorial/ha_numbering/ha_numbering.html)

## References

Burke DF, Smith DJ. A recommended numbering scheme for influenza A HA subtypes. PLoS One. 2014 Nov 12;9(11):e112302. doi: 10.1371/journal.pone.0112302. PMID: 25391151; PMCID: PMC4229193.



