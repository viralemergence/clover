# A harmonized mammal-virus association database

[![DOI](https://zenodo.org/badge/312029561.svg)](https://zenodo.org/badge/latestdoi/312029561)

You're looking for data for your research on wildlife virology, macroecology, or zoonotic risk prediction. But there are too many options! They're all slightly different, and some haven't been updated for years. How do you know which one to use? If everyone is using different data, could that be the reason they have different findings?

## What is 🍀?

CLOVER is a reconciled aggregate of four popular datasets that catalog host-virus associations: EcoHealth Alliance's [HP3](https://github.com/ecohealthalliance/HP3), UGA's [GMPD2](http://onlinelibrary.wiley.com/doi/10.1002/ecy.1799/suppinfo), U. Liverpool's [EID2](https://eid2.liverpool.ac.uk/), and an unnamed dataset curated by [Shaw](https://doi.org/10.6084/m9.figshare.8262779) _et al._ The four datasets have all been harmonized internally and then against the internal taxonomy used by NCBI on popular data portals like GenBank. Then, we've merged them to create the best available data source on the host-virus network (and save you the trouble of merging on your own). 

Buy 1️⃣ get 3️⃣ free!

## Who curates this?
Rory Gibb, Gregory Albery, Timothée Poisot, Colin Carlson, and Max Farrell

## Key files

The integrated CLOVER dataset of associations between mammal hosts and viruses (database and separate csv of field descriptions)

 ↳ clover → clover_0.1_mammalviruses → CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv

 ↳ clover → clover_0.1_mammalviruses → CLOVER_ColumnDescriptions.csv
 
Metadata and scripts for linking CLOVER to a [mammal phylogeny](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494) and information on host species' domestication status

 ↳ clover → clover_0.1_mammalviruses → phylogenies

 ↳ clover → clover_0.1_mammalviruses → domestic_status
 
A full reconciled dataset of all four databases that includes non-mammal hosts, and non-virus pathogens and parasites. This is standardised to the NCBI's taxonomic backbone, but higher taxonomy for a subset of fuzzy matched species is currently missing (3 database files for different pathogen subgroups, and a column descriptors csv)

 ↳ clover → clover_1.0_allpathogens
 
## Everything else you might want to know

### Data usage agreement

To use these data, please cite the DOI for the repository provided by Zenodo, and a citation of the preprint:

R Gibb, GF Albery, DJ Becker, L Brierley, R Connor, TA Dallas, EA Eskew, MJ Farrell, AL Rasmussen, SJ Ryan, AR Sweeny, CJ Carlson, T Poisot. Data proliferation, reconciliation, and synthesis in viral ecology. First posted January 16, 2021. bioRxiv DOI: 10.1101/2021.01.14.426572v1.

### Contact information

For inquiries contact Rory Gibb or Colin Carlson.
