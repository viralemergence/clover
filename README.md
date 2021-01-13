# Clover: A harmonized mammal-virus association database

[![DOI](https://zenodo.org/badge/312029561.svg)](https://zenodo.org/badge/latestdoi/312029561)

You're looking for data for your research on wildlife virology, macroecology, or zoonotic risk prediction. But there are too many options! They're all slightly different, and some haven't been updated for years. How do you know which one to use? If everyone is using different data, could that be the reason they have different findings?

## What is 🍀?

CLOVER is a reconciled aggregate of four popular datasets that catalog host-virus associations: EcoHealth Alliance's [HP3](https://github.com/ecohealthalliance/HP3), UGA's [GMPD2](http://onlinelibrary.wiley.com/doi/10.1002/ecy.1799/suppinfo), U. Liverpool's [EID2](https://eid2.liverpool.ac.uk/), and an unnamed dataset curated by [Shaw](https://doi.org/10.6084/m9.figshare.8262779) _et al._ The four datasets have all been harmonized internally and then against the internal taxonomy used by NCBI on popular data portals like GenBank. Then, we've merged them to create the best available data source on the host-virus network (and save you the trouble of merging on your own). 

Buy 1️⃣ get 3️⃣ free!

## Who curates this?
Rory Gibb, Gregory Albery, Timothée Poisot, Colin Carlson, and Max Farrell

## Key files

The integrated CLOVER dataset, currently limited to mammal viruses (database and separate csv of field descriptions)

 ↳ clover → Clover_v1.0_NCBIreconciled_20201218.csv

 ↳ clover → Clover_v1.0_ColumnDescriptions_20201218.csv
 
Metadata and scripts for linking CLOVER to a [mammal phylogeny](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494) and information on host species' domestication status

 ↳ clover ↳ phylogenies

 ↳ clover ↳ domestic_status
 
A full reconciled dataset of all four databases that includes non-virus pathogens and parasites (currently internally harmonised but not yet standardised to NCBI Taxonomy)

 ↳ output ↳ hostpathogen_harmonised → AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv
 
## Everything else you might want to know

### Data usage agreement

To use these data, please cite the DOI for the repository provided by Zenodo, and a manuscript citation (_forthcoming_).

### Contact information

For inquiries contact Rory Gibb or Colin Carlson.
