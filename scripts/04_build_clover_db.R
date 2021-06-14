

# ===================== Builds CLOVER database ====================

# 1. Resolves higher taxonomy with ref to NCBITaxonomy
# 2. Splits PubMed/Nucleotide fields
# 3. Adds ICTV flag 
# 4. Finalises column names and saves database with metadata


# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/repos/clover/")
library(dplyr); library(magrittr)
source("./scripts/00_highertax_func.R")


# ================ associations data ====================

# associations data 
assoc = read.csv("./output/clover_versions/CLOVER_Associations_Initial.csv", stringsAsFactors = FALSE) 





# --------------- 1. resolve higher taxonomy to NCBI: generate host and pathogen dictionaries --------------

# hosts
# hosts = unique(assoc$Host_Original)
# hhtax = hdict(hosts)
#write.csv(hhtax, "./output/highertax_dict/hosts_highertax.csv", row.names=FALSE)

# pathogens
# paths = unique(assoc$Pathogen_Harm2)
# phtax = vdict(paths)
#write.csv(phtax, "./output/highertax_dict/paths_highertax.csv", row.names=FALSE)

# combine associations with higher tax after manual checks
assoc = assoc %>%
  left_join(
    read.csv("./output/highertax_dict/hosts_highertax_manualcheck_r.csv") %>%
      dplyr::select(-ManualResolved_Note) %>% dplyr::rename(Host_Original="HostOriginal")
  ) %>%
  left_join(
    read.csv("./output/highertax_dict/paths_highertax_manualcheck_r.csv") %>%
      dplyr::select(-PathogenType) %>%
      dplyr::rename("Pathogen_Harm2"=1, "PathogenTaxID"=2, "PathogenNCBIResolved"=3,
                    "Pathogen"=4, "PathogenGenus"=5, "PathogenFamily"=6, "PathogenOrder"=7, "PathogenClass"=8)
  ) %>%
  dplyr::select(-Pathogen_Harmonised, -Pathogen_Harm2) %>%
  dplyr::rename("HostOriginal"=Host_Original,
                "PathogenOriginal"=Pathogen_Original)



# ------------------ 2. Split CitationID field ------------------

# more useful fields
assoc$NCBIAccession = assoc$CitationID
assoc$NCBIAccession[ (assoc$CitationIDType != "NCBI Nucleotide" & !is.na(assoc$CitationIDType)) ] = NA
assoc$PMID = assoc$CitationID
assoc$PMID[ (assoc$CitationIDType != "PMID" & !is.na(assoc$CitationIDType)) ] = NA
assoc = dplyr::select(assoc, -CitationID, -CitationIDType)



# ----------------- 3. Ratified by ICTV flag for viruses -----------------

# list of ICTV ratified viruses
ictv = read.csv("./data/ictv_ratifiedlist/ICTV Master Species List 2019.v1.csv")
assoc$ICTVRatified = ifelse(tolower(assoc$Pathogen) %in% tolower(ictv$Species), TRUE, FALSE)



# ---------------- 4. Manual fixes to some taxonomy issues -------------------

# bse/cjd 
# "torque teno virus" is redundant as species name, but still to be updated in NCBI; update to unclassified Alphatorquevirus
assoc %<>% dplyr::filter( !assoc$PathogenOriginal %in% c("bse agent", "cjd agent") )
assoc$PathogenGenus[ assoc$Pathogen == "Torque teno virus" ] = "Alphatorquevirus"
assoc$Pathogen[ assoc$Pathogen == "Torque teno virus" ] = NA





# ==================== Build final CLOVER dataset ======================

# rename columns
clover = assoc %>%
  dplyr::rename("DetectionMethod" = DetectionMethod_Harmonised,
                "DetectionMethodOriginal" = DetectionMethod_Original) %>%
  dplyr::mutate(PathogenType = Hmisc::capitalize(PathogenType))

# pathogen type fix
clover$PathogenType[ clover$PathogenType %in% c("Bacteria", "Bacteria/rickettsia") ] =  "Bacteria/Rickettsia"

# database version
database_version = data.frame(Database=c("EID2", "Shaw", "GMPD2", "HP3"),
                              DatabaseDOI = c("https://doi.org/10.1038/sdata.2015.49", "https://doi.org/10.1111/mec.15463", "https://doi.org/10.1002/ecy.1799", "https://doi.org/10.5281/zenodo.596810"),
                              DatabaseVersion = c("Wardeh et al. 2015 Sci Data", "Shaw et al. 2020 Mol Ecol", "Stephens et al. 2017 Ecology", "Olival et al. 2017 Nature"))
clover = left_join(clover, database_version)

# compile and reorder columns
clover = clover %>%
  dplyr::select(Host, HostClass, HostOrder, HostFamily, HostGenus, HostTaxID, HostNCBIResolved,
                Pathogen, PathogenType, PathogenClass, PathogenOrder, PathogenFamily, PathogenGenus, 
                PathogenTaxID, PathogenNCBIResolved, ICTVRatified,
                PublicationYear, ReleaseYear, ReferenceText, PMID, NCBIAccession,
                Database, DatabaseVersion, DatabaseDOI,
                DetectionMethod, Detection_NotSpecified, Detection_Serology, Detection_Genetic, Detection_Isolation,
                HostOriginal, PathogenOriginal, DetectionMethodOriginal) %>%
  distinct() %>%
  dplyr::arrange(PathogenType, HostClass, Host, Pathogen)

# tolower
clover = clover %>%
  dplyr::mutate_at(c("Host", "HostClass", "HostOrder", "HostFamily", "HostGenus",
                     "Pathogen", "PathogenType", "PathogenClass", "PathogenOrder", "PathogenFamily", "PathogenGenus"),
                   tolower)




# =========================== save database ============================


# --------------------- 1. Save mammal viruses CLOVER -------------------------

clovm = clover %>%
  dplyr::filter(HostClass == "mammalia" & PathogenType == "virus") %>%
  dplyr::arrange(HostOrder, Host, Pathogen) %>%
  dplyr::rename("Virus"=Pathogen, "VirusClass"=PathogenClass, "VirusOrder"=PathogenOrder, "VirusFamily"=PathogenFamily,
                "VirusGenus"=PathogenGenus, "VirusTaxID"=PathogenTaxID, "VirusNCBIResolved"=PathogenNCBIResolved,
                "VirusOriginal"=PathogenOriginal) %>%
  dplyr::select(-PathogenType)
write.csv(clovm, "./clover/clover_0.1_mammalviruses/CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv", row.names=FALSE)

meta = data.frame(ColName = colnames(clovm))
meta$Description = c("Host species",
                     "Host taxonomic class",
                     "Host taxonomic order", 
                     "Host taxonomic family",
                     "Host genus",
                     "Host NCBI taxid",
                     "Was host matched in NCBI Taxonomy database? TRUE/FALSE",
                     "Virus species",
                     "Virus taxonomic class",
                     "Virus taxonomic order", 
                     "Virus taxonomic family", 
                     "Virus genus",
                     "Virus NCBI taxid",
                     "Was virus matched in NCBI Taxonomy database? TRUE/FALSE",
                     "Is virus ratified by ICTV? TRUE/FALSE",
                     "Year association was published in scientific literature (if applicable and known)",
                     "Year association record was reported on NCBI Nucleotide database (if applicable)",
                     "Text description of source reference; either transferred verbatim from original database, accessed via PMID, or NCBI Nucleotide",
                     "PubMed ID if primary source is in PubMed",
                     "NCBI Nucleotide accession is primary source is NCBI Nucleotide",
                     "Source database",
                     "Version of source database that was accessed for CLOVER",
                     "DOI for source database version",
                     "Detection method (reconciled and harmonised to a standardised classification system)",
                     "True/false flag",
                     "True/false flag",
                     "True/false flag",
                     "True/false flag",
                     "Host species as listed in source database",
                     "Virus species as listed in source database",
                     "Detection method as described in source database")
write.csv(meta, "./clover/clover_0.1_mammalviruses/CLOVER_ColumnDescriptions.csv")




# ------------------ 2. Save full all-hosts-all-pathogens CLOVER -------------------

# save in separate pathogen types to fit into GitHub size
write.csv(clover[ clover$PathogenType == "virus", ], "./clover/clover_1.0_allpathogens/CLOVER_1.0_Viruses_AssociationsFlatFile.csv", row.names=FALSE)
write.csv(clover[ clover$PathogenType == "bacteria/rickettsia", ], "./clover/clover_1.0_allpathogens/CLOVER_1.0_Bacteria_AssociationsFlatFile.csv", row.names=FALSE)
write.csv(clover[ !clover$PathogenType %in% c("bacteria/rickettsia", "virus"), ], "./clover/clover_1.0_allpathogens/CLOVER_1.0_HelminthProtozoaFungi_AssociationsFlatFile.csv", row.names=FALSE)

# create descriptors
meta = data.frame(ColName = colnames(clover))
meta$Description = c("Host species",
                     "Host taxonomic class",
                     "Host taxonomic order", 
                     "Host taxonomic family",
                     "Host genus",
                     "Host NCBI taxid",
                     "Was host matched in NCBI Taxonomy database? TRUE/FALSE",
                     "Pathogen species",
                     "Pathogen type (bacteria/rickettsia, virus, helminth, protozoa, fungi, helminth, other)",
                     "Pathogen taxonomic class",
                     "Pathogen taxonomic order", 
                     "Pathogen taxonomic family", 
                     "Pathogen taxonomic genus",
                     "Pathogen NCBI taxid",
                     "Was pathogen matched in NCBI Taxonomy database? TRUE/FALSE",
                     "Is virus species ratified by ICTV? TRUE/FALSE",
                     "Year association was published in scientific literature (if applicable and known)",
                     "Year association record was reported on NCBI Nucleotide database (if applicable)",
                     "Text description of source reference; either transferred verbatim from original database, accessed via PMID, or NCBI Nucleotide",
                     "PubMed ID if primary source is in PubMed",
                     "NCBI Nucleotide accession is primary source is NCBI Nucleotide",
                     "Source database",
                     "Version of source database that was accessed for CLOVER",
                     "DOI for source database version",
                     "Detection method (reconciled and harmonised to a standardised classification system)",
                     "True/false flag",
                     "True/false flag",
                     "True/false flag",
                     "True/false flag",
                     "Host species as listed in source database",
                     "Pathogen species as listed in source database",
                     "Detection method as described in source database")
write.csv(meta, "./clover/clover_1.0_allpathogens/CLOVER_ColumnDescriptions.csv")



