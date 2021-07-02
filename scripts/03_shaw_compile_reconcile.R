



# ================== Combine GMPD2/HP3/EID2 associations with Shaw, reconcile names, and standardise detection methods ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/repos/clover/")
pacman::p_load("dplyr", "magrittr", "rentrez")



# ================= read datasets ================

# read datasets (associations)
assoc = read.csv("./output/hostpathogen_reconciled/EID2GMPD2HP3_Associations_CLOVERT_Harmonised_Jan2021.csv", stringsAsFactors = FALSE)

# shaw
shaw = read.csv("./data/source_databases/Shaw_Database_Rentrez_CLOVERT.csv", stringsAsFactors = FALSE)
shaw$Pathogen_Orig2 = shaw$Species
shaw$Host_Original = tolower(shaw$HostSpecies.new)
shaw$Host_AsReported = tolower(shaw$HostSpecies)
shaw$Pathogen_AsReported = NA

# initially harmonised pathogen names for Shaw and save
shawp = read.csv("./data/pathogennames_harmonised_rg/Shaw_allpathogens_harmonised_rg_gfa_20210407.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Pathogen_Original, Pathogen_Harmonised, Pathogen_Harm2)
shaw = left_join(shaw, shawp)
write.csv(shaw, "./output/hostpathogen_reconciled/Shaw_CLOVERT_Harmonised_Jan2021.csv", row.names=FALSE)




# ================= Combine into full associations dataframe with other datasets ====================

# subset to relevant columns
shaw_sub = shaw[ , c("Pathogen_Original", "Pathogen_Orig2", "Pathogen_AsReported", "Pathogen_Harmonised", "Pathogen_Harm2", "Type", "Host_Original", "Host_AsReported",
                     "PublicationYear", "CitationID", "ReferenceText", "Method")] %>%
  dplyr::rename("PathogenType" = Type, "DetectionMethod_Original" = Method) %>%
  dplyr::mutate(Database = "Shaw",
                PathogenType = tolower(PathogenType),
                CitationIDType = ifelse(!is.na(CitationID), "PMID", NA), 
                PublicationDate=NA, 
                ReleaseYear=NA, 
                ReleaseDate=NA)
assoc = rbind(assoc, shaw_sub)



# ------------ manual corrections to pathogen synonymy in complete database ---------------

# 1. "human polyomavirus" just listed as polyomavirus in Olival but 2 strains in Shaw; check Olival assoc and change to either
# 2. "human mastadenovirus a-f" and "human adenovirus a-f" are synonymous: 
# 3. "human enterovirus b and c" -> "enterovirus b - c"
# 3. Shaw lists 3 sub-strains of feline sarcoma virus, but just listed as one in other datasets: have aggregated to one for now

# correction of problems
assoc = assoc[ !is.na(assoc$Pathogen_Harmonised), ]
assoc[ assoc$Pathogen_Harmonised == "human enterovirus b", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "enterovirus b"
assoc[ assoc$Pathogen_Harmonised == "human enterovirus c", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "enterovirus c"
assoc[ assoc$Pathogen_Harmonised == "human mastadenovirus a", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human adenovirus a"
assoc[ assoc$Pathogen_Harmonised == "human mastadenovirus b", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human adenovirus b"
assoc[ assoc$Pathogen_Harmonised == "human mastadenovirus c", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human adenovirus c"
assoc[ assoc$Pathogen_Harmonised == "human mastadenovirus d", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human adenovirus d"
assoc[ assoc$Pathogen_Harmonised == "human mastadenovirus e", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human adenovirus e"
assoc[ assoc$Pathogen_Harmonised == "human mastadenovirus f", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human adenovirus f"
assoc[ assoc$Pathogen_Harmonised == "human polyomavirus", c("Pathogen_Harmonised", "Pathogen_Harm2") ] = "human polyomavirus 1"

# order and keep only distinct records
# most of the duplicates are from GMPD2 and are multiple records that originally had different metadata (e.g. prevalence) attached
table(assoc$Database[ which(duplicated(assoc))] )
assoc = assoc %>%
  distinct()




# ========================= harmonise detection methods ==================================

# any EID2 entries based on NCBI Nucleotide labelled as "PCR"
assoc$DetectionMethod_Original[ assoc$Database == "EID2" & assoc$CitationIDType == "NCBI Nucleotide" ] = "NCBI Nucleotide"

# detection methods
dm = assoc[ !duplicated(assoc$DetectionMethod_Original), c("DetectionMethod_Original"), drop=FALSE ]

# none reported
det0 = c("na", "CargoCarrier", "Other")
dm$Detection_NotSpecified = ifelse(dm$DetectionMethod_Original %in% det0 | is.na(dm$DetectionMethod_Original), TRUE, FALSE)

# serology/antibody-based methods
sero = c("bcELISA", "Antibodies", "Serology", "plaque reduction neutralization", 
         "ELISA", "SNT", "ELISA, SNT", "PAGE", "VNT", "PAGE", "hemmaglutination", "virus neutralization", 
         "Direct Fluorescent Antibody Testing", "CF", "VNT, ELISA", "agar gel immunodiffusion", 
         "IFA", "antigen", "Antigens", "Antigen", "Hemagglutination inhibition assay; Neutralization test",
         "PRNT", "NT", "neutralization test", "ELISA and plaque reduction", "PHA", "SouthernBlot", "Antibodies and Isolation")
dm$Detection_Serology = ifelse(dm$DetectionMethod_Original %in% sero, TRUE, FALSE)

# genetic detection methods
gen_det = c("PCR", "Cell culture, PCR", "PCR, Isolation", "DNA RFLP", 
            "EM, DNA (i.e. more than just PCR)", "Pyroseq", "RNA", "NCBI Nucleotide")
dm$Detection_Genetic = ifelse(dm$DetectionMethod_Original %in% gen_det, TRUE, FALSE)

# pathogen direct isolation/observation within/from tissue
# this could be more nuanced
iso_det = c("DirectBlood", "DirectFecal", "DirectOther", "Tissue", "Fecal", 
            "Isolation", "histopath, e-microscopy", "histopath; e microscopy", "isolation", 
            "Cell culture", "micropsy", "Antibodies and Isolation", "Pathology", 
            "Histopathology", "Histology", "PCR, Isolation", "Cell culture, PCR")
dm$Detection_Isolation = ifelse(dm$DetectionMethod_Original %in% iso_det, TRUE, FALSE)

# composite column of highest quality detection
dm$DetectionMethod_Harmonised = "Not specified"
dm$DetectionMethod_Harmonised[ dm$Detection_Serology == TRUE ] = "Antibodies"
dm$DetectionMethod_Harmonised[ dm$Detection_Genetic == TRUE ] = "PCR/Sequencing"
dm$DetectionMethod_Harmonised[ dm$Detection_Isolation == TRUE ] = "Isolation/Observation"

# combine and save for NCBITaxonomy resolving
assoc = left_join(assoc, dm, by="DetectionMethod_Original")

# remove release date (keep year)
assoc = assoc %>%
  dplyr::select(-PublicationDate, -ReleaseDate, -Pathogen_Original) %>%
  dplyr::rename("Pathogen_Original"=Pathogen_Orig2)

# save
write.csv(assoc, "./output/clover_versions/CLOVER_Associations_Initial.csv", row.names=FALSE)

