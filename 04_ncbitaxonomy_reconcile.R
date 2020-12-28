


# ================= Reconcile mammal host and virus taxonomy against NCBITaxonomy.jl lookup ==================

# Host and virus harmonised names are standardised against the NCBI Taxonomy backbone
# This script resolves names based on NCBITaxonomy.jl lookups to ensure consistent nomenclature
# Finally, any outstanding host species without clear matches in NCBI are cross-referenced against IUCN and Wilson & Reeder (~100 species)

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("dplyr", "magrittr")

# associations data
assoc = read.csv("./output/intermediate_versions/Clover_reconciledassociations_v1_20201120.csv", stringsAsFactors = FALSE)
assoc = assoc[ !assoc$Pathogen_Harmonised %in% c("bse agent", "pepper mild mottle"), ]

# remove nucleotide year (lookup issues)
#assoc$Year[ assoc$YearType == "Nucleotide" ] = NA



# =============== Resolve viruses against NCBITaxonomy lookup ======================

# viruses from NCBI taxonomy, with fuzzy matches manually checked
vir_checked = read.csv("./output/poisot_ncbitaxonomy/clover-viruses-poisot-rgchecked.csv", stringsAsFactors = FALSE)

# correct all data for species with fuzzy matches
viruses = vir_checked %>%
  dplyr::filter(rg_updated==TRUE & !is.na(rg_updated))
viruses$species = viruses$species_updated
viruses$order = viruses$order_updated
viruses$family = viruses$family_updated
viruses$class = viruses$class_updated
viruses$genus = viruses$genus_updated
viruses$ncbiexact = ifelse(viruses$unresolved_notes == "", TRUE, FALSE)
viruses = viruses[ , c("name", "class", "order", "family", "genus", "species", "ncbiexact") ]

# viruses without issues
vir_exact = vir_checked[ !vir_checked$name %in% viruses$name, ] %>%
  dplyr::select(name, class, order, family, genus, species) %>%
  dplyr::mutate(ncbiexact = TRUE)

# combine all viruses with taxonomic information
vir = rbind(viruses, vir_exact)
vir = vir[ vir$name != "bse agent", ]
#all(vir$name %in% assoc$Pathogen_Harmonised); all(assoc$Pathogen_Harmonised %in% vir$name)
names(vir) = c("Pathogen_Harmonised", "VirusClass", "VirusOrder", "VirusFamily", "VirusGenus", "Virus", "Virus_NCBIResolved")

# combine and rename cols
assoc = left_join(assoc, vir, by="Pathogen_Harmonised") %>%
  dplyr::select(-Pathogen_Harmonised) %>%
  dplyr::rename("Virus_Original" = Pathogen_Original)



# ========================= resolve hosts from NCBITaxonomy ========================

# 2a. update hosts from NBCITaxonomy lookup
hosts = read.csv("./output/poisot_ncbitaxonomy/conflicts-poisot-rg-family.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(type == "hosts") %>%
  dplyr::filter(rg_updated==TRUE & !is.na(rg_updated))

# replace corrected fuzzy match records (rg)
hosts$species = hosts$species_updated
hosts$order = hosts$order_updated
hosts$family = hosts$family_updated
hosts$class = hosts$class_updated
hosts$ncbiexact = ifelse(hosts$unresolved.notes == "", TRUE, FALSE)
hosts = hosts[ , c("name", "order", "family", "species", "ncbiexact") ]

# combine, and for all viruses where pathogen_harmonised is correct, assign this to virusspecies_harmonised column
hosts = hosts %>%
  dplyr::rename("Host_Harmonised"=1, "HostOrder_2"=2, "HostFamily_2" =3, "HostHarm_2" = 4, "Host_NCBIResolved" = 5)
assoc = left_join(assoc, hosts, by="Host_Harmonised")
assoc$HostHarm_2[ is.na(assoc$HostHarm_2) ] = assoc$Host_Harmonised[ is.na(assoc$HostHarm_2) ]
assoc$HostOrder_2[ is.na(assoc$HostOrder_2) ] = assoc$HostOrder[ is.na(assoc$HostOrder_2) ]
assoc$HostFamily_2[ is.na(assoc$HostFamily_2) ] = assoc$HostFamily[ is.na(assoc$HostFamily_2) ]
assoc$Host_NCBIResolved[ is.na(assoc$Host_NCBIResolved) ] = TRUE
assoc$HostHarm_2 = Hmisc::capitalize(assoc$HostHarm_2)

# replace and rename host harmonised
assoc = assoc %>%
  dplyr::mutate(Host_Harmonised = HostHarm_2,
                HostOrder = HostOrder_2,
                HostFamily = HostFamily_2) %>%
  dplyr::select(-HostHarm_2, -HostOrder_2, -HostFamily_2) %>%
  dplyr::rename("Host" = Host_Harmonised,
                "DetectionMethod" = DetectionMethod_Harmonised) %>%
  dplyr::mutate(Host_Original = Hmisc::capitalize(Host_Original),
                Virus_Original = Hmisc::capitalize(Virus_Original),
                PathogenType = "Virus")



# =========================== reorder columns and save ===========================

# Host = host species (first taxized and then resolved against NCBI)
# HostClass, HostOrder, HostFamily = taxonomic information (resolved against NCBI)
# Virus = virus species (manually harmonised then reconciled with ref to NCBI)
# VirusClass, VirusOrder, VirusFamily, VirusGenus = taxonomic info (accessed from NCBI)
# Year = year reported
# YearType = how was year obtained? (author = in source db; pubmed = from pubmed scrape; nucleotide = from nucleotide scrape)
# Database = source database
# DatabaseVersion = version of source database that was accessed for CLOVER
# DetectionMethod = detection method (reconciled and harmonised to a simple classification system)
# DetectionMethod_NotSpecified, Serology, Genetic, Isolation = TRUE/FALSE flags for detection method
# Host_Original = host species as listed in source database
# Virus_Original = virus species as listed in source database
# DetectionMethod_Original = detection method as reported in source database
# Host or Virus_NCBIResolved = does Host/Virus column have an exact match in NCBITaxonomy? TRUE/FALSE
# HostSynonyms = synonyms of host species (accessed from taxize)

# reorder columns
assoc = assoc %>%
  dplyr::select(Host, HostClass, HostOrder, HostFamily, Virus, VirusClass, VirusOrder, VirusFamily, VirusGenus, Year, YearType, Database, DatabaseVersion, Source,
                DetectionMethod, Detection_NotSpecified, Detection_Serology, Detection_Genetic, Detection_Isolation,
                Host_Original, Virus_Original, DetectionMethod_Original, Host_NCBIResolved, Virus_NCBIResolved, HostSynonyms) %>%
  distinct()

# # save corrected records for Tim
# vir_corrected = vir_checked %>%
#   dplyr::filter(rg_updated==TRUE & !is.na(rg_updated)) %>%
#   dplyr::select(-checked_prev)
# hosts_corrected =  read.csv("./output/poisot_ncbitaxonomy/conflicts-poisot-rg-family.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(type == "hosts") %>%
#   dplyr::filter(rg_updated==TRUE & !is.na(rg_updated)) %>%
#   dplyr::rename("unresolved_notes" = unresolved.notes) %>%
#   dplyr::mutate(class = "", 
#                 family = "", 
#                 genus = "")
# corrected = rbind(vir_corrected, hosts_corrected)
# write.csv(corrected, "./output/poisot_ncbitaxonomy/nbcitax-amendedrecords-rg.csv", row.names=FALSE)




# ================ final reconcilation stage: comparison of host names against IUCN records ===============

# a substantial number of host names in previous round didn't return a clear result in NCBITaxonomy
# cross-reference host lookups to IUCN range map data for mammals and manually check issues against NCBI
# 101 non-matches; flag issues to Tim and update records accordingly

# iucn
iucn = sf::st_read("C:/Users/roryj/Documents/PhD/202011_clover/data/iucn_range/MAMMALS/MAMMALS.shp")

# issues
issues = assoc[ !assoc$Host %in% iucn$binomial, ] %>% 
  dplyr::select(Host, HostClass, HostOrder, HostFamily, Host_Original, HostSynonyms, Host_NCBIResolved) %>%
  distinct()
#write.csv(issues, "./output/iucn_crossref/HostIssues_CLOVER.csv", row.names=FALSE)

# manual comparison of NCBITaxonomy, IUCN and full database to determine current host name
# the vast majority are homotypic synonyms in NCBI; manually corrected for IUCN and phylogeny harmonising and synonyms updated
issues_rg = read.csv("./output/iucn_crossref/HostIssues_CLOVER_rg.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Record_Updated == TRUE)

# correct issues in assoc database 1: domestic dog across all records ahnd "Artiodactyla"
assoc$Host[ assoc$Host_Original == "Canis lupus familiaris" ] = "Canis lupus familiaris"
assoc$HostOrder[ assoc$HostOrder == "Artiodactyla" ] = "Cetartiodactyla"

# resolve 
foo = assoc[ assoc$Host %in% issues_rg$Host, ] %>%
  left_join(issues_rg) %>%
  dplyr::mutate(Host = Host2,
                HostSynonyms = HostSynonyms2, 
                HostFamily = HostFamily2) %>%
  dplyr::select(-Host_IUCNMatch, -Record_Updated, -Description_Of_Issue, -Host2, -HostSynonyms2, -HostFamily2)

# combine, order and save
assoc_updated = rbind( assoc[ !assoc$Host %in% issues_rg$Host, ], foo) %>%
  dplyr::arrange(HostOrder, Host, VirusOrder, Virus)
write.csv(assoc_updated, "./output/Clover_v1.0_NBCIreconciled_20201218.csv", row.names=FALSE)




# ================= column descriptors CSV ==================

# create
meta = data.frame(ColName = colnames(assoc_updated))
meta$Description = c("Host species (first taxized and then resolved against NCBI Taxonomy, with outstanding issues cross-referenced to IUCN)",
                     "Host taxonomic class",
                     "Host taxonomic order", 
                     "Host taxonomic family",
                     "Virus species (manually harmonised across datasets, then resolved against NCBI Taxonomy)",
                     "Virus taxonomic class",
                     "Virus taxonomic order", 
                     "Virus taxonomic family", 
                     "Virus taxonomic genus",
                     "Year association was reported (if multiple year were included in source database record, the earliest is reported)",
                     "How was year obtained? Either scientific publication, listed in source database, or scraped from PubMed PMID or NCBI Nucleotide. If year is NA, it was not possible to resolve from source database",
                     "Source database",
                     "Version of source database that was accessed for CLOVER",
                     "Was association reported in a publication, or accessed from NCBI Nucleotide?",
                     "Detection method (reconciled and harmonised to a simple classification system)",
                     "True/false flag",
                     "True/false flag",
                     "True/false flag",
                     "True/false flag",
                     "Host species as listed in source database",
                     "Virus species as listed in source database",
                     "Detection method as described in source database",
                     "Does host species name have an exact match in NCBITaxonomy? True/false",
                     "Does virus species name have an exact match in NCBITaxonomy? True/false",
                     "Synonyms of host species, accessed from taxize")
write.csv(meta, "./output/Clover_v1.0_ColumnDescriptions_20201218.csv")
