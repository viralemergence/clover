


# ================= Reconcile host and virus taxonomy using NCBITaxonomy.jl lookup ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("RISmed", "dplyr", "magrittr")

# associations data
assoc = read.csv("./output/Clover_reconciledassociations_v1_20201120.csv", stringsAsFactors = FALSE)
assoc = assoc[ assoc$Pathogen_Harmonised != "bse agent", ]


# =============== NCBITaxonomy output from Tim ======================

# 1. update viruses from NBCITaxonomy lookup (only flagged issues)
viruses = read.csv("./output/poisot_ncbitaxonomy/poisot_ncbitaxonomy_rg.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(type == "viruses")

# replace corrected fuzzy match records (rg)
viruses$species[ viruses$rg_updated == TRUE & !is.na(viruses$rg_updated) ] = viruses$species_updated[ viruses$rg_updated == TRUE & !is.na(viruses$rg_updated) ]
viruses$order[ viruses$rg_updated == TRUE & !is.na(viruses$rg_updated) ] = viruses$order_updated[ viruses$rg_updated == TRUE & !is.na(viruses$rg_updated) ]

# combine, and for all viruses where pathogen_harmonised is correct, assign this to virusspecies_harmonised column
viruses = viruses[ , c("name", "order", "species")] %>%
  dplyr::rename("Pathogen_Harmonised"=1, "VirusOrder"=2, "VirusSpecies_Harmonised"=3)
assoc = left_join(assoc, viruses, by="Pathogen_Harmonised")
assoc$VirusSpecies_Harmonised[ is.na(assoc$VirusSpecies_Harmonised) ] = assoc$Pathogen_Harmonised[ is.na(assoc$VirusSpecies_Harmonised) ]
assoc$VirusSpecies_Harmonised = Hmisc::capitalize(assoc$VirusSpecies_Harmonised)

# replace and rename pathogen harmonised
assoc = assoc %>%
  dplyr::mutate(Pathogen_Harmonised = VirusSpecies_Harmonised) %>%
  dplyr::select(-VirusSpecies_Harmonised) %>%
  dplyr::rename("Virus" = Pathogen_Harmonised,
                "Virus_Original" = Pathogen_Original)

# 2. update hosts from NBCITaxonomy lookup
hosts = read.csv("./output/poisot_ncbitaxonomy/poisot_ncbitaxonomy_rg.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(type == "hosts")

# replace corrected fuzzy match records (rg)
hosts$species[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ] = hosts$species_updated[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ]
hosts$order[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ] = hosts$order_updated[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ]

# combine, and for all viruses where pathogen_harmonised is correct, assign this to virusspecies_harmonised column
hosts = hosts[ , c("name", "species")] %>%
  dplyr::rename("Host_Harmonised"=1, "HostHarm_2"=2)
assoc = left_join(assoc, hosts, by="Host_Harmonised")
assoc$HostHarm_2[ is.na(assoc$HostHarm_2) ] = assoc$Host_Harmonised[ is.na(assoc$HostHarm_2) ]
assoc$HostHarm_2 = Hmisc::capitalize(assoc$HostHarm_2)

# replace and rename host harmonised
assoc = assoc %>%
  dplyr::mutate(Host_Harmonised = HostHarm_2) %>%
  dplyr::select(-HostHarm_2) %>%
  dplyr::rename("Host" = Host_Harmonised,
                "DetectionMethod" = DetectionMethod_Harmonised) %>%
  dplyr::mutate(Host_Original = Hmisc::capitalize(Host_Original),
                Virus_Original = Hmisc::capitalize(Virus_Original),
                PathogenType = "Virus")


# =========================== reorder columns ===========================

# reorder columns
assoc = assoc %>%
  dplyr::select(Host, HostClass, HostOrder, HostFamily, Virus, VirusOrder, Year, YearType, Database, DatabaseVersion, 
                DetectionMethod, Detection_NotSpecified, Detection_Serology, Detection_Genetic, Detection_Isolation,
                Host_Original,  Virus_Original, DetectionMethod_Original) 

# save
write.csv(assoc, "./output/Clover_v1.0_NBCIreconciled_20201211.csv", row.names=FALSE)
