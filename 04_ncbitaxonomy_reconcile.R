


# ================= Reconcile host and virus taxonomy using NCBITaxonomy.jl lookup ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("RISmed", "dplyr", "magrittr")

# associations data
assoc = read.csv("./output/intermediate_versions/Clover_reconciledassociations_v1_20201120.csv", stringsAsFactors = FALSE)
assoc = assoc[ assoc$Pathogen_Harmonised != "bse agent", ]

# remove nucleotide year (lookup issues)
assoc$Year[ assoc$YearType == "Nucleotide" ] = NA



# =============== Resolve viruses using NCBITaxonomy lookup ======================

# viruses from NCBI taxonomy, with  fuzzy matches manually checked
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



# ============================ resolve hosts using NCBITaxonomy ========================

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


# =========================== reorder columns ===========================

# Host = host species (reconciled with ref to NCBI)
# HostClass, HostOrder, HostFamily = taxonomic information (resolved against NCBI)
# Virus = virus species (reconciled with ref to NCBI)
# VirusClass, VirusOrder, VirusFamily, VirusGenus = taxonomic info (resolved against NCBI)
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
  dplyr::select(Host, HostClass, HostOrder, HostFamily, Virus, VirusClass, VirusOrder, VirusFamily, VirusGenus, Year, YearType, Database, DatabaseVersion, 
                DetectionMethod, Detection_NotSpecified, Detection_Serology, Detection_Genetic, Detection_Isolation,
                Host_Original, Virus_Original, DetectionMethod_Original, Host_NCBIResolved, Virus_NCBIResolved, HostSynonyms) 

# save
write.csv(assoc, "./output/Clover_v1.0_NBCIreconciled_20201211.csv", row.names=FALSE)



# ================ visual inspection and correction of host record issues ====================

# corrections of host name errata: examine species names that do not match to IUCN species names
# manual comparison of NCBITaxonomy, IUCN and full database to determine current host name
# the vast majority are homotypic synonyms in NCBI; manually corrected for IUCN and phylogeny harmonising and synonyms updated
assoc$Host[ assoc$Host == "Dusicyon thous"] = "Cerdocyon thous"
assoc$Host[ assoc$Host_Original == "Canis lupus familiaris" ] = "Canis lupus familiaris"
assoc$Host[ assoc$Host_Original == "Lycalopex gymnocercus" ] = "Lycalopex gymnocercus"
assoc$HostSynonyms[ assoc$Host == "Galerella pulverulenta" ] = "Herpestes pulverulentus"
assoc$HostSynonyms[ assoc$Host == "Galerella sanguinea" ] = "Herpestes sanguineus"
assoc$Host[ assoc$Host == "Mustela vison" ] = "Neovison vison"
assoc$HostSynonyms[ assoc$Host == "Neovison vison" ] = "Mustela vison"
assoc$Host[ assoc$Host == "Puma yagouaroundi" ] = "Herpailurus yagouaroundi"
assoc$HostSynonyms[ assoc$Host == "Herpailurus yagouaroundi" ] = paste(assoc$HostSynonyms[ assoc$Host_Original == "Herpailurus yagouaroundi" ], "Puma yagouaroundi", sep=", ")
assoc$Host[ assoc$Host == "Uncia uncia" ] = "Panthera uncia"
assoc$Host[ assoc$Host == "Urva javanica" ] = "Herpestes javanicus"
assoc$HostSynonyms[ assoc$Host == "Herpestes javanicus" ] = "Urva javanica"
assoc$Host[ assoc$Host == "Alcelaphus lichtensteinii" ] = "Alcelaphus buselaphus"
assoc$Host[ assoc$Host == "Alces americanus" ] = "Alces alces"
assoc$Host[ assoc$Host == "Bos grunniens mutus" ] = "Bos mutus"
assoc$Host[ assoc$Host == "Capra hircus aegagrus" ] = "Capra aegagrus"
assoc$Host[ assoc$Host == "Ovis aries orientalis" ] = "Ovis orientalis"
assoc$Host[ assoc$Host == "Taurotragus oryx" ] = "Tragelaphus oryx"
assoc$HostSynonyms[ assoc$Host == "Tragelaphus oryx" ] = "Taurotragus oryx"
assoc$Host[ assoc$Host == "Aeorestes cinereus" ] = "Lasiurus cinereus"
assoc$Host[ assoc$Host == "Aeorestes egregius" ] = "Lasiurus egregius"
assoc$HostSynonyms[ assoc$Host == "Lasiurus cinereus" ] = "Aorestes cinereus"
assoc$HostSynonyms[ assoc$Host == "Lasiurus egregius" ] = "Aorestes egregius"
assoc$Host[ assoc$Host == "Artibeus cinereus" ] = "Dermanura cinerea"
assoc$HostSynonyms[ assoc$Host == "Dermanura cinerea" ] = "Artibeus cinereus" 
assoc$Host[ assoc$Host == "Artibeus phaeotis" ] = "Dermanura phaeotis"
assoc$HostSynonyms[ assoc$Host == "Dermanura phaeotis" ] = "Artibeus phaeotis" 
assoc$Host[ assoc$Host == "Artibeus toltecus" ] = "Dermanura tolteca"
assoc$HostSynonyms[ assoc$Host == "Dermanura tolteca" ] = "Artibeus toltecus" 
assoc$Host[ assoc$Host == "Dasypterus ega" ] = "Lasiurus ega"
assoc$HostSynonyms[ assoc$Host == "Lasiurus ega" ] = "Dasypterus ega" 
assoc$Host[ assoc$Host == "Dasypterus intermedius" ] = "Lasiurus intermedius"
assoc$HostSynonyms[ assoc$Host == "Lasiurus intermedius" ] = "Dasypterus intermedius" 
assoc$Host[ assoc$Host == "Dasypterus xanthinus" ] = "Lasiurus xanthinus"
assoc$HostSynonyms[ assoc$Host == "Lasiurus xanthinus" ] = "Dasypterus xanthinus" 
assoc$Host[ assoc$Host_Original == "Pipistrellus abramus" ] = "Pipistrellus abramus"
assoc$Host[ assoc$Host_Original == "Pteronotus davyi" ] = "Pteronotus davyi"
assoc$Host[ assoc$Host_Original == "Hipposideros abae" ] = "Hipposideros abae"
assoc$Host[ assoc$Host_Original == "Hipposideros crumeniferus" ] = "Hipposideros crumeniferus"
assoc$Host[ assoc$Host == "Macronycteris commersonii" ] = "Macronycteris commersoni"
assoc$Host[ assoc$Host == "Megaderma lyra" ] = "Lyroderma lyra"
assoc$HostSynonyms[ assoc$Host == "Lyroderma lyra"  ] = "Megaderma lyra"
assoc$HostSynonyms[ assoc$Host == "Miniopterus fuliginosus"  ] = "Miniopterus schreibersii"
assoc$Host[ assoc$Host == "Myonycteris angolensis" ] = "Lissonycteris angolensis"
assoc$HostSynonyms[ assoc$Host == "Myotis ricketti"  ] = paste(assoc$HostSynonyms[ assoc$Host == "Myotis ricketti"  ], "Myotis pilosus", sep=", ")
assoc$Host[ assoc$Host == "Perimyotis subflavus subflavus" ] = "Perimyotis subflavus"
assoc$Host[ assoc$Host == "Tadarida plicata" ] = "Chaerephon plicatus"
assoc$HostSynonyms[ assoc$Host ==  "Chaerephon plicatus"  ] = "Tadarida plicata"
assoc$Host[ assoc$Host == "Micoureus demerarae" ] = "Marmosa demerarae"
assoc$Host[ assoc$Host == "Dorcopsis veterum" ] = "Dorcopsis luctuosa"
assoc$Host_NCBIResolved[ assoc$Host == "Dorcopsis luctuosa" ] = FALSE
assoc$Host[ assoc$Host == "Notamacropus agilis" ] = "Macropus agilis"
assoc$Host[ assoc$Host == "Notamacropus dorsalis" ] = "Macropus dorsalis"
assoc$Host[ assoc$Host == "Notamacropus eugenii" ] = "Macropus eugenii"
assoc$Host[ assoc$Host == "Notamacropus parma" ] = "Macropus parma"
assoc$Host[ assoc$Host == "Notamacropus parryi" ] = "Macropus parryi"
assoc$Host[ assoc$Host == "Notamacropus rufogriseus" ] = "Macropus rufogriseus"
assoc$Host[ assoc$Host == "Osphranter robustus" ] = "Macropus robustus"
assoc$Host[ assoc$Host == "Osphranter rufus" ] = "Macropus rufus"
assoc$Host[ assoc$Host == "Macaca balantak" ] = "Macaca tonkeana"
assoc$HostSynonyms[ assoc$Host == "Macaca tonkeana" ] = "Macaca balantak"
assoc$Host[ assoc$Host == "Coendou rothschildi" ] = "Coendou quichua"
assoc$Host[ assoc$Host == "Dipodillus dasyurus" ] = "Gerbillus dasyurus"
assoc$HostSynonyms[ assoc$Host == "Gerbillus dasyurus" ] =  "Dipodillus dasyurus"
assoc$Host[ assoc$Host == "Erethizon dorsatus" ] = "Erethizon dorsatum"
assoc$Host[ assoc$Host == "Gerbilliscus kempii" ] = "Gerbilliscus kempi"
assoc$Host[ assoc$Host == "Lasiopodomys gregalis" ] = "Microtus gregalis"
assoc$HostSynonyms[ assoc$Host == "Microtus gregalis" ] = "Lasiopodomys gregalis"
assoc$Host[ assoc$Host == "Liomys adspersus" ] = "Heteromys adspersus"
assoc$HostSynonyms[ assoc$Host == "Heteromys adspersus" ] = "Liomys adspersus"
assoc$Host[ assoc$Host == "Liomys salvini" ] = "Heteromys salvini"
assoc$HostSynonyms[ assoc$Host == "Heteromys salvini" ] = "Liomys salvini"
assoc$Host[ assoc$Host == "Myotomys unisulcatus" ] = "Otomys unisulcatus"
assoc$HostSynonyms[ assoc$Host == "Otomys unisulcatus" ] = "Myotomys unisulcatus"
assoc$Host[ assoc$Host == "Rattus flavipectus" ] = "Rattus tanezumi"
assoc$HostSynonyms[ assoc$Host == "Rattus tanezumi" ] = "Rattus flavipectus"
assoc$Host[ assoc$Host == "Tamias amoenus" ] = "Neotamias amoenus"
assoc$HostSynonyms[ assoc$Host == "Neotamias amoenus" ] = "Tamias amoenus"
assoc$Host[ assoc$Host == "Tamias minimus" ] = "Neotamias minimus"
assoc$HostSynonyms[ assoc$Host == "Neotamias minimus" ] = "Tamias minimus"
assoc$Host[ assoc$Host == "Tamias quadrivittatus" ] = "Neotamias quadrivittatus"
assoc$HostSynonyms[ assoc$Host == "Neotamias quadrivittatus" ] = "Tamias quadrivittatus"
assoc$Host[ assoc$Host == "Tamias sibiricus" ] = "Eutamias sibiricus"
assoc$HostSynonyms[ assoc$Host == "Eutamias sibiricus" ] = "Tamias sibiricus"
assoc$Host[ assoc$Host == "Tamias umbrinus" ] = "Neotamias umbrinus"
assoc$HostSynonyms[ assoc$Host == "Neotamias umbrinus" ] = "Tamias umbrinus"
assoc$Host[ assoc$Host == "Tupaia chinensis" ] = "Tupaia belangeri"

# save
write.csv(assoc, "./output/Clover_v1.0_NBCIreconciled_20201211.csv", row.names=FALSE)












# xx = assoc[ !is.na(assoc$Year), ] %>%
#   group_by(Host, Virus) %>%
#   dplyr::summarise(MinYear = min(Year, na.rm=TRUE)) %>%
#   dplyr::group_by(MinYear) %>%
#   dplyr::summarise(UniqueHVAssoc = length(Host))
# pp = ggplot(xx) + 
#   geom_bar(aes(MinYear, UniqueHVAssoc), stat="identity", col="black", fill="skyblue") + 
#   theme_minimal() +
#   xlab("Year") + ylab("Unique host-virus associations reported") +
#   geom_vline(xintercept=2010, lty=2, size=1, col="darkred")
# ggsave(pp, file="barplot_timeseries.png", device="png", width=8, height=6, units="in", dpi=600)
# 
# err = assoc[ !assoc$Host %in% iucn$binomial, ]
# write.csv(err, "errata.csv", row.names=FALSE)


# ================ record of corrected issues ================

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







# # 2b. hosts (correct) from NBCI scrape; exact matches
# hosts_exact = read.csv("./output/poisot_ncbitaxonomy/clover-complete-poisot.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(type == "hosts") %>%
#   dplyr::filter(!name %in% hosts$name) %>%
#   dplyr::select(name, order, class, family, species) %>%
#   dplyr::mutate(ncbiexact = TRUE)
# 
# # combine all hosts
# hosts = rbind(hosts_exact, hosts)
# all(hosts$name %in% assoc$Host_Harmonised); all(assoc$Host_Harmonised %in% hosts$name)
# names(vir) = c("Pathogen_Harmonised", "VirusClass", "VirusOrder", "VirusFamily", "VirusGenus", "Virus", "Virus_NCBIResolved")
# 
# # replace corrected fuzzy match records (rg)
# hosts$species[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ] = hosts$species_updated[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ]
# hosts$order[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ] = hosts$order_updated[ hosts$rg_updated == TRUE & !is.na(hosts$rg_updated) ]
