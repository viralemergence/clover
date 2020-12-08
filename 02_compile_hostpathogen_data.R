

# ================== Compile EID2, GMPD2 and HP3 with harmonised pathogen names and host synonyms ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/")
pacman::p_load("dplyr", "magrittr")



# ================== Compile and process EID2 ======================

# date for nucleotide represents date of last modification rather than origin

# eid2 for mammals and birds cross-referenced to pubmed/nucleotide
# set year for nucleotide: if !is.na(Date2), use Date1 (date of last modification)
eid = read.csv( "./data/EID2_SpeciesInteractions_YearScrapeRG_20102020.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(YearType = ifelse(!is.na(Date_NCBI2), "Nucleotide2", "Nucleotide1"))
eid$YearType[ eid$Database == "pubmed"] = "PubMed"
eid$Year[ eid$YearType == "Nucleotide2"] = substr(eid$Date_NCBI2[ eid$YearType == "Nucleotide2"], 8, 12)
eid$Year[ eid$YearType == "Nucleotide1"] = substr(eid$Date_NCBI1[ eid$YearType == "Nucleotide1"], 8, 12)

# create parasite list, remove 'virus' from end of virus description
par_temp = strsplit(tolower(eid$Cargo), " ")
par_temp1 = lapply(par_temp, function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) })
par_temp1 = unlist(lapply(par_temp1, paste, collapse=" "))

# create parasite species list
par_temp2 = paste(par_temp1, "xxxx", " ")
par_temp2 = tolower(unlist(lapply(strsplit(par_temp2, " "), "[[", 2)))
par_temp2[ par_temp2 == "xxxx" ] = ""

# create updated df
eid2 = data.frame(HostType = eid$Carrier.classification,
                  Host = tolower(eid$Carrier),
                  HostGenus = tolower(unlist(lapply(strsplit(eid$Carrier, " "), "[[", 1))),
                  HostSpecies = tolower(unlist(lapply(strsplit(eid$Carrier, " "), "[[", 2))),
                  IsDomestic = ifelse(eid$Carrier.classification == "Domestic", 1, 0),
                  ParasiteType = tolower(eid$Cargo.classification),
                  Parasite = par_temp1,
                  ParasiteGenus = tolower(unlist(lapply(strsplit(par_temp1, " "), "[[", 1))),
                  ParasiteSpecies = par_temp2,
                  Database = "EID2",
                  Source = eid$Database,
                  Year = eid$Year,
                  NCBI_ID = eid$id,
                  YearType = eid$YearType,
                  stringsAsFactors=F)

# subset to only endoparasites (leaving n=6860 unique pathogens)
eid2 = eid2[ eid2$ParasiteType %in% c("bacteria", "helminth", "virus", "fungi", "protozoa", "others"), ]
eid2 = eid2[ order(eid2$HostType, eid2$Host), ]

# remove humans and non-animals from TaxClass (leaving n=4489 unique pathogens)
eid2 = eid2[ !eid2$HostType %in% c("Human", "Higher Plants", "Protozoa", "Bacteria", "Others", "Fungi", "Bryozoa", "Other Plant", "Green Algae", "Methanothermobacter"), ]

# EID2 pathogen names harmonised (to match to other databases)
assoc = read.csv("./data/harmonisednames_rg/EID2_pathogendata_allpathogens_harmonised_20180916.csv", stringsAsFactors = FALSE) %>%
  select(Pathogen_Original, PathogenName_Harmonised, PathogenType, HumanInfective_Any, DiseaseAgent, IsZoonotic, Disease_GIDEON, Route_GIDEON, UN_subregion, Countries_EID2Wertheim)
eid2 = left_join(eid2, assoc, by=c("Parasite"="Pathogen_Original")) 
eid2$ParasiteType = eid2$PathogenType
eid2 = eid2 %>% select(-PathogenType)
eid2$DetectionQuality = 0
eid2$DetectionMethod = "CargoCarrier"

# save
write.csv(eid2, "./output/hostpathogen_harmonised/EID2_MammalsBirds_Harmonised_Oct2020.csv", row.names=FALSE)



# ===================== Compile and process GMPD2 =========================

# main GMPD2 dataset
gmpd2 = read.csv("./data/GMPD2_main.csv", encoding = "latin1", stringsAsFactors = F)

# extract numerics (year) and trim out second year marker
years = regmatches(gmpd2$Citation, gregexpr("[[:digit:]]+", gmpd2$Citation))
gmpd2$Year = as.numeric(unlist(lapply(years, "[", 1)))
gmpd2 = gmpd2[ !is.na(gmpd2$Year), ]

# could add long-lat but for now not including geo info
# remove records with no binomial name for host
gmpd2 = data.frame(HostType = gmpd2$Group,
                   Host = tolower(gmpd2$HostCorrectedName),
                   HostOrder = tolower(gmpd2$HostOrder),
                   HostFamily = tolower(gmpd2$HostFamily),
                   Parasite = tolower(gmpd2$ParasiteCorrectedName),
                   ParasiteType = tolower(gmpd2$ParType),
                   ParasitePhylum = tolower(gmpd2$ParPhylum),
                   ParasiteClass = tolower(gmpd2$ParClass),
                   Year = gmpd2$Year,
                   YearType = "ScientificPublicationDate",
                   Database = "GMPD2",
                   HostsSampled = gmpd2$HostsSampled,
                   Seroprevalence = gmpd2$Prevalence,
                   Location = gmpd2$LocationName,
                   Longitude = gmpd2$Longitude,
                   Latitude = gmpd2$Latitude,
                   HostSex = gmpd2$HostSex,
                   HostAge = gmpd2$HostAge,
                   DetectionMethod = gmpd2$SamplingType,
                   DetectionQuality = ifelse(gmpd2$SamplingType %in% c("PCR", "DirectOther", "DirectBlood", "DirectFaecal", "Other", "Tissue", "Fecal"), 2, 1))
gmpd2 = gmpd2[ gmpd2$Host != "no binomial name", ]

# sort pathogen names: remove 'virus' from end, remove virus classification from start
path = as.vector(gmpd2$Parasite)
path = unlist(lapply(lapply(strsplit(path, " "), function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) }), paste, collapse=" "))
remove_virus_start = function(x){
  if(length(x) %in% c(0, 1)){ return(x)
  } else if(sum(grep("virus", x[1]))>0){ return(x[2:length(x)])
  } else(return(x))
}
path = unlist(lapply(lapply(strsplit(path, " "), remove_virus_start), paste, collapse=" "))
gmpd2$Parasite = path

# remove records with pathogen name 'abolished'
# and remove ectoparasites
gmpd2 = gmpd2[ !gmpd2$Parasite %in% c("abolished", "no binomial name"), ]
gmpd2 = gmpd2[ gmpd2$ParasiteType %in% c("helminth", "protozoa", "virus", "fungi", "bacteria"), ]

# harmonised pathogen names
assoc = read.csv("./data/harmonisednames_rg/GMPD2_pathogendata_allpathogens_harmonised_20180916.csv", stringsAsFactors = FALSE) %>%
  select(Pathogen_Original, PathogenName_Harmonised, PathogenType, HumanInfective_Any, DiseaseAgent, IsZoonotic, Disease_GIDEON, Route_GIDEON, UN_subregion, Countries_EID2Wertheim)

# get data on missing pathogens from eid2 (some pathogen metadata missing in updated pipeline - not sure why but still in harmonised EID2 file)
missing = gmpd2[ which(!gmpd2$Parasite %in% assoc$Pathogen), ]
missingrecs = read.csv("./data/harmonisednames_rg/EID2_pathogendata_allpathogens_harmonised_20180916.csv", stringsAsFactors = FALSE) %>%
  select(Pathogen_Original, PathogenName_Harmonised, PathogenType, HumanInfective_Any, DiseaseAgent, IsZoonotic, Disease_GIDEON, Route_GIDEON, UN_subregion, Countries_EID2Wertheim) %>%
  filter(Pathogen_Original %in% missing$Parasite)
assoc = rbind(assoc, missingrecs)

# combine and save
gmpd2 = left_join(gmpd2, assoc, by=c("Parasite" = "Pathogen_Original"))
write.csv(gmpd2, "./output/hostpathogen_harmonised/GMPD2_Mammals_Harmonised_Oct2020.csv", row.names=FALSE)



# ======================== Compile and process HP3 (Olival) ============================

# host-pathogen associations based on both serology and infection status
oliv_assoc = read.csv("./data/HP3_associations.csv", stringsAsFactors = F) %>%
  select(-WildDomInReference)
oliv_vir = read.csv("./data/HP3_viruses.csv", stringsAsFactors = F)
oliv_mam = read.csv("./data/HP3_hosts.csv", stringsAsFactors = F)

# standardise names
oliv_assoc$vVirusNameCorrected = tolower(unlist(lapply(strsplit(as.vector(oliv_assoc$vVirusNameCorrected), "_"), paste, collapse=" ")))
oliv_assoc$hHostNameFinal = tolower(unlist(lapply(strsplit(as.vector(oliv_assoc$hHostNameFinal), "_"), paste, collapse=" ")))
oliv_vir$vVirusNameCorrected = tolower(unlist(lapply(strsplit(as.vector(oliv_vir$vVirusNameCorrected), "_"), paste, collapse=" ")))
oliv_mam$hHostNameFinal = tolower(unlist(lapply(strsplit(as.vector(oliv_mam$hHostNameFinal), "_"), paste, collapse=" ")))

# subset virus/host data to relevant columns
oliv_vir = oliv_vir[ , c("vVirusNameCorrected", "vOrder", "vFamily", "vSubfamily", "vGenus", "vSSoDS", "vDNAoRNA", "vEnvelope", "vSegmentedTF", "vCytoReplicTF", "vGenomeAveLength", "ReverseZoonoses", "IsZoonotic", "IsZoonotic.stringent")]
oliv_mam = oliv_mam[ , c("hHostNameFinal", "hOrder", "hFamily", "hGenus", "hSpecies") ]

# combine
hp3 = left_join(oliv_assoc, oliv_mam)
hp3 = left_join(hp3, oliv_vir)

# remove 'virus' from pathogen names to standardise with other databases
hp3$vVirusNameCorrected = unlist(lapply(lapply(strsplit(hp3$vVirusNameCorrected, " "), function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) }), paste, collapse=" "))

# get numeric years and take earliest reported year
years = regmatches(hp3$Reference, gregexpr("[[:digit:]]+", hp3$Reference))
yy = unlist(lapply(years, function(x) min(as.numeric(x))))
yy[ yy == Inf ] = NA
hp3$Year = yy
hp3$YearType = "ScientificPublicationDate_Earliest"

# rename columns to standardise with other datasets
names(hp3) = c("Parasite", "Host", "DetectionMethod", "DetectionQuality", "Citation", "HostOrder", "HostFamily", "HostGenus", "HostSpecies", 
               "ParasiteOrder", "ParasiteFamily", "ParasiteSubfamily", "ParasiteGenus", "SSorDS", "DNAorRNA", "Envelope", "Segmented", 
               "CytoReplic", "GenomeAveLength", "ReverseZoonosis", "IsZoonotic_HP3", "IsZoonotic_Strict_HP3", "Year", "YearType")
hp3$ParasiteType = "virus"
hp3$Database = "HP3"

# harmoniesd pathogen names and definitions
assoc = read.csv("./data/harmonisednames_rg/Olival_pathogendata_allpathogens_harmonised_20180916.csv", stringsAsFactors = FALSE) %>%
  select(Pathogen_Original, PathogenName_Harmonised, PathogenType, HumanInfective_Olival,  Zoonotic_Olival_Strict, HumanInfective_External, IsZoonotic, DiseaseAgent, Disease_GIDEON, Route_GIDEON) %>%
  dplyr::mutate(HumanInfective_Any = ifelse(HumanInfective_Olival == 1 | HumanInfective_External == 1, 1, 0)) %>%
  select(-HumanInfective_Olival, -Zoonotic_Olival_Strict)
hp3 = left_join(hp3, assoc, by=c("Parasite"="Pathogen_Original"))

# standardise columns (tolower for all taxonomy)
hp3[ , c("HostOrder", "HostFamily", "HostGenus", "HostSpecies", "ParasiteOrder", "ParasiteFamily", "ParasiteSubfamily", "ParasiteGenus")] = apply(hp3[ , c("HostOrder", "HostFamily", "HostGenus", "HostSpecies", "ParasiteOrder", "ParasiteFamily", "ParasiteSubfamily", "ParasiteGenus")], 2, tolower)

# save
write.csv(hp3, "./output/hostpathogen_harmonised/HP3_MammalViruses_Harmonised_Oct2020.csv", row.names=FALSE)



# ===================== Combine and standardise definitions and host names ====================

# read in files
eid2 = read.csv("./output/hostpathogen_harmonised/EID2_MammalsBirds_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)
gmpd2 = read.csv("./output/hostpathogen_harmonised/GMPD2_Mammals_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)
hp3 = read.csv("./output/hostpathogen_harmonised/HP3_MammalViruses_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)

# combine keeping columns present in all 3 databases
cols = Reduce(intersect, list(names(eid2), names(gmpd2), names(hp3)))
assoc = rbind(eid2[ , cols], gmpd2[ , cols], hp3[ , cols])

# remove detection quality (finalise when combined with Shaw)
assoc = assoc[ , - which(names(assoc) == "DetectionQuality")]

# # harmonise detection methods
# dm = assoc$DetectionMethod
# harmDM = function(x){
#   d = dm[x]
#   dm_return = c()
#   if(any(c("Cargocarrier", "NA", "Other") %in% d)) dm_return = c(dm_return, "Cargocarrier")
#   if(any(c("Serology", "Antibodies", "bcELISA", "IFA", "CF", "antigen", "Antigen", "Antigens", "micropsy", "PHA", "hemmaglutination",
#            "PAGE", "Direct Fluorescent Antibody Testing", "agar gel immunodiffusion", "Antibodies and Isolation", "SouthernBlot") %in% d)) dm_return = c(dm_return, "Serology")
#   if(any(c("SNT", "VNT", "plaque reduction neutralization", "PRNT", "neutralization test", "NT", "ELISA and plaque reduction",
#            "Hemagglutination inhibition assay; Neutralization test", "virus neutralization") %in% d)) dm_return = c(dm_return, "VirusNeutralisationAssay")
#   if(any(c("PCR", "PCR/Direct") %in% d)) dm_return = c(dm_return, "PCR")
#   if(any(c("DNA RFLP", "EM, DNA (i.e. more than just PCR)", "Pyroseq", "RNA") %in% d)) dm_return = c(dm_return, "DirectDetection_GeneticOther")
#   if(any(c("DirectBlood", "DirectFecal", "DirectOther", "Tissue", "Fecal", "Isolation", "histopath, e-microscopy", "histopath; e microscopy", "isolation", "Cell culture", "Antibodies and Isolation") %in% d)) dm_return = c(dm_return, "DirectIsolationOrObservation")
#   if("Reservoir" %in% d) dm_return = c(dm_return, "Reservoir")
#   if(length(dm_return)>1) dm_return = paste(dm_return, collapse=", ")
#   if(is.null(dm_return)) dm_return = "Cargocarrier"
#   return(dm_return)
# }
# dmx = unlist(lapply(1:length(dm), harmDM))
# assoc$DetectionMethod_Harmonised = dmx
# 
# # harmonise detection quality
# assoc$DetectionQuality_Harmonised = 0
# assoc$DetectionQuality_Harmonised[ grep("Serology", assoc$DetectionMethod_Harmonised) ] = 1
# assoc$DetectionQuality_Harmonised[ grep("VirusNeutralisationAssay", assoc$DetectionMethod_Harmonised) ] = 1
# assoc$DetectionQuality_Harmonised[ grep("PCR", assoc$DetectionMethod_Harmonised) ] = 2
# assoc$DetectionQuality_Harmonised[ grep("DirectIsolationOrObservation", assoc$DetectionMethod_Harmonised) ] = 2
# assoc$DetectionQuality_Harmonised[ grep("DirectDetection_GeneticOther", assoc$DetectionMethod_Harmonised) ] = 2

# standardise host names (resolve synonyms using ITIS/CoL)
# data compiled during PREDICTs/pathogens paper (Nature, 2020)
all_syns = read.csv("./data/hostsynonyms_rg/hostdata_allsynonyms_itis_20180917.csv", stringsAsFactors=F, encoding="latin1")
all_syns1 = all_syns %>%
  group_by(Original) %>%
  dplyr::summarise(Host_Harmonised = tolower(head(Accepted_name, 1)),
                   HostClass = head(Selected_class, 1),
                   HostOrder = head(Selected_order, 1),
                   HostFamily = head(Selected_family, 1),
                   HostSynonyms = tolower(paste(unique(Synonyms), collapse=", ")))
all_syns1$HostSynonyms[ all_syns1$HostSynonyms == "na" ] = ""

# manual fix on species unsuccessful in scrape
all_syns2 = data.frame(Original = c("pica pica", "felis catus", "inia geoffrensis", "mesoplodon grayi", "physeter catodon", "physeter macrocephalus"),
                       Host_Harmonised =c("pica pica", "felis catus", "inia geoffrensis", "mesoplodon grayi", "physeter macrocephalus", "physeter macrocephalus"),
                       HostSynonyms = c("", "", "", "", "physeter catodon", "physeter catodon"),
                       HostClass = c("Aves", "Mammalia", "Mammalia", "Mammalia", "Mammalia", "Mammalia"), 
                       HostOrder = c("Passeriformes", "Carnivora", "Cetartiodactyla", "Cetartiodactyla", "Cetartiodactyla", "Cetartiodactyla"),
                       HostFamily = c("Corvidae", "Felidae", "Cetacea", "Cetacea", "Cetacea", "Cetacea"))
all_syns = rbind(all_syns1, all_syns2)

# combine
assoc = left_join(assoc, all_syns, by=c("Host"="Original"))

# order colummns
assoc = assoc %>%
  rename("Host_Original" = Host,
         "Pathogen_Original" = Parasite,
         "PathogenType" = ParasiteType,
         "DetectionMethod_Original" = DetectionMethod,
         "Pathogen_Harmonised" = PathogenName_Harmonised)
assoc = assoc[ , c("Database", "Pathogen_Original", "Pathogen_Harmonised", "PathogenType", 
                   "Host_Original", "Host_Harmonised", "HostClass", "HostOrder", "HostFamily", "HostSynonyms",
                   "Year", "YearType", "DetectionMethod_Original", "HumanInfective_Any", "DiseaseAgent", "IsZoonotic", "Disease_GIDEON", "Route_GIDEON")]

# issue with year lookup for nucleotide records: set all to NA for now
assoc$Year[ grep("Nucleotide", assoc$YearType) ] = NA

# save final dataset for analysis
write.csv(assoc, "./output/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", row.names=FALSE)

