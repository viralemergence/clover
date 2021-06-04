


# ================== Compile EID2, GMPD2 and HP3 with harmonised pathogen names and host synonyms ==================

# Process, reconcile structure and combine EID2, GMPD2 and HP3 into initial associations dataset
# Pathogen names initially reconciled across databases manually; both will be matched to the NCBI Taxonomy database using 'taxize'

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/repos/clover/")
pacman::p_load("dplyr", "magrittr")



# ================== Compile and process EID2 ======================

# eid2 cross-referenced to pubmed/nucleotide for publication details
eid = read.csv( "./data/source_databases/EID2_SpeciesInteractions_Rentrez_CLOVERT.csv", stringsAsFactors = FALSE)

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
                  Parasite_og = tolower(eid$Cargo),
                  Database = "EID2",
                  PublicationDate = eid$PublicationDate,
                  PublicationYear = eid$PublicationYear,
                  ReleaseDate = eid$ReleaseDate,
                  ReleaseYear = eid$ReleaseYear,
                  CitationID = eid$CitationID,
                  CitationIDType = eid$CitationIDType,
                  ReferenceText = eid$ReferenceText,
                  stringsAsFactors=F)

# subset to only endoparasites
eid2 = eid2[ eid2$ParasiteType %in% c("bacteria", "helminth", "virus", "fungi", "protozoa", "others"), ]
eid2 = eid2[ order(eid2$HostType, eid2$Host), ]

# remove non-animal hosts
# n.b. humans were previously removed; have harmonised viruses but need to update all other pathogen types for full database
eid2 = eid2[ !eid2$HostType %in% c("Higher Plants", "Protozoa", "Bacteria", "Others", "Fungi", "Bryozoa", "Other Plant", "Green Algae", "Methanothermobacter"), ]

# EID2 pathogen names harmonised (to match to other databases)
assoc = read.csv("./data/pathogennames_harmonised_rg/EID2_pathogendata_allpathogens_harmonised_20210407.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Pathogen_Original, PathogenName_Harmonised, PathogenName_Harm2, PathogenType)
eid2 = left_join(eid2, assoc, by=c("Parasite"="Pathogen_Original")) 
eid2$ParasiteType = eid2$PathogenType
eid2 = eid2 %>% dplyr::select(-PathogenType)
eid2$DetectionMethod = "CargoCarrier"

# save
write.csv(eid2, "./output/hostpathogen_reconciled/EID2_CLOVERT_Harmonised_Jan2021.csv", row.names=FALSE)



# ===================== Compile and process GMPD2 =========================

# main GMPD2 dataset
gmpd2 = read.csv("./data/source_databases/GMPD2_main.csv", encoding = "latin1", stringsAsFactors = F)

# extract numerics (year) and trim out second year marker
years = regmatches(gmpd2$Citation, gregexpr("[[:digit:]]+", gmpd2$Citation))
gmpd2$PublicationYear = as.numeric(unlist(lapply(years, "[", 1)))

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
                   PublicationYear = gmpd2$PublicationYear,
                   Database = "GMPD2",
                   ReferenceText = gmpd2$Citation,
                   HostsSampled = gmpd2$HostsSampled,
                   Seroprevalence = gmpd2$Prevalence,
                   Location = gmpd2$LocationName,
                   Longitude = gmpd2$Longitude,
                   Latitude = gmpd2$Latitude,
                   HostSex = gmpd2$HostSex,
                   HostAge = gmpd2$HostAge,
                   DetectionMethod = gmpd2$SamplingType)
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
gmpd2$Parasite_og = gmpd2$Parasite
gmpd2$Parasite = path

# remove records with pathogen name 'abolished'
# and remove ectoparasites
gmpd2 = gmpd2[ !gmpd2$Parasite %in% c("abolished", "no binomial name"), ]
gmpd2 = gmpd2[ gmpd2$ParasiteType %in% c("helminth", "protozoa", "virus", "fungi", "bacteria"), ]

# harmonised pathogen names
assoc = read.csv("./data/pathogennames_harmonised_rg/GMPD2_pathogendata_allpathogens_harmonised_20210407.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Pathogen_Original, PathogenName_Harmonised, PathogenName_Harm2, PathogenType)
  
# get data on missing pathogens from eid2 (some pathogen names metadata missing in updated pipeline - still in harmonised EID2 file)
missing = gmpd2[ which(!gmpd2$Parasite %in% assoc$Pathogen_Original), ]
missingrecs = read.csv("./data/pathogennames_harmonised_rg/EID2_pathogendata_allpathogens_harmonised_20210407.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Pathogen_Original, PathogenName_Harmonised, PathogenName_Harm2, PathogenType) %>%
  dplyr::filter(Pathogen_Original %in% missing$Parasite)
assoc = rbind(assoc, missingrecs)

# combine and save
gmpd2 = left_join(gmpd2, assoc, by=c("Parasite" = "Pathogen_Original"))
write.csv(gmpd2, "./output/hostpathogen_reconciled/GMPD2_CLOVERT_Harmonised_Jan2021.csv", row.names=FALSE)



# ======================== Compile and process HP3 (Olival) ============================

# host-pathogen associations based on both serology and infection status
oliv_assoc = read.csv("./data/source_databases/HP3_associations.csv", stringsAsFactors = F) %>%
  dplyr::select(-WildDomInReference)
oliv_vir = read.csv("./data/source_databases/HP3_viruses.csv", stringsAsFactors = F)
oliv_mam = read.csv("./data/source_databases/HP3_hosts.csv", stringsAsFactors = F)

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

# get numeric years and take earliest reported year (multiple years do not specify different detection methods, and many are secondary sources)
years = regmatches(hp3$Reference, gregexpr("[[:digit:]]+", hp3$Reference))
yy = unlist(lapply(years, function(x) min(as.numeric(x))))
yy[ yy == Inf ] = NA
yy[ grep("http|414", hp3$Reference) ] = NA # anything with http or strange date formatting
hp3$PublicationYear = yy

# rename columns to standardise with other datasets
names(hp3) = c("Parasite", "Host", "DetectionMethod", "DetectionQuality", "ReferenceText", "HostOrder", "HostFamily", "HostGenus", "HostSpecies", 
               "ParasiteOrder", "ParasiteFamily", "ParasiteSubfamily", "ParasiteGenus", "SSorDS", "DNAorRNA", "Envelope", "Segmented", 
               "CytoReplic", "GenomeAveLength", "ReverseZoonosis", "IsZoonotic_HP3", "IsZoonotic_Strict_HP3", "PublicationYear")
hp3$ParasiteType = "virus"
hp3$Database = "HP3"

# remove 'virus' from pathogen names to standardise with other databases
hp3$Parasite_og = hp3$Parasite
hp3$Parasite = unlist(lapply(lapply(strsplit(hp3$Parasite, " "), function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) }), paste, collapse=" "))

# harmonised pathogen names and definitions
assoc = read.csv("./data/pathogennames_harmonised_rg/Olival_pathogendata_allpathogens_harmonised_20210407.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Pathogen_Original, PathogenName_Harmonised, PathogenName_Harm2, PathogenType)
hp3 = left_join(hp3, assoc, by=c("Parasite"="Pathogen_Original"))

# save
write.csv(hp3, "./output/hostpathogen_reconciled/HP3_CLOVERT_Harmonised_Jan2021.csv", row.names=FALSE)



# ===================== Combine and standardise definitions and host names ====================

# read in files
eid2 = read.csv("./output/hostpathogen_reconciled/EID2_CLOVERT_Harmonised_Jan2021.csv", stringsAsFactors = FALSE) 
gmpd2 = read.csv("./output/hostpathogen_reconciled/GMPD2_CLOVERT_Harmonised_Jan2021.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(PublicationDate=NA, ReleaseDate=NA, ReleaseYear=NA, CitationID=NA, CitationIDType=NA)
hp3 = read.csv("./output/hostpathogen_reconciled/HP3_CLOVERT_Harmonised_Jan2021.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(PublicationDate=NA, ReleaseDate=NA, ReleaseYear=NA, CitationID=NA, CitationIDType=NA)

# combine keeping columns present in all 3 databases
cols = Reduce(intersect, list(names(eid2), names(gmpd2), names(hp3)))
assoc = rbind(eid2[ , cols], gmpd2[ , cols], hp3[ , cols])

# order colummns
assoc = assoc %>%
  rename("Host_Original" = Host,
         "Pathogen_Original" = Parasite,
         "Pathogen_Orig2"=Parasite_og,
         "PathogenType" = ParasiteType,
         "DetectionMethod_Original" = DetectionMethod,
         "Pathogen_Harmonised" = PathogenName_Harmonised,
         "Pathogen_Harm2"=PathogenName_Harm2)
assoc = assoc[ , c("Database", "ReferenceText", "CitationID", "CitationIDType", 
                   "Pathogen_Original", "Pathogen_Orig2", "Pathogen_Harmonised", "Pathogen_Harm2", "PathogenType", 
                   "Host_Original", 
                   #"Host_Harmonised", "HostClass", "HostOrder", "HostFamily", "HostSynonyms",
                   "PublicationYear", "PublicationDate", "ReleaseYear", "ReleaseDate", "DetectionMethod_Original")]

# save dataset for combining with Shaw in next script
write.csv(assoc, "./output/hostpathogen_reconciled/EID2GMPD2HP3_Associations_CLOVERT_Harmonised_Jan2021.csv", row.names=FALSE)

