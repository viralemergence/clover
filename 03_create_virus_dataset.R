



# ================== Combine GMPD2/HP3/EID2 associations with Shaw viruses and standardise detection methods ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("RISmed", "dplyr", "magrittr", "rentrez")



# ================= read datasets ================

# read datasets (associations; shaw mammals/viruses only)
assoc = read.csv("./output/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)

# shaw mammals/viruses and clip pathogen names
shaw = read.csv("./data/Shaw_Database_main.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(HostGroup %in% c("Human", "Ungulates", "Carnivora", "Cetacea", "Rodentia", "Mammalia", "Primates", "Chiroptera")) %>%
  dplyr::filter(Type == "Virus") %>%
  dplyr::mutate(Species = tolower(Species))
foo = strsplit(shaw$Species, " ")
foo = lapply(foo, function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) })
foo = unlist(lapply(foo, paste, collapse=" "))
shaw$Pathogen_Original = foo
shaw$PMID[ is.na(shaw$PMID) ] = ""
shaw$ReferenceComb = paste(shaw$Reference, shaw$Additional, sep="; ")
shaw$ReferenceComb = paste("PMID:", shaw$PMID, "_Citations:", shaw$ReferenceComb, sep="")


# # =============== access publication year for Shaw; first for refs, then scrape PubMed ================
# 
# # extract numerics (year) and trim out second year marker
# years = regmatches(shaw$Reference, gregexpr("[[:digit:]]+", shaw$Reference))
# shaw$YearRef = as.numeric(unlist(lapply(years, "[", 1)))
# 
# # create dataframe of PMIDs to lookup
# shaw$idx = 1:nrow(shaw)
# query_df = shaw %>%
#   dplyr::select(idx, PMID, Pathogen_Original) %>%
#   dplyr::mutate(Database = "pubmed") %>%
#   dplyr::filter(!is.na(PMID))
# 
# # function to lookup create and update dates from PubMed
# pubMedDateLookup = function(ids){
#   
#   # result
#   e = simpleError("lookup error")
#   
#   # lookup 5 attempts then break
#   e = simpleError("test error")
#   for(attempt in 1:5){
#     search = tryCatch(entrez_summary(db="pubmed", id=ids), error=function(e) e)
#     if(class(search)[1] != "simpleError"){ break }
#     Sys.sleep(0.5)
#   }
#   
#   # if throws error
#   if(class(search)[1] == "simpleError"){
#     resx = data.frame(id = ids)
#     resx$PubDate = NA; resx$PubYear = NA
#     resx$Lookup = "Fail"
#     return(resx)
#     
#   } else{
#     resx = data.frame(id = names(search))
#     pubdatenullfilter = function(x){ 
#       lx = x$pubdate
#       if(is.null(lx)){ lx = NA }
#       lx
#     }
#     resx$PubDate = unlist(lapply(search, pubdatenullfilter))
#     resx$PubYear = substr(resx$PubDate, 1, 4)
#     resx$Lookup = "Success"
#   }
#   
#   # return
#   return(resx)
# }
# 
# # batches
# batch = rep(1:1000, each=100); batch = batch[ 1:nrow(query_df) ]
# query_df$batch = batch
# 
# # create filenames
# output_loc = "./output/"
# #save_file = paste(output_loc, "Shaw_Pubmedscrape_28122020_rentrez.csv", sep="")
# 
# # append each new query to csv
# for(i in 1:n_distinct(query_df$batch)){
#   
#   # run query
#   cat(paste(i, "...", sep=""))
#       
#   e = simpleError("lookup error")
#   ids = query_df$PMID[ query_df$batch == unique(query_df$batch)[i] ]
#   lookupx = tryCatch(pubMedDateLookup(ids))
#   lookupx = rename(lookupx, "PMID" = id)
#   lookupx$PMID = as.numeric(as.vector(lookupx$PMID))
#   
#   # combine with eid2 records
#   lookupx = left_join(query_df[ query_df$batch == unique(query_df$batch)[i], ], lookupx,  by="PMID")
#   
#   # initialise file on first iteration, and then append
#   if(class(lookupx)[1] == "simpleError"){ next
#   } else if(i == 1){
#     write.csv(lookupx, save_file, row.names=FALSE)
#   } else{
#     write.table(lookupx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
#   }
#   
#   # sleep system to reduce overload
#   Sys.sleep(0.5)
# }


#' # =============== save host names for lookup via taxize ================
#' 
#' # hostnames, run through taxize script
#' hostnames = shaw %>%
#'   select(HostSpecies.new)
#' write.csv(hostnames, "./output/crossref_temp/Shaw_HostNames_formatching.csv", row.names=FALSE)
#' 


# =============== reconcile pathogen names ==================
# 
# # subset to just pathogen names and check whether name already in associations
# shawp = shaw %>%
#   dplyr::select(Pathogen_Original, Synonym) %>%
#   distinct() %>%
#   dplyr::mutate(InAssoc = ifelse(Pathogen_Original %in% assoc$Pathogen_Harmonised, TRUE, FALSE),
#                 Synonym = tolower(Synonym)) %>%
#   arrange(Pathogen_Original)
# shawp$Pathogen_Harmonised = ""
# shawp$Pathogen_Harmonised[ shawp$InAssoc == TRUE ] = shawp$Pathogen_Original[ shawp$InAssoc == TRUE ] 
# 
# # association names for cross-ref
# assocp = assoc %>%
#   dplyr::filter(PathogenType == "virus") %>%
#   dplyr::select(Pathogen_Harmonised) %>%
#   distinct() %>%
#   arrange(Pathogen_Harmonised)

# save both for cross-ref
# write.csv(shawp, "./output/crossref_temp/Shaw_virusnames_toharmonise.csv", row.names=FALSE)
# write.csv(assocp, "./output/crossref_temp/Assoc_virusnames_forreference.csv", row.names=FALSE)




# ====================== combine all datasets together and save Shaw metadata =======================

# shaw pathogen names reconciled
path = read.csv("./output/crossref_temp/Shaw_virusnames_harmonised_rg_gfa.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(Pathogen_Original, Pathogen_Harmonised)
shaw = left_join(shaw, path, by="Pathogen_Original")

# shaw years including additional column
years1 = regmatches(shaw$Reference, gregexpr("[[:digit:]]+", shaw$Reference))
additional = shaw$Additional
additional[ grep("DOI|doi|http", additional) ] = NA
years2 = regmatches(additional, gregexpr("[[:digit:]]+", additional))
for(i in 1:length(years1)){ years1[[i]] = c(years1[[i]], years2[[i]]) }
shaw$YearRef = as.numeric(as.vector(unlist(lapply(years1, function(x) min(x[ nchar(x) == 4])))))
shaw$YearRef[ shaw$YearRef >2019 & !is.na(shaw$YearRef)] = NA # couple of issues: remove year > 2019

# years scraped from pubmed
# years_pm = read.csv("./output/crossref_temp/Shaw_Pubmedscrape_12112020.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(Lookup_Successful == TRUE) %>%
#   dplyr::select(PMID, Year) %>%
#   dplyr::rename("YearPubMed" = Year) %>%
#   dplyr::filter(!duplicated(PMID))
# shaw = left_join(shaw, years_pm)

# years
years_pm = read.csv("./output/crossref_temp/Shaw_Pubmedscrape_28122020_rentrez.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(PMID, PubYear) %>%
  dplyr::rename("YearPubMed" = PubYear) %>%
  dplyr::filter(!duplicated(PMID)) %>%
  dplyr::mutate(PMID = as.character(PMID))
shaw = left_join(shaw, years_pm)

# how many have years? 19% in database; 90% from PubMed
sum(!is.na(shaw$YearRef))/nrow(shaw)
sum(!is.na(shaw$YearPubMed))/nrow(shaw)
# foo = shaw[ !is.na(shaw$YearPubMed), ]
# sum(foo$YearRef == foo$YearPubMed, na.rm=TRUE) / nrow(foo[ !is.na(foo$YearRef),]) 

# take earliest year of listed and PubMed scrape
shaw$Year = apply(shaw[ , c("YearRef", "YearPubMed")], 1, function(x) min(x[ !is.na(x) ]) )
shaw$Year[ shaw$Year == Inf ] = NA
shaw$YearType = ifelse(shaw$Year == shaw$YearPubMed, "PubMed publication", "Publication date from source database (earliest)")

# taxized hostnames
load("./output/crossref_temp/Shaw_names_matched_list.R")
sortSyns = function(x){
  if(length(x) > 0){ 
    return(  data.frame(Host_Original = x$Original[1],
                        Host_Harmonised = tolower(x$Accepted_name[1]),
                        # HostOrder = x$Selected_order,
                        # HostFamily = x$Selected_family,
                        # HostClass = x$Selected_class,
                        HostSynonyms = tolower(paste(x$Synonyms, collapse=", "))) )
    }
  }
syns = do.call(rbind.data.frame, lapply(synonyms_lookup, sortSyns)) %>%
  dplyr::filter(!duplicated(Host_Original))

# combine with shaw and save shaw harmonised with metadata
shaw$Host_Original = tolower(shaw$HostSpecies.new)
shaw = left_join(shaw, syns)
shaw$HostClass = "Mammalia"
write.csv(shaw, "./output/hostpathogen_harmonised/Shaw_MammalViruses_AllHarmonised_RG_GA.csv", row.names=FALSE)

# subset to relevant columns, combine with assoc, subset to mammals only
shaw_sub = shaw[ , c("Pathogen_Original", "Pathogen_Harmonised", "Type", "Host_Original", "Host_Harmonised", "HostClass",
                     "HostOrder", "HostFamily", "HostSynonyms", "Year", "YearType", "Method", "ReferenceComb")] %>%
  dplyr::rename("PathogenType" = Type, "DetectionMethod_Original" = Method, "Reference" = "ReferenceComb") %>%
  dplyr::mutate(Database = "Shaw",
                ReferenceType = "Publication and/or PMID",
                PathogenType = tolower(PathogenType))
assoc = assoc[ , which(names(assoc) %in% names(shaw_sub)) ]
assoc = assoc[ assoc$HostClass == "Mammalia" & assoc$PathogenType == "virus", ]
assoc = rbind(assoc, shaw_sub)

# note on database version
database_version = data.frame(Database=c("EID2", "Shaw", "GMPD2", "HP3"), DatabaseVersion = c("Wardeh et al. 2015 Sci Data", "Shaw et al. 2020 Mol Ecol", "Stephens et al. 2017 Ecology", "Olival et al. 2017 Nature"))
assoc = left_join(assoc, database_version)



# ------------ manual corrections to pathogen synonymy in complete database ---------------

# 1. "human polyomavirus" just listed as polyomavirus in Olival but 2 strains in Shaw; check Olival assoc and change to either
# 2. "human mastadenovirus a-f" and "human adenovirus a-f" are synonymous: 
# 3. "human enterovirus b and c" -> "enterovirus b - c"
# 3. Shaw lists 3 sub-strains of feline sarcoma virus, but just listed as one in other datasets: have aggregated to one for now

# correction of problems: remove pathogens not classified to genus/sp
assoc = assoc[ !is.na(assoc$Pathogen_Harmonised), ]
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human enterovirus b" ] = "enterovirus b"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human enterovirus c" ] = "enterovirus c"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human mastadenovirus a" ] = "human adenovirus a"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human mastadenovirus b" ] = "human adenovirus b"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human mastadenovirus c" ] = "human adenovirus c"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human mastadenovirus d" ] = "human adenovirus d"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human mastadenovirus e" ] = "human adenovirus e"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human mastadenovirus f" ] = "human adenovirus f"
assoc$Pathogen_Harmonised[ assoc$Pathogen_Harmonised == "human polyomavirus" ] = "human polyomavirus 1"

# hosts that do not harmonise to taxize for shaw; corrected for prev datasets
assoc$Host_Harmonised[ assoc$Host_Original == "felis catus" & assoc$Database == "Shaw" ] = assoc$Host_Harmonised[ assoc$Host_Original == "felis catus" & assoc$Database == "EID2" ][1] 
assoc$HostSynonyms[ assoc$Host_Original == "felis catus" & assoc$Database == "Shaw" ] = assoc$HostSynonyms[ assoc$Host_Original == "felis catus" & assoc$Database == "EID2" ][1] 
assoc$Host_Harmonised[ assoc$Host_Original == "physeter catodon" & assoc$Database == "Shaw" ] = "physeter macrocephalus"
assoc$HostSynonyms[ assoc$Host_Original == "physeter catodon" & assoc$Database == "Shaw" ] = "physeter catodon"
assoc$Host_Harmonised[ assoc$Host_Original == "balaena mysticetus" & assoc$Database == "Shaw" ] = "balaena mysticetus"
assoc$HostSynonyms[ assoc$Host_Original == "balaena mysticetus" & assoc$Database == "Shaw" ] = ""

# fix NAs in synonyms
assoc$HostSynonyms[ assoc$HostSynonyms == "na" ] = ""

# set artiodactyls
assoc$HostOrder[ assoc$HostOrder %in% c("Artiodactyla", "Cetacea") ] = "Cetartiodactyla"

# order, keep only distinct records, and save
assoc = assoc %>%
  dplyr::arrange(HostOrder, Host_Harmonised, Pathogen_Harmonised, Database) %>%
  dplyr::select(Host_Harmonised, Host_Original, HostClass, HostOrder, HostFamily, Pathogen_Harmonised, Pathogen_Original, PathogenType, Database, DatabaseVersion, Reference, ReferenceType, Year, YearType,
                DetectionMethod_Original, HostSynonyms) %>%
  distinct()



# ============= EID2 reports multiple instances of same association using nuccore/PMID; collapse down to one record per year ==============

# combine multiple records from EID2
ei = assoc[ assoc$Database == "EID2", ]

# subset
combx = ei %>%
  group_by(Host_Harmonised, Pathogen_Harmonised, Year, ReferenceType) %>%
  dplyr::summarise(Reference = paste(Reference, collapse="; "))

# combine with eid
ei = ei %>%
  dplyr::select(-Reference) %>%
  distinct() %>%
  left_join(combx)
assoc = rbind(ei, assoc[ assoc$Database != "EID2", ])

# write 
#write.csv(assoc, "./output/Clover_reconciledassociations_v1_20201120.csv", row.names=FALSE)



# ========================= harmonise detection methods ==================================

# assoc
#assoc = read.csv("./output/Clover_reconciledassociations_v1_20201120.csv", stringsAsFactors = FALSE)

# any EID2 entries based on NCBI Nucleotide labelled as "PCR"
assoc$DetectionMethod_Original[ assoc$Database == "EID2" & assoc$Source == "NCBI Nucleotide" ] = "PCR"

# detection methods
dm = assoc[ !duplicated(assoc$DetectionMethod_Original), c("DetectionMethod_Original"), drop=FALSE ]

# none reported
det0 = c("na", "CargoCarrier", "Other")
dm$Detection_NotSpecified = ifelse(dm$DetectionMethod_Original %in% det0 | is.na(dm$DetectionMethod_Original), TRUE, FALSE)

# serology/antibody-based methods
sero = c("bcELISA", "Antibodies", "Serology", "plaque reduction neutralization", "ELISA", "SNT", "ELISA, SNT", "PAGE", "VNT", "PAGE", "hemmaglutination", "virus neutralization", 
         "Direct Fluorescent Antibody Testing", "CF", "VNT, ELISA", "agar gel immunodiffusion", "IFA", "antigen", "Antigens", "Antigen", "Hemagglutination inhibition assay; Neutralization test",
         "PRNT", "NT", "neutralization test", "ELISA and plaque reduction", "PHA", "SouthernBlot", "Antibodies and Isolation")
dm$Detection_Serology = ifelse(dm$DetectionMethod_Original %in% sero, TRUE, FALSE)

# genetic detection methods
gen_det = c("PCR", "Cell culture, PCR", "PCR, Isolation", "DNA RFLP", "EM, DNA (i.e. more than just PCR)", "Pyroseq", "RNA")
dm$Detection_Genetic = ifelse(dm$DetectionMethod_Original %in% gen_det, TRUE, FALSE)

# pathogen direct isolation/observation within/from tissue
# this could be more nuanced
iso_det = c("DirectBlood", "DirectFecal", "DirectOther", "Tissue", "Fecal", "Isolation", "histopath, e-microscopy", "histopath; e microscopy", "isolation", "Cell culture", "micropsy", "Antibodies and Isolation", "Pathology", "Histopathology", "Histology", "PCR, Isolation", "Cell culture, PCR")
dm$Detection_Isolation = ifelse(dm$DetectionMethod_Original %in% iso_det, TRUE, FALSE)

# composite column of highest quality detection
dm$DetectionMethod_Harmonised = "Not specified"
dm$DetectionMethod_Harmonised[ dm$Detection_Serology == TRUE ] = "Antibodies"
dm$DetectionMethod_Harmonised[ dm$Detection_Genetic == TRUE ] = "PCR/Sequencing"
dm$DetectionMethod_Harmonised[ dm$Detection_Isolation == TRUE ] = "Isolation/Observation"

# create detection quality column (simple for now, can update)
# dm$DetectionQuality = 0
# dm$DetectionQuality[ dm$DetectionMethod_Harmonised == "Antibodies" ] = 1
# dm$DetectionQuality[ dm$DetectionMethod_Harmonised %in% c("Isolation/Observation", "PCR/Sequencing") ] = 2

# combine and save for NCBITaxonomy resolving
assoc = left_join(assoc, dm, by="DetectionMethod_Original")
write.csv(assoc, "./output/intermediate_versions/Clover_reconciledassociations_v1_20201120.csv", row.names=FALSE)








############################ OLD CODE #########################



# =============== access publication year for Shaw; first for refs, then scrape PubMed ================

#' # extract numerics (year) and trim out second year marker
#' years = regmatches(shaw$Reference, gregexpr("[[:digit:]]+", shaw$Reference))
#' shaw$YearRef = as.numeric(unlist(lapply(years, "[", 1)))
#' 
#' # create dataframe of PMIDs to lookup
#' shaw$idx = 1:nrow(shaw)
#' query_df = shaw %>%
#'   select(idx, PMID, Pathogen_Original) %>%
#'   dplyr::mutate(Database = "pubmed") %>%
#'   dplyr::filter(!is.na(PMID)) 
#' 
#' # function to scrape pubmed results
#' #' @param x 1:nrow(query_df), i.e. the nth row of query_df
#' 
#' searchNCBI = function(x){
#'   
#'   # print identifier
#'   print(sprintf("Processing: %s", query_df$PMID[x]))
#'   
#'   # run scrape (5 attempts and exit if successful)
#'   e = simpleError("test error")
#'   for(attempt in 1:5){
#'     search = tryCatch(RISmed::EUtilsGet(query_df$PMID[x], type="efetch", db=query_df$Database[x]), error=function(e) e)
#'     if(class(search)[1] != "simpleError"){ break }
#'     Sys.sleep(0.5)
#'   }
#'   
#'   # PubMed
#'   if(query_df$Database[x] == "pubmed"){
#'     
#'     # if not returned successfully
#'     if(class(search)[1] != "Medline"){
#'       res = data.frame(
#'         idx = query_df$idx[x],
#'         PMID = query_df$PMID[x],
#'         Database = query_df$Database[x],
#'         Lookup_Successful = FALSE,
#'         PMID_lookup = NA,
#'         Year = NA,
#'         Journal = NA,
#'         Title = NA
#'       ) 
#'       Sys.sleep(0.5)
#'       return(res)
#'     }
#'     
#'     # if no records exist
#'     if(length(search@PMID) == 0){
#'       res = data.frame(
#'         idx = query_df$idx[x],
#'         PMID = query_df$PMID[x],
#'         Database = query_df$Database[x],
#'         Lookup_Successful = TRUE,
#'         PMID_lookup = "no record",
#'         Year = NA,
#'         Journal = NA,
#'         Title = NA
#'       ) 
#'     } else{
#'       # otherwise return records
#'       res = data.frame(
#'         idx = query_df$idx[x],
#'         PMID = query_df$PMID[x],
#'         Database = query_df$Database[x],
#'         Lookup_Successful = TRUE,
#'         PMID_lookup = search@PMID,
#'         Year = search@YearPubmed,
#'         Journal = search@MedlineTA,
#'         Title = search@ArticleTitle
#'       )
#'     }
#'   }
#'   
#'   # sleep for 0.5 seconds to prevent over-requesting on API (might be able to get away with shorter)
#'   Sys.sleep(0.5)
#'   return(res)
#' }
#' 
#' # create filenames for running scrape
#' output_loc = "./output/"
#' save_file = paste(output_loc, "Shaw_Pubmedscrape_12112020_v2.csv", sep="")
#' 
#' # append each new query to csv
#' for(i in 1:nrow(query_df)){
#' 
#'   # run query
#'   cat(paste(i, "...", sep=""))
#'   e = simpleError("test error")
#'   resx = tryCatch(searchNCBI(i))
#'   
#'   # initialise file on first iteration, and then append
#'   if(class(resx)[1] == "simpleError"){ next
#'   } else if(i == 1){
#'     write.csv(resx, save_file, row.names=FALSE)
#'   } else{
#'     write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
#'   }
#' }
#' 