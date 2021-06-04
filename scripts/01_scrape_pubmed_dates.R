

# ============= Scrape NCBI databases (PubMed/Nucleotide) for EID2 and Shaw publication details inc. year using rentrez ============

# splits EID2 into 1 record per accession, accesses date/publication metadata for each accession from rentrez, and saves for later use
# splits Shaw into 1 record per cited reference, accesses date/publication metadata for each accession from rentrez, and saves for later use

# dependencies and basedir
#setwd("C:/Users/roryj/Documents/PhD/202011_clover/repos/clover/")
pacman::p_load("dplyr", "magrittr", "rentrez")




# ================ process EID2 =================

# eid2 for vertebrates
eid2 = read.csv("./data/source_databases/EID2_SpeciesInteractions_Wardeh2015.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Carrier.classification %in% c("Human", "Mammal", "Domestic", "Primate", "Rodent", "Aves", "Fish", "Reptile", "Amphibian")) 

# EID2 data with PMIDs (stored in 'Publications' column)
# 'database' field corresponds to exact name of db in NCBI databases (used to lookup relevant db)
pm = eid2[ eid2$Publications != "", -which(names(eid2) == "Sequences")] %>%
  tidyr::separate_rows(Publications, sep=";") %>%
  dplyr::rename("id" = Publications)
pm$Database = "pubmed"

# data with NCBI nucleotide IDs
# set database to nuccore
nuc = eid2[ eid2$Sequences != "", -which(names(eid2) == "Publications")] %>%
  tidyr::separate_rows(Sequences, sep=";") %>%
  dplyr::rename("id" = Sequences)
nuc$Database = "nuccore"





# ----------- run date lookup in pubmed ----------------

# function to lookup create and update dates from PubMed
pubMedDateLookup = function(ids){
  
  # result
  e = simpleError("lookup error")
  
  # lookup 5 attempts then break
  for(attempt in 1:5){
    search = tryCatch(rentrez::entrez_summary(db="pubmed", id=ids), error=function(e) e)
    if(class(search)[1] != "simpleError"){ break }
    Sys.sleep(0.5)
  }
  
  # if throws error
  if(class(search)[1] == "simpleError"){
    resx = data.frame(id = ids)
    resx$PubDate = NA; resx$PubYear = NA
    resx$Lookup = "Fail"
    return(resx)
    
  } else{
    
    resx = data.frame(id = names(search)) 
    
    # publication year
    pubdatenullfilter = function(x){
      lx = x$pubdate
      if(is.null(lx)){ lx = NA }
      lx
      }
    resx$PubDate = unlist(lapply(search, pubdatenullfilter))
    resx$PubYear = substr(resx$PubDate, 1, 4)
    
    # citation info
    lookupCitation = function(x){
      if(is.null(x$pubdate)){ return(NA)
      } else{
          cite = paste(x$sortfirstauthor, substr(x$pubdate, 1, 4), x$source, x$volume, sep=", ")
          return(cite)
        }
    }
    resx$Citation = unlist(lapply(search, lookupCitation))
    
    resx$Lookup = "Success"
  }
  
  # return
  return(resx)
}

# batches
batch = rep(1:1000, each=225); batch = batch[ 1:nrow(pm) ]
pm$batch = batch

# create filenames
output_loc = "./output/intermediate/"
save_file = paste(output_loc, "EID2_PubMed_DateScrape_26012021.csv", sep="")

# append each new query to csv
for(i in 1:n_distinct(pm$batch)){
  
  # run query
  cat(paste(i, "...", sep=""))
  e = simpleError("lookup error")
  ids = unique(pm$id[ pm$batch == unique(pm$batch)[i] ])
  lookupx = tryCatch(pubMedDateLookup(ids))
  
  # combine with eid2 records
  lookupx = left_join(pm[ pm$batch == unique(pm$batch)[i], ], lookupx,  by="id")
  
  # initialise file on first iteration, and then append
  if(class(lookupx)[1] == "simpleError"){ 
    print(paste("failed on batch: ", unique(pm$batch)[i], sep=""))
    next
    
  } else if(i == 1){
    write.csv(lookupx, save_file, row.names=FALSE)
    
  } else{
    write.table(lookupx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
  
  # sleep system to reduce overload
  Sys.sleep(0.75)
}



# --------------- run date lookup in NCBI Nucleotide ------------------

# function to lookup create and update dates from Nucleotide
nuccoreDateLookup = function(ids){
  
  # result
  resx = data.frame(id = ids)
  e = simpleError("lookup error")
  
  # lookup
  for(attempt in 1:5){
    search = tryCatch(entrez_summary(db="nuccore", id=ids, version="2.0"), error=function(e) e)
    if(class(search)[1] != "simpleError"){ break }
    Sys.sleep(0.5)
  }
  
  # if throws error
  if(class(search)[1] == "simpleError"){
    resx$CreateDate = NA; resx$CreateYear = NA
    resx$UpdateDate = NA; resx$UpdateYear = NA
    resx$Lookup = "Fail"
    return(resx)
  } else{
    resx$CreateDate = unlist(lapply(search, function(x) x$createdate))
    resx$CreateYear = substr(resx$CreateDate, 1, 4)
    resx$UpdateDate = unlist(lapply(search, function(x) x$updatedate))
    resx$UpdateYear = substr(resx$UpdateDate, 1, 4)
    resx$Lookup = "Success"
  }
  
  # return
  return(resx)
}

# batches
batch = rep(1:1000, each=225); batch = batch[ 1:nrow(nuc) ]
nuc$batch = batch

# create filenames
output_loc = "./output/intermediate/"
save_file = paste(output_loc, "EID2_Nuccore_DateScrape_26012021.csv", sep="")

# append each new query to csv
for(i in 1:n_distinct(nuc$batch)){
  
  # run query
  cat(paste(i, "...", sep=""))
  e = simpleError("lookup error")
  ids = unique(nuc$id[ nuc$batch == unique(nuc$batch)[i] ])
  lookupx = tryCatch(nuccoreDateLookup(ids))
  
  # combine with eid2 records
  lookupx = left_join(nuc[ nuc$batch == unique(nuc$batch)[i], ], lookupx,  by="id")
  
  # initialise file on first iteration, and then append
  if(class(lookupx)[1] == "simpleError"){ next
  } else if(i == 1){
    write.csv(lookupx, save_file, row.names=FALSE)
  } else{
    write.table(lookupx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
  
  # sleep system to reduce overload
  Sys.sleep(0.5)
}



# -------------- combine and save -------------------

# pubmed and nuccore year scrapes
scr1 = read.csv("./output/intermediate/EID2_PubMed_DateScrape_26012021.csv", stringsAsFactors = FALSE)
scr2 = read.csv("./output/intermediate/EID2_Nuccore_DateScrape_26012021.csv", stringsAsFactors = FALSE) 

# harmonise column names
scr1 = scr1 %>%
  dplyr::rename("PublicationDate" = PubDate, "PublicationYear" = PubYear, "ReferenceText"=Citation) %>%
  #dplyr::mutate("YearType" = "PubMed publication") %>%
  dplyr::mutate(ReleaseDate=NA, ReleaseYear=NA) %>%
  dplyr::select(-batch, -Lookup)
scr2 = scr2 %>%
  dplyr::rename("ReleaseDate" = CreateDate, "ReleaseYear" = CreateYear) %>%
  #dplyr::mutate("YearType" = "Nucleotide create date") %>%
  dplyr::mutate(ReferenceText = "NCBI Nucleotide", PublicationDate=NA, PublicationYear=NA) %>%
  dplyr::select(-batch, -Lookup, -UpdateDate, -UpdateYear)
scr = rbind(scr1, scr2)

# 
scr = scr %>%
  dplyr::rename("CitationID"=id, "CitationIDType"=Database) %>%
  dplyr::mutate(CitationIDType = ifelse(CitationIDType == "pubmed", "PMID", "NCBI Nucleotide"))

# remove duplicate records and save as partitioned EID2 in data folder
scr = distinct(scr)
write.csv(scr, "./data/source_databases/EID2_SpeciesInteractions_Rentrez_CLOVERT.csv", row.names=FALSE)






# ======================= Process Shaw database and extract metadata =======================

# shaw mammals/viruses and clip pathogen names
shaw = read.csv("./data/source_databases/Shaw_Database_main.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(Species = tolower(Species))

# clip "virus" off pathogen names (to match other source datasets)
foo = strsplit(shaw$Species, " ")
foo = lapply(foo, function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) })
foo = unlist(lapply(foo, paste, collapse=" "))
shaw$Pathogen_Original = foo



# ------------- split citations across multiple lines ----------------

result = data.frame()
for(i in 1:nrow(shaw)){
  
  # ith record
  print(i)
  shawi = shaw[ i, ]
  
  # reference and PMID
  cit_ref = shawi$Reference
  cit_pmid = shawi$PMID[ !is.na(shawi$PMID) ]
  
  # additional citations, and correct erroneous splits
  cit_add = strsplit(shawi$Additional, ";")
  process_erroneous_splits = function(x){
    if(any(x == 2)) print("corrected")
    x = x[ x != 2 ]
    x[ x == "http://dx.doi.org/10.1647/1082-6742(2001)015[0204:CTIIAR]2.0.CO" ] = "http://dx.doi.org/10.1647/1082-6742(2001)015[0204:CTIIAR]2.0.CO;2"
    return(x)
  }
  cit_add = lapply(cit_add, process_erroneous_splits)
  cit_add = unlist(cit_add)
  
  # split into PMID and non-PMID additionals
  cit_add_pmid = cit_add[ !is.na(as.numeric(cit_add)) ]
  cit_add_text = cit_add[ is.na(as.numeric(cit_add)) ]
  
  # create expanded shaw with explicit PMID and text citation fields
  shawi_r = shawi %>%
    dplyr::mutate(Citation_PMID = "",
                  Citation_Text = "")
  
  # PMID 
  shawi_pmid = shawi_r[ rep(1, length(cit_pmid)), ] 
  if(nrow(shawi_pmid) > 0){
    shawi_pmid = shawi_pmid %>%
      dplyr::mutate(Citation_PMID = cit_pmid,
                    Citation_Text = NA)
  }
  
  # "Additional" that are PMIDs
  shawi_pmid_add = shawi_r[ rep(1, length(cit_add_pmid)), ] 
  if(nrow(shawi_pmid_add) > 0){
    shawi_pmid_add = shawi_pmid_add %>%
      dplyr::mutate(Citation_PMID = cit_add_pmid,
                    Citation_Text = NA)
  }
  
  # If reference provided but PMID is not, create field
  if(nrow(shawi_pmid)==0 & cit_ref != ""){
    shawi_text_add = shawi_r %>%
      dplyr::mutate(Citation_PMID = NA,
                    Citation_Text = cit_ref)
  
  # Otherwise add ref text from additional field
  } else{
    
    shawi_text_add = shawi_r[ rep(1, length(cit_add_text)), ] 
    if(nrow(shawi_text_add) > 0){
      shawi_text_add = shawi_text_add %>%
        dplyr::mutate(Citation_PMID = NA,
                      Citation_Text = cit_add_text)
    }
  }
  
  # combine
  resi = rbind(shawi_pmid, shawi_pmid_add, shawi_text_add)
  
  # if there are no fields for either PMID or Additional, populate with "Reference" (handful like this)
  if(nrow(resi) == 0){
    resi = shawi_r %>%
      dplyr::mutate(Citation_PMID = NA,
                    Citation_Text = cit_ref)
  }
  
  # set any "" to NA in citation_text and add identifier
  resi$Citation_Text[ resi$Citation_Text == "" ] = NA
  
  # add to df
  result = rbind(result, resi)
}

shaw$idx = shaw$X
write.csv(result, "./output/intermediate/Shaw_citationsplit_14022021.csv", row.names=FALSE)
#shaw = read.csv("./output/intermediate/Shaw_citationsplit_14022021.csv")



# --------------- lookup metadata for PMIDs via rentrez --------------------

# create query dataframe
query_df = data.frame(id = as.character(unique(shaw$Citation_PMID))) %>%
  dplyr::filter(!is.na(id) & id != "NA")
query_df$identifier = 1:nrow(query_df)

# batches
batch = rep(1:1000, each=250); batch = batch[ 1:nrow(query_df) ]
query_df$batch = batch

# create filenames
output_loc = "./output/intermediate/"
save_file = paste(output_loc, "Shaw_PubMed_DateScrape_260012021.csv", sep="")

# append each new query to csv
for(i in 1:n_distinct(query_df$batch)){
  
  # run query
  cat(paste(i, "...", sep=""))
  e = simpleError("lookup error")
  ids = unique(query_df$id[ query_df$batch == unique(query_df$batch)[i] ])
  lookupx = tryCatch(pubMedDateLookup(ids))
  
  # combine with eid2 records
  lookupx = left_join(query_df[ query_df$batch == unique(query_df$batch)[i], ], lookupx,  by="id")
  
  # initialise file on first iteration, and then append
  if(class(lookupx)[1] == "simpleError"){ next
  } else if(i == 1){
    write.csv(lookupx, save_file, row.names=FALSE)
  } else{
    write.table(lookupx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
  
  # sleep system to reduce overload
  Sys.sleep(0.5)
}



# ------------------ write text citations for cross-checking ------------------
 
# all text citations
text_cites = data.frame(idx = shaw$idx, Citation_Text = shaw$Citation_Text) %>%
  dplyr::filter(!is.na(Citation_Text))
#write.csv(text_cites, "./dev/output/intermediate/Shaw_text_cites.csv", row.names=FALSE)

# doi cites: future task to extract publication dates from DOI
# doi_cites = text_cites[ grep("doi|DOI", text_cites$Citation_Text), ]



# ------------------ combine all citation cross-ref info into full df ------------------

# lookup info
pm = read.csv("./output/intermediate/Shaw_PubMed_DateScrape_260012021.csv") %>%
  dplyr::select(id, PubYear, Citation) %>%
  dplyr::filter(!is.na(PubYear)) %>%
  dplyr::rename("CitationID"=id, "ReferenceText"=Citation, "PublicationYear"=PubYear) %>%
  dplyr::mutate(CitationID = as.numeric(CitationID), CitationIDType = "PMID")
tc = read.csv("./output/intermediate/Shaw_text_cites_rg2.csv") %>%
  dplyr::filter(!is.na(Citation_Text)) %>%
  dplyr::select(-to_update)

# combine with shaw
shaw = dplyr::rename(shaw, "CitationID" = Citation_PMID)
shaw = left_join(shaw, pm, by="CitationID")
shaw = left_join(shaw, tc, by=c("idx", "Citation_Text"))

# standardise columns
shaw$PublicationYear[ !is.na(shaw$YearRef) ] = shaw$YearRef[ !is.na(shaw$YearRef) ]
shaw$ReferenceText[ !is.na(shaw$Citation_Text) ] = shaw$Citation_Text[ !is.na(shaw$Citation_Text) ]
shaw = dplyr::select(shaw, -YearRef, -Citation_Text, -Reference, -PMID, -Additional, -X)
shaw = shaw %>% dplyr::rename("AssocID"=idx)

# save Shaw for later harmonisation
write.csv(shaw, "./data/source_databases/Shaw_Database_Rentrez_CLOVERT.csv", row.names=FALSE)




# ------------------ Save Shaw host names for taxize reconciliation ------------------

# taxized hostnames from virus scrape (during CLOVER mammals build)
load("./data/hostsynonyms_taxize/Shaw_names_matched_list.R")
sortSyns = function(x){
  if(length(x) > 0){ 
    return(  data.frame(Host_Original = x$Original[1],
                        Host_Harmonised = tolower(x$Accepted_name[1]),
                        HostOrder = x$Selected_order,
                        HostFamily = x$Selected_family,
                        HostClass = x$Selected_class,
                        HostSynonyms = tolower(paste(x$Synonyms, collapse=", "))) )
  }
}
syns = do.call(rbind.data.frame, lapply(synonyms_lookup, sortSyns)) %>%
  dplyr::filter(!duplicated(Host_Original))

# host names that are not already looked up
# hostnames = shaw %>%
#   dplyr::select(HostSpecies.new, HostOrder, HostGroup) %>%
#   dplyr::filter(! tolower(HostSpecies.new) %in% syns$Host_Original) %>%
#   dplyr::filter(!duplicated(HostSpecies.new))
# write.csv(hostnames, "./output/crossref_temp/Shaw_HostNames_CLOVERT_fortaxize.csv", row.names=FALSE)

load("./data/hostsynonyms_taxize/Shaw_names_matched_list_CLOVERT.R")
sortSyns = function(x){
  if(length(x) > 0){ 
    return(  data.frame(Host_Original = x$Original[1],
                        Host_Harmonised = tolower(x$Accepted_name[1]),
                        HostOrder = x$Selected_order,
                        HostFamily = x$Selected_family,
                        HostClass = x$Selected_class,
                        HostSynonyms = tolower(paste(x$Synonyms, collapse=", "))) )
  }
}
syns_clovert = do.call(rbind.data.frame, lapply(synonyms_lookup, sortSyns)) %>%
  dplyr::filter(!duplicated(Host_Original))

# combine together
syns_all = rbind(syns, syns_clovert)
syns_all$HostSynonyms[ syns_all$HostSynonyms == "na" ] = ""
write.csv(syns_all, "./data/hostsynonyms_taxize/Shaw_HostNames_taxize_CLOVERT_26012021.csv", row.names=FALSE)



