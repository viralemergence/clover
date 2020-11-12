



# ================== Initial looking at combining Shaw with harmonised data ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("RISmed", "dplyr", "magrittr")



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



# =============== access publication year for Shaw; first for refs, then scrape PubMed ================

# extract numerics (year) and trim out second year marker
years = regmatches(shaw$Reference, gregexpr("[[:digit:]]+", shaw$Reference))
shaw$YearRef = as.numeric(unlist(lapply(years, "[", 1)))

# create dataframe of PMIDs to lookup
shaw$idx = 1:nrow(shaw)
query_df = shaw %>%
  select(idx, PMID, Pathogen_Original) %>%
  dplyr::mutate(Database = "pubmed") %>%
  dplyr::filter(!is.na(PMID)) 

# function to scrape pubmed results
#' @param x 1:nrow(query_df), i.e. the nth row of query_df

searchNCBI = function(x){
  
  # print identifier
  print(sprintf("Processing: %s", query_df$PMID[x]))
  
  # run scrape (5 attempts and exit if successful)
  e = simpleError("test error")
  for(attempt in 1:5){
    search = tryCatch(RISmed::EUtilsGet(query_df$PMID[x], type="efetch", db=query_df$Database[x]), error=function(e) e)
    if(class(search)[1] != "simpleError"){ break }
    Sys.sleep(0.5)
  }
  
  # PubMed
  if(query_df$Database[x] == "pubmed"){
    
    # if not returned successfully
    if(class(search)[1] != "Medline"){
      res = data.frame(
        idx = query_df$idx[x],
        PMID = query_df$PMID[x],
        Database = query_df$Database[x],
        Lookup_Successful = FALSE,
        PMID_lookup = NA,
        Year = NA,
        Journal = NA,
        Title = NA
      ) 
      Sys.sleep(0.5)
      return(res)
    }
    
    # if no records exist
    if(length(search@PMID) == 0){
      res = data.frame(
        idx = query_df$idx[x],
        PMID = query_df$PMID[x],
        Database = query_df$Database[x],
        Lookup_Successful = TRUE,
        PMID_lookup = "no record",
        Year = NA,
        Journal = NA,
        Title = NA
      ) 
    } else{
      # otherwise return records
      res = data.frame(
        idx = query_df$idx[x],
        PMID = query_df$PMID[x],
        Database = query_df$Database[x],
        Lookup_Successful = TRUE,
        PMID_lookup = search@PMID,
        Year = search@YearPubmed,
        Journal = search@MedlineTA,
        Title = search@ArticleTitle
      )
    }
  }
  
  # sleep for 0.5 seconds to prevent over-requesting on API (might be able to get away with shorter)
  Sys.sleep(0.5)
  return(res)
}

# create filenames for running scrape
output_loc = "./output/"
save_file = paste(output_loc, "Shaw_Pubmedscrape_12112020_v2.csv", sep="")

# append each new query to csv
for(i in 1:nrow(query_df)){

  # run query
  cat(paste(i, "...", sep=""))
  e = simpleError("test error")
  resx = tryCatch(searchNCBI(i))
  
  # initialise file on first iteration, and then append
  if(class(resx)[1] == "simpleError"){ next
  } else if(i == 1){
    write.csv(resx, save_file, row.names=FALSE)
  } else{
    write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
}





# =============== harmonise names ==================

# subset to just pathogen names and check whether name already in associations
shawp = shaw %>%
  dplyr::select(Pathogen_Original, Synonym) %>%
  distinct() %>%
  dplyr::mutate(InAssoc = ifelse(Pathogen_Original %in% assoc$Pathogen_Harmonised, TRUE, FALSE),
                Synonym = tolower(Synonym)) %>%
  arrange(Pathogen_Original)
shawp$Pathogen_Harmonised = ""
shawp$Pathogen_Harmonised[ shawp$InAssoc == TRUE ] = shawp$Pathogen_Original[ shawp$InAssoc == TRUE ] 

# association names for cross-ref
assocp = assoc %>%
  dplyr::filter(PathogenType == "virus") %>%
  dplyr::select(Pathogen_Harmonised) %>%
  distinct() %>%
  arrange(Pathogen_Harmonised)

# save both for cross-ref
write.csv(shawp, "./output/crossref_temp/Shaw_virusnames_toharmonise.csv", row.names=FALSE)
write.csv(assocp, "./output/crossref_temp/Assoc_virusnames_forreference.csv", row.names=FALSE)



  