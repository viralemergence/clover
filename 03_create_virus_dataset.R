



# ================== Combine assoc viruses with Shaw viruses (reconcile names) ==================

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



# =============== save host names for harmonising via taxize ================

# hostnames, run through taxize script
hostnames = shaw %>%
  select(HostSpecies.new)
write.csv(hostnames, "./output/crossref_temp/Shaw_HostNames_formatching.csv", row.names=FALSE)



# =============== reconcile pathogen names ==================

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

# scraped from pubmed
years_pm = read.csv("./output/crossref_temp/Shaw_Pubmedscrape_12112020.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Lookup_Successful == TRUE) %>%
  dplyr::select(PMID, Year) %>%
  dplyr::rename("YearPubMed" = Year) %>%
  dplyr::filter(!duplicated(PMID))
shaw = left_join(shaw, years_pm)

# how many have years? 19% in database; 90% from PubMed
sum(!is.na(shaw$YearRef))/nrow(shaw)
sum(!is.na(shaw$YearPubMed))/nrow(shaw)
foo = shaw[ !is.na(shaw$YearRef), ]
sum(foo$YearRef == foo$YearPubMed, na.rm=TRUE) / nrow(foo[ !is.na(foo$YearRef),]) 

# take earliest year of listed and PubMed scrape
shaw$Year = apply(shaw[ , c("YearRef", "YearPubMed")], 1, function(x) min(x[ !is.na(x) ]) )
shaw$Year[ shaw$Year == Inf ] = NA
shaw$YearType = ifelse(shaw$Year == shaw$YearPubMed, "PubMed", "Author")

# resolved hostnames
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
                     "HostOrder", "HostFamily", "Year", "YearType", "Method")] %>%
  dplyr::rename("PathogenType" = Type, "DetectionMethod_Original" = Method) %>%
  dplyr::mutate(DetectionMethod_Harmonised = NA,
                DetectionQuality_Original = NA,
                DetectionQuality_Harmonised = NA,
                Database = "Shaw",
                PathogenType = tolower(PathogenType))
assoc = assoc[ , which(names(assoc) %in% names(shaw_sub)) ]
assoc = assoc[ assoc$HostClass == "Mammalia" & assoc$PathogenType == "virus", ]
assoc = rbind(assoc, shaw_sub)


# ------------ correct issues with pathogen synonymy in complete database ---------------

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

# save
write.csv(assoc, "./output/Clover_reconciledassociations_v1_20201120.csv", row.names=FALSE)

# some initial viz (exclude nucleotide from year)
p1 = ggplot(assoc[assoc$YearType != "Nucleotide" & assoc$Year <= 2017 & assoc$Host_Harmonised != "homo sapiens", ] ) + 
  geom_histogram(aes(Year), binwidth=1, fill="grey80", col="grey20") + theme_classic() + 
  ggtitle("Total associations reported by year")

a2 = assoc[ assoc$Year <= 2017 & assoc$Host_Harmonised != "homo sapiens", ] %>%
  dplyr::group_by(Host_Harmonised, Pathogen_Harmonised) %>%
  dplyr::summarise(HostOrder = head(HostOrder, 1),
                   Year = min(Year, na.rm=TRUE))
p2 = ggplot(a2) + geom_histogram(aes(Year, fill=HostOrder), binwidth=5) + theme_classic() + 
  ggtitle("Unique host-virus associations (first year reported)") +
  geom_histogram(aes(Year), binwidth=5, col="grey20", fill=NA)

a3 = assoc[ assoc$Year <= 2017 & assoc$Host_Harmonised != "homo sapiens", ] %>%
  group_by(Host_Harmonised) %>%
  dplyr::summarise(HostOrder = head(HostOrder, 1),
                   ViralRichness = n_distinct(Pathogen_Harmonised))
p3 = ggplot(a3[ !is.na(a3$HostOrder), ]) + 
  geom_histogram(aes(ViralRichness, fill=HostOrder), binwidth=5) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("Num host species") + 
  ggtitle("Viral richness by host")

a4 = assoc[ assoc$Year <= 2017 & assoc$Host_Harmonised != "homo sapiens", ] %>%
  group_by(Pathogen_Harmonised) %>%
  dplyr::summarise(HostRichness = n_distinct(Host_Harmonised))
p4 = ggplot(a4) + 
  geom_histogram(aes(HostRichness), binwidth=2, fill="grey80", col="grey20") +
  theme_classic() +
  ylab("Num virus species") + 
  ggtitle("Host richness by virus")

px = gridExtra::grid.arrange(grobs=list(p1, p2, p3, p4), nrow=2, ncol=2, widths=c(1, 1.3))
ggsave(px, file="./output/clover_v1reconciled_initialviz.png", device="png", dpi=300, width=10, height=8, units="in")
