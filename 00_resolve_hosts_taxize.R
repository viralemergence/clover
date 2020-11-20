
# taxize/GBIF
library(taxize)
library(rgbif)
library(doParallel)
library(plyr)

# function to find and resolve taxonomic synonyms based on Encyclopedia of Life
findSyns2 = function(x){
  
  # get specific species name
  taxname = hosts_vec[x]
  
  # print progress
  print(paste("Processing:", taxname, sep=" "))
  
  # phyla to consider
  phyla = c("Chordata","Arthropoda","Gastropoda", "Mollusca")
  
  # (1) resolve misspellings
  taxname_resolved = gnr_resolve(taxname, with_canonical_ranks = TRUE)$matched_name2[1]
  if(!is.null(taxname_resolved)){ if(length(strsplit(taxname_resolved, " ", fixed=TRUE)[[1]]) == 2 ){ taxa = taxname_resolved }}
  if(!is.null(taxname_resolved)){ if(length(strsplit(taxname_resolved, " ", fixed=TRUE)[[1]]) > 2 ){ taxa = paste(strsplit(taxname_resolved, " ", fixed=TRUE)[[1]][1:2], collapse=" ")} }
  
  # if taxa == NA, return list with nothing defined 
  if(is.na(taxa)){   if(class(syns)[1] == 'simpleError'){ return(data.frame(Original=taxname, Submitted=taxname_resolved, Accepted_name=NA, Selected_family=NA, Selected_order=NA, Selected_class=NA, Synonyms=NA))} }
  
  # (2) remove sub-species categorisations and set 'genus' and 'species' variables
  genus = NULL
  if(length(strsplit(taxa, " ", fixed=TRUE)[[1]]) %in% c(2,3)){ genus = strsplit(taxa," ",fixed=TRUE)[[1]][1]; species = strsplit(taxa," ",fixed=TRUE)[[1]][2] }
  if(length(strsplit(taxa, "_", fixed=TRUE)[[1]]) %in% c(2,3)){ genus = strsplit(taxa,"_",fixed=TRUE)[[1]][1]; species = strsplit(taxa,"_",fixed=TRUE)[[1]][2] }
  if(length(strsplit(taxa, " ", fixed=TRUE)[[1]]) >3 | length(strsplit(taxa, "_" , fixed=TRUE)[[1]][1]) > 3){ return("name error") }
  if(is.null(genus)){ genus = taxa; species = NA }
  
  # (3) use genus to lookup family, order, class
  syns = tryCatch( name_lookup(genus)$data, error = function(e) e)
  if(class(syns)[1] == 'simpleError'){ return(data.frame(Original=taxname, Submitted=taxa, Accepted_name=NA, Selected_family=NA, Selected_order=NA, Selected_class=NA, Synonyms=NA))}
  
  # for cases where the lookup does not find a phylum within the specified range
  if(all(! syns$phylum %in% phyla)){
    fam1 = syns$family[ !is.na(syns$family) & !is.na(syns$phylum) ]
    order1 = syns$order[ !is.na(syns$family) & !is.na(syns$phylum) ]
    class1 = syns$class[ !is.na(syns$family) & !is.na(syns$phylum) ]
    datfam = data.frame(fam1=fam1, order=1:length(fam1), order1=order1, class1=class1)
    # select highest frequency fam/class/order combo
    fam2 = as.data.frame( table(datfam[ , c(1,3,4)]) )
    family2 = as.vector(fam2[ fam2$Freq==max(fam2$Freq, na.rm=TRUE), "fam1"] ) 
    order2 = as.vector(fam2[ fam2$Freq==max(fam2$Freq, na.rm=TRUE), "order1"] )
    class2 = as.vector(fam2[ fam2$Freq==max(fam2$Freq, na.rm=TRUE), "class1"] )
    if(length(fam2) > 1){
      datfam2 = datfam[datfam$fam1 %in% family2, ]
      family2 = as.vector(datfam2[datfam2$order == min(datfam2$order, na.rm=TRUE), "fam1"])
      order2 = as.vector(datfam2[datfam2$order == min(datfam2$order, na.rm=TRUE), "order1"])
      class2 = as.vector(datfam2[datfam2$order == min(datfam2$order, na.rm=TRUE), "class1"])
    }
  } else {	# for everything else
    fam1 = syns$family[ !is.na(syns$family) & !is.na(syns$phylum) & (syns$phylum %in% phyla) ]
    order1 = syns$order[ !is.na(syns$family) & !is.na(syns$phylum) & (syns$phylum %in% phyla) ]
    class1 = syns$class[ !is.na(syns$family) & !is.na(syns$phylum) & (syns$phylum %in% phyla) ]
    datfam = data.frame(fam1=fam1, order=1:length(fam1), order1 = order1, class1=class1)
    # select highest frequency fam/class/order combo
    fam2 = as.data.frame( table(datfam[ , c(1,3,4)]) )
    family2 = as.vector(fam2[ fam2$Freq==max(fam2$Freq, na.rm=TRUE), "fam1"] ) 
    order2 = as.vector(fam2[ fam2$Freq==max(fam2$Freq, na.rm=TRUE), "order1"] )
    class2 = as.vector(fam2[ fam2$Freq==max(fam2$Freq, na.rm=TRUE), "class1"] )
    # select highest in list if more than one max
    if(length(family2) > 1){
      datfam2 = datfam[datfam$fam1 %in% family2, ]
      family2 = as.vector(datfam2[datfam2$order == min(datfam2$order, na.rm=TRUE), "fam1"])
      order2 = as.vector(datfam2[datfam2$order == min(datfam2$order, na.rm=TRUE), "order1"])
      class2 = as.vector(datfam2[datfam2$order == min(datfam2$order, na.rm=TRUE), "class1"])
    } 
  }
  
  # # (4) search for species records in EOL
  # eol = tryCatch(as.data.frame(eol_search(taxa)), error = function(e) e)
  # if(class(eol)[1] == 'simpleError'){ return(list(Original=taxname, Submitted=taxa, Taxa_synonyms="EOL lookup failed", Selected_family="EOL lookup failed"))}
  # 
  # # (5) get EOL species records and identify current synonyms
  # if(nrow(eol)>5) eol = eol[ 1:5, ]
  # getEOLrecs = function(x){ as.data.frame(eol_pages(eol$pageid[x])$scinames) }
  # recs = tryCatch(do.call(rbind.data.frame, lapply(1:nrow(eol), getEOLrecs)), error=function(e) e)
  # #recs = tryCatch(eol_pages(eol$pageid[1])$scinames, error=function(e) e)
  # if(class(recs)[1] == 'simpleError'){ return(list(Original=taxname, Submitted=taxa, Taxa_synonyms="EOL processing failed", Selected_family="EOL processing failed"))}
  # synonyms = unique(recs$canonicalform)
  
  # (4) search for species synonyms in ITIS
  syns = tryCatch(suppressMessages(synonyms(taxa, db='itis')), error=function(e) e)
  if(class(syns)[1] == 'simpleError'){ return(data.frame(Original=taxname, Submitted=taxa, Accepted_name="failed", Selected_family=family2, Selected_order=order2, Selected_class=class2, Synonyms="failed"))}
  syns = as.data.frame(syns[[1]])
  
  # get info
  original = taxa
  accepted_name = taxa # save accepted name as original searched name
  if("acc_name" %in% names(syns)){ accepted_name = syns$acc_name } # unless search shows that this is not the accepted name
  # if("message" %in% names(syns) & syns$message == "no syns found"){ synonyms = NA 
  # } else if("syn_name" %in% names(syns)){ synonyms = unique(syns$syn_name) 
  # } else{ synonyms = NA
  # } # add synonyms
  if("syn_name" %in% names(syns)){ synonyms = unique(syns$syn_name) 
  } else{ synonyms = NA }
  
  # combine into list and add synonyms 
  result = data.frame(Original=taxname, 
                      Submitted=taxa,
                      Accepted_name=accepted_name,
                      Selected_family=family2,
                      Selected_order=order2,
                      Selected_class=class2)
  result = do.call("rbind", replicate(length(synonyms), result[1, ], simplify = FALSE))
  result$Synonyms = synonyms
  return(result)
}

# nest function within a tryCatch call in case of any errors
findSyns3 = function(x){
  result = tryCatch(findSyns2(x), error=function(e) NULL)
  return(result)
}

# # read in host data
hostdata = read.csv("./output/crossref_temp/Shaw_HostNames_formatching.csv", stringsAsFactors=F)
hosts_vec = tolower(unique(hostdata$HostSpecies.new))

# ### run lookup: non parallelised
start = Sys.time()
synonyms_lookup = lapply(1:length(hosts_vec), findSyns3)
end = Sys.time()
end-start

# save R object and also save as dataframe
save(synonyms_lookup, file="./output/crossref_temp/Shaw_names_matched_list.R")
syns_df = do.call(rbind.data.frame, synonyms_lookup)
write.csv(syns_df, "./output/crossref_temp/Shaw_HostNames_matched.csv", row.names=F)

# 
# ### run full lookup in parallel
# # detect number of available cores and subtract 2 (use 6)
# cores = detectCores()-2
# 
# # create core cluster and register it with the 'foreach' package
# cl = makeCluster(cores)
# registerDoParallel(cl)
# 
# # run synonym queries
# start = Sys.time()
# foreach(i=1:length(hosts_vec))  %dopar% {
#   library(taxize)
#   require(rgbif)
#   syns = findSyns3(i)
#   if(is.null(syns)){ next 
#     } else{ write.csv(syns, paste("./output/host_parasite_processed/synonyms/Syns", i, as.character(syns$Original[1]), ".csv", sep="_"), row.names=F) }
# }
# end = Sys.time() 
# 
# # end cluster
# stopCluster(cl)
# 


# ### ----- load all synonyms and match to host dataset ----------
# 
# # data
# files = list.files("./output/host_parasite_processed/synonyms/", pattern=".csv", full.names = T)
# all_syns = do.call(rbind.data.frame, lapply(files, read.csv))
# write.csv(all_syns, "./output/host_parasite_processed/hostdata_allsynonyms_itis_20180917.csv", row.names=F)
# 
# # summarise by original
# all_syns1 = ddply(all_syns, .(Original), plyr::summarise,
#                  AcceptedBinomial=head(Accepted_name, 1),
#                  Synonyms=paste(unique(Synonyms), collapse=", "))
# all_syns1$AcceptedBinomial = tolower(all_syns1$AcceptedBinomial)
# all_syns1$Synonyms = tolower(all_syns1$Synonyms)
# all_syns1$Synonyms[ all_syns1$Synonyms == "na" ] = ""
# 
# # create an 'all synonyms' column
# allSyns = function(x){
#   r = all_syns1[x, ]
#   orig = as.character(r$Original)
#   syns = c(r$AcceptedBinomial, unlist(strsplit(r$Synonyms, ", ")))
#   syns = syns[ syns != orig ]
#   return(paste(syns, collapse=", "))
# }
# all_syns1$LookupSynonyms = unlist(lapply(1:nrow(all_syns1), allSyns))
# 
# 
# 
# # ------- add synonyms info to host dataset ----------
# 
# # add synonyms info
# names(hostdata)[5] = "Synonyms_old"
# hostdata$AcceptedBinomial = all_syns1$AcceptedBinomial[ match(hostdata$Host, all_syns1$Original) ]
# hostdata$AcceptedSynonyms = all_syns1$Synonyms[ match(hostdata$Host, all_syns1$Original) ]
# hostdata$LookupSynonyms = all_syns1$LookupSynonyms[ match(hostdata$Host, all_syns1$Original) ]
# 
# # sort out binomial for non-matches
# hostdata$AcceptedBinomial[ is.na(hostdata$AcceptedBinomial) ] = hostdata$Host[ is.na(hostdata$AcceptedBinomial) ]
# 
# # save
# write.csv(hostdata, "./output/host_parasite_processed/hostdata_alldatasets_harmonised_20180917.csv", row.names=F)
# 

