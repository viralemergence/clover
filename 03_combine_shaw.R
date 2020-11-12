



# ================== Initial looking at combining Shaw with harmonised data ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("dplyr", "magrittr")



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

# subset to just pathogen names and check whether name already in associations
shawp = shaw %>%
  dplyr::select(Pathogen_Original) %>%
  distinct() %>%
  dplyr::mutate(InAssoc = ifelse(Pathogen_Original %in% assoc$Pathogen_Harmonised, TRUE, FALSE)) %>%
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



  