



# ================== Initial looking at combining Shaw with harmonised data ==================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/")
pacman::p_load("dplyr", "magrittr")



# ================= read datasets ================

# read datasets
shaw = read.csv("./data/Shaw_Database_main.csv", stringsAsFactors = FALSE)
assoc = read.csv("./output/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)

# create parasite list, remove 'virus' from end of virus description
par_temp = strsplit(tolower(shaw$Species), " ")
par_temp1 = lapply(par_temp, function(x){ if(x[length(x)] == "virus"){ return(x[1:length(x)-1]) } else(return(x)) })
par_temp1 = unlist(lapply(par_temp1, paste, collapse=" "))
shaw$Pathogen_Original = par_temp1

# all pathogen names - keep mammals/viruses only
shawp = shaw %>%
  dplyr::filter(HostGroup %in% c("Human", "Ungulates", "Carnivora", "Cetacea", "Rodentia", "Mammalia", "Primates", "Chiroptera")) %>%
  dplyr::select(Pathogen_Original, Type) %>%
  distinct() %>%
  dplyr::mutate(Pathogen_Original = tolower(Pathogen_Original)) %>%
  dplyr::filter(Type == "Virus")
shawp$id_shaw = 1:nrow(shawp)
assocp = assoc %>%
  dplyr::select(Pathogen_Harmonised) %>%
  distinct()
assocp$id_assoc = 1:nrow(assocp)

# quick look
# shaw also in assoc?
# ~500 viruses to harmonise/check
shawp$Pathogen_Original[ shawp$Pathogen_Original %in% assocp$Pathogen_Harmonised ]
shawp$Pathogen_Original[ !shawp$Pathogen_Original %in% assocp$Pathogen_Harmonised ]

  