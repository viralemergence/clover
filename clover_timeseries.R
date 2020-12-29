

# ===================== view temporal trends in detection =======================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202011_clover/clover/")
pacman::p_load("dplyr", "magrittr", "ggplot2")

# data
clover = read.csv("./output/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE)



# plot 1: overall temporal trends in associations
cl1 = clover %>%
  dplyr::group_by(Year, Database) %>%
  dplyr::summarise(NReports = length(Year))
p1 = ggplot(cl1) + 
  geom_bar(aes(Year, NReports, fill=Database), stat="identity") + 
  theme_minimal() + 
  scale_fill_viridis_d(begin=0.05, end=0.95) + 
  xlab("Year") + 
  ylab("Number of associations") + 
  ggtitle("Total host-virus associations by source database") +
  scale_x_continuous(limits = c(1920, 2020), breaks=seq(1920, 2020, by=20)) + 
  theme(plot.title =  element_text(size=13.5, hjust=0.5, vjust=-4),
        legend.position=c(0.08, 0.8), legend.title=element_blank(),
        axis.text = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12))

# plot 2: temporal trends in novel associations
cl2 = clover %>%
  dplyr::mutate(HostVirus = paste(Host, Virus, sep="_")) %>%
  dplyr::group_by(Database, HostVirus) %>%
  dplyr::summarise(Year = min(Year, na.rm=TRUE)) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(NReports = n_distinct(HostVirus))
cl2$Year[ cl2$Year == Inf ] = NA
p2 = ggplot(cl2) + 
  geom_bar(aes(Year, NReports), stat="identity", fill="grey50", col=NA) + 
  theme_minimal() + 
  xlab("Year") + 
  ylab("Number of associations") + 
  ggtitle("Novel host-virus associations (previously unreported)") +
  scale_x_continuous(limits = c(1920, 2020), breaks=seq(1920, 2020, by=20)) + 
  theme(plot.title =  element_text(size=13.5, hjust=0.5, vjust=-4),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12))

# combine
p_comb = gridExtra::grid.arrange(p1, p2, ncol=1)
ggsave(p_comb, file="./output/figures/MS_TemporalTrends.png", device="png", units="in", scale=0.95, width=9, height=6, dpi=300)


# 
# # plot 2: temporal trends in novel associations by database
# cl2 = clover %>%
#   dplyr::mutate(HostVirus = paste(Host, Virus, sep="_")) %>%
#   dplyr::group_by(Database, HostVirus) %>%
#   dplyr::summarise(Year = min(Year, na.rm=TRUE)) %>%
#   dplyr::group_by(Database, Year) %>%
#   dplyr::summarise(NReports = n_distinct(HostVirus))
# cl2$Year[ cl2$Year == Inf ] = NA
# ggplot(cl2) + 
#   geom_bar(aes(Year, NReports, fill=Database), stat="identity") + 
#   theme_minimal() + 
#   scale_fill_viridis_d(begin=0.0, end=0.95)
