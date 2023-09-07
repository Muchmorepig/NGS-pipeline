# ARF2_Clu_select_GO

# prepare -----------------------------------------------------------------

# load package
library(dplyr)
library(ggplot2)
library(BuenColors)
library(patchwork)

list_file <- list.files("rawdata_processed/ARF2_Clu_GO",full.names = T)
source("~/scripts/normal_script/ggplot2_plot/GO_plot/GOplot_barplot.R")

list_plot <- list()
for (Cl in list_file){
  Cl_name <- basename(Cl) %>% gsub("_GO\\.txt","",.)
  
  Cl_info <- read.table(Cl,sep = "\t",header = T)
  
  p <- GOplot_barplot(data = Cl_info,title_name = Cl_name)
  
  list_plot[[Cl_name]] <- p
}


wrap_plots(list_plot,nrow = 1) -> p_list


ggsave(p_list,filename = "plot/WT0h_ARF2_Clu_GO.pdf",width = 40)
