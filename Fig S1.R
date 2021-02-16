install.packages("remotes")
remotes::install_github("madsalbertsen/ampvis2@*release")
library(ampvis2)

myotutable <- read.csv("Fig S1 arc_otu_info.csv", check.names = FALSE)
mymetadata <- read.csv("Fig S1 arc_meta_info.csv", check.names = FALSE)
arc <- amp_load(otutable = myotutable,
              metadata = mymetadata)
(pa_pcoa=amp_ordinate(arc, 
             type = "pcoa",
             distmeasure = "bray",
             species_label_size = 11,
             sample_color_by = "MPP",
             sample_colorframe = TRUE,
             sample_colorframe_label = "NAME") + theme(legend.position = "NAME")
)
# Venn diagram grouped by WWTP
(pa_venn=amp_venn(arc,
         cut_a=0.001,
         cut_f=80,
         group_by = "MPP"))

###Bacteria
myotutable <- read.csv("Fig S1 bac_otu_info.csv", check.names = FALSE)
mymetadata <- read.csv("Fig S1 bac_meta_info.csv", check.names = FALSE)
bac <- amp_load(otutable = myotutable,
              metadata = mymetadata)
(pb_pcoa=amp_ordinate(bac, 
             type = "pcoa",
             distmeasure = "bray",
             species_label_size = 11,
             sample_color_by = "MPP",
             sample_colorframe = TRUE,
             sample_colorframe_label = "NAME") + theme(legend.position = "NAME")
)
# Venn diagram grouped by WWTP
(pb_venn=amp_venn(bac,
         cut_a=0.001,
         cut_f=80,
         group_by = "MPP")
)

library(cowplot)
plot_grid(pa_pcoa, pb_pcoa,pa_venn,pb_venn, labels = c("A", "B","C","D"))
ggsave(file="Fig S1.pdf",dpi=300,width = 7.5, units="in",limitsize = FALSE)

