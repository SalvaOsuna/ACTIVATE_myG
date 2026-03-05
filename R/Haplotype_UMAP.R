#Load relevant libraries
 #install.packages(c("crosshap", "umap", "gganimate"))
library(crosshap)
library(umap)
library(gganimate)
library(ggplot2)
library(dplyr)
library(tibble)

#Read in example data (LD Matrix, VCF, Haplotype object and chosen epsilon)
LD <- crosshap::LD
vcf <- crosshap::vcf
HapObject <- crosshap::HapObject
epsilon <- 0.6

#Capturing Marker Groups with UMAP
#Run UMAP on LD matrix to get x & y coordinates for plotting
umap_in <- umap::umap(LD, min_dist = 2, spread = 2.5, n_neighbors = 30)

#Add UMAP coordinates to SNPs with Marker Group information
MGumap_data <- umap_in$layout %>% tibble::as_tibble() %>%
  cbind(rownames(umap_in$data)) %>% tibble::as_tibble() %>%
  dplyr::rename(UMAP1 = V1, UMAP2 = V2, ID = 'rownames(umap_in$data)') %>%
  #HapObject$Haplotypes_MGmin30_E0.6$Varfile contains variant information at epsilon = 0.6
  dplyr::left_join(HapObject$Haplotypes_MGmin30_E0.6$Varfile %>% 
                     dplyr::select(ID, MGs), by = "ID") %>%
  dplyr::mutate(MGs = gsub('0', NA, MGs))

#Plot SNPs by UMAP coordinates, coloured by Marker Group 
MGumap <- ggplot2::ggplot(MGumap_data, ggplot2::aes(x = UMAP1, y = UMAP2)) +
  ggplot2::geom_point(alpha = 0.4, ggplot2::aes(colour = MGs), size = 1, na.rm = F)+
  ggplot2::theme_minimal() +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 5, alpha = 0.7), title = "Marker group"))

MGumap
