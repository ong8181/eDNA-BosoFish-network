####
#### Boso Peninsula project
#### Load fish images
#### 2022.11.11 Ushio (R4.2.1)
####


# ------------------------------------- #
# Load libraries
# ------------------------------------- #
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.11
library(ggimage); packageVersion("ggimage") # 0.3.1, 2022.11.11
library(magick); packageVersion("magick") # 2.7.3, 2021.12.8
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1


# ------------------------------------- #
# Load fish images
# ------------------------------------- #
return_gg <- function (img) {
  return(ggdraw() + draw_image(img))
}

## 1-10
Acanthopagrus_schelegeli <- "0_Illustrations/fish_img/Acanthopagrus_schelegeli.png" %>% image_read() %>% return_gg()
Atherion_elymus <- "0_Illustrations/fish_img/Atherion_elymus.png" %>% image_read() %>% return_gg()
Bathygobius_fuscus <- "0_Illustrations/fish_img/Bathygobius_fuscus.png" %>% image_read() %>% return_gg()
Canthigaster_rivulata <- "0_Illustrations/fish_img/Canthigaster_rivulata.png" %>% image_read() %>% return_gg()
Chaenogobius_annularis <- "0_Illustrations/fish_img/Chaenogobius_annularis.png" %>% image_read() %>% return_gg()
Chaenogobius_gulosus <- "0_Illustrations/fish_img/Chaenogobius_gulosus.png" %>% image_read() %>% return_gg()
Cheilodactylus_zonatus <- "0_Illustrations/fish_img/Cheilodactylus_zonatus.png" %>% image_read() %>% return_gg()
Dictyosoma_rubrimaculatum <- "0_Illustrations/fish_img/Dictyosoma_rubrimaculatum.png" %>% image_read() %>% return_gg()
Ditrema_temminckii_temminckii <- "0_Illustrations/fish_img/Ditrema_temminckii_temminckii.png" %>% image_read() %>% return_gg()
## 11-20
Enchelycore_pardalis <- "0_Illustrations/fish_img/Enchelycore_pardalis.png" %>% image_read() %>% return_gg()
Engraulis_japonicus <- "0_Illustrations/fish_img/Engraulis_japonicus.png" %>% image_read() %>% return_gg()
Enneapterygius_etheostomus <- "0_Illustrations/fish_img/Enneapterygius_etheostomus.png" %>% image_read() %>% return_gg()
Entomacrodus_stellifer_stellifer <- "0_Illustrations/fish_img/Entomacrodus_stellifer_stellifer.png" %>% image_read() %>% return_gg()
Etrumeus_teres <- "0_Illustrations/fish_img/Etrumeus_teres.png" %>% image_read() %>% return_gg()
Eviota_abax <- "0_Illustrations/fish_img/Eviota_abax.png" %>% image_read() %>% return_gg()
Girella_lenonina <- "0_Illustrations/fish_img/Girella_lenonina.png" %>% image_read() %>% return_gg()
Girella_punctata <- "0_Illustrations/fish_img/Girella_punctata.png" %>% image_read() %>% return_gg()
Gymnothorax_kidako <- "0_Illustrations/fish_img/Gymnothorax_kidako.png" %>% image_read() %>% return_gg()
## 21-30
Halichoeres_tenuispinis <- "0_Illustrations/fish_img/Halichoeres_tenuispinis.png" %>% image_read() %>% return_gg()
Hypoatherina_tsurugae <- "0_Illustrations/fish_img/Hypoatherina_tsurugae.png" %>% image_read() %>% return_gg()
Iso_flosmaris <- "0_Illustrations/fish_img/Iso_flosmaris.png" %>% image_read() %>% return_gg()
Istiblennius_enosimae <- "0_Illustrations/fish_img/Istiblennius_enosimae.png" %>% image_read() %>% return_gg()
Kyphosus_bigibbus <- "0_Illustrations/fish_img/Kyphosus_bigibbus.png" %>% image_read() %>% return_gg()
Lateolabrax_japonicus <- "0_Illustrations/fish_img/Lateolabrax_japonicus.png" %>% image_read() %>% return_gg()
Lateolabrax_latus <- "0_Illustrations/fish_img/Lateolabrax_latus.png" %>% image_read() %>% return_gg()
Microcanthus_strigatus <- "0_Illustrations/fish_img/Microcanthus_strigatus.png" %>% image_read() %>% return_gg()
Mugil_cephalus <- "0_Illustrations/fish_img/Mugil_cephalus.png" %>% image_read() %>% return_gg()
Neoclinus_bryope <- "0_Illustrations/fish_img/Neoclinus_bryope.png" %>% image_read() %>% return_gg()
## 31-40
Omobranchus_elegans <- "0_Illustrations/fish_img/Omobranchus_elegans.png" %>% image_read() %>% return_gg()
Oplegnathus_fasciatus <- "0_Illustrations/fish_img/Oplegnathus_fasciatus.png" %>% image_read() %>% return_gg()
Ostracion_immaclatus <- "0_Illustrations/fish_img/Ostracion_immaclatus.png" %>% image_read() %>% return_gg()
Parablennius_yatabei <- "0_Illustrations/fish_img/Parablennius_yatabei.png" %>% image_read() %>% return_gg()
Paralichthys_olivaceus <- "0_Illustrations/fish_img/Paralichthys_olivaceus.png" %>% image_read() %>% return_gg()
Pempheris_schwenkii <- "0_Illustrations/fish_img/Pempheris_schwenkii.png" %>% image_read() %>% return_gg()
Prionurus_scalprum <- "0_Illustrations/fish_img/Prionurus_scalprum.png" %>% image_read() %>% return_gg()
Pseudoblennius_marmoratus <- "0_Illustrations/fish_img/Pseudoblennius_marmoratus.png" %>% image_read() %>% return_gg()
Pseudoblennius_percoides <- "0_Illustrations/fish_img/Pseudoblennius_percoides.png" %>% image_read() %>% return_gg()
Pseudolabrus_eoethinus <- "0_Illustrations/fish_img/Pseudolabrus_eoethinus.png" %>% image_read() %>% return_gg()
Pseudolabrus_sieboldi<- "0_Illustrations/fish_img/Pseudolabrus_sieboldi.png" %>% image_read() %>% return_gg()
## 41-50
Pterogobius_elapoides <- "0_Illustrations/fish_img/Pterogobius_elapoides.png" %>% image_read() %>% return_gg()
Sardinops_melanostictus <- "0_Illustrations/fish_img/Sardinops_melanostictus.png" %>% image_read() %>% return_gg()
Scomber_spp <- "0_Illustrations/fish_img/Scomber_spp.png" %>% image_read() %>% return_gg()
Sebastes_spp <- "0_Illustrations/fish_img/Sebastes_spp.png" %>% image_read() %>% return_gg()
Siganus_fuscescens <- "0_Illustrations/fish_img/Siganus_fuscescens.png" %>% image_read() %>% return_gg()
Spratelloides_gracilis <- "0_Illustrations/fish_img/Spratelloides_gracilis.png" %>% image_read() %>% return_gg()
Stethojulis_interrupta_terina <- "0_Illustrations/fish_img/Stethojulis_interrupta_terina.png" %>% image_read() %>% return_gg()
Takifugu_niphobles <- "0_Illustrations/fish_img/Takifugu_niphobles.png" %>% image_read() %>% return_gg()
Thalassoma_cupido <- "0_Illustrations/fish_img/Thalassoma_cupido.png" %>% image_read() %>% return_gg()
Trachurus_japonicus <- "0_Illustrations/fish_img/Trachurus_japonicus.png" %>% image_read() %>% return_gg()
Zoarchias_neglectus <- "0_Illustrations/fish_img/Zoarchias_neglectus.png" %>% image_read() %>% return_gg()



# ------------------------------------- #
# Wrap fish images
# ------------------------------------- #
ggfish <- list(Acanthopagrus_schelegeli,
               Atherion_elymus,
               Bathygobius_fuscus,
               Canthigaster_rivulata,
               Chaenogobius_annularis,
               Chaenogobius_gulosus,
               Cheilodactylus_zonatus,
               Dictyosoma_rubrimaculatum,
               Ditrema_temminckii_temminckii,
               Enchelycore_pardalis,
               Engraulis_japonicus,
               Enneapterygius_etheostomus,
               Entomacrodus_stellifer_stellifer,
               Etrumeus_teres,
               Eviota_abax,
               Girella_lenonina,
               Girella_punctata,
               Gymnothorax_kidako,
               Halichoeres_tenuispinis,
               Hypoatherina_tsurugae,
               Iso_flosmaris,
               Istiblennius_enosimae,
               Kyphosus_bigibbus,
               Lateolabrax_japonicus,
               Lateolabrax_latus,
               Microcanthus_strigatus,
               Mugil_cephalus,
               Neoclinus_bryope,
               Omobranchus_elegans,
               Oplegnathus_fasciatus,
               Ostracion_immaclatus,
               Parablennius_yatabei,
               Paralichthys_olivaceus,
               Pempheris_schwenkii,
               Prionurus_scalprum,
               Pseudoblennius_marmoratus,
               Pseudoblennius_percoides,
               Pseudolabrus_eoethinus,
               Pseudolabrus_sieboldi,
               Pterogobius_elapoides,
               Sardinops_melanostictus,
               Scomber_spp,
               Sebastes_spp,
               Siganus_fuscescens,
               Spratelloides_gracilis,
               Stethojulis_interrupta_terina,
               Takifugu_niphobles,
               Thalassoma_cupido,
               Trachurus_japonicus,
               Zoarchias_neglectus)

# Add names
names(ggfish) <- c("Acanthopagrus_schelegeli",
                   "Atherion_elymus",
                   "Bathygobius_fuscus",
                   "Canthigaster_rivulata",
                   "Chaenogobius_annularis",
                   "Chaenogobius_gulosus",
                   "Cheilodactylus_zonatus",
                   "Dictyosoma_rubrimaculatum",
                   "Ditrema_temminckii_temminckii",
                   "Enchelycore_pardalis",
                   "Engraulis_japonicus",
                   "Enneapterygius_etheostomus",
                   "Entomacrodus_stellifer_stellifer",
                   "Etrumeus_teres",
                   "Eviota_abax",
                   "Girella_lenonina",
                   "Girella_punctata",
                   "Gymnothorax_kidako",
                   "Halichoeres_tenuispinis",
                   "Hypoatherina_tsurugae",
                   "Iso_flosmaris",
                   "Istiblennius_enosimae",
                   "Kyphosus_bigibbus",
                   "Lateolabrax_japonicus",
                   "Lateolabrax_latus",
                   "Microcanthus_strigatus",
                   "Mugil_cephalus",
                   "Neoclinus_bryope",
                   "Omobranchus_elegans",
                   "Oplegnathus_fasciatus",
                   "Ostracion_immaclatus",
                   "Parablennius_yatabei",
                   "Paralichthys_olivaceus",
                   "Pempheris_schwenkii",
                   "Prionurus_scalprum",
                   "Pseudoblennius_marmoratus",
                   "Pseudoblennius_percoides",
                   "Pseudolabrus_eoethinus",
                   "Pseudolabrus_sieboldi",
                   "Pterogobius_elapoides",
                   "Sardinops_melanostictus",
                   "Scomber_spp",
                   "Sebastes_spp",
                   "Siganus_fuscescens",
                   "Spratelloides_gracilis",
                   "Stethojulis_interrupta_terina",
                   "Takifugu_niphobles",
                   "Thalassoma_cupido",
                   "Trachurus_japonicus",
                   "Zoarchias_neglectus")

# Multi-language version
#ggfish_all <- list(ggfish_en, ggfish_jp)
#names(ggfish_all) <- c("ggfish_en", "ggfish_jp")
#names(ggfish_all$ggfish_en)

