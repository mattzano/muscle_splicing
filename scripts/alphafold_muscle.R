##alphafold
library(r3dmol)
#install.packages("bio3d")
library(bio3d)
library(tidyverse)



RColorBrewer::display.brewer.all(colorblindFriendly = T)
RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")

r3dmol() %>%
  m_add_model(data = here::here("data/tcr_modelling/final_models_hdgfl2_correct/ranked_0.pdb"), format = "pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_surface()) %>% 
  #m_add_style(style = m_style_cartoon(color = "white"),
  #            m_sel(ss = c("h","s"))) %>% 
  m_add_style(style = m_style_cartoon(color = "#ABD9E9", thickness = 0.4),
              m_sel(chain = "C")) %>% 
  m_add_style(style = m_style_cartoon(color = "#313695", thickness = 0.01),
              m_sel(chain = c("A"))) %>% 
  m_add_style(style = m_style_cartoon(color = "#D73027", thickness = 0.01),
              m_sel(chain = c("E"))) %>%
  m_add_style(style = m_style_cartoon(color = "#FEE090", thickness = 0.01),
              m_sel(chain = c("D")))
####compare this to random generated TCRs - use validated from shawn


r3dmol() %>%
  m_add_model(data = here::here("data/tcr_modelling/final_models_hdgfl2_cmvpept/ranked_0.pdb"), format = "pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_surface()) %>% 
  #m_add_style(style = m_style_cartoon(color = "white"),
  #            m_sel(ss = c("h","s"))) %>% 
  m_add_style(style = m_style_cartoon(color = "#ABD9E9", thickness = 0.4),
              m_sel(chain = "C")) %>% 
  m_add_style(style = m_style_cartoon(color = "#313695", thickness = 0.01),
              m_sel(chain = c("A"))) %>% 
  m_add_style(style = m_style_cartoon(color = "#D73027", thickness = 0.01),
              m_sel(chain = c("E"))) %>%
  m_add_style(style = m_style_cartoon(color = "#FEE090", thickness = 0.01),
              m_sel(chain = c("D")))
