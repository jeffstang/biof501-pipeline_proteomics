library(dplyr)
library(tidyr)
library(tibble)

# Load the RDA downloaded from: https://github.com/jinworks/CellChat/blob/main/data/CellChatDB.mouse.rda
load(file = "CellChatDB.mouse.rda")

# Process the interactions table
cellchat_mouse_LRI <- CellChatDB.mouse$interaction %>% 
  # I selected specifically on the gene symbols
  select(ligand.symbol, receptor.symbol) %>%
  # Each value in ligand.symbol and receptor.symbol had multiple genes separated by a comma
  # I wanted to decouple these values and flatten it all
  # to one ligand and its multiple receptor associated interactions
  mutate(receptor.symbol = strsplit(receptor.symbol, ",\\s*")) %>%
  unnest(receptor.symbol) %>%
  mutate(ligand.symbol = strsplit(ligand.symbol, ",\\s*")) %>%
  unnest(ligand.symbol) %>%
  group_by(ligand.symbol) %>%
  distinct() %>%
  # collapse receptors into a list
  summarise(Receptors = list(receptor.symbol), .groups = "drop") %>%
  deframe()

# Save this transformed list into an RDA object for downstream analysis
save(cellchat_mouse_LRI, file = "cellchatv2_mouseLRI.rda")