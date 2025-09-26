# this script is to extract clades from different infraorders using Wiens squamate tree
# then collect each trait data across these infraorders based on SquamBase dataset. 

# Note some species can have the same synonym due to them being treated as a same species in the past
# so I ignore those species when getting clades. 

#### LOAD LIBRARIES ####
library(diversitree)
library(ape)
library(openxlsx)

#### FILE SETTING ####
fp      = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
in_fp   = paste0(fp,"data/Empirical/squamate_phylogeny/")
out_fp  = paste0(in_fp,"trait_data/")
tree_fp = paste0(in_fp,"clade_data/")

setwd(in_fp)

# load squamate trait data from SquamBase 
squam_trait_data = read.xlsx("trait_data/Supplementary_Table_S1_-_squamBase1.xlsx")
col_names = names(squam_trait_data)
squam_trait_data[, 1] <- gsub(" ", "_", squam_trait_data[, 1]) # add "_" between genus and species name 
# add "_" between genus and species name for synonyms too (same speces can have multipe synonyms)
squam_trait_data[, 6] <- ifelse(
  is.na(squam_trait_data[, 6]),
  NA,
  sapply(strsplit(squam_trait_data[, 6], ","), function(parts) {
    parts <- trimws(parts)           # remove spaces before/after species
    parts <- gsub(" ", "_", parts)   # replace spaces inside species names
    paste(parts, collapse = ", ")    # join back with comma + space
  })
)
#

# subset according to each infraorder
gekkota_trait         = squam_trait_data[squam_trait_data$infraorder=="Gekkota",] #Gekkota
alethinophidia_trait  = squam_trait_data[squam_trait_data$infraorder=="Alethinophidia",] #Alethinophidia
acrodontia_trait      = squam_trait_data[squam_trait_data$infraorder=="Acrodontia",] #Acrodontia
laterata_trait        = squam_trait_data[squam_trait_data$infraorder=="Laterata",] #Laterata
scolecophidia_trait   = squam_trait_data[squam_trait_data$infraorder=="Scolecophidia",] #Scolecophidia
anguimoprha_trait     = squam_trait_data[squam_trait_data$infraorder=="Anguimoprha",] #Anguimoprha
pleurodonta_trait     = squam_trait_data[squam_trait_data$infraorder=="Pleurodonta",] #Pleurodonta
dibamidae_trait       = squam_trait_data[squam_trait_data$infraorder=="Dibamidae",] #Dibamidae
scincoidea_trait      = squam_trait_data[squam_trait_data$infraorder=="Scincoidea",] #Scincoidea
amphisbaenia_trait    = squam_trait_data[squam_trait_data$infraorder=="Amphisbaenia (Laterata)",] #Amphisbaenia

# convert to binary trait: Carnivorous -> no_plant, Omnivorous & Herbivorous -> plant 
gekkota_trait$Diet[gekkota_trait$Diet == "Carnivorous" | gekkota_trait$Diet == "*Carnivorous" ] <- "no_plant"
gekkota_trait$Diet[gekkota_trait$Diet == "Omnivorous" | gekkota_trait$Diet == "omnivorous" | gekkota_trait$Diet == "Herbivorous" ] <- "plant"

alethinophidia_trait$Diet[alethinophidia_trait$Diet == "Carnivorous" | alethinophidia_trait$Diet == "*Carnivorous" ] <- "no_plant"
alethinophidia_trait$Diet[alethinophidia_trait$Diet == "Omnivorous" | alethinophidia_trait$Diet == "omnivorous" | alethinophidia_trait$Diet == "Herbivorous" ] <- "plant"

acrodontia_trait$Diet[acrodontia_trait$Diet == "Carnivorous" | acrodontia_trait$Diet == "*Carnivorous" ] <- "no_plant"
acrodontia_trait$Diet[acrodontia_trait$Diet == "Omnivorous" | acrodontia_trait$Diet == "omnivorous" | acrodontia_trait$Diet == "Herbivorous" ] <- "plant"

laterata_trait$Diet[laterata_trait$Diet == "Carnivorous" | laterata_trait$Diet == "*Carnivorous" ] <- "no_plant"
laterata_trait$Diet[laterata_trait$Diet == "Omnivorous" | laterata_trait$Diet == "omnivorous" | laterata_trait$Diet == "Herbivorous" ] <- "plant"

scolecophidia_trait$Diet[scolecophidia_trait$Diet == "Carnivorous" | scolecophidia_trait$Diet == "*Carnivorous" ] <- "no_plant"
scolecophidia_trait$Diet[scolecophidia_trait$Diet == "Omnivorous" | scolecophidia_trait$Diet == "omnivorous" | scolecophidia_trait$Diet == "Herbivorous" ] <- "plant"

anguimoprha_trait$Diet[anguimoprha_trait$Diet == "Carnivorous" | anguimoprha_trait$Diet == "*Carnivorous" ] <- "no_plant"
anguimoprha_trait$Diet[anguimoprha_trait$Diet == "Omnivorous" | anguimoprha_trait$Diet == "omnivorous" | anguimoprha_trait$Diet == "Herbivorous" ] <- "plant"

pleurodonta_trait$Diet[pleurodonta_trait$Diet == "Carnivorous" | pleurodonta_trait$Diet == "*Carnivorous" ] <- "no_plant"
pleurodonta_trait$Diet[pleurodonta_trait$Diet == "Omnivorous" | pleurodonta_trait$Diet == "omnivorous" | pleurodonta_trait$Diet == "Herbivorous" ] <- "plant"

dibamidae_trait$Diet[dibamidae_trait$Diet == "Carnivorous" | dibamidae_trait$Diet == "*Carnivorous" ] <- "no_plant"
dibamidae_trait$Diet[dibamidae_trait$Diet == "Omnivorous" | dibamidae_trait$Diet == "omnivorous" | dibamidae_trait$Diet == "Herbivorous" ] <- "plant"

scincoidea_trait$Diet[scincoidea_trait$Diet == "Carnivorous" | scincoidea_trait$Diet == "*Carnivorous" ] <- "no_plant"
scincoidea_trait$Diet[scincoidea_trait$Diet == "Omnivorous" | scincoidea_trait$Diet == "omnivorous" | scincoidea_trait$Diet == "Herbivorous" ] <- "plant"

amphisbaenia_trait$Diet[amphisbaenia_trait$Diet == "Carnivorous" | amphisbaenia_trait$Diet == "*Carnivorous" ] <- "no_plant"
amphisbaenia_trait$Diet[amphisbaenia_trait$Diet == "Omnivorous" | amphisbaenia_trait$Diet == "omnivorous" | amphisbaenia_trait$Diet == "Herbivorous" ] <- "plant"

# load squamate tree from Wiens et al
squam_tree = read.tree("squamate.txt")
squam_tree = force.ultrametric(squam_tree) #force to be ultrametric due to measurement error

#### GEKKOTA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
gekkota_diet = gekkota_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
gekkota_diet = gekkota_diet[!is.na(gekkota_diet$Diet),]
# get species names belong to gekkota
gekkota_spec = gekkota_diet$`Species.name.(Binomial)`
# get list of missing gekkota species from the tree 
gekkota_missing       = setdiff(gekkota_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(gekkota_missing, function(spec) { #loop over each missing species from tree
 syn <- gekkota_trait$Synonyms[gekkota_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- gekkota_missing
syn_map <- unlist(syn_map)

# Replace in gekkota_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    gekkota_spec[gekkota_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
gekkota_spec <- gekkota_spec[gekkota_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
gekkota_spec <- unique(gekkota_spec)

# Remove those absent species from character data 
gekkota_diet <- gekkota_diet[gekkota_diet$`Species.name.(Binomial)` %in% gekkota_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to gekkota_diet
if (nrow(gekkota_diet) != length(gekkota_spec)){
  which_syn = setdiff(gekkota_spec,gekkota_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = gekkota_trait$Diet[gekkota_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      gekkota_diet = rbind(gekkota_diet,c(which_syn[i],"Gekkota",diet_added))
    }
  }
}

# get the clade
gekkota_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, gekkota_diet$`Species.name.(Binomial)`))

#### ALETHINOPHIDIA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
alethinophidia_diet = alethinophidia_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
alethinophidia_diet = alethinophidia_diet[!is.na(alethinophidia_diet$Diet),]
# get species names belong to alethinophidia
alethinophidia_spec = alethinophidia_diet$`Species.name.(Binomial)`
# get list of missing alethinophidia species from the tree 
alethinophidia_missing       = setdiff(alethinophidia_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(alethinophidia_missing, function(spec) { #loop over each missing species from tree
  syn <- alethinophidia_trait$Synonyms[alethinophidia_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- alethinophidia_missing
syn_map <- unlist(syn_map)

# Replace in alethinophidia_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    alethinophidia_spec[alethinophidia_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
alethinophidia_spec <- alethinophidia_spec[alethinophidia_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
alethinophidia_spec <- unique(alethinophidia_spec)

# Remove those absent species from character data 
alethinophidia_diet <- alethinophidia_diet[alethinophidia_diet$`Species.name.(Binomial)` %in% alethinophidia_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to alethinophidia_diet
if (nrow(alethinophidia_diet) != length(alethinophidia_spec)){
  which_syn = setdiff(alethinophidia_spec,alethinophidia_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    # Some species have the same synonym from the database due to it's being considered a different species now,
    # So, I'll just ignore those species (from the synonym)
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = alethinophidia_trait$Diet[alethinophidia_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      alethinophidia_diet = rbind(alethinophidia_diet,c(which_syn[i],"Alethinophidia",diet_added))
    }
  }
}

# get the clade
alethinophidia_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, alethinophidia_diet$`Species.name.(Binomial)`))


#### ACRODONTIA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
acrodontia_diet = acrodontia_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
acrodontia_diet = acrodontia_diet[!is.na(acrodontia_diet$Diet),]
# get species names belong to acrodontia
acrodontia_spec = acrodontia_diet$`Species.name.(Binomial)`
# get list of missing acrodontia species from the tree 
acrodontia_missing       = setdiff(acrodontia_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(acrodontia_missing, function(spec) { #loop over each missing species from tree
  syn <- acrodontia_trait$Synonyms[acrodontia_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- acrodontia_missing
syn_map <- unlist(syn_map)

# Replace in acrodontia_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    acrodontia_spec[acrodontia_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
acrodontia_spec <- acrodontia_spec[acrodontia_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
acrodontia_spec <- unique(acrodontia_spec)

# Remove those absent species from character data 
acrodontia_diet <- acrodontia_diet[acrodontia_diet$`Species.name.(Binomial)` %in% acrodontia_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to acrodontia_diet
if (nrow(acrodontia_diet) != length(acrodontia_spec)){
  which_syn = setdiff(acrodontia_spec,acrodontia_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = acrodontia_trait$Diet[acrodontia_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      acrodontia_diet = rbind(acrodontia_diet,c(which_syn[i],"Acrodontia",diet_added))
    }
  }
}

# get the clade
acrodontia_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, acrodontia_diet$`Species.name.(Binomial)`))

#### LATERATA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
laterata_diet = laterata_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
laterata_diet = laterata_diet[!is.na(laterata_diet$Diet),]
# get species names belong to laterata
laterata_spec = laterata_diet$`Species.name.(Binomial)`
# get list of missing laterata species from the tree 
laterata_missing       = setdiff(laterata_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(laterata_missing, function(spec) { #loop over each missing species from tree
  syn <- laterata_trait$Synonyms[laterata_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- laterata_missing
syn_map <- unlist(syn_map)

# Replace in laterata_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    laterata_spec[laterata_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
laterata_spec <- laterata_spec[laterata_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
laterata_spec <- unique(laterata_spec)

# Remove those absent species from character data 
laterata_diet <- laterata_diet[laterata_diet$`Species.name.(Binomial)` %in% laterata_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to laterata_diet
if (nrow(laterata_diet) != length(laterata_spec)){
  which_syn = setdiff(laterata_spec,laterata_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = laterata_trait$Diet[laterata_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      laterata_diet = rbind(laterata_diet,c(which_syn[i],"Laterata",diet_added))
    }
  }
}

# get the clade
laterata_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, laterata_diet$`Species.name.(Binomial)`))

#### SCOLECOPHIDIA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
scolecophidia_diet = scolecophidia_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
scolecophidia_diet = scolecophidia_diet[!is.na(scolecophidia_diet$Diet),]
# get species names belong to scolecophidia
scolecophidia_spec = scolecophidia_diet$`Species.name.(Binomial)`
# get list of missing scolecophidia species from the tree 
scolecophidia_missing       = setdiff(scolecophidia_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(scolecophidia_missing, function(spec) { #loop over each missing species from tree
  syn <- scolecophidia_trait$Synonyms[scolecophidia_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- scolecophidia_missing
syn_map <- unlist(syn_map)

# Replace in scolecophidia_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    scolecophidia_spec[scolecophidia_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
scolecophidia_spec <- scolecophidia_spec[scolecophidia_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
scolecophidia_spec <- unique(scolecophidia_spec)

# Remove those absent species from character data 
scolecophidia_diet <- scolecophidia_diet[scolecophidia_diet$`Species.name.(Binomial)` %in% scolecophidia_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to scolecophidia_diet
if (nrow(scolecophidia_diet) != length(scolecophidia_spec)){
  which_syn = setdiff(scolecophidia_spec,scolecophidia_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = scolecophidia_trait$Diet[scolecophidia_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      scolecophidia_diet = rbind(scolecophidia_diet,c(which_syn[i],"Scolecophidia",diet_added))
    }
  }
}

# get the clade
scolecophidia_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, scolecophidia_diet$`Species.name.(Binomial)`))

#### ANGUIMORPHA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
anguimoprha_diet = anguimoprha_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
anguimoprha_diet = anguimoprha_diet[!is.na(anguimoprha_diet$Diet),]
# get species names belong to anguimoprha
anguimoprha_spec = anguimoprha_diet$`Species.name.(Binomial)`
# get list of missing anguimoprha species from the tree 
anguimoprha_missing       = setdiff(anguimoprha_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(anguimoprha_missing, function(spec) { #loop over each missing species from tree
  syn <- anguimoprha_trait$Synonyms[anguimoprha_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- anguimoprha_missing
syn_map <- unlist(syn_map)

# Replace in anguimoprha_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    anguimoprha_spec[anguimoprha_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
anguimoprha_spec <- anguimoprha_spec[anguimoprha_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
anguimoprha_spec <- unique(anguimoprha_spec)

# Remove those absent species from character data 
anguimoprha_diet <- anguimoprha_diet[anguimoprha_diet$`Species.name.(Binomial)` %in% anguimoprha_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to anguimoprha_diet
if (nrow(anguimoprha_diet) != length(anguimoprha_spec)){
  which_syn = setdiff(anguimoprha_spec,anguimoprha_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = anguimoprha_trait$Diet[anguimoprha_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      anguimoprha_diet = rbind(anguimoprha_diet,c(which_syn[i],"Anguimoprha",diet_added))
    }
  }
}

# get the clade
anguimoprha_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, anguimoprha_diet$`Species.name.(Binomial)`))

#### PLEURODONTA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
pleurodonta_diet = pleurodonta_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
pleurodonta_diet = pleurodonta_diet[!is.na(pleurodonta_diet$Diet),]
# get species names belong to pleurodonta
pleurodonta_spec = pleurodonta_diet$`Species.name.(Binomial)`
# get list of missing pleurodonta species from the tree 
pleurodonta_missing       = setdiff(pleurodonta_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(pleurodonta_missing, function(spec) { #loop over each missing species from tree
  syn <- pleurodonta_trait$Synonyms[pleurodonta_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- pleurodonta_missing
syn_map <- unlist(syn_map)

# Replace in pleurodonta_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    pleurodonta_spec[pleurodonta_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
pleurodonta_spec <- pleurodonta_spec[pleurodonta_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
pleurodonta_spec <- unique(pleurodonta_spec)

# Remove those absent species from character data 
pleurodonta_diet <- pleurodonta_diet[pleurodonta_diet$`Species.name.(Binomial)` %in% pleurodonta_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to pleurodonta_diet
if (nrow(pleurodonta_diet) != length(pleurodonta_spec)){
  which_syn = setdiff(pleurodonta_spec,pleurodonta_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = pleurodonta_trait$Diet[pleurodonta_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      pleurodonta_diet = rbind(pleurodonta_diet,c(which_syn[i],"Pleurodonta",diet_added))
    }
  }
}

# get the clade
pleurodonta_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, pleurodonta_diet$`Species.name.(Binomial)`))

#### DIBAMIDAE ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
dibamidae_diet = dibamidae_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
dibamidae_diet = dibamidae_diet[!is.na(dibamidae_diet$Diet),]
# get species names belong to dibamidae
dibamidae_spec = dibamidae_diet$`Species.name.(Binomial)`
# get list of missing dibamidae species from the tree 
dibamidae_missing       = setdiff(dibamidae_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(dibamidae_missing, function(spec) { #loop over each missing species from tree
  syn <- dibamidae_trait$Synonyms[dibamidae_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- dibamidae_missing
syn_map <- unlist(syn_map)

# Replace in dibamidae_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    dibamidae_spec[dibamidae_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
dibamidae_spec <- dibamidae_spec[dibamidae_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
dibamidae_spec <- unique(dibamidae_spec)

# Remove those absent species from character data 
dibamidae_diet <- dibamidae_diet[dibamidae_diet$`Species.name.(Binomial)` %in% dibamidae_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to dibamidae_diet
if (nrow(dibamidae_diet) != length(dibamidae_spec)){
  which_syn = setdiff(dibamidae_spec,dibamidae_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = dibamidae_trait$Diet[dibamidae_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      dibamidae_diet = rbind(dibamidae_diet,c(which_syn[i],"Dibamidae",diet_added))
    }
  }
}

# get the clade
dibamidae_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, dibamidae_diet$`Species.name.(Binomial)`))

#### SCINCOIDEA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
scincoidea_diet = scincoidea_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
scincoidea_diet = scincoidea_diet[!is.na(scincoidea_diet$Diet),]
# get species names belong to scincoidea
scincoidea_spec = scincoidea_diet$`Species.name.(Binomial)`
# get list of missing scincoidea species from the tree 
scincoidea_missing       = setdiff(scincoidea_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(scincoidea_missing, function(spec) { #loop over each missing species from tree
  syn <- scincoidea_trait$Synonyms[scincoidea_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- scincoidea_missing
syn_map <- unlist(syn_map)

# Replace in scincoidea_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    scincoidea_spec[scincoidea_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
scincoidea_spec <- scincoidea_spec[scincoidea_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
scincoidea_spec <- unique(scincoidea_spec)

# Remove those absent species from character data 
scincoidea_diet <- scincoidea_diet[scincoidea_diet$`Species.name.(Binomial)` %in% scincoidea_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to scincoidea_diet
if (nrow(scincoidea_diet) != length(scincoidea_spec)){
  which_syn = setdiff(scincoidea_spec,scincoidea_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = scincoidea_trait$Diet[scincoidea_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      scincoidea_diet = rbind(scincoidea_diet,c(which_syn[i],"Scincoidea",diet_added))
    }
  }
}

# get the clade
scincoidea_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, scincoidea_diet$`Species.name.(Binomial)`))

#### AMPHISBAENIA ####
###### TRAIT: DIET ######
# get trait data based on "Diet"
amphisbaenia_diet = amphisbaenia_trait[c("Species.name.(Binomial)","infraorder","Diet")]
# remove all NA in Diet
amphisbaenia_diet = amphisbaenia_diet[!is.na(amphisbaenia_diet$Diet),]
# get species names belong to amphisbaenia
amphisbaenia_spec = amphisbaenia_diet$`Species.name.(Binomial)`
# get list of missing amphisbaenia species from the tree 
amphisbaenia_missing       = setdiff(amphisbaenia_spec,squam_tree$tip.label)

# Replace missing species with synonyms and check if they exist in the tree
syn_map <- lapply(amphisbaenia_missing, function(spec) { #loop over each missing species from tree
  syn <- amphisbaenia_trait$Synonyms[amphisbaenia_trait$`Species.name.(Binomial)` == spec]
  
  if(length(syn) > 0 && !is.na(syn) && syn != "") {
    syn_list <- unlist(strsplit(syn, ",\\s*"))
    found_on_tree <- syn_list[syn_list %in% squam_tree$tip.label] #check if the synonyms exist on the tree
    if(length(found_on_tree) > 0) return(found_on_tree[1]) #get the first synonym if multiple syn exist
  }
  return(NA)  # assign NA if no synonym found (also not on the tree)
})
names(syn_map) <- amphisbaenia_missing
syn_map <- unlist(syn_map)

# Replace in amphisbaenia_spec to contain tips on the tree (based on synonym and no synonym)
for(spec in names(syn_map)) {
  if(!is.na(syn_map[spec])) {
    amphisbaenia_spec[amphisbaenia_spec == spec] <- syn_map[spec]
  }
}

# Remove any species not in the tree (those with NA i.e. no synonym)
amphisbaenia_spec <- amphisbaenia_spec[amphisbaenia_spec %in% squam_tree$tip.label]

# remove same synonym use for multiple species 
amphisbaenia_spec <- unique(amphisbaenia_spec)

# Remove those absent species from character data 
amphisbaenia_diet <- amphisbaenia_diet[amphisbaenia_diet$`Species.name.(Binomial)` %in% amphisbaenia_spec, ]
# if there is a mismatch because it is using the synonym, add that synonym to amphisbaenia_diet
if (nrow(amphisbaenia_diet) != length(amphisbaenia_spec)){
  which_syn = setdiff(amphisbaenia_spec,amphisbaenia_diet$`Species.name.(Binomial)`)
  for (i in 1:length(which_syn)){
    if (length(names(which(syn_map==which_syn[i])))==1){
      diet_added = amphisbaenia_trait$Diet[amphisbaenia_trait$`Species.name.(Binomial)`==names(which(syn_map==which_syn[i]))]
      amphisbaenia_diet = rbind(amphisbaenia_diet,c(which_syn[i],"Amphisbaenia",diet_added))
    }
  }
}

# get the clade
amphisbaenia_tree <- drop.tip(squam_tree, setdiff(squam_tree$tip.label, amphisbaenia_diet$`Species.name.(Binomial)`))

#SAVE THE TREE AND TRAIT DATA (DIET)
#
write.table(gekkota_diet,paste0(out_fp,"gekkota_diet.csv"),sep = ";",row.names = F)
write.table(alethinophidia_diet,paste0(out_fp,"alethinophidia_diet.csv"),sep = ";",row.names = F)
write.table(acrodontia_diet,paste0(out_fp,"acrodontia_diet.csv"),sep = ";",row.names = F)
write.table(laterata_diet,paste0(out_fp,"laterata_diet.csv"),sep = ";",row.names = F)
write.table(scolecophidia_diet,paste0(out_fp,"scolecophidia_diet.csv"),sep = ";",row.names = F)
write.table(anguimoprha_diet,paste0(out_fp,"anguimoprha_diet.csv"),sep = ";",row.names = F)
write.table(pleurodonta_diet,paste0(out_fp,"pleurodonta_diet.csv"),sep = ";",row.names = F)
write.table(dibamidae_diet,paste0(out_fp,"dibamidae_diet.csv"),sep = ";",row.names = F)
write.table(scincoidea_diet,paste0(out_fp,"scincoidea_diet.csv"),sep = ";",row.names = F)
write.table(amphisbaenia_diet,paste0(out_fp,"amphisbaenia_diet.csv"),sep = ";",row.names = F)
#
write.tree(gekkota_tree,paste0(tree_fp,"gekkota.tre"))
write.tree(alethinophidia_tree,paste0(tree_fp,"alethinophidia.tre"))
write.tree(acrodontia_tree,paste0(tree_fp,"acrodontia.tre"))
write.tree(laterata_tree,paste0(tree_fp,"laterata.tre"))
write.tree(scolecophidia_tree,paste0(tree_fp,"scolecophidia.tre"))
write.tree(anguimoprha_tree,paste0(tree_fp,"anguimoprha.tre"))
write.tree(pleurodonta_tree,paste0(tree_fp,"pleurodonta.tre"))
write.tree(dibamidae_tree,paste0(tree_fp,"dibamidae.tre"))
write.tree(scincoidea_tree,paste0(tree_fp,"scincoidea.tre"))
write.tree(amphisbaenia_tree,paste0(tree_fp,"amphisbaenia.tre"))
