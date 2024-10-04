#!/usr/bin/env Rscript

# Author: Alison Menezes
# Description: Generate phylogenetic trees from .nwk files, root them, modify tip labels, and save figures in SVG format.

# Load necessary libraries
library(ggtree)
library(ggplot2)
library(TDbook)
library(tidyverse)
library(RColorBrewer)
library(geiger)
library(phytools)
library(ape)

# Set working directory
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

### FIGURES ###

### SMF ###
# #load phylogeny
# smf_tree <- read.tree("SMF.nwk")
# smf_tree <- midpoint.root(smf_tree)
# #rename tip labels
# new_labels <- str_replace(smf_tree$tip.label, "_", " ")
# new_tiplabels <- str_replace(new_labels, "!", " | ")
# smf_tree$tip.label <- new_tiplabels
# #plot tree
# ggtree(smf_tree, size=.5)+
#   geom_tiplab(size=6)+ xlim(NA, 1.5) + ylim(0, 92)
# #save figure
# ##ggsave("SMF.svg", width = 18, height = 18, dpi = 1200)

### SPCH
#load phylogeny
spch_tree <- read.tree("SPCH.nwk")
spch_tree <- midpoint.root(spch_tree)
spch_tree <- root(spch_tree, outgroup = c("Arabidopsis_thaliana!SPCH", "Eucalyptus_grandis!SPCH"), resolve.root = TRUE)
#rename tip labels
new_labels <- str_replace(spch_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
spch_tree$tip.label <- new_tiplabels
#plot tree
ggtree(spch_tree, size=.5)+
  geom_tiplab(size=5)+ xlim(NA, 1.5) + ylim(0, 45)+
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 80))
#save figure
#ggsave("SPCH.svg", width = 18, height = 18, dpi = 1200)

### MUTE
#load phylogeny
mute_tree <- read.tree("old/MUTE.nwk")
mute_tree <- midpoint.root(mute_tree)
#rename tip labels
new_labels <- str_replace(mute_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
mute_tree$tip.label <- new_tiplabels
#plot tree
ggtree(mute_tree, size=.5)+
  geom_tiplab(size=6)+ xlim(NA, 1.5) + ylim(0, 23)
#save figure
#ggsave("MUTE.svg", width = 18, height = 18, dpi = 1200)

### FAMA
#load phylogeny
fama_tree <- read.tree("FAMA.nwk")
fama_tree <- midpoint.root(fama_tree)
#rename tip labels
new_labels <- str_replace(fama_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
fama_tree$tip.label <- new_tiplabels
#plot tree
ggtree(fama_tree, size=.5)+
  geom_tiplab(size=10)+ xlim(NA, .7) + ylim(0, 38)
#save figure
#ggsave("FAMA.svg", width = 20, height = 20, dpi = 1200)

### SCRM ###
#load phylogeny
scrm_tree <- read.tree("SCRM.nwk")
scrm_tree <- midpoint.root(scrm_tree)
#rename tip labels
new_labels <- str_replace(scrm_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
scrm_tree$tip.label <- new_tiplabels
#plot tree
ggtree(scrm_tree, size=.5)+
  geom_tiplab(size=7)+ xlim(NA, 1.5) + ylim(0, 65)
#save figure
#ggsave("SCRM.svg", width = 20, height = 20, dpi = 1200)

### EPFL9 ###
#load phylogeny
epfl9_tree <- read.tree("EPFL9.nwk")
epfl9_tree <- midpoint.root(epfl9_tree)
#rename tip labels
new_labels <- str_replace(epfl9_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
epfl9_tree$tip.label <- new_tiplabels
#plot tree
ggtree(epfl9_tree, size= 1)+
  geom_tiplab(size=10)+ xlim(NA, 4) #+ ylim(0, 75) 
#save figure
ggsave("EPFL9.svg", width = 18, height = 18, dpi = 1200)

### SUPPLEMENTARY FIGURES ###

### ABI ###
#load phylogeny
abi_tree <- read.tree("ABI1.nwk")
abi_tree <- midpoint.root(abi_tree)
#rename tip labels
new_labels <- str_replace(abi_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
abi_tree$tip.label <- new_tiplabels
#plot tree
ggtree(abi_tree, size=.5)+
  geom_tiplab(size=8)+ xlim(NA, 1.5) + ylim(0, 55)
#save figure
#ggsave("ABI.svg", width = 18, height = 18, dpi = 1200)

### NAP1 ###
#load phylogeny
nap_tree <- read.tree("NAP1.nwk")
nap_tree <- midpoint.root(nap_tree)
#rename tip labels
new_labels <- str_replace(nap_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
nap_tree$tip.label <- new_tiplabels
#plot tree
ggtree(nap_tree, size=.5)+
  geom_tiplab(size=8)+ xlim(NA, .6) + ylim(0, 38)
#save figure
#ggsave("NAP1.svg", width = 18, height = 18, dpi = 1200)

### PAN1 ###
#load phylogeny
pan1_tree <- read.tree("PAN1.nwk")
pan1_tree <- midpoint.root(pan1_tree)
#rename tip labels
new_labels <- str_replace(pan1_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
pan1_tree$tip.label <- new_tiplabels
#plot tree
ggtree(pan1_tree, size=.5)+
  geom_tiplab(size=10)+ xlim(NA, 1.5) + ylim(0, 37)
#save figure
#ggsave("PAN1.svg", width = 18, height = 18, dpi = 1200)

### PAN2 ###
#load phylogeny
pan2_tree <- read.tree("PAN2.nwk")
pan2_tree <- midpoint.root(pan2_tree)
#rename tip labels
new_labels <- str_replace(pan2_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
pan2_tree$tip.label <- new_tiplabels
#plot tree
ggtree(pan2_tree, size=.5)+
  geom_tiplab(size=10)+ xlim(NA, 1) + ylim(0, 25)
#save figure
#ggsave("PAN2.svg", width = 18, height = 18, dpi = 1200)

### PIR1212 ###
#load phylogeny
pir_tree <- read.tree("PIR121.nwk")
pir_tree <- midpoint.root(pir_tree)
#rename tip labels
new_labels <- str_replace(pir_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
pir_tree$tip.label <- new_tiplabels
#plot tree
ggtree(pir_tree, size=.5)+
  geom_tiplab(size=10)+ xlim(NA, .5) + ylim(0, 32)
#save figure
#ggsave("PIR.svg", width = 18, height = 18, dpi = 1200)

### POLAR ###
#load phylogeny
polar_tree <- read.tree("POLAR.nwk")
polar_tree <- midpoint.root(polar_tree)
#rename tip labels
new_labels <- str_replace(polar_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
polar_tree$tip.label <- new_tiplabels
#plot tree
ggtree(polar_tree, size=.5)+
  geom_tiplab(size=9)+ xlim(NA, 2.5) + ylim(0, 58)
#save figure
#ggsave("POLAR.svg", width = 20, height = 20, dpi = 1200)

### SCAR ###
#load phylogeny
scar_tree <- read.tree("SCAR.nwk")
scar_tree <- midpoint.root(scar_tree)
#rename tip labels
new_labels <- str_replace(scar_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
scar_tree$tip.label <- new_tiplabels
#plot tree
ggtree(scar_tree, size=.5)+
  geom_tiplab(size=10)+ xlim(NA, 3) + ylim(0, 40)
#save figure
#ggsave("SCAR.svg", width = 18, height = 18, dpi = 1200)

### SCARECROW ###
#load phylogeny
scr_tree <- read.tree("SCR.nwk")
scr_tree <- midpoint.root(scr_tree)
#rename tip labels
new_labels <- str_replace(scr_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
scr_tree$tip.label <- new_tiplabels
#plot tree
ggtree(scr_tree, size=.5)+
  geom_tiplab(size=8)+ xlim(NA, 1.5) + ylim(0, 40)
#save figure
#ggsave("SCR.svg", width = 18, height = 18, dpi = 1200)

### SHORTROOT ###
#load phylogeny
shr_tree <- read.tree("SHR.nwk")
shr_tree <- midpoint.root(shr_tree)
#rename tip labels
new_labels <- str_replace(shr_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
shr_tree$tip.label <- new_tiplabels
#plot tree
ggtree(shr_tree, size=.5)+
  geom_tiplab(size=8)+ xlim(NA, 1) + ylim(0, 40)
#save figure
#ggsave("SHR.svg", width = 18, height = 18, dpi = 1200)

### YODA ###
#load phylogeny
yda_tree <- read.tree("YODA.nwk")
yda_tree <- midpoint.root(yda_tree)
#rename tip labels
new_labels <- str_replace(yda_tree$tip.label, "_", " ")
new_tiplabels <- str_replace(new_labels, "!", " | ")
yda_tree$tip.label <- new_tiplabels
#plot tree
ggtree(yda_tree, size=.5)+
  geom_tiplab(size=8)+ xlim(NA, 1) + ylim(0, 63)
#save figure
ggsave("YODA.svg", width = 18, height = 18, dpi = 1200)
