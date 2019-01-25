# Keystone 2019 Day 3: 2019-01-16

## 8:00am - Stephen Quake - Stanford and CZI Biohub
Single Cell Genomics - A tool for drug discovery

Characterizing IgE producing B cells in PMBCs of allergic individuals
IgE is one of the rarest antibody isotypes, usually the culprit for allergies.
Primary purpose is for antihelmeth activity.

Enriching this rare population with antibody sorts/FACS is challenging.  Did
scRAN-seq + VDJ analysis.  From the VDJ analysis, they construct the antigen
receptors.
-   IgE B cells are predominantly plasmablasts (PBs).
-   IgE antibody diversity mimics other class switched isotypes.
-   Isolated 6 antibodies from their data, and found that they are highly
    mutated, and tested their binding affinity to peanut antigens.  They bound
    very highly and were cross reactive for multiple peanut antigens.  Appears
    the mutated light chain is needed for this high binding affinity.
-   They can further design antibodies starting from these naturally found
    antibodies that have picomolar cross-reactive binding affinity to these
    peanut allergies

The Tabula Muris


## 8:30am - Bart Deplancke - Swiss Federal Institute of Technology
Low-input single- and intra-cell RNA sequencing technologies

asap.epfl.ch

Developing low input sequencing tech so that he can study single samples (not
pooled data) and support clinical samples that have low cell counts.

Drop-seq and InDrop require a base minimum of 10k and 2k input cells. Drop-seq
has 90% cell loss, 10X has 50%, and InDrop has 25%.

Developed DeterminIStic CO-encapsulation (DisCo).  Designed to remove the
single/double Poisson processes involved in cell-bead capture.  Uses machine
vision to have cells/beads line up at gates which only open when both are present.
Co-encapsulation efficiency is ~80%. However, it's slow:  400-500 cells/hr.
Bead processing also leads to substantial loss of beads; you lose <60% of beads
during washing.  Had to reinvent a new bead recovery mechanism which reduces
bead loss to only 10%.  Also allows for highly diluted samples.

Due to machine vision, they can do size/shape sorts on the chip.  In the future,
this could be used to detect certain types of fluorescence, etc.

Dream is to link high-resolution image of the cell to the transcriptome (via
image recognition of the bead).

Live-seq: scRNA-seq while keeping cells alive
Using a force-microscopy head, you can carefully puncture individual cells,
apply counter pressure and suck up small amounts of fluid/molecules from inside
the cell.  With 0.1 - 3.2 pL capture, 82% of cells remain viable.

Using a modified Smart-Seq2, they can sequence sub-picoliter volumes of mRNA.
So far only able to capture ~1000 genes per cell, and only ~14 cells.


## 9:30am - Georg Seelig - UWash
Split-Pool barcoding

Split-seq ~11 million barcodes.  Validates with 1800 cells.  With only 4000
cells, you get 50% chance of collisions.  They say 1% doublet rate, meaning 10k
doublets.

Split-seq does better than 10X for nuclei but slightly worse for cells.  Both
are worse than smart-seq2.   However, all technologies capture sufficient
information to identify cell types (of GABAergic neurons at least).  Comparison
was done at Allen Brain Institute.

Sample multiplexing is natural with split-seq, as first found is with RT primer,
which can be used as a sample identifier. So put each sample into a different
well to start.

Technology also lets you fix samples then freeze, and then sequence later, so
sample prep seems much more relaxed.

This is being commercialized: www.splitbiosciences.com

Poster 3047: pairing split-seq with TCR/BCR profiling.  Does not rely on PCR or
pull-down with biotinalated beads.


## 10:00am - Alex K Shalek - MIT/Broad
Identifying and rationaly modulating cellular drivers of enhanced immunity

Studying when your immune system works and how vaccines/protective measures help
when your immune system doesn't work.

Looking at granulomas from mTB infection (Travis Hughes, MD/PhD student).
Very hard to get healthy controls and finding granulomas is difficult.
Collaborating with Joanne Flynn at Pitt for Macaque model of mTB.

Using Seq-Well protocol to study these data.  Collected 40 granulomas from 6
animals, and there is substantial variability in bacterial burden among
granulomas from the same animal.  They can then correlate cell type numbers with
bacterial burden.  Macrophages, plasma cells, and mast cells are associated with
higher bacterial burden. Finds Th17 cells in addition to T-effector cells that
are associated with bacterial control.

The TB vaccine isn't 100% effecitve.  He wants to know why.  Turns out the
vaccine is administered in different ways and animal models suggest that IV
innoculation is the most effective.  Compares scRNA-seq where the conidition is
the vaccine route.

Shalek lab, multi-species, multi-dataset database that is queriable.  To be
public "soon".

## 10:30am--11:15am - Morning short talks

### 10:30am - Gary Hon - UT Southwestern Medical School
Reprogramming via TFs.

Reprogram-seq: mouse embryonic fibroblasts, hit them with a TF retroviral
library and sequence them.

### 10:45am - Bre D'Andreth - MIT
One-pot transfection method for rapid characterization and optimization of
genetic systems.

### 11:00am - Luca Rappez - EMBL
SpaceM, Spatial single-cell profiling of intracellular metabolomes in situ

Poster 3041

[preprint](https://www.biorxiv.org/content/early/2019/01/02/510222)


## Garry Nolan couldn't make his talk.

## 5:00pm - Xiaowei Zhuang - Harvard
Single-cell transcriptome and chromosome imaging

Ad for MERFISH [protocols](http://zhuang.harvard.edu/merfish.html)
-   sequential imaging with combinatorial labeling
-   yields binary expression matrix with 2^N rows and N columns
-   uses FISH to add probes in a error-correcting manner through a seperate
    readout
-   detection efficiency is >90%.

Spatial and functional atlas of mouse hypothalamus preoptic region using MERFISH
and 10X
-   9 major cell classes, >60 neuronal cell populations from 10X
-   yields marker genes for MERFISH imaging
-   imaged 155 genes with MERFISH, 85 from previous knowledge, 70 from 10X
-   were able to image 1.8mm x 1.8mm x 0.6mm per animal across multiple mice
    yielding 1.1 cells imaged.
-   strong correlation between MERFISH and scRNA-seq

Exposed virgin/parent male/female combination mice to parental conditions
-   use cFos as an early marker for activation.
-   looking for genes that are co-expressed with cFos in 70+ neuronal cell types
-   identified essentially a single "parenting cluster" that's activated across
    all combinations except virgin males, who show aggression instead. Parent
    males show aggression but it is outweighed (either by being downstream or
    directly inhibited) by parenting cluster.

Tracing 3D conformation of chromatin
-   multiplexed FISH identifies genomic loci and links them together (how?
    through ordering of specific FISH probes?)
-   Can get down to 30kb resolution.
-   Validate with HI-C data, but since it's single cell data, they can look at
    single cell contact maps.


## 5:30pm - Rick Horwitz - Allen Institute for Cell Science
Conjoining single cell imaging, genomics, and computaion to create a stem cell
state space

3 principles of Allen Institute for Cell Science
-   High dimensional cell atlas
-   Principles of cell organization
-   Gene regulatory circuits and morphogenesis

Top down approach -> combining genomics, organelles, phenotypes across scales.
They do it all in stem cells.  "Dynamic, data-driven virtual cell"

Cell Collection: FP-tagged hiPCS lines
Imaging uses automated high replicate microscropy with 4 channels
-   membrane
-   gene edited structure
-   nucleus
-   bright field

Segementation: Allen Cell structure segmenter
[github](https://github.com/AllenInstitute/aics-segmentation)
-   works across a variety cell morphologies
-   has iterative ML process

Analyses of cell variation: ANOCVA
-   all data is available
-   including the segmented images

UNet trained model to identify moitosis structure in bright-field.  Using this,
don't need fluorecence and can produce dynamic, 3-D images of multiple cells
containing multiple subcellular structures.

Trained a model to generate a continuous landscape of cellular and sub-cellular
morphologies.

Next steps involve integrating cell states from cellular neighborhoods,
integrating transcriptomics -> FISH -> imaging.


## 6:15pm - Remco Loos - Celgene Institute for Translational Research Europe
Understanding heterogeneity in AML through bulk and single-cell gene expression
profiling

intra-patient profiling: scRNA-seq
inter-patient profiling: bulk RNA-seq


## 6:40pm - Fan Zhang - Harvard Med
Defining inflammatory cell states in rhumatoid arthritis

[preprint](https://www.biorxiv.org/content/early/2018/06/20/351130)


## Day 3 reflection

### Personal Chats

Why don't we always align to pre-mRNA references? Aren't nuclear mRNAs present
in the lysed cell?

Harmony guy is open to someone implementing in python, either by
reimplementation or by just writing a frontend to the C++ backend.

Ttyh2 Mature Oligodendro marker
Pdgfra Immature Oligodendro marker
