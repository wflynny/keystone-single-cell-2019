# Keystone 2019 Day 4: 2019-01-17

## 8:00am - Thomas Hofer - Deutsches Krebsforschungszentrum, Germany
Modeling Population Dynamics of Lineage Pathways through Barcoding Analysis

Using population models to study fate mapping in differentiation
n = exp(f(proliferation, differentiation, migration, death) * t)

Number of HSCs in development is controversial topic.  To find out, they used
noninvasive genetic barcoding in HSCs using inducible Cre-polylox system.

Prob(barcode generation probability) matches what they see.

HSCs appear to come from a small number of precursors and proliferate greatly.


## 8:45am - Barbara Wold - California Institute of Technology, USA
Following Cellular Reprogramming at the Single Cell Level

scRNA-seq of forelimb development.


## 9:40am - Wolf Reik - Babraham Institute, UK
Single Cell Epigenome Landscape of Development and Aging

Epigenetics confers memory

scNMT-seq - both transcriptome, chromatin accessibility, and methylome

Can they resolve the timing of epigenetic and transcriptomic events

MOFA - multiomic factor analysis

[preprint](https://www.biorxiv.org/content/early/2019/01/13/519207)


## 10:15am - Fernando D. Camargo - Boston Children's Hospital, USA
Cellular Barcoding in Mammalian Systems

Use transposons to tag cells.  All cells have the same donor locus of the
transposon, but when doxycyclin is added, it reinserts somewhere completely
random.  They consider the adjacent DNA to the reinserted transponon to be the
tag.

HSCs are primarily megakaryocyte progenitors whereas MPPs, not, HSCs are
pregenetiors of most other lineages.

## 10:50am--11:45am - Morning short talks

### 10:50am - Aziz Al'Khafaji - University of Texas, USA
Short Talk: Studying Therapeutic Resistance in Chronic Lymphocytic Leukemia
using Functionalized Lineage Tracing

Method name: COLBERT

### 11:05am - Bushra Raj - Harvard University, USA
Short Talk: Reconstruction of Brain Specification and Lineage Trees with
Single-Cell Profiling and CRISPR Recorders

URD - new pseudotime algo
[paper](http://science.sciencemag.org/content/360/6392/eaar3131)

scGESTALT - lineage tracing methodology



## 2:30pm - Keyue Ma - UT Austin
High-throughput single-cell linking of antigen specificity to TCR sequences

No fast/cost-effective methods to generate tagging peptides
New tech: TetTCR-seq, [paper](https://www.nature.com/articles/nbt.4282)

## 3:00pm - Chris McGinnis - UCSF
MULTI-seq: scalable sample multiplexing for single-cell RNA sequencing
[preprint](https://www.biorxiv.org/content/early/2018/08/08/387241)

Wants to move sample multiplexing prior to droplet microfluidics devices
(instead of at the end when the i7 indices are added).

Lipid-modified-oligos (LMOs)
-   Scalable: modular rapid quenchable cheap
-   Compatible with nuclei since it's targeting the lipid membrane
-   Compatible with all species/whatever

So far, doesn't work too well on fixed samples yet. Click-tags, however, work
super well on fixed samples.

In conjunction with doubletfinder or scrublet, you can essentially remove all
doublets since cells are labeled beforehand.

96 samples cost $5k to prep vs. $120k with vanilla 10X.

Not too much exchange within ~6hrs in terms of lipid exchange.

Can't make reagents at scale yet, but get in contact:
chris.mcginnis@ucsf.edu


## 3:30pm - Arnav Moudgil - WashU
Simultaneous TF-binding assays and RNA-seq at single cell resolution

Method: Calling Cards

Transposons are carrying reporters and are linked to TFs. When TFs bind, the
transposon inserts itself (and thus the reporter) nearby.  For scRNA-seq, they
need to circularize with biotin pull-down.

They can make the reporter fluorescent to do this in vivo.


## 3:45pm - Manu Singh - Garvan Institute of Medical Research, Australia
Linking high-throughput single cell RNA_seq with long read TCR/BCR transcripts.

Method: RAGE-seq
[preprint](https://www.biorxiv.org/content/early/2018/09/24/424945)

After 10X cDNA library generation, they do hybridization capture and Oxford
nanopore sequencing in addition to short read sequencing.

They do de novo assembly of nanopore reads on a per-cell basis.

Validated usind Jurkat and Ramos cell lines with known antigen receptors, and
they get 98% accuracy.


## 4:00pm - Chengzhe Tian - University of Colorado Boulder, USA
A "Global Tracker" for Hard-to-Track Cancer Cells Reveals Substantial
Heterogeneity in the Dynamics of Single-Cell Drug Responses

One parameter set for segmentation across all frames.


## 4:15pm - David M. Suter - École Polytechnique Fédérale de Lausanne, Switzerland
Single Live Cell Monitoring of Protein Turnover Reveals Intercellular
Variability and Cell Cycle Dependence of Degradation Rates

Using tandem slow/fast maturing fluorescent proteins so that fluorescent signal
changes color over time.

Looking at protein synthesis and degradation rates, they can show that protein
accumulation in cells after cell division is predominantly due to increases in
synthesis rates and not decreases in degradation rates.  The opposite is true
before and during cell division.

Synthesis and degradations rates within individual cells are highly correlated
which allows cells to buffer variability in gene expression.

Degradation rates aren't really gene-specific


## 4:30pm - Emily B. Fabyanic - University of Pennsylvania, USA
Probing the Mammalian Hydroxymethylome at Single-Cell Resolution Using a DNA
Deaminase

Method: ACE-seq
Leveraging APOBEC3A which deaminases cytosine and methylctyosine with much more
efficiency than hydroxy-methylcytosine, so that they can see epigentic signatures
due to hmC vs mC and C, which in the neurons is the terminal methylated state of
cytosine (for the most part).

bACE-seq (bisulfite-assisted ACE-seq) confers protection to hmC.


## 5:00pm - Hongkui Zeng - Allen Institute for Brain Science
Cell Type Classification in the Mouse Brain

## 5:45pm - Barbara Treutlein - Max Planck Institute, Germany
Single-Cell Transcriptomics Uncovers Convergence of Cell Identities during
Axolotl Limb Regeneration


## 6:20 - Farimah Mapar, MIT
Short Talk: Design and Analysis of Staged Mutual Inhibition: Using Single Neuron
Computation to Implement Bi-Stable Neuronal Toggle Switch

This was the most incomprehensible talk I've ever witnessed.  And she's standing
in the way of food and beer and dance party!


## Day 4 reflection

### Personal Chats

Can we define doublets in the way shown in the poster of Decibel Therapeutics?
How many fake doublets do we generate?
```
doubs = adata.copy()
sizes = adata.obs.groupby(cluster_key).size()

for cluster in adata.obs[cluster_key].unique():

```
