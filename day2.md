# Keystone 2019 Day 2: 2019-01-15

## 8:00am - Lucas Pelkmans - University of Zurich
Understanding emergent phenomena across biological scales

Across scales, studying similar phenomena: emergent self-organizing properties.

Tracking mRNA transcripts across thousands of cells on the scale of millions of
sub-cellular objects.

He's concerned with preventing photoactivated crosslinking of antibodies so that
you can effectively elute your antibodies without destroying your sample.
Crosslinking occurs because you flurophore emission sometimes releases singlet
oxygens which activate crosslinking reactions.  His solution is to add a soluble
O- acceptor that prevents crosslinking.  This allows repeated staining and
elution without loss of signal.

Method name: 4i

Each pixel in their 4i data contains 40 intensity measurements of different
probes.  They uses these vectors to cluster each pixels into units and then maps
these units back into the spatial representation of the cell.  This effective
coarse graining allows you to understand proximity of units to others, classify
entire cells by their makeup of units, etc.

Perhaps you don't even need segmentation to see patterns.  This has implications
for Santhoshes data.

## 8:30am - Ellen Rothenberg - CalTech
Making Tcells from multipotent progenitors and understand commitment decisions
in differentiation fates.

They have a single gene, Bc11b, that is tightly correlated with T cell
commitment that is part of a feed-forward gene regulatory network.

## 9:30am - Caroline Uhler - MIT
Mathematical and computational tools to study high throughput measurements of
single cells.

Causal inference and gene regulatory networks. Sets of nonlinear equations to
describe a causal graph modeling a regulatory network.

## 10:00am - Tim Schroeder - ETH Zurich
Be cautious in interpreting dynamics from snapshot data.

[The (Timm's) Tracking Tool](https://www.bsse.ethz.ch/csd/software/ttt-and-qtfy.html)

## 10:45am--12:00pm - Morning short talks

###Rene Maehr - UMass Medical School

Poster 3003

Human embryo thymus development using scATAC- and scRNA-seq.

### Min Xue - UC Riverside

Poster 4???

Building FRET-based reporting system based on structural biology to
understanding cell-to-cell signaling.

### Masahiro Ueda - RIKEN
Automated live single-cell imaging analysis of intracellular signal transduction
with single-molecule sensitivity

Uses Dicty as a model system.

PTEN and PIP3 are expressed in a mutually exclusive manner on the membrane.
This forms a somewhat molecular compass.


## 5:00pm - Cristina Lo Celso - Imperial College London
Healthy and malignant haematopoiesis in the bone marrow: dynamic cells in an
evolving environment

Studying stress in HSCs using confocal and 2-photon microscopy.  Specifically
-   how HCSs competes with Leukaemia,

    T-ALL and AML (their models for Leukaemia) cells are mobile (moreso after
    chemotherapy) and highly express CXCR4. Inhibiting CXCR4 mobilizes these
    Leukaemia cells in the bone marrow.  However, these two models have different
    interactions with the bone marrow microenvironment and CXCR4.

-   how HCSs and the bone marrow microenvironment changes during infection
    Malaria infection.  Ran 10X on infected and control mice.  Not really
    doing a controlled analysis of case vs. control.  Seems to making
    conclusions from relative cell numbers in datasets that haven't been
    aggregated and clustered together.

-   how HCSs behave early after engraftment.


## 5:30pm - Prisca Liberali - FMI, Basel
Self-organization in intestinal organoid development

"Self-organization, symmetry breaking, and collective behavior" - molecular and
cellular mechanisms that drive emergent behavior in cell populations.

Calls changes in population's macroscropic properties "symmetry breaking".  What
is the reaction coordinate and critical point?  How to define?

Using 4i to image organoids. Segment with CNNs at nuclei, cell, and organoid
level.

Symmetry breaking isn't universal in organoids. Some never phase change.
Organoid fate appears to be time sensitive; if organoid doesn't generate stem
(planet) cell by some time threshold, it never will.  Using light sheet
microscopy imaging the developing organoid from 1 cell for a week, they can use
dynamic time warping to align pseudotime and real time.  This time threshold
seems to be around 48 hours.

They identified a gene, Yap1, which is needed for intesinal regeneration that
is expressed variably during organoid development.  Both knocking it out and
making its expression isotropic inhibit symmetry breaking. Suggests variable
expression of this gene is required for symmetry breaking.

They hit ~500,000 organoids with a Novartis drug library (2+ compounds per
condition), imaged with on-the-fly organoid identification and segmentation, and
clustered organoid morphology to define a landscape of development organoid
phenotypes.  This can be used to generate a drug-gene-phenotype mapping.


## 6:00pm - Sinem Saka - Harvard
Highly multiplexed in situ protein imaging with signal amplification

[preprint](https://www.biorxiv.org/content/early/2018/12/28/507566)

Signal amplification by exchange reaction (SABER)
Uses primer exchange reaction which generations long concatemers with many
binding sites for imaging probes.

You stain with tagged antibodies all at once, then apply SABER, which
multiplexes the anitobody tags in one go. Through branching, they can get large
dynamic range of signal detection.

Can also apply this for
-   RNA/DNA detection where FISH probes are ligated to SABER concatemers.
-   expansion microscopy.

Poster 3054


## 7:00pm - Leo Kunz - ETH Zurich
Visualization and quantification of cytokines at the single molecule level in
situ

Poster 2045

## Day 2 reflection
### Personal Chats

A grad student at Mount Sinai suggested that capturing neutrophils with
scRNA-seq is challenging because as soon as they are lysed, they spill all their
mRNA?  Shouldn't that still be captured by the GEM?  He claims they will always
have very few counts, often below filtering thresholds.

Need to use veloctyo on Colin's data to see if there is actual flow between
immature astrocytes to chemosensing vs mature astrocytes.  Show pseudotime
coloring as `cmo.phase`.

Create perceptually uniform monochromatic colormaps starting at light gray and
ending at red/blue/green. Essentially gray to black through a color.

I wonder if Jeremy's interest in the topological spaces of phylogenetic trees
could be useful in understand lineage tracing/pseudotime.  Given a tree, a
knockout or other alteration to the system can be viewed as an operator to move
you from one area in the topology to another (a different tree).  What types of
transformations have the largest "distance" or smallest "distance"?  Is there a
similar topological space that spans cyclic graphs (could be used to study
regulatory networks)? What about using this to understand relationships between
different trees learned from different patients/replicates?

Discussed with 10X about:
1.  No tdTomato reads.  They took my contact info and will discuss with customer
    and technical support.
2.  Why did they switch to 91bp for read 2 in v3? They didn't know.  They likely
    think that the 100bp kits have more than 25bp leeway.
