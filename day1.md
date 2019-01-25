# Keystone 2019 Day 1: 2019-01-14

## 8:00am - Keynote speaker - Scott Fraser - USC
Multimodal multiplex biomedical imaging

Problem: single cell analyses necessitate that cellular dynamics are lost and
normal context is violated.  Results in more heterogeneity than is likely there
in vivo.

Solution: studdy embodied systems, embrace dynamics, variance.

Tried combining HCR (hybridization chain reaction) with IMC, but HCR signal was
lost.  Instead, constructed MUSE (multimodal universal signal enhancement) where
hairpins contain multiple offshoots where you can add dyes, heavy metals, etc.
This allows multiple signals without destruction, meaning the scan rate of the
microscope is 10,000x faster than IMC. 5 gigavoxels takes ~10 hours, and this is
the low end.  Additionally, number of colors is limited.

Tried: Light sheet microscopy (2 photon SPIM). Very low phototoxicity.  Great
signal to noise ratio.  100x faster than MUSE.  Tunable excitation depths. Can
use multiple channels if you know what the base spectra are via fourier space
unmixing.  However, the problem is that timesteps between multiple channel
capture loses dynamic information; synchronicity is lost between 3D channels.

Tried light field microscopy to overcome the synchronicity problem.  Selective
volume imaging + LF microscropy (LF-SVIM).  ASOSVIM.

Have to use at least 2 channels for nucleus and cell boundary detection. Cell
segmentation is still a big problem for this technology.

Lattice Light sheet at Genelia?  Lower field of view but excellent performance.

## 9:00am - Fabian Theis - Munich
Modeling differentiation and stimulation response in single cell genmoics

Showing PAGA data that's all available on github.
Showing scvelo data.  Arguing that rna-velocity is deterministic steady state
whereas scvelo is a stochastic model. Uses scvelo to change PAGA connectivities
to directed graph structure.

Autoencoders to learn lower dimensional representation of count matrices and to
denoise data.  Stressing that finding the correct latent space is of utmost
importance.
[https://github.com/theislab/scGen](https://github.com/theislab/scGen)

scGen can be used to correct batch effects as well.  Purportedly better than
bbknn.  But what about performance?

How does he handle orthologs (homologs?) across species?

## 10:00am - Rahul Satija - NYGC
Comprehensive single cell integration

STARmap - transcriptomics in situ (Garry Nolan's lab).

Leveraging aspects of CCA and mutual nearest neighbor matching for batch correction.

Determining technical vs biological is difficult and even harder is a natural
biological effects vs sample preparation biological effects.

Imputing CITEseq data in datasets that are scRNA-seq only by CCAing in other
datasets that have corresponding CITEseq data.  Planning to use HCA as reference
datasets to map extra modalities into.

Finding nearest neighbors between scATACseq and scRNAseq
-   Assumes correspondence between gene and modules of chromatin accessibility.
    Essentially, that one ATACseq peak isn't indicative of transcription, but
    multiple peaks within a gene module implies transcription.

Integrating sequencing and imaging data using STARmap.  STARmap has limited gene
capture (~1000 genes).  Use scRNAseq to impute full gene expression in STARmap
data.

Can try to do the opposite: impute spatial location of scRNAseq clusters using
STARmap datasets.

Guy admits that previous versions of multiCCA sucked.  Claims previous version
is better.

## 10:30am - Dana Pe'er - Sloan Kettering Institute
Mapping the emergence of organ identities in time and space

Cell states are points on phenotypic manifolds
Cell fates should also live on that manifold

Discussing Palantir.  Treating cell fate as a probability.  Trying to separate
out uncertainty due to true stochasticity or missing data.
-  Strong assumption: cells can't go back in time and differentiation is a
   markov process.
-  Entropy of cell fate probability is a measure of a cell's
   "differentiability" or "differential potential".
[https://github.com/dpeerlab/Palantir](https://github.com/dpeerlab/Palantir)

Manifold alignment techniques assume you're looking at the same
thing/timepoints.  Use a MNN correction to align different timepoints that
integrates between sample correlations and within sample correlations.
[https://github.com/dpeerlab/Harmony](https://github.com/dpeerlab/Harmony)

**Interesting point**: real trajectories are usually made up of small steps
(pseudotime or transcription changes) whereas trajectories with large/variable
step sizes are likely not real.  How does this affect development of
chemosensing astrocytes in Colin's data?

AP pseudo-space?

## 11:15am--12:00pm - Morning short talks

### 11:15am - Chuangqi Wang - Worchester Polytechnic Institute (Phd student)
Studying heterogeneity in cell protusion.

### 11:30am - Hanchuan Peng - Allen Institute for Brain Science
3D neuron imaging at singlecell level.

Vaa3D image software


## 1:30pm--4:00pm - Afternoon rapid fire talks

###  Robin Browaeys - NicheNet - Belgium
Inference of intercellular communication

Input: expression data.
Method: uses 20 open databased on gene regulatory networks, ligand/receptor
interactions, Uses network propagation methods (like what? belief propagation?)
Output: predicts genes effected by ligand/receptor expression.

Case study: tumor microenvironment. Which CAFs regulate the expression of p-EMT
in malignant cells?

Poster 1028.

### Andrew Butler - Comprehensive integration of sc data - NYGC
Goals: create larger reference datasets and identify common cell types across
datasets.

Problem: transferring cell type annotations from one type from another.  Placing
new data onto a reference.

Poster 1030

### George C Linderman - Zero-preserving imputation - Yale
Trying to disentangle technical vs biological zeros in expression matrices.

ALRA method. Rank k-svd approximation + adaptative thresholding

Assumption for thresholding: Values corresponding to biological zeros are
symmetric around zero.  So if you threshold at the abs value of the largest
negative value, this is should correspond to the largest incorrectly imputed
biological zero.  Everything else should be technical zeros.

Poster 2056.

### Ilya Korsunsky - Batch effect correction with Harmony - Harvard

[preprint](https://www.biorxiv.org/content/early/2018/11/05/461954)

Linear mixture models to account for batch effects. Iterative method that first
estimates shared cell types, then uses linear regression to regress out batches.

Poster 2039

### Gan Yu - Integrating differential expression and gene set analysis - UMich

Method name: iDEA

Thinks current methodology which feeds DE genes into GSEA isn't correct as both
methods are related.  Instead, proposes iterative analysis.  Seems to provide
marginal improvement over existing methods in finding DEs and enriched gene
sets.

Poster 4061

### Gregoire Altan-Bonnet - hierarchical agglomaterative clustering in CyTOF- NCI

Method name: HAL-x or ALPS

Poster 1008

### Shiyi Yang - Removing contamination from ambient RNA in scRNA-seq - Boston University
Bayesnian mixture model

[github](https://github.com/compbiomed/celda/blob/master/R/decon.R)

Assumptions: each cell contains two distributions of counts: 1 for true
expression and another for counts by contaminations.

Method name: DecontX

Poster 4054

### Gregory Schwartz - Too Many Cells - U Penn

[github](https://github.com/GregorySchwartz/too-many-cells)

Hierarchical clustering algorithm and associated visualization tools.

Poster 3057

### Rongxin Fang - Single cell regulomes - UCSD/Salk

[github](https://github.com/r3fang/snATAC)

SNAPATAC - comprehensive solution for analysis of scATAC-seq

Poster 1056

### Junil Kim - Reconstruction of gene regulatory network from pseudotime - U Copenhagen

New algorithm for large scale causal network reconstruction.

Poster 2033

### Corey Williams - Trajectory mapping using FLOW-MAP - UVA

[github](https://github.com/zunderlab/FLOWMAP)

Understand molecular processes that determine cell state and cell fate.

Use time course data to constrain trajectory graph construction. Prohibit knn
connections between non-adjacent timepoints.

Poster 4043

### Robrecht Cannoodt - Dynbenchmark: benchmark of pseudotime methos - Ghent University, Belgium

[github](https://github.com/dynverse/dynbenchmark)

Poster 1033

### Wouter Saelens - Dyno: Inferring and intepreting trajectories - Ghent University, Belgium

[github](https://github.com/dynverse/dyno)

Poster 3053

## 4:00pm--4:30pm - Panel discussion
Biggest challenges:
-   Alon Klein - Integration of reference datasets into future analyses.
-   Dana Pe'er - Mechanistic models utilizing time and spatial information
    and batch effects
-   Fabian Theis - Generative models and annotations.
-   Rahul Satija - Inference of various sorts, and quantifying uncertainty in
    methods development.

There aren't platforms for sharing and making accessible yet, but HCA hopefully
fill this need.

Dana thinks that methods are so context specific that finding the "best" method
might be an ill-posed problem.  Instead, people need to validate that using a
method on their data is a robust result.

Fabian thinks that it is critical we make our code/figures/validation as
accessible as possible.

Alon thinks that it's important to identify metrics that identify if your method
works well.  And possibly to show/describe when/how your method breaks and is
limited.


## Day 1 reflection

### Personal chats
MD/PhD students Andrew and Joel from Mount Sinai.  Guy from MARS-seq
implementing maximum likelihood method to embed cells into k-means clustering
generated using a subset of a dataset.

Lior Patcher lab doing non-antibody tags. Jase Gehring
[preprint](https://www.biorxiv.org/content/early/2018/05/05/315333)


### Other thoughts
Can we do leave one out cross validation of bbknn data.

Can I talk to 10X people about missing mTomato/mCherry transcripts.
