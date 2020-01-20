# SoupX

An R package for the estimation and removal of cell free mRNA contamination in droplet based single cell RNA-seq data.

The problem this package attempts to solve is that all droplet based single cell RNA-seq experiments also capture ambient mRNAs present in the input solution along with cell specific mRNAs of interest.  This contamination is ubiquitous and can vary hugely between experiments (2% - 50%), although around 10% seems reasonably common.

There's no way to know in advance what the contamination is in an experiment, although solid tumours and low-viability cells tend to produce higher contamination fractions.  As the source of the contaminating mRNAs is lysed cells in the input solution, the profile of the contamination is experiment specific and produces a batch effect. 

Even if you decide you don't want to use the SoupX correction methods for whatever reason, you should at least want to know how contaminated your data are.

**NOTE:** The package focus has changed from v1.0.0 onwards.  The basic syntax is mostly the same as before, but there are three key differences in approach.  Firstly, you can no longer load multiple channels simultaneously and the concept of a "SoupChannelList" has been removed.  This is emphasises that the contamination rate should be calculated independently for each channel.  Secondly, the heuristic tools for guiding the choice of genes to use in estimating the conamination have been reworked.  In particular, they will no longer explicitly return gene names in any fashion.  This is to prevent users from using these genes blindly **which they absolutely should not do**.  The `plotMarkerDistribution` function will still suggest genes, but they will only be plotted, not returned directly in any way.  Finally, the default estimation has been switched to estimate one global contamination fraction across all cells in a channel.  In most cases there is little evidence of a strong per-cell difference in contamination rate and even where this does exist, there is seldom the power to calculate it.

## Installation

The package can be installed by running

```R
devtools::install_github("constantAmateur/SoupX")
```

If you encounter errors saying `multtest` is unavalibale, please install this manually from bioconductor with:

```R
BiocManager::install('multtest')
```

## Documentation

The methodology implemented in this package is explained in detail in [this paper](https://doi.org/10.1101/303727).  

A detailed vignette is provided with the package and can be viewed [here](https://cdn.rawgit.com/constantAmateur/SoupX/master/doc/pbmcTutorial.html).  

## Frequently Asked Questions

### I can't find a good set of genes to estimate the contamination fraction.

Generally the gene sets that work best are sets of genes highly specific to a cell type that is present in your data at low frequency.  Think HB genes and erythrocytes, IG genes and B-cells, TPSB2/TPSAB1 and Mast cells, etc.  Before trying anything more esoteric, it is usually a good idea to at least try out the most commonly successful gene sets, particularly HB genes.  If this fails, the `plotMarkerDistribution` function can be used to get further inspiration as described in the vignette.  If all of this yields nothing, we suggest either leaving your data uncorrected or trying a range of corrections to see what effect this has on your downstream analysis.  In our experience most experiments have somewhere between 2-10% contamination.

### `estimateNonExpressingCells` can't find any cells to use to estimate contamiantion.

At this point we assume that you have chosen a set (or sets) of genes to use to estimate the contamination.  The default behaviour (with 10X data) is to look for cells with strong evidence of endogenous expression of these gene sets in all cells, then exclude any cluster with a cell that has strong evidence of endogenous expression.  This conservative behaviour is designed to stop the over-estimation of the contamination fraction, but can sometimes make estimation difficult.  If all clusters have at least one cell that "looks bad" you have 3 options.
1. Recluster the data to produce more clusters with fewer cells per cluster.  This is the preferred option, but requires more work on the users part.
2. Make the criteria for declaring a cell to be genuinely expressing a gene set less strict.  This seldom works, as usually when a cell is over the threshold, it's over by a lot.  But in some cases tweaking the values `maximumContamination` and/or `pCut` can yield usable results.
3. Set `clusters=FALSE` to force `estimateNonExpressingCells` to consider each cell independently.  If you are going to do this, it is worth making the criteria for excluding a cell more permissive by decreasing `maximumContamination` as much as is reasonable.

## Changelog

### v1.0.0

Review of method, with focus on simplification of code.  Functions that were being used to "automate" selection of genes for contamination estimation have been removed as they were being misused.  Clustering is now used to guide selection of cells where a set of genes is not expressed.  Default now set to use global estimation of rho.  A hierarchical bayes routine has been added to share information between cells when the user does use cell specific estimation.  See NOTE for further details.

### v0.3.0

Now passes R CMD check without warnings or errors.  Added extra vignette on estimating contamination correctly.  Changed the arguments for the interpolateCellContamination function and made montonically decreasing lowess the default interpolation method.  A number of other plotting improvements.

### v0.2.3

Added lowess smoothing to interpolation and made it the default.  Modified various functions to allow single channel processing in a more natural way.  Some minor bug fixes.

### v0.2.2

Integrated estimateSoup into class construction to save memory when loading many channels.
Added function to use tf-idf to quickly estimate markers.
Some minor bug fixes and documentation updates.

### v0.2.1

Update documentation and modify plot functions to return source data.frame.

### v0.2.0

A fairly major overhaul of the data structures used by the package.  Not compatible with previous versions.

### v0.1.1

Some bug fixes to plotting routines.

# License

```
Copyright (c) 2018 Genome Research Ltd. 
Author: Matthew Young <my4@sanger.ac.uk> 
 
This program is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License version 3 
as published by the Free Software Foundation. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details <http://www.gnu.org/licenses/>. 
```
