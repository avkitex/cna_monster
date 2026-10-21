# Allele-specific copy number calling from capture sequencing: a concept for a Sequenza-based caller with FACETS SNP loci, diploid-aware model selection, and capture-profile (AUC) depth ratios

**Author:** Nikita Kotlov
**First public disclosure:** 21 October 2026
**Document type:** method description (design specification). No implementation or benchmark results are reported here.

## Summary

This document describes a method for allele-specific copy number (ASCN), purity and ploidy estimation from matched tumor/normal hybrid-capture sequencing. The method builds on the Sequenza likelihood model [1] and changes three components:

1. **Input loci (Version 1).** Read counts at a fixed set of common polymorphic loci, heterozygous and homozygous in the normal, as defined by FACETS [2], replace the full-pileup input of Sequenza.
2. **Model selection (Versions 1 and 2).** A ploidy prior, an explicit test of the 2× copy-number degeneracy, a whole-genome doubling (WGD) evidence test, and a post-fit refinement of purity from B-allele frequency (BAF) using several segment-weighting schemes and per-chromosome-arm aggregation.
3. **Depth ratio measurement (Version 2).** Per-target tumor/normal ratios of the area under the capture coverage profile, computed on capture-kit-specific filtered region sets, replace single-position depth ratios. BAF continues to be taken from heterozygous SNP loci.

## Contents

1. [Background](#1-background)
2. [Version 1](#2-version-1-facets-loci-into-sequenza-diploid-aware-fit-baf-based-post-processing)
3. [Version 2](#3-version-2-capture-kit-aware-depth-ratio-from-coverage-profiles)
4. [Relation to prior art](#4-relation-to-prior-art)
5. [Implementation outline](#5-implementation-outline)
6. [Evaluation plan](#6-evaluation-plan)
7. [References](#references)

---

## 1. Background

ASCN calling from tumor/normal exome or panel sequencing requires joint estimation of segment boundaries, tumor purity and tumor ploidy. Existing open-source tools differ in input signal and fitting strategy:

| Tool | Input signal | Purity/ploidy model | Component modified in this method |
|---|---|---|---|
| Sequenza [1] | Per-position depth ratio and BAF from a full pileup (`sequenza-utils bam2seqz`) | Likelihood grid over cellularity × ploidy, joint on depth ratio and BAF | Input positions; depth ratio definition; solution selection |
| FACETS [2] | Read counts at a fixed set of polymorphic loci (`snp-pileup`), heterozygous and homozygous in the normal | logR and log-odds-ratio, EM fit, diploid logR (dipLogR) selection | Loci set reused as input |
| CNVkit [3] | Mean depth per capture target and off-target bin | No joint purity/ploidy fit | — |
| PureCN [4] | Interval-level coverage on kit-specific intervals, panel of normals, SNV allelic fractions | Likelihood over purity/ploidy | — |
| ASCAT [5], ABSOLUTE [6] | logR and BAF (SNP arrays or sequencing) | Goodness-of-fit; karyotype priors | — |

The method is based on two properties of the data:

1. **Depth ratio and BAF have different error sources.** BAF at a heterozygous site is a ratio of allele counts within one sample and does not depend on capture efficiency. Depth ratio at a single position depends on capture efficiency, GC content, fragment length and position within the capture target. Depth at neighboring positions within one target is correlated, so per-position depth ratios are not independent observations.
2. **The 2× degeneracy.** A solution at purity *p* and ploidy *ψ* and a solution with all integer copy numbers doubled, at an adjusted purity *p′*, can produce near-identical expected depth ratios and BAF. The two are distinguishable only through segment states that have no counterpart under the other solution (for example, odd total copy number in the doubled frame, or LOH at two different total copy numbers). Section 2.2 tests for such states explicitly. Given a state assignment, expected BAF does not depend on ploidy, which Section 2.3 uses to re-estimate purity.

Version 1 changes the input loci and the solution selection. Version 2 additionally changes how depth ratio is measured.

---

## 2. Version 1: FACETS loci into Sequenza, diploid-aware fit, BAF-based post-processing

### 2.1 Input: replace the full pileup with FACETS polymorphic loci

- Run `snp-pileup` (FACETS) on tumor and normal BAMs against the common-SNP VCF used by FACETS (dbSNP common variants), keeping **both** sites heterozygous and sites homozygous in the normal.
- Classify sites by normal genotype: heterozygous if normal BAF is within an interval around 0.5 with minimum depth (FACETS-style thresholds, configurable); homozygous otherwise.
- Convert directly to Sequenza's `.seqz` format, one row per site:
  `chromosome, position, base.ref, depth.normal, depth.tumor, depth.ratio, Af, Bf, zygosity.normal, GC.percent, good.reads, AB.normal, AB.tumor, tumor.strand`
  - heterozygous sites: `zygosity.normal = het`, contribute depth ratio **and** BAF;
  - homozygous sites: `zygosity.normal = hom`, contribute depth ratio only;
  - `GC.percent` from the standard Sequenza GC wiggle (or computed over the read window around the site).
- Then run the standard Sequenza flow (`sequenza.extract` → segmentation → `sequenza.fit`) on this reduced, genotype-aware position set.

Properties of this input: positions are restricted to known polymorphic sites, the heterozygous/homozygous assignment is explicit, and the position set is identical across samples. The number of positions processed is smaller than in a full pileup.

### 2.2 Solution selection: diploid baseline unless rejected by the data

Sequenza returns a log-posterior surface over (cellularity, ploidy). The modification:

1. Add a log-prior on ploidy centered at 2 (configurable strength), or equivalently a penalty for solutions with ploidy > ~2.7.
2. **Degeneracy check.** For the best non-diploid solution, map every segment to its integer state (A, B). If all segment states are even (A and B both divisible by 2, or reducible to the diploid solution by halving), the diploid solution explains the data equally well and is selected (parsimony).
3. **Whole-genome doubling (WGD) evidence test.** Accept the higher-ploidy solution only if it is supported by states that cannot be represented in the diploid frame, in particular:
   - LOH segments at total copy number 2 (A = 2, B = 0, copy-neutral LOH after doubling) coexisting with LOH at total copy number 1;
   - odd total copy number states at high BAF resolution;
   - fraction of the genome with LOH above a configurable threshold.
4. Both solutions are reported with their log-likelihood difference; a flag is set when the difference is below a threshold.

### 2.3 Post-processing: BAF-driven refinement of purity and ploidy

For a segment with clonal tumor state (A, B) at purity *p*, normal 1:1, and tumor ploidy *ψ*:

$$
\mathrm{BAF}_{\text{minor}} = \frac{pB + (1-p)}{p(A+B) + 2(1-p)}, \qquad
R = \frac{p(A+B) + 2(1-p)}{p\,\psi + 2(1-p)}
$$

For a given state assignment, BAF does not depend on *ψ*; purity is therefore estimated from BAF, and *ψ* is then obtained from the depth ratio equation. Steps:

1. **Per-segment purity from BAF.** For each segment with an assigned state (A, B), the BAF equation is solved for *p* using the mirrored BAF of its heterozygous sites. For a fixed *p*, the deviation of expected BAF from 0.5 is largest for clonal LOH states (B = 0), so these segments carry the most information on *p*.
2. **Alternative segment weighting.** Purity and ploidy are re-estimated under several weighting schemes; the spread across schemes is reported as an uncertainty measure:
   - equal weight per segment;
   - weight by genomic length;
   - weight by number of informative heterozygous sites;
   - inverse-variance weight (binomial variance of BAF given depth and site count).
   Disagreement between schemes (for example, long segments and SNP-dense segments implying different purities) is reported as an indicator of a subclonal event or an incorrect state assignment.
3. **Per-chromosome-arm analysis.** Segments are aggregated by arm (p/q, using cytoband centromere coordinates for the reference build). Arm-level copy number changes are common in solid tumors [9] and provide large blocks with many heterozygous sites. For each arm:
   - estimate purity from its BAF/depth ratio pair;
   - test whether arms are consistent with a single purity under integer states;
   - arms that disagree with the consensus are flagged as candidate subclonal events and down-weighted in the global estimate.
4. **Consensus.** Final purity is a robust combination (median or trimmed mean) of arm- and segment-level estimates; ploidy follows from the depth ratio equation. If the consensus differs from the Sequenza grid optimum by more than a set tolerance, the grid is re-searched in a window around the consensus.

### 2.4 Outputs

- Selected purity and ploidy, with the alternative (doubled or halved) solution and the log-likelihood difference.
- Spread of purity estimates across weighting schemes and arms.
- Per-segment and per-arm states, with flags for arms inconsistent with the consensus purity.

---

## 3. Version 2: capture-kit-aware depth ratio from coverage profiles

### 3.1 Principle

In hybrid-capture sequencing, captured fragments overlap the bait, so read depth across a target forms a unimodal ("bell-shaped") profile: maximal over the bait and decreasing towards and beyond the target edges over a distance determined by fragment length. Depth at a single base samples one point of this profile; the integral of the profile uses all reads assigned to the target.

Version 2 replaces single-position depth ratios with **per-target ratios of the area under the coverage profile (AUC) between tumor and matched normal**, computed on regions defined for the specific capture kit. BAF from heterozygous sites is retained as in Version 1.

### 3.2 Kit-specific region sets

- Maintain a registry of capture designs, one region set per kit and genome build, for example Agilent SureSelect Human All Exon V6, V7, V8, and other vendors' exome and panel designs. Each entry holds the vendor target BED, the bait (probe) BED where available, and the derived filtered "bell" BED (3.4).
- The kit is either specified by the user or detected automatically by correlating the normal sample's coverage with each registered design.

### 3.3 Signal per target

For each filtered target *t*, with window *W_t* = target ± flank (flank set from the library's insert-size distribution):

$$
\mathrm{AUC}_t^{T} = \sum_{x \in W_t} d_T(x), \qquad
\mathrm{AUC}_t^{N} = \sum_{x \in W_t} d_N(x)
$$

$$
R_t = \frac{\mathrm{AUC}_t^{T} / L_T}{\mathrm{AUC}_t^{N} / L_N}
$$

where *d(x)* is per-base depth (or, equivalently and faster, fragment counts assigned to the window) and *L* is a library-size normalizer (total on-target reads, or median target AUC). *R_t* is then GC-corrected using fragment-level GC of the target window and median-centered.

Alternative estimator: **shape-template fit.** The normal profile is used as a template and the tumor scale factor is estimated by least squares, $s_t = \arg\min_s \lVert d_T - s\,d_N \rVert^2$. Bases with low normal depth (target edges) receive proportionally less weight, and the fit residual identifies targets whose tumor profile shape deviates from the normal (local artefacts, breakpoints within a target).

### 3.4 Filtering the region BED

Targets are retained only if they pass the following filters (all thresholds configurable):

- mappability below threshold, or overlap with segmental duplications / blacklisted regions;
- extreme GC content of the target window;
- targets shorter than a minimum length, or with normal AUC below a minimum count;
- high variance of normalized AUC in a panel of normals sequenced with the same kit (where available);
- adjacent targets closer than approximately one fragment length have overlapping profiles and are merged into one region with one AUC;
- sex chromosomes handled separately (PAR, expected copy number by sex).

The filtered BED is versioned per kit so results are reproducible.

### 3.5 Feeding Sequenza

- Each target is represented by one record in the `.seqz` stream: position = profile maximum (or target midpoint), `depth.ratio = R_t`, `zygosity.normal = hom`; the record contributes to depth ratio only.
- Heterozygous SNP sites (FACETS loci, as in Version 1) are kept as separate records and contribute **BAF only**. Their per-position depth ratio is not used, so that reads within a target are not counted twice.
- Windowing and segmentation weight each target by its normal read count rather than one weight per record.
- Solution selection and BAF post-processing are as in Version 1 (Sections 2.2–2.3).

### 3.6 Differences from Version 1 input

- One depth ratio observation per target instead of one per polymorphic position.
- Depth ratio is defined for targets without polymorphic loci.
- Region definitions are specific to each capture kit and genome build.

---

## 4. Relation to prior art

- Using a fixed set of polymorphic loci including homozygous sites is FACETS' input model [2]; the new element is using it as the input to Sequenza's joint likelihood and segmentation.
- Mean depth per capture target (which equals AUC divided by target length) is used by CNVkit [3] and PureCN [4]. The new elements here are: windows defined by the empirical capture profile (target + fragment-length flank, merged bells), the optional shape-template ratio, kit-specific filtered region sets, and feeding target-level ratios into an allele-specific joint purity/ploidy fit while BAF comes from separate SNP records.
- Preference for a diploid baseline exists in FACETS (dipLogR) [2] and karyotype priors exist in ABSOLUTE [6]. The new elements are the explicit degeneracy test on integer states, the WGD evidence test, and the post-fit BAF refinement with multiple segment-weighting schemes and per-arm consistency.

---

## 5. Implementation outline

```
cnv-caller/
  kits/                      # per-kit, per-build region sets
    agilent_v6/hg38/{targets.bed, baits.bed, bells.filtered.bed}
    agilent_v7/hg38/...
    agilent_v8/hg38/...
  src/
    loci/        snp-pileup wrapper, het/hom classification
    bells/       per-target AUC / template-fit, GC correction, BED filtering
    seqz/        writers for Version 1 and Version 2 .seqz streams
    fit/         patched sequenza.fit: ploidy prior, degeneracy and WGD tests
    post/        BAF-based purity refinement, arm-level analysis, consensus
    report/      solutions, likelihood gap, uncertainty, QC plots
  benchmarks/
```

- Sequenza and FACETS are distributed under GPL licenses; derived code is released under a compatible license.
- The unmodified Sequenza path is retained as a baseline mode, so that each modification can be evaluated separately.

---

## 6. Evaluation plan

Comparison against Sequenza, FACETS, PureCN, CNVkit and ASCAT on the following data and metrics.

**Datasets with orthogonal truth**
- SEQC2 somatic reference (HCC1395 tumor / HCC1395BL normal), sequenced as WGS and WES across centers, including tumor/normal titrations with known purity [7, 8].
- In silico admixtures of tumor and matched normal reads at defined purities (e.g. 10–90%), per capture kit.
- Exomes with matched SNP-array ASCN calls (e.g. TCGA, with ASCAT/ABSOLUTE purity and ploidy from SNP6 arrays [9]).

**Metrics**
- Purity: absolute error vs known titration / array-based estimate.
- Ploidy: error and rate of 2× (diploid/tetraploid) swaps; WGD call accuracy.
- Segment-level total and allele-specific copy number concordance (fraction of genome correct).
- Breakpoint precision/recall.
- Stability: same sample across kits (Agilent V6/V7/V8) and across sequencing depth downsampling.
- Runtime.

**Ablations**
- Version 1: FACETS loci input only; + solution selection (2.2); + BAF post-processing (2.3).
- Version 2: AUC ratio vs template-fit ratio vs single-position depth ratio.
- Each region filter (3.4) enabled vs disabled.

---

## References

1. Favero F, et al. Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data. *Ann Oncol.* 2015;26(1):64–70.
2. Shen R, Seshan VE. FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. *Nucleic Acids Res.* 2016;44(16):e131.
3. Talevich E, et al. CNVkit: genome-wide copy number detection and visualization from targeted DNA sequencing. *PLoS Comput Biol.* 2016;12(4):e1004873.
4. Riester M, et al. PureCN: copy number calling and SNV classification using targeted short read sequencing. *Source Code Biol Med.* 2016;11:13.
5. Van Loo P, et al. Allele-specific copy number analysis of tumors. *PNAS.* 2010;107(39):16910–16915.
6. Carter SL, et al. Absolute quantification of somatic DNA alterations in human cancer. *Nat Biotechnol.* 2012;30(5):413–421.
7. Fang LT, et al. Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. *Nat Biotechnol.* 2021;39:1151–1160.
8. Xiao W, et al. Toward best practice in cancer mutation detection with whole-genome and whole-exome sequencing. *Nat Biotechnol.* 2021;39:1141–1150.
9. Taylor AM, et al. Genomic and functional approaches to understanding cancer aneuploidy. *Cancer Cell.* 2018;33(4):676–689.

---

*This document is a public disclosure of the concepts described above by the author, dated as shown.*
