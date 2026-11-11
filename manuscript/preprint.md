# Most candidate integration sites in host sequencing data are library artefacts, and a Poisson test removes them

Iraaj Gangavaram

School of Biological Sciences, University of East Anglia, Norwich, United Kingdom

ORCID: 0009-0002-4615-2846

Correspondence: [email]

---

## Abstract

Host sequencing data routinely contains reads that a taxonomic classifier assigns to bacteria or viruses. Where several such reads share a position in the host genome, that position is often reported as a candidate integration site. Most of these positions are artefacts. During library preparation a host fragment and a foreign fragment are joined into one molecule, and the resulting read pair is indistinguishable from a genuine junction. On simulated data built from real human and bacterial genomes, taking every position with at least three supporting pairs returned 24,861 candidate sites of which 390 were real.

hgt-site tests each candidate position against the number of chance chimeras expected there, estimated from the genome-wide rate for that taxon, and reports the sites that exceed it. On the same datasets it returned 364 sites, all of which were real, while recovering all 360 insertions. The method works for bacterial donors, viral donors and supplied vector sequences, places junctions to the base pair where reads cross them, and runs in seconds on an aligned sample. It is distributed through PyPI and Bioconda.

All results presented here are from simulated reads built from real genomes. The tool has not been evaluated on a sequencing dataset with an independently confirmed integration site, and results should be treated accordingly.

## Introduction

A host sequencing library almost always contains foreign DNA. Some of it is genuine: a prophage, an integrated virus, an inserted transgene. Much more of it is not. Reagents and kits carry bacterial DNA, low-biomass samples are substantially reagent in origin, and classifiers misassign host reads in repetitive regions.

The hardest of these to deal with is the chimeric molecule. Library preparation occasionally joins a host fragment to a foreign fragment, by ligation or during amplification. Sequenced, that molecule gives a read pair with one mate on the host and one mate on something foreign, consistent with each other and with a junction. Nothing about an individual read pair distinguishes it from a real insertion.

The usual response is to require several supporting reads at a position. This does not work, because chimeras cluster too. They accumulate where coverage is high, and in a whole-genome library there are enough of them that thousands of positions pass any threshold low enough to retain real sites. Raising the threshold loses true sites faster than false ones.

What distinguishes a real insertion is not how many foreign reads sit at a position, but how many more sit there than chance alone would put there. That quantity can be estimated from the data itself, and that is what hgt-site does.

Existing tools approach the problem differently. LGTSeek identifies lateral gene transfer in host data through read classification and subsequent filtering, and has had no commit since 2016. VIRUSBreakend assembles the junction with GRIDSS and is the established choice for viral integration in cancer genomes. ViFi uses hidden Markov models over viral families. nf-core/hgtseq aligns reads, classifies the unmapped ones and reports the positions of their mates. None of these models the rate at which library preparation produces chimeras, which is the dominant source of false positives in this setting.

## Implementation

hgt-site takes an alignment of paired-end reads to the host genome, together with per-read taxonomic calls for the reads that did not align. Output from Kraken2, KrakenUniq, Kaiju, Centrifuge, DIAMOND and BLAST is read directly, as is a two-column file of read name and taxid. Where the foreign sequence is known in advance, as with a vector or a particular virus, its FASTA can be supplied instead and the unmapped reads are matched against it by exact 31-mer matching in both orientations, which removes the need for a classifier database.

**Collecting anchored pairs.** The alignment is read once. Pairs where both mates map are discarded immediately, which in a host library is well over 99% of records and is most of the reason the tool is fast. Pairs where one mate maps with mapping quality at least 20 and the other does not are kept as anchors. Pairs with both mates unmapped are counted as background, since they indicate how much of each taxon is present without indicating where.

**Clustering.** Anchors are grouped by the end facing the junction, so that reads from both sides of one insertion fall in the same cluster. Cluster width is capped to stop scattered chimeric reads chaining into one long false cluster. The limits are set from fragment sizes measured in the library itself, which matters for long-fragment libraries where fixed limits split a single insertion into two reported sites.

**Testing against chance.** For each cluster, the number of pairs from its dominant taxon is compared with a Poisson expectation: the number of pairs of that taxon that would fall in a window of that size if the taxon's anchored pairs were distributed at random across the host genome. P-values are adjusted by the Benjamini-Hochberg procedure over every genome window that could have produced a cluster, for every taxon present, rather than over the observed clusters alone, which would make the correction too lenient.

The assumption that chimeras are distributed evenly across the genome fails wherever host coverage is uneven, since a chimera's host half is drawn from the host DNA present in the library. A region at six copies contributes six times the chimeras. An optional mode counts host coverage in 10 kb bins during the same pass and scales the expectation by the local rate over the genome mean. The adjustment only ever raises the expectation and never lowers it, because the local estimate carries a few per cent of sampling noise and allowing it to lower the bar passes marginal candidates.

**Locating the junction.** Read pairs place a junction only to within a fragment length. Reads that cross the junction align to the host up to that point and are soft-clipped beyond it, so several clips ending at the same base give the position to the base pair. Both sides are reported, since they differ where the insertion duplicated or deleted host sequence.

**Reporting.** Every candidate position appears in the output with either PASS or the reasons it failed: insufficient support, support from several taxa, low mapping quality, low sequence complexity, or no enrichment over chance. Keeping the failures visible is deliberate, since a rejected position in a region of interest is often more informative than the list of accepted ones. Output is written as a table, optionally as BED, as VCF 4.2 with the positional uncertainty recorded in CIPOS and CIEND, and as a MultiQC section.

## Results

All benchmarks use simulated reads built from real genome sequence: 10 Mb of human chromosome 21 as the host, Escherichia coli as the donor, and four other bacterial species as contaminants. Insertions are placed at known positions, chance chimeras are generated at a controlled rate, and the whole dataset is rebuilt from a seed.

**Separating insertions from chimeras.** Across six combinations of coverage and chimera rate, with three seeds each, hgt-site reported 364 sites of which all 364 were real, and recovered all 360 insertions. On the same data, taking every position with at least three supporting pairs returned 24,861 sites of which 390 were real. LGTSeek's classification step behaved like the support-only threshold, with precision falling to 0.007 at the highest chimera rate. Repeating the comparison with an independently published read simulator gave the same conclusion.

**Viral donors.** Substituting a viral genome for the bacterial donor gave a precision of 0.989 over eighteen runs, with 357 of 360 insertions recovered. The virus used does not integrate and stands in for size and base composition only.

**Samples containing free donor.** A transduced culture or a virus-positive tissue sample contains free donor alongside the integrated copy, and that free donor forms chimeras carrying the same taxon as the genuine insertion, so classification offers no help in separating them. With six chance chimeras in ten drawn from the donor, precision remained 1.000 and recall fell from 1.000 to 0.933 as the chimera load rose. The sites that were missed had few supporting pairs and failed the enrichment test, which is the test behaving conservatively rather than incorrectly: abundant free donor genuinely raises the rate at which that taxon appears anywhere.

**Uneven host coverage.** With a 500 kb region present at six copies, the genome-wide background reported 69 sites across three seeds of which 10 were false, every one inside the amplified region. The local-coverage background reported 59 sites with none false and lost no insertion. Where coverage was even, the two produced identical results.

**Vectors and transgenes.** Constructs of 1 kb, 3 kb and 8 kb, called from the construct sequence with no classifier, gave 100 sites across five runs, all real, with all 100 insertions found and every junction exact. Because only donor-matching reads are classified in this mode, chimeras built from the contaminants are invisible and cannot produce false sites, so that result overstates the difficulty. Repeating it with nine chance chimeras in ten carrying the construct itself gave 6,187 candidate positions, of which 302 had enough support to be tested, all 302 failed, and the 20 real insertions were the only sites reported.

**Junction accuracy.** Over 8 runs and 157 insertions, 71% of junctions were exact and 98% were within 2 bp, against 41% and 57% from read pairs alone.

**Speed.** Calling takes about 3 seconds on a sample of 2.1 million reads using 2 cores.

**Where it does not work.** On simulated RNA-seq, where coverage follows expression across four orders of magnitude, precision reached 0.773 at best against 1.000 on DNA. With reads sampled from spliced transcripts and aligned with a DNA aligner, 72 of 133 false sites fell on an exon boundary, because a read crossing a junction is soft-clipped at exactly the position the tool inspects. hgt-site should not be used on RNA-seq.

## Limitations

The central limitation is that every result above comes from simulated reads. Simulated data is cleaner than real data in ways that matter here: the classifier database holds six genomes rather than tens of thousands, chimeras are placed at random rather than in the hotspots real protocols produce, insertions are clean with no rearrangement at the junction, and the host is 10 Mb of one chromosome rather than a genome with its full complement of repeats and segmental duplications. A tool that performs well on this benchmark has demonstrated that the statistics work. It has not demonstrated that it works on a patient sample.

The method assumes an unenriched library. Amplicon protocols that deliberately enrich around junctions, as used in integration site analysis for gene therapy, violate that assumption outright and should be analysed with tools built for them.

Long reads are not supported. A single long read spanning a junction is stronger evidence than any paired-read method can produce.

The composition of the foreign reads at each site is reported but is deliberately not used as a filter, since measurement showed it does not separate chance chimeras from real insertions: a chimera contains genuinely foreign sequence and therefore looks exactly as foreign.

A site reported by hgt-site should be treated as a lead to confirm by PCR across the junction or by inspection in a genome browser, not as a finding.

## Availability

hgt-site is written in Python, requires Python 3.10 or later and has one runtime dependency. It is available from PyPI and from Bioconda, with a container built automatically from the Bioconda recipe. Source, documentation, benchmark scripts and raw benchmark results are at https://github.com/iraajgangavaram/hgt-site. Releases are archived on Zenodo under the concept DOI 10.5281/zenodo.23842712. The software is released under the MIT licence.

Reports of runs on real sequencing data are actively wanted, including runs that found nothing, and an issue template is provided for them. Real data is the one thing the project cannot supply for itself.

## References

To be completed: LGTSeek (Riley et al.), VIRUSBreakend (Cameron et al. 2021, Bioinformatics), ViFi (Nguyen et al.), nf-core/hgtseq (Carpanzano et al.), Kraken2 (Wood et al. 2019), Benjamini and Hochberg (1995), and the independent read simulator used for the replication.
