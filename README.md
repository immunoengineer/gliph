# gliph

* GROUPING OF LYMPHOCYTE INTERACTIONS BY PARATOPE HOTSPOTS - GLIPH

Questions? Comments? Contact
   jakeg@stanford.edu / jake@hackbio.com
   huang7@stanford.edu

Name
----

  gliph - Grouping of Lymphocyte Interactions by Paratope Hotspots

Synopsis
--------

  gliph-group-discovery.pl [options ] --tcr TCR_TABLE

  gliph-group-scoring.pl --convergence_file TCR_TABLE-convergence-groups.txt \
                         --clone_annotations TCR_TABLE \
                         --hla_file HLA_TABLE \
                         --motif_pval_file TCR_TABLE.minp.ove10.txt

Description
-----------

GLIPH clusters TCRs that are predicted to bind the same MHC-restricted
peptide antigen. When multiple donors have contributed to the clusters,
and HLA genotypes for those donors are available, GLIPH additionally 
can provide predictions of which HLA-allele is presenting the antigen.

Typically the user will pass in a sequence set of hundreds to thousands
of TCR sequences. This dataset will be analyzed for very similar TCRs,
or TCRs that share CDR3 motifs that appear enriched in this set relative
to their expected frequencies in an unselected naive reference TCR set. 

GLIPH returns significant motif lists, significant TCR convergence groups,
and for each group, a collection of scores for that group indicating
enrichment for motif, V-gene, CDR3 length, shared HLA among contributors,
and proliferation count. When HLA data is available, it predicts the
likely HLA that the set of TCRs recognizes. 

Installation
------------

First, unpack gliph:

tar -xzvf gliph-1.0.tgz

You are done. The gliph commands are found in 
             
gliph/bin/gliph-group-discovery.pl
gliph/bin/gliph-group-scoring.pl

Optionally, for added convenience, you could add the path to your gliph/bin
directory to your system $PATH.  

Options
-------

Required Data Inputs
~~~~~~~~~~~~~~~~~~~~

The user must provide a table of TCR sequences. 

  --tcr TCR_TABLE

The format of the table is tab delimited, expecting the following columns in this 
order. Only TCRb is required for the primary component of the algorithm to function, 
but patient identity is required for HLA prediction. 

Example:

CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	Patient Counts
CAADTSSGANVLTF	TRBV30	TRBJ2-6	CALSDEDTGRRALTF	TRAV19		TRAJ5	09/0217 1
CAATGGDRAYEQYF	TRBV2	TRBJ2-7	CAASSGANSKLTF	TRAV13-1	TRAJ56	03/0492 2 
CAATQQGETQYF	TRBV2	TRBJ2-5	CAASYGGSARQLTF	TRAV13-1	TRAJ22	02/0259 1
CACVSNTEAFF	TRBV28	TRBJ1-1	CAGDLNGAGSYQLTF	TRAV25		TRAJ28	PBMC863 1
CAGGKGNSPLHF	TRBV2	TRBJ1-6	CVVLRGGSQGNLIF	TRAV12-1	TRAJ42	02/0207 1
CAGQILAGSDTQYF	TRBV6-4	TRBJ2-3	CATASGNTPLVF	TRAV17		TRAJ29	09/0018 1
CAGRTGVSTDTQYF	TRBV5-1	TRBJ2-3	CAVTPGGGADGLTF	TRAV41		TRAJ45	02/0259 1
CAGYTGRANYGYTF	TRBV2	TRBJ1-2	CVVNGGFGNVLHC	TRAV12-1	TRAJ35	01/0873 3


Optional Data Inputs
~~~~~~~~~~~~~~~~~~~~

The user may additional supply a table of HLA genotyping for each subject.

  --hla HLA_TABLE

The format of the table is tab delimited, with each row beginning with the identity
of a subject, and then two or more following column providing HLA identification. 
The number of total columns (HLA defined genotypes) is flexible. 

Example:

09/0217	DPA1*01:03	DPA1*02:02	DPB1*04:01	DPB1*14:01	DQA1*01:02
09/0125	DPA1*02:02	DPA1*02:02	DPB1*05:01	DPB1*05:01	DQA1*06:01
03/0345	DPA1*02:01	DPA1*02:01	DPB1*17:01	DPB1*01:01	DQA1*01:03
03/0492	DPA1*01:03	DPA1*02:01	DPB1*03:01	DPB1*11:01	DQA1*01:02
02/0259	DPA1*01:03	DPA1*01:03	DPB1*104:01	DPB1*02:01	DQA1*02:01


Optional Arguments
~~~~~~~~~~~~~~~~~~~

  --refdb DB             optional alternative reference database

  --gccutoff=1           global covergence distance cutoff. This is the maximum
			 CDR3 Hamming mutation distance between two clones sharing
			 the same V, same J, and same CDR3 length in order for them
		         to be considered to be likely binding the same antigen. 
			
                         This number will change depending on sample depth, as
			 with more reads, the odds of finding a similar sequence 
			 increases even in a naive repertoire. 

			 This number will also change depending on the species
			 evaluated and even the choice of reference database (memory
			 TCRs will be more likely to have similar TCRs than naive 
			 TCR repertoires). Thus, by default this is calculated at 
			 runtime if not specified.

  --simdepth=1000        simulated resampling depth for non-parametric convergence 
			 significance tests. This defines the number of random
			 repeat samplings into the reference distribution that
			 GLIPH performs when analyzing 
				1) global similarity cutoff
				2) local similarity motif enrichment
				3) V-gene enrichment
				4) CDR3 length enrichment
				5) clonal proliferation enrichment
			 A higher number will take longer to run but will produce 
			 more reproducible results. 

  --lcminp=0.01          local convergence minimum probability score cutoff. 
			 The score reports the probability that a random sample of the
			 same size as the sample set would but into the reference set
			 (i.e. naive repertoire) would generate an enrichment of the
			 given motif at least as high as has been observed in the 
			 sample set. It is set to 0.01 by default.

  --lcminove=10          local convergence minimum observed vs expected fold change.
			 This is a cutoff for the minimum fold enrichment over a
			 reference distribution that a given motif should have in 
			 the sample set in order to be considered for further evaluation.
			 It is set to 10 by default.

  --kmer_mindepth=3      minimum observations of kmer for it to be evaluated. This is
			 the minimum number of times a kmer should be observed in 
			 the sample set in order for it to be considered for further
			 evaluation. The number can be set higher to provide less
			 motif-based clusters with higher confidence. This could be
			 recommended if the sample set is greater than 5000 reads. Lowering
			 the value to 2 will identify more groups but likely at a cost
			 of an increase False Discovery Rate. 

  --global=1             Search for global TCR similarity (Default 1)

  --local=1              Search for local TCR similarity (Default 1)

  --make_depth_fig=0     Perform repeat random samplings at the sample set depth in 
			 order to visualize convergence

  --discontinuous=0      Allow discontinuous motifs (Default 0)

  --positional_motifs=0  Restrict motif clustering to a shared position that is fixed
			 from the N-terminal end of CDR3

  --cdr3len_stratify=0   Stratify by shared cdr3length distribution (Default 0)

  --vgene_stratify=0     Stratify by shared V-gene frequency distribution (Default 0)
 
  --public_tcrs=0        Reward motifs in public TCRs (Default 0)

Controlling Output
~~~~~~~~~~~~~~~~~~

GLIPH produces multiple output files. 

  --output FILE          Place command output into the named file.

  --verbose              Have *tmo* print messages describing
                         

Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Convergence group  - a set of multiple TCRs from one or more individuals that bind
	  	     the same antigen in a similar manner through similar TCR 
		     contacts. GLIPH predicts convergence groups by asking "what
		     is the probability that this cluster of similar TCRs could
		     have appeared without the selection of common antigen?"

Global convergence - A pair of TCRs that share the same length CDR3 and differ
		     by less than a certain number of amino acids in those CDR3s. 
		     Example: in a set of 300 random TCRs, finding two TCRs that
		     only differ by one amino acid in their CDR3 would be highly
		     unlikely. 

Local convergence  - A pair of TCRs that share in their CDR3 regions an amino acid
		     motif that appears enriched in their sample set. Optionally
		     this common motif could be positionally constrained. 
		     Example: in the malaria set, an enriched QRW motif was found
		     in 23 unique TCRs from 12 individuals in a conserved position
		     in their CDR3. 

Reference set	   - A large database of TCR sequences that are not expected to
		     be enriched for the specificities found in the sample set. 
		     Example: by default, GLIPH uses as a reference database 
		     over 200,000 nonredundant naive CD4 and CD8 TCRb sequences 
		     from 12 healthy controls. 

Sample set         - The input collection of TCRs under evaluation that are
		     potentially enriched for a specificity not present in the
		     reference set. 

Function of Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GLIPH first makes the input sample set non-redundant at the CDR3 amino acid level.
The identity of the V-genes, the J-genes, the clonal frequencies, and the associated
HLA information for each TCR is hidden during the CDR3 analysis in order to act
as independent confirmation variables of the resulting clusters. 

GLIPH next counts the number of unique CDR3 sequences in the sample dataset. 

GLIPH next calculates the minimum hamming distance from each CDR3 in the sample 
dataset to the next closest same-length sequence in the sample set. This result is 
the minimum distance distribution.
 
GLIPH next loads the reference dataset. 

GLIPH next repeatedly resamples from the reference database at the nonredundant
sample set depth, calculating the minimum hamming distance distribution in each case.
This is repeated simdepth times (1000 by default). Probability of encountering
next closest members at distant 0, 1, 2... is tabulated across all random samples.  

GLIPH next compares the sample set minimum hamming distance distribution to that of
the reference set, to assign a probability to each hamming distance that sequences
of that distance are likely to appear randomly. 

Example: after analysis, it is found that out of 1000 simulations

   hamming distance in-test-set  in-ref-set    probability
   1                 14          0             0.001
   2                 37          8+/-3         0.02
   3                 194         147+/-35      0.24

GLIPH then selects a maximum hamming distance cutoff for grouping similar TCRs. 

NOTE: For the standard reference set this analysis has been precomputed for a range
of sampling depths, enabling a lookup table rather than simulation as a runtime
accelerator. 
 
GLIPH then analyzes all possible 3mer, 4mer, 5mer and discontinuous 3mer and 4mer
motifs for their frequency in the sample set. This is done excluding the first three
and last three residues in the CDR, as they are not observed to be in contact with
antigen in known crystal structures. This is done for every non-redundant CDR3 in
the reference set. 

Example: for a scan of the following CDR3, the following 3mers would be collected 

  CAS SFGSGGHYE TYF
  xxx|SFGSGGHYE|xxx
      SFG
       FGS
        GSG
         SGG
          GGH
           GHY
            HYE
      SFGS
       FGSG
        GSGG
         SGGH
          GGHY
           GHYE
      S.GS
       F.SG
      (etc)

The frequences off all motifs across the sample set are calculated. 

In order to establish whether these motifs are just naturally abundant or specifically
enriched by antigen, these frequencies are compared to a repeat random sampling of
the reference set at the same depth as the non-redundant database. The frequency
of every motif is collected with each random sampling in order to create a distribution.

If each motif could only be observed in a given sequence once, then the distribution of sampling
motif frequency means u become normally distributed and this result is equivalent
to calculating the frequencies of all motifs in the reference database, and then calculating
one-sided confidence intervals for expected frequencies of any given motif in the reference
database at any given sampling depth:
   
               _
CI(99.9% os) = y + t        (s / sqrt(n))
                 -  0.0005
                                                        _
where n is the sample set non-redundant CDR3 sample size, y is the motif mean frequency, and
s is the SD estimate at that sampling depth for the motif, as
       ________________
      |  _ _       _ 2
s =   |  \   (y1 - y)
      |  /__ 
      | --------------
     \|   (n - 1)

This approximation provides a runtime acceleration.

GLIPH then compares the probability that each motif in the sample set would have appeared
at the observed frequency by random change (calculated as the number of times the motif
was found at at least that depth in a random sampling, divided by the total number of
random samplings, or alternatively though the CI). 

GLIPH next generates a network of all CDR3s as nodes, and all edges as either global or
local interactions. The clusters can be optionally filtered with the following arguements

  --positional_motifs=0  Restrict motif edges to CDR3s where the motif is in a shared 
			 position that is fixed from the N-terminal end of CDR3

  --public_tcrs=0        Allow CDR3s from multiple subjects to increase the size of
			 their clusters. Rewards public TCRs but risks contamination
			 clouding results. 

  --global_vgene=0       Restricts global relationships to TCRs of common V-gene
			 Recommended but invalidates V-gene probability score

GLIPH then scores each individual cluster (candidate convergence group) by evaluating a
set of features that are independet of the CDR3 observations, assigning a probability 
p to each feature, and then combining those probabilities into a single score by conflation. 

For tests i through N, testing Pi(X=C) probability that cluster X is convergent, the combined
conflation score is given as

                      __N
P(X=C) =              ||  P (X=C)
                        i  i
           ----------------------------
            __N             __N
            ||  P (X=C)  +  ||  P (X!=C)
              i  i            i  i

The individual Pi(X=C) tests include 

1) global similarity probability
2) local motif probability
3) network size
4) enrichment of V-gene within cluster
5) enrichment of CDR3 length (spectratype) within cluster
6) enrichment of clonal expansion within cluster
7) enrichment of common HLA among donor TCR contributors in cluster 

Individual score components are calculated as follows:

== 3) Calculating network size p ==

For each discrete cluster, the probability p of a given cluster topology can be 
obtained from the number of members of the cluster by comparison to a lookup table 
calculated from repeat random sampling and GLIPH clustering of naive TCR sequences 
at a range of different sampling depths n from 25 to 5000, each performed 1000 times 
each. 

Example: at sampling depth n=500, clusters of size 5 have a probabilty p=0.002 of occuring
in naive TCR sample sets.

== 4) enrichment of V-gene within cluster ==

GLIPH hides V-gene usage for all clones prior to CDR3 analysis, and does not explore the
V-gene template endcoded amino acids in the local (motif) search. 

To evaluate whether there is an enrichment of common V-segment within the cluster over
the degree that would be expected from an unbiased sampling of TCRs of cluster size n
from the total dataset, we perform repeat random sampling at sample size n from the
total dataset, each time obtaining the V-genes for each clone and calculating a 
 V-gene Simpson's Diversity Index D for each sample, where D is interpreted as the probability
that any two members within the cluster would share a V-gene and is calculated as

      __V
     \   v(v-1)
     /__1
D = -----------
      __V
     \   n(n-1)
     /__1

where V is the number of all V-genes, v is the total counts of a given v-gene, and n 
is the sampling size (cluster size). 
 
The probability of the observed D for candidate convergence group is obtained as the 
one-tailed probability of observing a score at least that high in the D score distribution
from random sampled clusters of same size n.

== 4) enrichment of CDR3 within cluster ==

For local (motif) scoring, no stratification based on CDR3 length is performed. 

To evaluate whether there is an enrichment of common CDR3 length within the cluster over
the degree that would be expected from an unbiased sampling of TCRs of cluster size n
from the total dataset, we perform repeat random sampling at sample size n from the
total dataset, each time obtaining the CDR3 lengths for each clone and calculating a
 V-gene Simpson's Diversity Index D for each sample, where D is interpreted as the probability
that any two members within the cluster would share a CDR3length and is calculated as

      __C
     \   c(c-1)
     /__1
D = -----------
      __C
     \   n(n-1)
     /__1

where C is the number of all CDR3 lengths, c is the total counts of a given CDR3 length, 
and n is the sampling size (cluster size).

The probability of the observed D for candidate convergence group is obtained as the
one-tailed probability of observing a score at least that high in the D score distribution
from random sampled clusters of same size n.


== 5) enrichment of clonal expansion within cluster ==

The first step in GLIPH is to make all CDR3s non-redundant, hiding their clone frequences.
To evaluate whether there is an enrichment of expanded clones within a given convergence
group, the number of counts for each clone are recovered, and a convergence group expansion
coefficient e is calculated as

e = total clones / total unique clones

To evaluate whether there is an enrichment of expanded clones within the cluster over
the degree that would be expected from an unbiased sampling of TCRs of cluster size n
from the total dataset, we perform repeat random sampling at sample size n from the
total dataset, each time obtaining the clone counts for each clone and calculating e
for each random sample to establish a distribution. 

The probability of the observed e for candidate convergence group is obtained as the
one-tailed probability of observing a score at least that high in the e score distribution
from random sampled clusters of same size n.

== 6) enrichment of common HLA among donor TCR contributors in cluster ==

TCRs that recognize a common antigen should be constrained by the same HLA. Thus, their
contributing donors should contain that HLA allele in their genotype. 

When HLA genotype information is available, GLIPH uses combinatorial resampling without
replacement to estimate the probability that the collection of TCRs in the convergent
group recognizes any given HLA. For each HLA allele that appears at least twice in the
candidate convergence group, GLIPH obtains

a the number of subjects in the convergence group that harbor that allele
A the number of subjects in the study that harbor that allele 
n the number of subjects in the convergence group
N the number of subjects in the study

GLIPH uses this information to calculate the probability that a given HLA allele is
present by chance.
 
Examples
--------

To run GLIPH on a TCR_TABLE file mytcrtable.txt, run:
gliph-group-discovery.pl --tcr mytcrtable.txt

Tun run GLIPH on a list of CDR3s mycdr3list.txt, run:
gliph-group-discovery.pl --tcr mycdr3list.txt

To run GLIPH with an alternative mouse reference DB mouseDB.fa, run:
gliph-group-discovery.pl --tcr mytcrtable --refdb=mouseDB.fa
                                       
To run GLIPH slower with a more thorough simdepth, run
gliph-group-discovery.pl --tcr mytcrtable.txt --simdepth=10000

To run GLIPH slower with a more thorough simdepth and an altered lcminp, run
gliph-group-discovery.pl --tcr mytcrtable.txt --simdepth=10000 --lcminp=0.001

To score GLIPH clusters, run 
  gliph-group-scoring.pl --convergence_file TCR_TABLE-convergence-groups.txt \
                         --clone_annotations TCR_TABLE \
                         --hla_file HLA_TABLE \
                         --motif_pval_file TCR_TABLE.minp.ove10.txt
