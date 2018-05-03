# reconstruct-accuracy
Assessment of 3D genome reconstruction accuracy.

1) Using multiplex FISH imaging data:

The "gold standard" derives from multiplex FISH imaging, in particular as described by Wang et al., Science  05 Aug 2016: Vol. 353, Issue 6299, pp. 598-602, DOI: 10.1126/science.aaf8084.
Reconstructions of individual chromosome 3D coordinates, as based on HiC or other conformational assays, are referenced to the consensus (averaged over replicates) multiplex FISH structure, with the distance (sum-of-squared deviations after Procrustes alignment) of the multiplex FISH replicates to the consensus providing a referent distribution.  Procrustes alignment is effected using the procrustes function from the R package vegan.
Bias resulting from data re-use (replicates used in obtaining consensus structure AND for referent distribution therefrom) is modest in view of large number (>100) of replicates.

2) Using genome architecture mapping (GAM) data:

The GAM assay (Beagrie et al., Nature Mar 23 2017: 543(7646):519-524. doi: 10.1038/nature21411.) provides a proximity-ligation free methodology of identifying chromatin interactions through use of numerous thin cryo-sections or nuclear profiles.  The planarity of these profiles can be used to evaluate 3D reconstructions based on GAM identified interactions via principal component analyses that contrast planar explained variances with within-profile permutation based measures thereof.
