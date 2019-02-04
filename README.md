# gene_sequence_clusteringV2
Version 2 of the gene clustering. Includes more accurate and efficient cross validation code
This has a newer version of the genetic sequences clustering. 
A little more about this: I am trying to detect novel signals in genetic data (currently in promoter regions)
The well known TATA box only accounts for 20% of the promoter regions in genes, what about the remaining 80%?
Current methods primarily focus in scanning for existing architectures, or if find new ones that are very overrepresented
My clustering method aims to cluster these sequences to be able to better understand the diversity in these genes.

Index:

crossValV2: The cross validation code that approximates the correct number of clusters for the sequences
orderByCluster: Clusters the given data into the # of clusters given by crossVal2
realRefinedGen: The code I use to generate my synthetic data (with varying alpha to vary noise and cluster distributions in data)

Some datasets I am working with
newFlyData: genetic sequences from the drusophila melanogaster genome that I am working with.(Data obtained from: 	Ni T.,  Corcoran D.L.,  Rach E.A.,  Song S.,  Spana E.P.,  Gao Y.,  Ohler U.,  Zhu J.. A paired-end sequencing strategy to map the complex landscape of transcription initiation, Nat. Methods , 2010, vol. 7 (pg. 521-527)
tbDATA
syntheticDATA: An example of the kind of synthetic data I generate (with the widely used Dirichlet Distribution) in order to 
test my algorithm for accuracy (since I know the truth).
