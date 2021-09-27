Sign prediction based on methods from: https://www.nature.com/nmeth/journal/v11/n1/full/nmeth.2733.html
Uses only z-scores
    is positive regulation if it's above +1.5
    is negative regulation below -1.5
    is 0 in-between +1.5 and -1.5
In a given RNAi screen, if both interacting proteins have nonzero values, then the relationship is classified as either a
positive correlation (both  +1  or  both  −1)  or  a  negative  correlation  (one  is  +1  and another is −1).
For each interacting pair, we computed the total number  of  positive  and  negative  correlations.
Creates gene-phenotype matrix based on data from: http://www.genomernai.org/DownloadAllExperimentsForm
    needed header data: PMID, organism, score type