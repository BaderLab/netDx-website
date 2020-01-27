# Clique filtering

Clique filtering is a technique used to remove "random-like" networks when working with binary similarity and very sparse data. Data types include copy number variations, which are few per patient.

*Motivation section TBA*

## Output format of `cliqueFilterNets.R`
A data.frame with per-network stats on clique-filtering:

* NETWORK: network name
* orig_pp: Num (+,+) interactions.
* orig_rest: Num (+,-) and (-,-) interactions
* ENR: Enrichment or bias of (+,+) interactions relative to other interactions. Specifically defined as `(orig_pp-orig_rest)/(orig_pp+orig_rest)`. Ranges between -1 (all non-(+,+)) to +1 (all (+,+)).
* TOTAL_INT: orig_pp+orig_rest (log-10 transformed)
* numPerm: num permutations done in clique filtering
* shuf_mu: mean ENR of permuted nets (i.e. mean null ENR)
* shuf_sigma: standard deviation of ENR of permuted nets (i.e. s.d. of null ENR)
* Z: Z-score of ENR in real network, relative to null distribution
* pctl: Percentile of ENR in real network, relative to null distribution. Also the p-value
* Q: Benjamini-Hochberg corrected pvalue.
