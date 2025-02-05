# Frequently asked questions

This section will be expanded as we find new questions and potential pitfalls. 

## Nanopore data

### Detecting contaminations

Generally speaking, detecting contaminations from read data is based on either the presence of variable sites when there shouldn't be any (bacteria are haploid, after all) or mixed taxonomic signals from e.g. Kmer analyses. 

Nanopore data poses a particular challenge here since the read data is comparatively noisy. While GABI does try to perform contamination checks on Nanopore data, the results are to be interpreted with great care. 
Specifically, low levels of intra-species contaminations are unlikely to show up in Nanopore data since the (potentially) small number of genetic differences are drowned by the noise. To this end, we run ConfindR 
with rather stringent settings to prevent the noise from triggering warnings (at least 10 reads must support a contaminating SNP, and only reads >= Q15 are used). Unfortunately, this still isn't a guarantee for a 
totally robust inference, depending on read depth and quality. Note that ConfindR is run *after* read filtering, so if your data is of particularily high quality and depth, it may make sense to specifically filter 
for those reads using `--ont_min_q` and `ont_min_length` to limit the "junk" that makes it into CondindR. 