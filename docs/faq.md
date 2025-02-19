# Frequently asked questions

This section will be expanded as we find new questions and potential pitfalls. 

## Which parameters should I pay attention to?

Generally, GABI runs fine with all-default settings. Parameterization is most typically needed for ONT data.

### ONT data

`--medaka_model`  The basecalling model used; only needed if your basecaller does not encode it in the sequence headers

`--ont_min_q` Minimum quality ONT reads to keep

`--ont_min_length` Minimum ONT read length to keep

The latter two options are more meant to nudge the dataset towards "longer and better". GABI will perform downsampling of the reads anyway (unless deactivated by the user); but it will normally not select for the "best" reads during that process. 

`--skip_homopolish` Skip polishing of homopolymer errors. Since this uses sequence information from related assemblies, you may or may not wish to skip this. 

## Technologies

### Do any ONT reads work with GABI?

Well...yes. But to be perfectly transparent, we are not testing with R9 reads or earlier and will not dedicate time to issues related to their support. With the latest R10 chemistry, results with GABI generally look good. We do recommend you use the super-accurate base calling (SUP) model 5.0.0 or later however, since our tests show this to be a prerequisite for consistently good assemblies. 

## Detecting contaminations

Generally speaking, detecting contaminations from read data is based on either the presence of variable sites when there shouldn't be any (bacteria are haploid, after all) or mixed taxonomic signals from e.g. Kmer analyses. 

### Nanopore data

Nanopore data poses a particular challenge for the detection of contamination from SNP data since the read data is comparatively noisy. While GABI does try to perform such contamination checks on Nanopore data, the results are to be interpreted with a big grain of salt. Essentially, low levels of intra-species contaminations are unlikely to show up in Nanopore data since the (potentially) small number of genetic differences are drowned by the noise. To this end, we run ConfindR with rather stringent settings to prevent the noise from triggering warnings (a contaminating SNP must be supported by at least 5 reads, which would correspond to 10% when sequencing to a recommended depth of 50X; and only reads >= Q20 are used). Unfortunately, this still isn't a guarantee for a totally robust inference, depending on read depth and quality. In fact, as hinted at earlier, this strategy  will obscure cases of true contamination when read coverage or levels of contamination are low and converesely trigger contamination warnings for no reason when coverage is really high.
