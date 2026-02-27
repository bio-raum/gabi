# Common issues

## Assembly and qc

### The pipeline did not close my assembly, even though I have used both Nanopore and Illumina reads. 

:   Many reasons can contribute to incomplete assemblies - from the starting material being of poor quality, insufficient sequencing depth, biases in your read data (i.e. loss of certain genomic regions during DNA extraction/preparation) or accidental mix-ups in sample assignment of a subset of your reads. 

    You can check results from ConfindR to see if all your reads are from the same strain and were not accidentally mixed up or in fact contaminated. 
    
    That said, if you did everything right, it could be that the assembly algorithm we employ in this pipeline simply wasn't up to the task. Please let us know if you suspect that to be the case! We try to use state-of-the-art methods, but are always happy to learn new things. 

### Quast reports many differences and an incomplete assembly

:   GABI tries to find the best matching reference genome from RefSeq Bacteria for each sample against which Quast then benchmarks the respective assembly. 

    That said, it is not guaranteed that the best match in RefSeq Bacteria is actually very closely related to your specimen, nor that the assembly is in fact "reference grade". Especially when using long reads for assembly, you may find that Quast reports many mismatches and errors. Our best guess is that the given RefSeq genome is more fragmented than your long-read or hybrid de-novo assembly, which then leads to such incongruencies. So most likely nothing wrong with your assembly but GABI simply being unable to match a reference grade genome to your sample. 

## Pacbio

### Can I use Pacbio subreads (CLR) with GABI?

:   Yes, but the HiFi format has been the defacto standard for Pacbio sequencing for a few years now. If you still have subread data, consider transforming it to CCS/HiFi using available [tools](https://ccs.how/). Alternatively, you may wish to provide complementary short-read data to help GABI improve your assembly. Subread-only assemblies will likely contain numerous errors. 

## Crashes

### I am behind a proxy and the pipeline fails

:   If you are running the pipeline behind a proxy in combination with a container manager (Apptainer, Docker, etc), you may notice pipeline failures with an error message concerning hosts being unreachable. This is because the tools running inside the container do not know about your local proxy settings.

    First, make sure that your proxy settings are correctly configured and stored in the default environment variable:

    ```bash
    echo $HTTPS_PROXY
    ```

    This should return your proxy information. If not, you can set this variable yourself from the command line - else ask to your local IT admin. 

    Next, you will want to configure nextflow to forward these settings to the container environments using a local or site-specific config file. This is done by adding the `envWhitelist` argument:

    ```Nextflow

    apptainer {
        enabled = true
        cacheDir = "$HOME/nextflow_envs_cache"
        envWhitelist = "HTTP_PROXY,HTTPS_PROXY"
    }
    ```

### The pipeline fails because a Conda/Mamba environment could not be solved

:   Conda, in our eyes, is not well suited for production purposes - for excactly this reason. Basically, each time you update your conda installation, there is a risk of certain packages no longer working. 

    And in some cases - and for reasons not entirely transparent - certain systems simply will have problems solving certain environments. Long story short, please consider using a container engine like [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/) or [Apptainer](https://apptainer.org/). 

### The pipeline crashes with an out-of-memory error in one of the processes. 

:   This could simply be an issue with your executing machine not having enough RAM to run some of the tools we put into this pipeline. 

    The exact amount of RAM needed is difficult to predict and can depend on factors like read length and/or sequencing depth - but we suspect that at least 32GB RAM should be available to avoid RAM-related issues. 

    It is also possible that you have not set a memory limit for your compute environment via a site-specifig [config file](https://github.com/bio-raum/nf-configs/) - in which case GABI will use the built-in default (64GB Ram); this perhaps exceed the limits of your system.  

### My ONT assembly crashes with an obscure error

:   Please check if the option `--onthq` is set to `true`. It's possible that this setting is not appropriate for your data, which can lead the assembler to exit on an empty Fasta file halfway through the assembly process; you can disable this option by setting `--onthq false` and resume the pipeline (`-resume`).

### The pipeline immediately fails with a "no such file" error

:   Most likely you saw something like this:

    ```bash
    ERROR ~ No such file or directory: 
    ```

    This is most likely happening because you passed the `reference_base` option from a custom config file via the "-c" argument. There is currently a [known bug](https://github.com/nextflow-io/nextflow/issues/2662) in Nextflow which prevents the correct passing of parameters from a custom config file to the workflow. Please use the command line argument `--reference_base` instead or consider contributing a site-specific [config file](https://github.com/bio-raum/nf-configs). 

## Performance

### Why is the pipeline so slow?

:   We assume you mean the overall start-up time - the performance of the individual processes is dictated by the capabilities of your hardware and the complexity/depth of your data. 

    Otherwise, if you run this pipeline without a site-specific config file, the pipeline will not know where to cache the various containers or conda environments. In such cases, it will install/download these dependencies into the respective work directory of your pipeline run, every time you run the pipeline. And yes, that can be slow - especially when using Conda. Consider adding your own config file to make use of the caching functionality.

### Sourmash `search` is very slow

:   We use sourmash to identify the best matching reference genome for each assembly. This database is currently over 10GB in size and highly contigious assemblies can produce very long run times (30mins+). 

    If you do not care about the best reference genome, but are happy to just find a closely related one so GABI knows which species this is, use the `--fast_ref` option. 


