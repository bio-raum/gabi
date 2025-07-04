# Installation

If you are new to our pipeline ecosystem, we recommend you first check out our general setup guide [here](https://github.com/bio-raum/nf-configs/blob/main/doc/installation.md). That said, the instructions below will probably be sufficient for most users. 

## Installing nextflow

Nextflow is a highly portable pipeline engine. Please see the official [installation guide](https://www.nextflow.io/docs/latest/getstarted.html#installation) to learn how to set it up.

This pipeline expects Nextflow version 24.04.4, available [here](https://github.com/nextflow-io/nextflow/releases/tag/v24.04.4). Other, more recent versions, will probably work also. 

## Software provisioning

This pipeline is set up to work with a range of software provisioning technologies - no need to manually install packages. 

You can choose one of the following options:

[Docker](https://docs.docker.com/engine/install/)

[Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/)

[Podman](https://podman.io/docs/installation)

[Apptainer](https://apptainer.org/)

There is also the option to use Conda, but we **strongly discourage** this because Conda environments are not guaranteed to be reproducible. 

[Conda](https://github.com/conda-forge/miniforge)

The pipeline comes with simple pre-set profiles for all of these as described [here](usage.md); if you plan to use this pipeline regularly, consider adding your own custom profile to our [central repository](https://github.com/bio-raum/nf-configs) to better leverage your available resources.

Also note that Nextflow supports additional software provisioning [frameworks](https://www.nextflow.io/docs/latest/container.html). These may very well work also, but have not been tested by us and would need to be configured as part of your [site-specific](#site-specific-config-file) config file. 

## Installing the references

This pipeline requires locally stored references, matched to the pipeline version you plan on using (-r). To build these, do:

```bash
nextflow run bio-raum/gabi -r main -profile singularity \\
--build_references \\
--run_name build_refs \\
--reference_base /path/to/references
```

where `/path/to/references` could be something like `/data/pipelines/references` or whatever is most appropriate on your system. On a distributed compute environment, this directory needs to live on a shared file system. If you already use a site-specific [config](https://github.com/bio-raum/nf-configs) file, the `--reference_base` option does not need to be set. 

If you do not have singularity on your system, you can also specify docker, podman or conda for software provisioning - see the [usage information](usage.md).

Please note that the build process will create a pipeline-specific subfolder (`gabi`) that must not be given as part of the `--reference_base` argument. GABI is part of a collection of pipelines that use a shared reference directory and it will choose/create the appropriate subfolder automatically. 

Finally, depending on your internet connection, the installation process can take a little while - some of the reference databases are "fairly" large (8-10GB). However, once installed you are all set and ready to go. 

## Site-specific config file

If you run on anything other than a local system, this pipeline requires a site-specific configuration file to be able to talk to your cluster or compute infrastructure. Nextflow supports a wide range of such infrastructures, including Slurm, LSF and SGE - but also Kubernetes and AWS. For more information, see [here](https://www.nextflow.io/docs/latest/executor.html).

Site-specific config-files for our pipeline ecosystem are stored centrally on [github](https://github.com/bio-raum/nf-configs). Please talk to us if you want to add your system. 

### Custom config

If you absolutely do not want to add your system to this repository, you can manually pass a compatible configuration to nextflow using the `-c`  command line option:

```bash
nextflow -c my.config run bio-raum/gabi -profile myprofile -r 1.0.1 --input samples.csv --run_name my_run_name --reference_base /path/to/references
```

A basic example using Singularity may look as follows:

```GROOVY
process {
  resourceLimits = [ cpus: 16, memory: 64.GB, time: 72.h ]
}

singularity {
  enabled = true
  cacheDir = "/path/to/singularity_cache"
}
``` 
This would be for a single computer, with 16 cores and 64GB Ram, using Singularity. Containers are cached to the specified location to be re-used on subsequent pipeline runs.  

Or with the Conda/Mamba package manager:

```GROOVY
process {
  resourceLimits = [ cpus: 16, memory: 64.GB, time: 72.h ]
}

conda {
  enabled = true
  useMamba = true
  cacheDir = "/path/to/conda_cache"
}
```

