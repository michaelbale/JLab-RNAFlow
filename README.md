# JLab-ATACFlow #


## Assay for Transposase-Accessible Chromatin using Sequencing (ATACseq) ##


ATACseq is a sequencing technique taking advantage of hyperactive Tn5 transposase to "tagment" so-called Transposase-Accessible chromatin by inserting sequencing adapters into open regions of chromatin. Tagmented DNA can be amplified and sequenced for downstream analysis of differential regions of open chromatin as well as transcription factor binding sites.

## This pipeline

This pipeline was developed by Michael J. Bale, MSc at Weill Cornell Medical College for principal use in processing ATACseq data generated by the lab of Steven Z. Josefowicz, PhD. Although some attempts at generalizing for wider use was made, some extra hacking might be required to get this pipeline to run on your specific machine and/or server. Use at your own risk; but I'm happy to help if I can!

=================================================


      A T A C s e q   P I P E L I N E v0.5
      
=================================================

Author: Michael J. Bale (mib4004@med.cornell.edu)

The typical command for running the pipeline is as follows:

`nextflow run michaelbale/JLab-ATACFlow --input 'data/*_R{1,2}.fastq.gz' --genome mouse [--peaks [--min-replicates n]] -profile singularity`

Mandatory arguments:
```
--input						  Path to input data (must be surrounded with quotes)

--genome					  Name of Genomes reference (current supported: mouse10 -- mm10, mouse39 -- mm39, human -- hg38)
							  For mouse, unless beginning new project, likely better to keep to mm10 (03/04/2022).

-profile                      Name of package manager system (available: docker, singularity, conda);
                              for WCM default -- singularity is recommended, but conda works while docker 
                              does not. For minimal - use conda.

--peaks			  			  Specifies whether or not to call peaks using HMMRATAC (default: will not call peaks)
							  N.B. do not do if running locally/without an HPC -- won't work well. Can still stop
							  at BAM file and call using <insert favorite peak-caller here>

--min-replicates			  Requires --peaks; minimum number of replicates required for overlapping of 
							  individual peak calls to consensus peak calls using ChIP-r.


Options:

--executorConfig              Path to config file that contains specifics for execution. 
                              Will default to WCM SCU-specific parameters. N.B. for
                              single-threaded running use --executorConfig conf/minimal.config

--singleSample                Specifies that the input is a single sample and will not generate a PCA graph

--PCATitle                    Title to be included in PCA graph; must be surrounded with quotes

--catLanes                    Tells CnRFlow to take input files and concatenate lanes into single fastq files

--name                        Project Name; cannot have whitespace characters

--addBEDFilesProfile          Path to csv file with info on additional BED files for generating
                              Sunset-style profile plots; csv format: rName,BEDPath

--addBEDFilesRefPoint         Path to csv file with info on additional BED files for generating
                              Torndao-style region profile plots; csv format: pName,BEDPath,PlusMinus,pLabel

--workDir                     Name of folder to output concatenated fastq files to (not used unless --catLanes)

--outdir                      Name of folder to output all results to (Default: results)

--genomeAssets                Home directory of where genome-specific files are. 
                              Defaults to /athena/josefowiczlab/scratch/szj2001/JLab-Flow_Genomes
```

## Citations


## Requirements:

nextflow v20 or higher

java 8+

Some combination of singularity, docker, or conda. If using -profile debug, all programs specified in the environment.yml file will have to be installed and accessible to your path.


# JLab-RNAFlow