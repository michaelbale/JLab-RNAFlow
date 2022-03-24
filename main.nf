//Genome specific
params.genome = ''
params.genomes = []
params.kallIndex = params.genome ? params.genomes[ params.genome ].kallIndex ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.genomeInfo = params.genome ? params.genomes[ params.genome ].genomeInfo ?: false : false
params.help = false
params.citations = false



version = 0.5


def helpMessage() {
	log.info """
		=================================================
            R N A  s e q  P I P E L I N E v${version}
        =================================================
		Author: Michael J. Bale (mib4004@med.cornell.edu)
    Usage:
	The typical command for running this pipeline is as folls:
	nextflow run michaelbale/JLab-RNAFlow --input 'project/*{R1,R2}*' --name myProject --genome mouse -profile singularity
--input						  Path to input data (must be surrounded with quotes)

--genome					  Name of Genomes reference (current supported: mouse10 -- mm10, mouse39 -- mm39, human -- hg38)
							  For mouse, unless beginning new project, likely better to keep to mm10 (03/04/2022).

-profile                      Name of package manager system (available: docker, singularity, conda);
                              for WCM default -- singularity is recommended, but conda works while docker 
                              does not. For minimal - use conda.

Options:

--executorConfig              Path to config file that contains specifics for execution. 
                              Will default to WCM SCU-specific parameters. N.B. for
                              single-threaded running use --executorConfig conf/minimal.config

--catLanes                    Tells CnRFlow to take input files and concatenate lanes into single fastq files

--boot						  Number of inferential replicates for kallisto to run; default: 100

--name                        Project Name; cannot have whitespace characters

--workDir                     Name of folder to output concatenated fastq files to (not used unless --catLanes)

--outdir                      Name of folder to output all results to (Default: results)

--genomeAssets                Home directory of where genome-specific files are. 
                              Defaults to /athena/josefowiczlab/scratch/szj2001/JLab-Flow_Genomes
	""".stripIndent()
}


def citationMessage() {
    log.info """
	Please cite the following tools if publishing a paper utilizing this pipeline:
	
	""".stripIndent()
}





/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (params.citations) {
	citationMessage()
	exit 0
}




log.info """\
		=================================================
            R N A S E Q    P I P E L I N E v${version}
        =================================================
		Author: Michael J. Bale (mib4004@med.cornell.edu)
		
		Project ID: ${params.name}
        Genome: ${params.genome}
        Reads: ${params.input}
        Publish Directory: ${params.outdir}
        """
         .stripIndent()
		 
		 
if(params.catLanes) {

    getSampleID = {
	    (it =~ /(.+)_S\d+_L\d{3}/)[0][1]
	}
    
	Channel
	  .fromFilePairs(params.input, flat: true)
	  .map { prefix, r1, r2 -> tuple(getSampleID(prefix), r1, r2) }
	  .groupTuple()
	  .set {inFq_ch}
	

      process catLanes {
	    tag "Concatenating lanes into $params.workDir"
		publishDir "$params.workDir/$sampleID", mode: 'copy', pattern: "*.gz"
		label 'small_mem'
		
		input:
		tuple val(sampleID), path(R1), path(R2) from inFq_ch
				
		output:
		tuple val(sampleID), path("${sampleID}_*_init.fq.gz") into reads_ch
		
		script:
		"""
		zcat $R1 > ${sampleID}_R1_init.fq
		zcat $R2 > ${sampleID}_R2_init.fq
		
		gzip ${sampleID}_R1_init.fq
		gzip ${sampleID}_R2_init.fq
		"""
	  }
} else {
    Channel
      .fromFilePairs(params.input)
      .set {reads_ch}
}


process trim {
	tag "Trimmomatic on ${pair_id}"
	label 'med_mem'

	input:
	tuple val(pair_id), path(reads) from reads_ch
	val(libPrep) from params.libPrep
	
	output:
	path("${pair_id}_trim.log") into trimmomaticLogs_ch
	tuple pair_id, path("${pair_id}*.fastq.gz") into trimmedReads_ch, tReadsFqc_ch

	script:
	"""
	trimmomatic PE \
	  -threads $task.cpus \
	  ${reads[0]} \
	  ${reads[1]} \
	  -baseout ${pair_id}_trim \
	  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 ILLUMINACLIP:${libPrep}:2:30:10 2> ${pair_id}_trim.log
	
	mv ${pair_id}_trim_1P ${pair_id}_trim_R1.fastq
	mv ${pair_id}_trim_2P ${pair_id}_trim_R2.fastq
	gzip ${pair_id}_trim_R1.fastq
	gzip ${pair_id}_trim_R2.fastq
	"""
}

process kallisto {
    tag "Pseudoalignment and quant with Kallisto"
	label 'big_mem'
	publishDir "$params.outdir/kallisto", mode: 'copy', pattern: '*.h5'
	publishDir "$params.outdir/pseudobam", mode: 'copy', pattern: '*.bam$'
	
	input:
	tuple val(pair_id), path(reads) from trimmedReads_ch
	path (kallIndex) from params.kallIndex
	path (gtf) from params.gtf
	path (chrInfo) from params.genomeInfo
	val(nBoot) from params.boot
	
	output:
	path("${pair_id}_abundance.h5") into kallH5_ch
	path("${pair_id}_kall.log") into kallLog_ch
	tuple pair_id, path("${pair_id}_pseudoalignments.bam.bai"), path("${pair_id}_pseudoalignments.bam") into genomebam_ch
	
	script:
	def threads = task.cpus - 4
	if(threads > 0)
	"""
	   kallisto quant \
		 -i $kallIndex \
		 -o . \
		 -b $nBoot \
		 --genomebam \
		 -g $gtf \
		 -c $chrInfo \
		 -t $threads \
		 ${reads[0]} ${reads[1]} 2> "${pair_id}_kall.log"
	
	mv abundance.h5 "${pair_id}_abundance.h5"
	mv pseudoalignments.bam "${pair_id}_pseudoalignments.bam"
	mv pseudoalignments.bam.bai "${pair_id}_pseudoalignments.bam.bai"
	"""
	else
	"""
	kallisto quant \
	  -i $kallIndex \
	  -o . \
	  -b $nBoot \
	  --genomebam \
	  -g $gtf \
	  -c chrInfo \
		 ${reads[0]} ${reads[1]} 2> "${pair_id}_kall.log"
	
	mv abundance.h5 "${pair_id}_abundance.h5"
	mv pseudoalignments.bam "${pair_id}_pseudoalignments.bam"
	mv pseudoalignments.bam.bai "${pair_id}_pseudoalignments.bam.bai"
	"""
}

process fastqc {
	
	tag "FASTQC on ${sample_id}"
	label 'small_mem'
	
	input:
	tuple val(sample_id), path(reads) from tReadsFqc_ch

	output:
	path("fastqc_${sample_id}_logs") into fastqc_ch

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""  
} 




process makeBigwig{

	tag "Creating ${sampleID} bigwig"
	publishDir "$params.outdir/bigwig", mode: 'move'
	label 'big_mem'

	input:
	tuple val(sampleID), file(bamI), file(bam) from genomebam_ch
	
	
	output:
	file("${sampleID}_CPMnorm.bw")

	//TODO: add -p $task.cpus
	script:
	"""
	bamCoverage -p $task.cpus \
	  --bam ${finalBam} \
	  -o ${sampleID}_CPMnorm.bw \
	  -bs 10 --smoothLength 50 \
	  --normalizeUsing RPKM \
	  --ignoreForNormalization chrX chrY  \
	  --skipNonCoveredRegions 
	"""
}


process multiqc {
	publishDir "$params.outdir/results", mode:'move'
	label 'small_mem'

	input:
	path('*') from fastqc_ch
	  .mix(trimmomaticLogs_ch)
	  .mix(kallLog_ch)
	  .collect()
	
	output:
	path('multiqc_report.html')

	script:
	"""
	multiqc .
	"""
}








  
