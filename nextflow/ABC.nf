#!/usr/bin/env nextflow
// Copyright (C) 2020 Kade Pettie

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Nextflow pipeline to make Activity-By-Contact enhancer-gene predictions
// using ATAC-seq + HiChIP data normalized for comparison between groups.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

///////// DOWNLOAD DATA OPTION ////////

process clone_abc {
  label 'git'
  storeDir params.git_dir

  when:
  params.mode =~ /(all)/

  input:
  val(gl) from Channel.of( params.download.abc_git )

  output:
  path("ABC-Enhancer-Gene-Prediction/**") into BIN1

  script:
  """
  git clone -b ${params.download.abc_branch} $gl
  """

}

process clone_genrich {
  label 'git'
  storeDir params.git_dir

  when:
  params.mode =~ /(all)/

  input:
  val(gl) from Channel.of( params.download.genrich_git )

  output:
  path("Genrich/**") into BIN2

  script:
  """
  git clone $gl; \
  cd Genrich; \
  make; \
  cd ..
  """

}


///////// PEAK CALLING INPUT //////////

// ATAC overrides DHS
(BAM_ATAC,
  BAM_DHS) = ( params.atac_bamdir==""
    ? [ Channel.empty(),
        Channel.fromPath( params.dhs_bamdir + '/' + params.dhs_glob + '.bam' )
          .map{it -> ['dhs', it]}]
    : [ Channel.fromPath( params.atac_bamdir + '/' + params.atac_glob + '.bam' )
          .map{it -> ['atac', it]},
        Channel.empty() ]
  )

BAM_ATAC
  .mix(BAM_DHS)
  .tap{ BSI }
  .tap{ GCRE_BAM }
  .combine( Channel.fromPath( params.chrsize ) )
  .combine( BIN1.combine(BIN2) )
  .map{ it -> [ it[0], it[1], it[2] ] } // add channels with code so it gets cloned and stored but don't link it in each workdir
  .set{ BAMS }

if (params.candidates_dir=="") {
  BAMS
    .set{ BAMS_CALL_PEAKS }
} else {
  BAMS_CALL_PEAKS = Channel.empty()
}

//////// cCRE INPUT ////////

Channel.fromPath( params.chrsize )
  .combine( Channel.fromPath(params.blocklist) )
  .combine( Channel.fromPath(params.promoters) )
  .set{ CRE }


/////// PROCESSES ////////

process call_peaks {
  label 'genrich'
  conda params.abc_code + "/macs.yml"
  publishDir "${params.outdir}/peaks/", pattern: "*.narrowPeak.sorted"
  publishDir "${params.outdir}/peaks/", pattern: "*.bedgraph"

  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
  maxRetries 5

  when:
  params.mode =~ /(all)/

  input:
  tuple seqtype, path(bam), path(sizes) from BAMS_CALL_PEAKS

  output:
  tuple gn, path("*.narrowPeak.sorted") into PEAKS
  path("*.bedgraph")

  script:
  // TODO: use macs2 if seqtype=='dhs'
  // TODO: ensure bam is name-sorted for genrich (?)
  gn = bam.baseName
  """
  ${params.genrich_bin}/Genrich -t $bam \
  -o ${gn}.narrowPeak \
  -v \
  -f ${gn}.pileup.qvals.bedgraph \
  -k ${gn}.pileup.pvals.bedgraph \
  -y \
  -j \
  -d 151 \
  -p 0.1; \
  bedtools sort -faidx $sizes \
  -i ${gn}.narrowPeak \
  > ${gn}.narrowPeak.sorted
  """

}

BSI
  .map{ it -> [it[0], it[1].baseName, it[1]] }
  .set{ BSIGN }

process sort_index {
  label 'samtools'
  publishDir "${params.outdir}/bams/", pattern: "*.psort.bam"
  publishDir "${params.outdir}/bams/", pattern: "*.psort.bam.bai"

  when:
  params.mode =~ /(all)/

  input:
  tuple seqtype, gn, path(bam) from BSIGN

  output:
  tuple gn, seqtype, path("*.bam"), path("*.bam.bai") into MBAM, CBAM

  script:
  outbam = bam.baseName + '.psort.bam'
  """
  samtools sort \
  --threads ${params.samcores} \
  -o $outbam \
  $bam; \
  samtools index \
  -@ ${params.samcores} \
  $outbam
  """

}

CRE
  .combine( PEAKS
              .combine(MBAM, by: 0) )
  .set{ MCRE }

process call_ccres {
  label 'ccres'
  publishDir "${params.outdir}/cCREs/", pattern: "*.bed"
  conda params.abc_code + "/abcenv.yml"

  when:
  params.mode =~ /(all)/

  input:
  tuple path(sizes), path(excl), path(incl), gn, path(peaks), seqtype, path(bam), path(idx) from MCRE

  output:
  tuple gn, seqtype, path("*.candidateRegions.bed") into GCRE_FROM_PEAK
  path("*.Counts.bed")

  script:
  // should peakExtendFromSummit be different for dhs vs atac?
  """
  python ${params.abc_code}/src/makeCandidateRegions.py \
  --narrowPeak $peaks \
  --bam $bam \
  --outDir ./ \
  --chrom_sizes $sizes \
  --regions_blocklist $excl \
  --regions_includelist $incl \
  --peakExtendFromSummit 250 \
  --nStrongestPeaks 150000
  """

}

// Mix GCRE from any provided candidate regions
// to skip peak calling/makeCandidateRegions
(GCRE_SUPP,
  GCRE_NEW) = ( params.candidates_dir==""
    ? [ Channel.empty(),
        GCRE_FROM_PEAK ]
    : [ GCRE_BAM.map{ it -> [ it[1].baseName, it[0] ] }
          .combine( Channel.fromPath( params.candidates_dir + '/' + params.candidates_glob ) ),
        Channel.empty() ]
  )

GCRE_NEW
  .mix(GCRE_SUPP)
  .combine(CBAM, by: [0,1]) // match un-merged bams with peaks called on corresponding merged bams
  .combine( Channel.fromPath(params.genes) )
  .combine( Channel.fromPath(params.expression) )
  .combine( Channel.fromPath(params.chrsize) )
  .combine( Channel.fromPath(params.chrsizebed) )
  .combine( Channel.fromPath(params.ubiquitous) )
  .set{ CATAC }

process atac_activity {
  label 'ccres'
  publishDir "${params.outdir}/neighborhoods/${sn}/"
  conda params.abc_code + "/abcenv.yml"

  when:
  params.mode =~ /(all)/

  input:
  tuple gn, seqtype, path(ccres), path(bam), path(idx), path(g), path(e), path(s), path(sb), path(u) from CATAC

  output:
  tuple sn, path("EnhancerList.txt"), path("GeneList.txt") into PDX

  script:
  sn = bam.simpleName
  if (seqtype=='atac') {
    """
    python ${params.abc_code}/src/run.neighborhoods.py \
    --candidate_enhancer_regions $ccres \
    --genes $g \
    --ATAC $bam \
    --expression_table $e \
    --chrom_sizes $s \
    --ubiquitously_expressed_genes $u \
    --gene_name_annotations symbol \
    --primary_gene_identifier symbol \
    --cellType $sn \
    --outdir ./
    """
  } else {
    """
    python ${params.abc_code}/src/run.neighborhoods.py \
    --candidate_enhancer_regions $ccres \
    --genes $g \
    --DHS $bam \
    --expression_table $e \
    --chrom_sizes $s \
    --ubiquitously_expressed_genes $u \
    --gene_name_annotations symbol \
    --primary_gene_identifier symbol \
    --cellType $sn \
    --outdir ./
    """
  }


}

HIC = ( params.hic_dir==""
        ? Channel.empty()
        : Channel.fromPath(params.hic_dir)
            .map{ it -> ['HiC', it] }
      )

PWRLAW = ( params.power_law
           ? Channel.of( ['powerlaw', ""] )
           : Channel.empty()
          )

PDX
  .combine( Channel.fromPath(params.chrsize) )
  .combine( HIC.mix(PWRLAW) )
  .set{ PRED }

process hichip_predict {
  label 'ccres'
  publishDir "${params.outdir}/predictions/${ctype}/${sn}/"
  conda params.abc_code + "/abcenv.yml"

  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
  maxRetries 5

  when:
  params.mode =~ /(all)/

  input:
  tuple sn, path(enh), path(gene), path(s), ctype, c from PRED

  output:
  tuple ctype, sn, path("EnhancerPredictionsAllPutative.txt.gz") into ANAL
  path("*NonExpressedGenes.txt.gz")
  path("*.bedpe")
  path("*.txt")

  script:
  if (ctype=="HiC") {
    if (params.hichip) {
      """
      python ${params.abc_code}/src/predict.py \
      --enhancers $enh \
      --genes $gene \
      --HiCdir $c \
      --hic_type bedpe \
      --hichip \
      --chrom_sizes $s \
      --hic_resolution 5000 \
      --window 5000000 \
      --score_column ABC.Score \
      --threshold .02 \
      --cellType $sn \
      --outdir ./ \
      --make_all_putative
      """
    } else {
      """
      python ${params.abc_code}/src/predict.py \
      --enhancers $enh \
      --genes $gene \
      --HiCdir $c \
      --hic_type ${params.hic_type} \
      --chrom_sizes $s \
      --hic_resolution 5000 \
      --window 5000000 \
      --scale_hic_using_powerlaw \
      --score_column ABC.Score \
      --threshold .02 \
      --cellType $sn \
      --outdir ./ \
      --make_all_putative
      """
    }
  } else {
    """
    python ${params.abc_code}/src/predict.py \
    --enhancers $enh \
    --genes $gene \
    --chrom_sizes $s \
    --hic_resolution 5000 \
    --window 5000000 \
    --scale_hic_using_powerlaw \
    --score_column powerlaw.Score \
    --threshold .02 \
    --cellType $sn \
    --outdir ./ \
    --make_all_putative
    """
  }


}
