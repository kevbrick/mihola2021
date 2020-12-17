nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "Mihola 2021 PIPELINE (Version 1.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/mihola2021.nf \\"
  log.info " --bam           <path to bwa bam files>"
  log.info " --outdir        <string: default = output> \\"
  log.info " "
  exit 1
  }

// Params:
params.name           = 'mihola_2021'
params.outdir         = "${launchDir}/output"
params.accessorydir   = "${launchDir}/accessoryFiles"

params.genome_fa      = "${params.accessorydir}/genome.fa"
params.genome_fai     = "${params.accessorydir}/genome.fa.fai"

process makeWinFiles {

  memory '4.GB'
  time '1.hour'

  container  "docker://kevbrick/mihola2021:1.0"

  output:
  path("win*bed", emit:bed)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles

  grep chr17 accessoryFiles/genome.fa.fai >genome.fai
  bedtools makewindows -g genome.fai -w 147         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 147 )\'  >win147.bed
  bedtools makewindows -g genome.fai -w 500 -s 50   | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500 )\'  >win500s50.bed
  bedtools makewindows -g genome.fai -w 1000 -s 150 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 1000 )\' >win1000s150.bed
  """
  }

process makeFigures {

  publishDir "${params.outdir}/coverage", mode: 'copy', overwrite: true, pattern: '*bedgraph'
  publishDir "${params.outdir}/coverage", mode: 'copy', overwrite: true, pattern: '*bam'
  publishDir "${params.outdir}/coverage", mode: 'copy', overwrite: true, pattern: '*gtf'
  publishDir "${params.outdir}/coverage", mode: 'copy', overwrite: true, pattern: '*tab'
  publishDir "${params.outdir}/fig", mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir "${params.outdir}/fig", mode: 'copy', overwrite: true, pattern: '*png'

  memory '4.GB'
  time '1.hour'

  container  "docker://kevbrick/mihola2021:1.0"

  input:
  path(winbeds)
  path(bams)

  output:
  path("*bedgraph", emit:bg)
  path("*bam", emit:bam)
  path("region.gtf", emit:gtf)
  path("*tab", emit:tab)
  path("*png", emit:png)
  path("*pdf", emit:pdf)

  script:
  """
  ln -s ${params.accessorydir} accessoryFiles

  echo -e "chr17\t17529000\t18040000"  >region.bed

  ## Get genes
  wget http://hgdownload.soe.ucsc.edu/goldenPath/rn5/bigZips/genes/rn5.ensGene.gtf.gz
  gunzip rn5.ensGene.gtf.gz
  intersectBed -a rn5.ensGene.gtf -b region.bed >region.gtf

  ## Get HS
  cp accessoryFiles/allMergedRatHotspots.finalTable.tab .
  cut -f1-3,5 allMergedRatHotspots.finalTable.tab |grep -v from >koHS.bedgraph
  intersectBed -a region.gtf  -b koHS.bedgraph -c |perl -lane '\$_ =~ /^.+\\s+(\\d+)\\s*\$/; \$hs = \$1; print join("\\t",@F[2..4],\$F[6],\$hs)' >region.gtf.forR.tab

  ## Parse region for each BAM
  intersectBed -a win1000s150.bed -b region.bed -wa -u >region.w1000s150.bed
  intersectBed -a win500s50.bed   -b region.bed -wa -u >region.w500s50.bed
  intersectBed -a win147.bed      -b region.bed -wa -u >region.ws147.bed

  echo -e "cs\\tfrom\\tto\\tcoverage\\tsample" >allCoverage.Rdata.tab

  for bam in *.bam; do
    regBam="region_"\$bam
    samtools view -hb \$bam chr17:17519000-18050000 >\$regBam
    samtools index \$regBam

    n1=\${regBam/region_Mihola_et_al_2021_/}
    name=\${n1/.rn5.bam/}

    bg1000=\${regBam/.bam/w1000s150.bedgraph}
    bg500=\${regBam/.bam/w500s50.bedgraph}
    bg147=\${regBam/.bam/ws147.bedgraph}

    intersectBed -a region.w1000s150.bed -b \$regBam  -c |sort -k1,1 -k2n,2n >\$bg1000
    intersectBed -a region.w500s50.bed   -b \$regBam  -c |sort -k1,1 -k2n,2n >\$bg500
    intersectBed -a region.ws147.bed     -b \$regBam -c  |sort -k1,1 -k2n,2n >\$bg147

    cat \$bg1000 |perl -lane 'chomp; print \$_."\\t'\$name'"' >>allCoverage.Rdata.tab

  done

  ## Make coverage fig
  R --no-save <accessoryFiles/scripts/drawCoverageFigRN5MiholaEtAl2021.R

  ## Make reviewer figs
  R --no-save <accessoryFiles/scripts/drawReviewerFigsMiholaEtAl2021.R
  """
  }

// OK ... let's start
workflow {
  bams = Channel.fromPath("${params.accessorydir}/bam/*bam*")

  wins = makeWinFiles()

  fig5 = makeFigures(wins.bed.collect(),
                   bams.collect())
  }
