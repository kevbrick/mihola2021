params.accessorydir   = "${launchDir}/accessoryFiles"

profiles {
  local {
    process.maxForks = 1
    process.executor='local'
    env{
      TMPDIR='./'
    }
  }

  slurm {
    process.executor='slurm'
    process.scratch = '/lscratch/$SLURM_JOBID'
    process.clusterOptions = ' --gres=lscratch:800 --partition=norm'
    env{
      TMPDIR='/lscratch/$SLURM_JOBID'
    }
  }

  none {
    // Add custom configs here
  }
}

singularity.enabled = true
singularity.autoMounts = true
singularity.envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,NXF_GENOMES'
singularity.runOptions=" -B \$NXF_GENOMES -B ${params.accessorydir} -B ${launchDir} "
process.container="docker://kevbrick/mihola2021:1.0"

params.outdir = "${launchDir}/output"

report {
  enabled = true
  file = "${params.outdir}/nxfReports/mihola_2021_report.html"
}

timeline {
  enabled = true
  file = "${params.outdir}/nxfReports/mihola_2021_timeline.html"
}

trace {
  enabled = true
  file = "${params.outdir}/nxfReports/mihola_2021_trace.txt"
}

manifest {
  description = '2020: Kevin Brick'
}

