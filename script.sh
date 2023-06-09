#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=10G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=20G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=120:00:00   # job requires up to 24 hours of runtime
#$ -r n               # if job crashes, it should be restarted

date
hostname

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                         # e.g. "did my job exceed its memory request?

INPUT=/wynton/scratch/finterly_run_results_newer/results_controls
OUTPUT=/wynton/scratch/finterly_vcf/vcfs_results_controls

NXF_VER=22.11.0-edge nextflow run Making_gVCFs.nf -profile sge,apptainer --inputdir $INPUT --outdir $OUTPUT 

exit 0