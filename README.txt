This was compiled with gcc/6.2.0

To submit job arrays of simulations/DE evaluations use:
    sbatch --array=1-20 job.slurm

# Update to R/3.6.1
module add R/3.6.1-gcb03
module add gcc/7.1.0-fasrc01

2020-08-30 and after

There's something funny with batch submitting evaluations on existing data sets: when sweeps over proportion
of DE genes and correlations between true and observed total counts are c(0.1, 0.2, ..., 0.8, 0.9) all cases
involving 0.3 and 0.7 are skipped???

I'm not sure why this happens but as there's no harm in re-running conditions that have already been evaluated
(they'll be skipped), this can be worked around by shuffling the order of evaluated conditions -- e.g.
c(0.3, 0.7, 0.1, ..., 0.9) These jobs seem to work if run any other way, so maybe it's an issue of timing/
some kind race condition?
