\This was compiled with gcc/6.2.0

To submit job arrays of simulations/DE evaluations use:
    sbatch --array=1-20 job.slurm

# Update to R/3.6.1
Use R/3.6.1-gcb03 and
    gcc/7.1.0-fasrc01

# 2020-08-30 and after
# There's something funny with batch submitting evaluations on existing data sets:
# when sweeps over proportion of DE genes and correlations between true and observed
# total counts are c(0.1, 0.2, ..., 0.8, 0.9) all cases involving 0.3 and 0.7 are 
# skipped???
# I'm not sure why this happens but as there's no harm in re-running conditions
# that have already been evaluated (they'll be skipped), this can be worked around 
# by shuffling the order of evaluated conditions -- e.g. c(0.3, 0.7, 0.1, ..., 0.9)
# These jobs seem to work if run any other way, so maybe it's an issue of timing/
# some kind race condition?

Stand-alone files in the root dir:
  generate_empirical_DA.R
    -- Use GTEx data to built a distribution of quantile-quantile empirical differential
       expression across tissue types
  pull_results.*
    -- these file parse `simulated_data/metadata.tsv` and pull the results of run
       (tp, fp, etc.) for use by the predictive modeling scripts in `exploratory`
  unfinished_jobs.*
    -- identify missing conditions in the "metadata" file
  write_metadata.*
    -- re-write or update the "metadata" file; in theory this shoudln't need to be
       manually updated but just in case something goes wrong...
