# --------------------------------------------------------------------------------------------------------
#   GENERATING DATA
# --------------------------------------------------------------------------------------------------------

1) Run `generate_settings.R --p=XXX` first to find missing datasets/parameterizations. This script sweeps 
through and enumerates all parameter settings, saving them to a file.

2) Run `batch_gen.pl`, setting P and TOTAL JOB NUMBER in the script. This will put the results in the 
`temp` folder.

    This calls `generate_data.R --input=XXX --output=XXX --start=XXX --end=XXX`, which generates the new 
    data sets in chunks and makes baseline differential abundance calls on the absolute values, writing
    results to a file tagged with row numbers.

3) Run `store_gen_results.R` to write results (in `temp`) to the database.

# --------------------------------------------------------------------------------------------------------
#   EVALUATING JOBS
# --------------------------------------------------------------------------------------------------------

1) Run `pull_eval_uuids.R` to enumerate jobs as:
   > Rscript pull_eval_uuids.R --p=100

2) Run `batch_eval.pl`, setting P and JOB START and END LINE NUBMERS in the script.
   This will put the results in the `temp` folder.

3) Run `store_eval_results.R`
   This will add ANY results that exist in files in the `temp` directory to the DB!

4) Run `render_rates.R` to calculate TPR and FPR for all results. Note: This takes about 2 hr. to run if
all TPR / FPR fields in the database must be updated!

# --------------------------------------------------------------------------------------------------------
#   CHARACTERIZING DATA SETS
# --------------------------------------------------------------------------------------------------------

1) Run `pull_char_uuids.R` to enumerate jobs as:
   > Rscript pull_char_uuids.R --p=100 --file=input

2) Run `batch_char.pl`, setting P and JOB START and END LINE NUBMERS in the script.
   This will put the results in the `temp` folder.

3) Run `store_char_results.R`
   This will add ANY results that exist in files in the `temp` directory to the DB!

# --------------------------------------------------------------------------------------------------------
#   OLDER SETUP STUFF
# --------------------------------------------------------------------------------------------------------

This was compiled with gcc/6.2.0

To submit job arrays of simulations/DE evaluations use:
    sbatch --array=1-20 job.slurm

# Update to R/3.6.1
module add R/3.6.1-gcb03
module add gcc/7.3.0-gcb01
module add libpng/1.5.21-fasrc01
# module add jpeg/6b-fasrc02

---------- 2021-03-16

INSTALLING DESeq2 + dependencies

Currently I'm having an issue with these dependencies of DESeq2
	Warning messages:
	1: In install.packages(...) :
	  installation of package ‘XML’ had non-zero exit status
	2: In install.packages(...) :
	  installation of package ‘png’ had non-zero exit status
	3: In install.packages(...) :
	  installation of package ‘jpeg’ had non-zero exit status
	4: In install.packages(...) :
	  installation of package ‘annotate’ had non-zero exit status
	5: In install.packages(...) :
	  installation of package ‘latticeExtra’ had non-zero exit status
	6: In install.packages(...) :
	  installation of package ‘genefilter’ had non-zero exit status
	7: In install.packages(...) :
	  installation of package ‘geneplotter’ had non-zero exit status
	8: In install.packages(...) :
	  installation of package ‘Hmisc’ had non-zero exit status
	9: In install.packages(...) :
	  installation of package ‘DESeq2’ had non-zero exit status

To install the XML package I had to:
- Use the updated version of libxml2 I have installed in /data/mukherjeelab/Mongrel/libxml2 (those path/link includes
  are already set up in /home/ker48/.bashrc)
- Download an older version of XML (3.9-4) and install it from source
- Include a bunch of flags in my /home/ker48/.R/Makevars BY HAND; these are specified in the "Configuration information"
  in the XML install.packages output but apparently don't get applied automatically; these are currently commented out
- Disable Intel/MKL for C++ compilation (apparently those libraries don't play nice)

To install the png packages I had to:
- Use module add libpng/1.5.21-fasrc01; I got strange compile-time issues with libpng/1.6.23-fasrc01 ("/usr/bin/ld: 
  skipping incompatible libpng16.so when searching for -lpng16) that implied it had been compiled for the wrong
  CPU architecture (???)

To install annotate I had to:
- REENABLE specifically the Intel/MKL components

To install jpeg I had to:
- module add jpeg/6b-fasrc02

The rest of the dependencies and DESeq2 installed fine after that.

INSTALLING Seurat + dependencies

R package RANN

	This requires the ANN library compiled (source here: http://www.cs.umd.edu/~mount/ANN/)
        Download and extract this and simply run `make` + OS choice; it will build to the extracted directory
        Minimally you need to add these include directories to the .R/Makevars:
		-I/data/mukherjeelab/Mongrel/ann_1.1.2/include
		-I/data/mukherjeelab/Mongrel/ann_1.1.2/include/ANN
	(I also added the /bin and /lib directories to the PATH and LD_LIBRARY_PATH respectively but I don't
	know if this was necessary.)

After this most dependencies fall into place except `spatstat`, which needs to be installed as < v2. I used
spatstat_1.64-1.

Then installed Seurat_3.2.3 (since Seurat 4 is not compatible with R 3.6.1).

---------- 2020-08-30

There's something funny with batch submitting evaluations on existing data sets: when sweeps over proportion
of DE genes and correlations between true and observed total counts are c(0.1, 0.2, ..., 0.8, 0.9) all cases
involving 0.3 and 0.7 are skipped???

I'm not sure why this happens but as there's no harm in re-running conditions that have already been evaluated
(they'll be skipped), this can be worked around by shuffling the order of evaluated conditions -- e.g.
c(0.3, 0.7, 0.1, ..., 0.9) These jobs seem to work if run any other way, so maybe it's an issue of timing/
some kind race condition?

