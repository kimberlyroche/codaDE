use strict;
use warnings;
use POSIX;

# bulk RNA-seq
# my @features = qw(20000);
# my @evaluate_alr = qw(TRUE);
# my @filter_abundance = qw(1);
# my $label = "bulkRNAseq";
# my $label = "bulkRNAseq_ALR";

# single-cell RNA-seq
# I've been using p = 100, 200, 500, 1000, 2000, 5000, 10000, 20000
# my @features = qw(100 200 500 1000 2000 5000);
my @features = qw(10000 20000);
my @evaluate_alr = qw(FALSE);
my @filter_abundance = qw(0);
my $NB_for_DE = "TRUE"; # TRUE: NB, FALSE: log-LM + permutation

my $filename = "job.slurm";
my $f = 0;
my $ea = "";
my $fa = 0;
my $walltime = 0;
for my $i (0 .. $#features) {
  for my $j (0 .. $#evaluate_alr) {
    for my $k (0 .. $#filter_abundance) {
      $f = $features[$i];
      $ea = $evaluate_alr[$j];
      $fa = $filter_abundance[$k];
      $walltime = ceil(0.002*$f);

      open(my $fh, '>', $filename);
      print $fh '#!/bin/bash'."\n";
      print $fh '#SBATCH -J DA'.$f."\n";
      print $fh '#SBATCH --mem=16GB'."\n";
      print $fh '#SBATCH --get-user-env'."\n";
      print $fh '#SBATCH --time='.$walltime.':00:00'."\n";
      print $fh '#'."\n\n";
  
      print $fh 'module add R/3.6.1-gcb03'."\n";
      print $fh 'module add gcc/7.1.0-fasrc01'."\n\n";

      print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

      # FALSE refers to rarefication
      print $fh 'srun Rscript run.R --p='.$f.' --n=250 --k=1 --NB_for_DE='.$NB_for_DE.' --filter_abundance='.$fa."\n\n";

      close $fh;

      my $call_str = "sbatch --array=1-20 $filename";
      print("Calling: ".$call_str."\n");
      `$call_str`;

      # the world's laziest delay
      my $lazy = 0;
      while($lazy < 1000000) { $lazy++; }
    }
  }
}


