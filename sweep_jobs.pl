use strict;
use warnings;
use POSIX;

my @features = qw(20000);
my @evaluate_alr = qw(TRUE);
my @filter_abundance = qw(0);

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
      $walltime = ceil(0.004*$f);

      open(my $fh, '>', $filename);
      print $fh '#!/bin/bash'."\n";
      print $fh '#SBATCH -J DA'.$f."\n";
      print $fh '#SBATCH --mem=16GB'."\n";
      print $fh '#SBATCH --get-user-env'."\n";
      print $fh '#SBATCH --time='.$walltime.':00:00'."\n";
      print $fh '#'."\n\n";
  
      print $fh 'module add R/3.5.1-gcb01'."\n";
      print $fh 'module add gcc/6.2.0-fasrc01'."\n\n";

      print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

      # FALSE refers to rarefication
      print $fh 'srun Rscript run.R '.$f.' 250 ALR_nofilter '.$ea.' '.$fa.' FALSE'."\n\n";

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


