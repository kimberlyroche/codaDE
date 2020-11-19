use strict;
use warnings;
use POSIX;

my @evaluate_alr = qw(FALSE);
my @filter_abundance = qw(5);

# sweep over number of features
#my @features = qw(1000 5000 10000 20000);
#my @de_sweep = qw(0.5);
#my @corr_sweep = qw(0);
#my $method = "edgeR";

# sweep over DE (discarding library size info)
#my @features = qw(20000);
#my @de_sweep = qw(0.2 0.4 0.6 0.8);
#my @corr_sweep = qw(0);
#my $method = "edgeR";

# sweep over DE (retaining library size info)
#my @features = qw(20000);
#my @de_sweep = qw(0.2 0.4 0.6 0.8);
#my @corr_sweep = qw(0.5);
#my $method = "NB";

# sweep over DE (using TMM normalization)
my @features = qw(2000);
my @de_sweep = qw(0.2 0.4 0.6 0.8);
my @corr_sweep = qw(0);
my $method = "edgeR_TMM";

my $filename = "job.slurm";
my $f = 0;
my $ea = "";
my $fa = 0;
my $de = 0;
my $sf_corr = 0;
my $walltime = 0;
for my $i (0 .. $#features) {
  for my $j (0 .. $#evaluate_alr) {
    for my $k (0 .. $#filter_abundance) {
      for my $de_idx (0 .. $#de_sweep) {
        for my $corr_idx (0 .. $#corr_sweep) {
          $f = $features[$i];
          $ea = $evaluate_alr[$j];
          $fa = $filter_abundance[$k];
          $de = $de_sweep[$de_idx];
          $sf_corr = $corr_sweep[$corr_idx];
          $walltime = ceil(0.001*$f);

          open(my $fh, '>', $filename);
          print $fh '#!/bin/bash'."\n";
          print $fh '#SBATCH -J DA'.$f."\n";
          print $fh '#SBATCH --mem=16GB'."\n";
          print $fh '#SBATCH --get-user-env'."\n";
          print $fh '#SBATCH --time='.$walltime.':00:00'."\n";
          print $fh '#'."\n\n";
  
          print $fh 'module add R/3.6.1-gcb03'."\n";
          print $fh 'module add gcc/7.3.0-gcb01'."\n\n";

          print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

          print $fh 'srun Rscript run.R --p='.$f.' --n=100 --k=1 --prop_da='.$de.' --sf_corr='.$sf_corr.' --method='.$method.' --filter_abundance='.$fa."\n\n";

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
  }
}


