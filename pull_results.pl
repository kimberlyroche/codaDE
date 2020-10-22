use strict;
use warnings;

my $label = "analysis1";
my @chunks = qw(0 200 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400 2600);
my $filename = "pull_results.slurm";
for my $ch (0 .. $#chunks) {
  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J results_pull'."\n";
  print $fh '#SBATCH --mem=16GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time=1:00:00'."\n";
  print $fh '#'."\n\n";
  print $fh 'module add R/3.6.1-gcb03'."\n";
  print $fh 'module add gcc/7.1.0-fasrc01'."\n\n";
  print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";
  print $fh 'srun Rscript pull_results.R '."$label $chunks[$ch]\n\n";
  close $fh;
  my $call_str = "sbatch pull_results.slurm";
  `$call_str`;
}
