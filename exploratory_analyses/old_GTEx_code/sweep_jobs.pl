use strict;
use warnings;
use POSIX;

my $filename = "job_baseline.slurm";
for my $i (1..49) {
  # Iterate tissue type indices 1 through 49
  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J base'.$i."\n";
  print $fh '#SBATCH --mem=20GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time=00:30:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'module add R/3.6.1-gcb03'."\n";
  print $fh 'module add gcc/7.1.0-fasrc01'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/codaDE/exploratory_analyses/GTEx_DE'."\n\n";

  print $fh 'srun Rscript build_baseline_model_1.R '.$i."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }
}

