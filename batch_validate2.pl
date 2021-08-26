use strict;
use warnings;
use POSIX;
use List::Util qw(min);

# Note: This is for likely non-compositional data sets.

#my @datasets =  qw(Muraro ESCA Gruen Morton Athanasiadou Ferreira Kimmerling Hashimshony);
#my @RAM =       qw(12     24   12    12     24           64       24         12);
#my @hours =     qw(2      2    2     2      2            2        2          2);
#my @threshold = qw(1      1    1     1      1            1        1          1);

my @datasets =  qw(Muraro Gruen Morton Athanasiadou Ferreira Kimmerling Hashimshony);
my @RAM =       qw(12     12    12     24           64       24         12);
my @hours =     qw(2      2     2      2            2        2          2);
my @threshold = qw(1      1     1      1            1        1          1);

my $filename = "job.slurm";
my $i = 0;
my $end = $#datasets;

#$i = 4;
#$end = 4;

while($i <= $end) {
  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J '.$datasets[$i]."\n";
  print $fh '#SBATCH --mem='.$RAM[$i].'GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time='.$hours[$i].':00:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

  print $fh 'srun Rscript validate_noprediction.R --dataset='.$datasets[$i].' --baseline=self --threshold='.$threshold[$i]."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }
  $i++;
}
