use strict;
use warnings;
use POSIX;
use List::Util qw(min);

#my @datasets =  qw(VieiraSilva Barlow Song Monaco Hagai Klein Owens Yu);
#my @RAM =       qw(12          12     24   32     32    64    64    64);
#my @hours =     qw(2           2      2    2      2     2     2     2);
#my @threshold = qw(1           1      1    2      3     2     1     1);

my @datasets = qw(Song Song Song Song Song Monaco Monaco Monaco Monaco Monaco Hagai Hagai Hagai Hagai Hagai Klein Klein Klein Klein Klein);
my @RAM = qw(24 24 24 24 24 32 32 32 32 32 32 32 32 32 32 64 64 64 64 64);
my @hours = qw(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
my @threshold = qw(1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5);

my $baseline = "self";
my $partial_flag = "nopartial";

my $filename = "job.slurm";
my $i = 0;
my $end = $#datasets;

while($i <= $end) {
  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J '.$threshold[$i].'_'.$datasets[$i]."\n";
  print $fh '#SBATCH --mem='.$RAM[$i].'GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time='.$hours[$i].':00:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

  print $fh 'srun Rscript validate.R --dataset='.$datasets[$i].' --baseline='.$baseline.' --threshold='.$threshold[$i].' --model_folder='.$baseline.'_'.$partial_flag."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }
  $i++;
}
