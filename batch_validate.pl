use strict;
use warnings;
use POSIX;
use List::Util qw(min);

my @datasets =  qw(VieiraSilva Barlow Song Monaco Hagai Klein Owens Yu);
my @RAM =       qw(12          12     24   32     32    64    64    64);
my @hours =     qw(2           2      2    2      2     2     2     2);
my @threshold = qw(1           1      1    2      1     1     1     1);

my $model_dir = "self_nopartial";
my $norm = "FALSE";
my $classify = "TRUE";
my $alpha = 0.9;
my $permodel = "TRUE";
my $self_baseline = "TRUE";
my $partials = "FALSE";
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

  print $fh 'srun Rscript validate.R --dataset='.$datasets[$i].' --threshold='.$threshold[$i].' --norm='.$norm.' --classify='.$classify.' --alpha='.$alpha.' --permodel='.$permodel.' --selfbaseline='.$self_baseline.' --partials='.$partials."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }
  $i++;
}
