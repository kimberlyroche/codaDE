use strict;
use warnings;
use POSIX;
use List::Util qw(min);

# For P = 5000 settings, a conservative estimate is ~2 min. per setting.

my $p = 5000;
my $input = "input_eval_".$p.".txt";
my $output = "output_eval_".$p;
my $start = 1;
my $end = 298;
my $chunks = 20;

my $filename = "job.slurm";
my $n = $end - $start + 1;
my $chunk_sz = ceil($n / $chunks);

my $i = $start;
my $j = min($i + $chunk_sz - 1, $end);

while($i <= $end) {
  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J eval_'.$i.'-'.$j."\n";
  print $fh '#SBATCH --mem=32GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time=1:00:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

  print $fh 'srun Rscript evaluate_methods.R --input='.$input.
                                           ' --output='.$output.
                                           ' --start='.$i.
                                           ' --end='.$j."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }

  $i = $j + 1;
  $j = min($i + $chunk_sz - 1, $end);
}
