use strict;
use warnings;
use POSIX;
use List::Util qw(min);

my $p = 100;
my $input = "input_".$p.".txt";
my $ln = 375;
my $chunks = 5;

my $filename = "job.slurm";

my $chunk_sz = floor($ln / $chunks);
my $start = 1;
my $end = -1;
while($start <= $ln) {
  $end = min($ln, $start + $chunk_sz - 1);

  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J gen_'.$start.'-'.$end."\n";
  print $fh '#SBATCH --mem=16GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time=2:00:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/codaDE'."\n\n";

  print $fh 'srun Rscript generate_data.R --input='.$input.
                                        ' --output=res_'.$p.
                                        ' --start='.$start.
                                        ' --end='.$end."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }

  $start = $start + $chunk_sz;
}
