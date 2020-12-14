use strict;
use warnings;
use List::Util qw(sum);

my @props = qw();
my $filename = "";
my $folder = "slurm_output";

opendir(DIR, $folder) or die "cannot open directory";
my @docs = grep(/\.out$/, readdir(DIR));
foreach my $file (@docs) {
  open(RES, $folder."/".$file) or die "could not open $file\n";
  while(my $line = <RES>) {
    chomp $line;
    if($line =~ /^Using tissues:/) {
      #print($line."\n");
    }
    if($line =~ /^Proportion DE genes: (.*?)$/) {
      print($1."\n");
      push(@props, $1);
    }
  }
  close RES;
}

# print("Average: ".(sum(@props)/@props)."\n");

