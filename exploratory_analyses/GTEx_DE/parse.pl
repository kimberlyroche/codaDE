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
  while(<RES>){
    open my $in, "<:encoding(utf8)", $folder."/".$file or die "$file: $!";
    while (my $line = <$in>) {
      chomp $line;
      if($_ =~ /^Proportion DE genes: (.*?)$/) {
        push(@props, $1);
      }
    }
    close $in;
  }
}

print("Average: ".(sum(@props)/@props)."\n");
