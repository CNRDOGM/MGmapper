#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
    $command .= "$ARGV[$i] ";
    $i++;
}
use Getopt::Std;
use strict;
use Cwd;
#
# Configuration
#

my $maxreadentries = 250; # In 1000
my $verbose=0;
*LOG=*STDERR;

#
# Process command line
#
sub Usage {
  print ("Usage: $0 [-h] [-s readsize] -f <name> -r <name> -a <name> -b <name> -A <name> -B <name> -i <name> -I -l <name> -v\n");
  print ("Description:\n");
  print ("$0 - find reads in commen given 2 input fastq files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -s  : read size (memory use and performance option, the bigger the better) default 250 which uses 1 GB RAM\n");
  print ("  -f  : input forward fastq\n");
  print ("  -r  : input reverse fastq\n");
  print ("  -a  : output forward fastq\n");
  print ("  -b  : output reverse fastq\n");
  print ("  -I  : output interleaved fastq [off]\n");
  print ("  -i  : output interleaved fastq [STDOUT]\n");
  print ("  -A  : output singletons from forward fastq, optional\n");
  print ("  -B  : output singletons from reverse fastq, optional\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : verbose mode [off]\n");
  print ("\n");
  exit;
} # Usage
getopts('hf:r:a:b:s:A:B:Ii:l:v')||Usage();
#
# Usage
#
&Usage($Getopt::Std::opt_h) if defined($Getopt::Std::opt_h);
&Usage() unless defined($Getopt::Std::opt_f) and defined($Getopt::Std::opt_r);
&Usage() unless (defined($Getopt::Std::opt_a) and defined($Getopt::Std::opt_b)) or defined($Getopt::Std::opt_I);
&Usage() if defined($Getopt::Std::opt_s) and $Getopt::Std::opt_s !~ m/^\d+$/;

# calulate read size (performance)
# The larger, the faster, the more memory. 250 is around 1 GB
$maxreadentries = $Getopt::Std::opt_s if defined $Getopt::Std::opt_s;
$maxreadentries *= 1000;

#
# Open forward fastq file
#
if (($Getopt::Std::opt_f=~/\.gz$/) or ($Getopt::Std::opt_f=~/\.Z$/)){
   open(F,'-|', 'gunzip', '-c', $Getopt::Std::opt_f) or die ("can't open file $Getopt::Std::opt_f: $!");
} else {
   open(F,'<', $Getopt::Std::opt_f) or die ("can't open file $Getopt::Std::opt_f: $!");
}

#
# Open reverse fastq file
#
if (($Getopt::Std::opt_r=~/\.gz$/) or ($Getopt::Std::opt_r=~/\.Z$/)){
   open(R,'-|', 'gunzip', '-c', $Getopt::Std::opt_r) or die ("can't open file $Getopt::Std::opt_r: $!");
} else {
   open(R,'<', $Getopt::Std::opt_r) or die ("can't open file $Getopt::Std::opt_r: $!");
}


# No logging
#if (defined($Getopt::Std::opt_l)){
#    open(LOG,">$Getopt::Std::opt_l");
#}


# Output files
if (! defined($Getopt::Std::opt_I)){
    if (($Getopt::Std::opt_a=~/\.gz$/) or ($Getopt::Std::opt_a=~/\.Z$/)){
	open(OUTF,'|-', "gzip -c > $Getopt::Std::opt_a") or die ("can't write file $Getopt::Std::opt_a: $!");
    } else {
	open(OUTF, '>', $Getopt::Std::opt_a) or die ("can't write file $Getopt::Std::opt_a: $!");
    }
    if (($Getopt::Std::opt_b=~/\.gz$/) or ($Getopt::Std::opt_b=~/\.Z$/)){
	open(OUTR,'|-', "gzip -c > $Getopt::Std::opt_b") or die ("can't write file $Getopt::Std::opt_b: $!");
    } else {
	open(OUTR, '>', $Getopt::Std::opt_b) or die ("can't write file $Getopt::Std::opt_b: $!");
    }
}
else{
    if (($Getopt::Std::opt_i=~/\.gz$/) or ($Getopt::Std::opt_i=~/\.Z$/)){
	open(OUTI,'|-', "gzip -c > $Getopt::Std::opt_i") or die ("can't write file $Getopt::Std::opt_i: $!");
    }
    elsif (defined($Getopt::Std::opt_i)){
	open(OUTI, '>', $Getopt::Std::opt_i) or die ("can't write file $Getopt::Std::opt_i: $!");
    }
    else{
	*OUTI = *STDOUT;
    }
}
# optional singleton output
unless (defined $Getopt::Std::opt_A) {
} elsif (($Getopt::Std::opt_A=~/\.gz$/) or ($Getopt::Std::opt_A=~/\.Z$/)){
   open(OUTSF,'|-', "gzip -c > $Getopt::Std::opt_A") or die ("can't write file $Getopt::Std::opt_A: $!");
} else {
   open(OUTSF, '>', $Getopt::Std::opt_A) or die ("can't write file $Getopt::Std::opt_A: $!");
}
unless (defined $Getopt::Std::opt_B) {
} elsif (($Getopt::Std::opt_B=~/\.gz$/) or ($Getopt::Std::opt_B=~/\.Z$/)){
   open(OUTSR,'|-', "gzip -c > $Getopt::Std::opt_B") or die ("can't write file $Getopt::Std::opt_B: $!");
} else {
   open(OUTSR, '>', $Getopt::Std::opt_B) or die ("can't write file $Getopt::Std::opt_B: $!");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
my $thisDir=cwd();
if ($verbose){
    print LOG "## Local date and time $datestring - Start program\n";
    print LOG "# $command\n";
    print LOG "# working dir: $thisDir\n\n";
}

# Initialise variables, make perl aware of the sizes
my $done = 0;
=pod
my $buf_interleaved;
my $buf_f;
my $buf_r;
if (defined($Getopt::Std::opt_I)){
    $buf_interleaved = 'x' x $buffersize*2;
}
else{
my $buf_f = 'x' x $buffersize;
my $buf_r = 'x' x $buffersize;
}
my @f;
$f[$arraysize] = 1;
my @r;
$r[$arraysize] = 1;
my @fkey;
$fkey[$arraysize] = 1;
my @rkey;
$rkey[$arraysize] = 1;
my %hash = map($_ => $_, 0..$arraysize);
$buf_f = '';
$buf_r = '';
@f = ();
@r = ();
@fkey = ();
@rkey = ();
=cut
my ($buf_f, $buf_r, $buf_interleaved) = ('', '', '');
my (@f, @r, @rkey, @fkey, %hash);
until ($done) {
    # Read a big block from both filehandles
    my $countf = $maxreadentries;
    while (--$countf and defined(my $line = <F>)) {
	next unless $line =~ m/^@(\S+)[\/|#|\t ]/;
        my $k = $1;
        $k =~ s/\/[12]$//;
      push(@fkey, $k);
      $line .= <F> . <F> . <F>;
      push(@f, $line);
   }
   my $countr = $maxreadentries;
   while (--$countr and defined(my $line = <R>)) {
      next unless $line =~ m/^@(\S+)[\/|#|\t ]/;
        my $k = $1;
        $k =~ s/\/[12]$//;
      push(@rkey, $k);
      $line .= <R> . <R> . <R>;
      push(@r, $line);
   }
   $done = 1 if $countr and $countf;
   # Find commons quick
   my ($indexf, $indexr) = (-1, -1);
   my (@fspos, $rspos);
   %hash = ();
   for (my $i = $#rkey; $i >= 0; $i--) {
      $hash{$rkey[$i]} = $i;
   }
   for (my $i = 0; $i <= $#fkey; $i++) {
      next unless exists $hash{$fkey[$i]};
      if (defined($Getopt::Std::opt_I)){
	  $buf_interleaved .= $f[$indexf = $i];
	  $buf_interleaved .= $r[$indexr = $hash{$fkey[$i]}];
      }
      else{
	  $buf_f .= $f[$indexf = $i];
	  $buf_r .= $r[$indexr = $hash{$fkey[$i]}];
      }
  }
   # Find singletons
   if (defined $Getopt::Std::opt_A) {
      my $target = $done ? $#fkey : $indexf-1;
      my $buffer = '';
      for (my $i = 0; $i <= $target; $i++) {
         $buffer .= $f[$i] unless exists $hash{$fkey[$i]};
      }
      print OUTSF $buffer;
   }
   if (defined $Getopt::Std::opt_B) {
      %hash = ();
      for (my $i = $#fkey; $i >= 0; $i--) {
         $hash{$fkey[$i]} = $i;
      }
      my $target = $done ? $#rkey : $indexr-1;
      my $buffer = '';
      for (my $i = 0; $i <= $target; $i++) {
         $buffer .= $f[$i] unless exists $hash{$rkey[$i]};
      }
      print OUTSR $buffer;
   }
   # Clean arrays
   if ($indexf >= $#fkey) {
      @f = ();
      @fkey = ();
   } else {
      splice(@f, 0, $indexf+1);
      splice(@fkey, 0, $indexf+1);
   }
   if ($indexr >= $#rkey) {
      @r = ();
      @rkey = ();
   } else {
      splice(@r, 0, $indexr+1);
      splice(@rkey, 0, $indexr+1);
   }
   # Output
      if (defined($Getopt::Std::opt_I)){
	  print OUTI $buf_interleaved;
      }
      else{
	  print OUTF $buf_f;
	  print OUTR $buf_r;
      }
      $buf_interleaved = '';
      $buf_f = '';
      $buf_r = '';
}

close(F);
close(R);
close(OUTF);
close(OUTR);
close(OUTSF) if defined $Getopt::Std::opt_A;
close(OUTSR) if defined $Getopt::Std::opt_B;
close(OUTI) if defined $Getopt::Std::opt_i;
