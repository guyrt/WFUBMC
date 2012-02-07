#!/usr/bin/perl
#
# run_stat_tests
#
# Runs SNPLASH and calls R to evaluate output.  One or more engines
# may be specified on the command line.
#
# David R. McWilliams <dmcwilli@wfubmc.edu>
#
# 11-Aug-2011 Start

# Return codes in this script follow the unix convention that
# not-zero is an error.

use strict ;
use warnings ;

my $usage = <<END;
usage:
run_stat_tests all|adtree|bagging|dandelion|dprime|intertwolog|qsnpgwa|snpgwa

Run the script from /trunk/t/ or edit the script to change the paths.
END

die $usage
  unless scalar @ARGV > 0 ;

my $snplash = "../snplash" ;
# my $snplash = "/home/drm/proj/snplash/WFUBMC/snplash" ;

my $sim2000_geno    = "./Sim2000/sim2000.geno" ;
my $sim2000_pheno   = "./Sim2000/sim2000.covphen" ;
my $sim2000_map     = "./Sim2000/sim2000.map" ;

my $hisp_geno  = "./QsnpgwaTest/hisp/hisp.geno" ;
my $hisp_pheno = "./QsnpgwaTest/hisp/hisp.phen" ;
my $hisp_map   = "./QsnpgwaTest/hisp/hisp.map" ;

my %engines = () ;
initialize(\%engines) ;

my @to_do = () ;
if (grep(/all/i, @ARGV)) {
  @to_do = keys %engines ;
} else {
  push @to_do, lc($_) foreach @ARGV ;
}

open LOG, ">", "run_stat_tests.log" or
  die "Could not open log file: $!" ;
print LOG "Starting snplash evaluation ", scalar localtime(), ".\n" ;
print LOG "Engines: ", join(", ", @to_do), " supplied to $0.\n" ;

my $failed = 0 ;

TEST:
foreach my $engine (sort @to_do) {
  if (defined $engines{$engine}) {
    if($engines{$engine}{run}->()) {
      $failed++ ;
      next TEST ;
    }
    if($engines{$engine}{test}->()) {
      $failed++ ;
      next TEST ;
    }
  } else {
    print LOG "Do not recognize engine \"$engine\". Skipping.\n" ;
    next TEST ;
  }
  print LOG "Done with $engine, ", scalar localtime(), ".\n" ;
} # end foreach engine

print "There were $failed tests; check log and output.\n"
  if $failed ;

print LOG "Finished snplash evaluation ", scalar localtime(), ".\n" ;

exit 0 ;

# end main

##############################################################################
#                                                                            #
#                             Functions                                      #
#                                                                            #
##############################################################################

##############################################################################
#
# intitialize
#
# Initialize the dispatch table of engines and tests.
#
# Accepts: Hash reference for the table
# Returns: Success
#
# 11-Aug-2011

sub initialize {
  my $href = shift ;

  $href->{adtree}{run}  = \&run_adtree ;
  $href->{adtree}{test} = \&test_adtree ;

  $href->{bagging}{run}  = \&run_bagging ;
  $href->{bagging}{test} = \&test_bagging ;

  $href->{dandelion}{run}  = \&run_dandelion ;
  $href->{dandelion}{test} = \&test_dandelion ;

  $href->{dprime}{run}  = \&run_dprime ;
  $href->{dprime}{test} = \&test_dprime ;

  $href->{intertwolog}{run}  = \&run_intertwolog ;
  $href->{intertwolog}{test} = \&test_intertwolog ;

  $href->{qsnpgwa}{run}  = \&run_qsnpgwa ;
  $href->{qsnpgwa}{test} = \&test_qsnpgwa ;

  $href->{snpgwa}{run}  = \&run_snpgwa ;
  $href->{snpgwa}{test} = \&test_snpgwa ;

  return 1 ;
} # end initialize


sub run_adtree {
  print LOG "No test for --engine adtree yet.\n" ;
  return 0 ;
} # end run_adtree

sub test_adtree {
  print LOG "No test for --engine adtree yet.\n" ;
  return 0 ;
} # end test_adtree

sub run_bagging {
  print LOG "No test for --engine bagging yet.\n" ;
  return 0 ;
} # end run_bagging

sub test_bagging {
  print LOG "No test for --engine bagging yet.\n" ;
  return 0 ;
} # end test_bagging

sub run_dandelion {
  print LOG "No test for --engine dandelion yet.\n" ;
  return 0 ;
} # end run_dandelion

sub test_dandelion {
  print LOG "No test for --engine dandelion yet.\n" ;
  return 0 ;
} # end test_dandelion

##############################################################################
#
# parse_dprime
#
# Parse the dprime output to make it conform to the R output.
#
# Accepts: Scalar with name of file to parse
# Returns: Success
#
# 01-Sep-2011
#
sub parse_dprime {
  my $fh = shift ;

  unless (open DP, "<", $fh) {
    print LOG "run_stat_tests::parse_dprime: Could not open $fh for input: $!" ;
    return 1 ;
  }

  # This is the file name used in test_dprime.Snw
  my $fout = "sim2000.dprime.just100.out.txt" ;
  unless (open OUT, ">", $fout) {
    print LOG "run_stat_tests::parse_dprime: Could not open $fout for output: $!" ;
    return 1 ;
  }

  my $skip = 18 ;
  while ($skip) {
    <DP> ;
    $skip-- ;
  }

  my $next_snp = 1 ;
  my $max_snp  = 1900 ;

  # Print the first appearance of a snp and the next 99.
 SNP100:
  while ( defined(my $line = <DP>) && ($next_snp <= $max_snp) ) {
    chomp $line ;
    my ($m1, $m2, $d, $dpmulti, $rsqbi, $deltabi) = split(/\s+/, $line) ;

    my $idx = 0 ;
    if ($m1 =~ m/^SNP(\d+)/) {
      $idx = $1 ;
    }
    if ($idx == $next_snp) {
      print OUT join("\t", $m1, $m2, $d, $dpmulti, $rsqbi, $deltabi), "\n" ;
      for (my $i=0; $i<99 and $line = <DP>; $i++) {
        chomp $line ;
        print OUT join("\t", split(/\s+/, $line)), "\n" ;
      }
      $next_snp++ ;
    }
  } # end while

  return 0 ;
}  # end parse_dprime

sub run_dprime {
  my $outfile = "sim2000.dprime.out.txt" ;
  my $cmd = join(" ", $snplash,
                 "-engine", "dprime",
                 "-geno",   $sim2000_geno,
                 "-phen",   $sim2000_pheno,
                 "-map",    $sim2000_map,
                 "-trait",  "ds",
                 "-out",    $outfile) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::run_dprime: Procedure run_dprime failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e $outfile) {
    print LOG "run_stat_tests::run_dprime: Output from snplash/dprime not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::run_dprime: dprime ran ok.\n" ;

  if (parse_dprime($outfile)) {
    print LOG "run_stat_tests::run_dprime: Failed to parse $outfile.\n" ;
    return 1 ;
  } else {
    print LOG  "run_stat_tests::run_dprime: Parsed $outfile.\n" ;
    return 0 ;
  }
} # end run_dprime

sub test_dprime {
  my $file = "test_dprime.Snw" ;

  my $cmd = join(" ", "R CMD Sweave", $file) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::test_dprime: Procedure test_dprime failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e "test_dprime.tex") {
    print LOG "run_stat_tests::test_dprime: test_qnspgwa.tex not produced.\n" ;
    return 1 ;
  }
  $cmd = "pdflatex test_dprime.tex" ;
  system($cmd) ;

  $val = $? ;
  if ($val) {
    print LOG "run_stat_tests::test_dprime: pdflatex failed; final report not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::test_dprime: R tests ran ok.\n" ;

  unlink("test_dprime.aux") ;
  unlink("test_dprime.log") ;
  unlink("test_dprime.tex") ;

  return 0 ;
} # end test_dprime

sub run_intertwolog {
  my $outfile = "sim2000.intertwolog.out.txt" ;
  my $cmd = join(" ", $snplash,
                 "-engine", "intertwolog",
                 "-geno",   $sim2000_geno,
                 "-phen",   $sim2000_pheno,
                 "-map",    $sim2000_map,
                 "-cov",    "cov1,cov2",
                 "-trait",  "ds",
                 "-out",    $outfile) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::run_intertwolog: Procedure run_intertwolog failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e $outfile) {
    print LOG "run_stat_tests::run_intertwolog: Output from snplash/intertwolog not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::run_intertwolog: intertwolog ran ok.\n" ;
  return 0 ;

} # end run_intertwolog

sub test_intertwolog {
  my $file = "test_intertwolog.Snw" ;

  my $cmd = join(" ", "R CMD Sweave", $file) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::test_intertwolog: Procedure test_intertwolog failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e "test_intertwolog.tex") {
    print LOG "run_stat_tests::test_intertwolog: test_qnspgwa.tex not produced.\n" ;
    return 1 ;
  }
  $cmd = "pdflatex test_intertwolog.tex" ;
  system($cmd) ;

  $val = $? ;
  if ($val) {
    print LOG "run_stat_tests::test_intertwolog: pdflatex failed; final report not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::test_intertwolog: R tests ran ok.\n" ;

  unlink("test_intertwolog.aux") ;
  unlink("test_intertwolog.log") ;
  unlink("test_intertwolog.tex") ;

  return 0 ;
} # end test_intertwolog

sub run_qsnpgwa {
  my $outfile = "hisp.qsnpgwa.out.txt" ;
  my $cmd = join(" ", $snplash,
                 "-engine", "qsnpgwa",
                 "-geno",   $hisp_geno,
                 "-phen",   $hisp_pheno,
                 "-map",    $hisp_map,
                 "-trait",  "response",
                 "-out",    $outfile) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::run_qsnpgwa: Procedure run_qsnpgwa failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e $outfile) {
    print LOG "run_stat_tests::run_qsnpgwa: Output from snplash/qsnpgwa not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::run_qsnpgwa: qsnpgwa ran ok.\n" ;
  return 0 ;

} # end run_qsnpgwa


sub test_qsnpgwa {
  my $file = "test_qsnpgwa.Snw" ;

  my $cmd = join(" ", "R CMD Sweave", $file) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::test_qsnpgwa: Procedure test_qsnpgwa failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e "test_qsnpgwa.tex") {
    print LOG "run_stat_tests::test_qsnpgwa: test_qnspgwa.tex not produced.\n" ;
    return 1 ;
  }
  $cmd = "pdflatex test_qsnpgwa.tex" ;
  system($cmd) ;

  $val = $? ;
  if ($val) {
    print LOG "run_stat_tests::test_qsnpgwa: pdflatex failed; final report not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::test_qsnpgwa: R tests ran ok.\n" ;

  unlink("test_qsnpgwa.aux") ;
  unlink("test_qsnpgwa.log") ;
  unlink("test_qsnpgwa.tex") ;

  return 0 ;
} # end test_qsnpgwa

sub run_snpgwa {
  my $outfile = "sim2000.snpgwa.out.txt" ;
  my $cmd = join(" ", $snplash,
                 "-engine", "snpgwa",
                 "-geno",   $sim2000_geno,
                 "-phen",   $sim2000_pheno,
                 "-trait",  "ds",
                 "-map",    $sim2000_map,
                 "-cov",    "cov1,cov2",
                 "--val",
                 "-out",    $outfile) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::run_snpgwa: Procedure run_snpgwa failed; system returned: $val.\n" ;
    return 1 ;
  }

  unless (-e $outfile) {
    print LOG "run_stat_tests::run_snpgwa: Output from snplash/snpgwa not produced.\n" ;
    return 1 ;
  }

  print LOG "run_stat_tests::run_snpgwa: snpgwa ran ok.\n" ;
  return 0 ;

} # end run snpgwa

sub test_snpgwa {
  my $file = "test_snpgwa.Snw" ;

  my $cmd = join(" ", "R CMD Sweave", $file) ;

  system($cmd) ;
  my $val = $? ;

  print $val, "\n" ;
  if ($val) {
    print LOG "run_stat_tests::test_snpgwa: Procedure test_snpgwa failed; system returned: $val.\n" ;
    return 1 ;
  }
  unless (-e "test_snpgwa.tex") {
    print LOG "run_stat_tests::test_snpgwa: test_snpgwa.tex not produced.\n" ;
    return 1 ;
  }

  $cmd = "pdflatex test_snpgwa.tex" ;
  system($cmd) ;

  $val = $? ;
  if ($val) {
    print LOG "run_stat_tests::test_snpgwa: pdflatex failed; final report not produced.\n" ;
    return 1
  }

  print LOG "run_stat_tests::test_snpgwa: R tests ran ok.\n" ;

  unlink("test_snpgwa.aux") ;
  unlink("test_snpgwa.log") ;
  unlink("test_snpgwa.tex") ;

  return 0 ;
} # end test snpgwa


# end functions

# end run_stat_tests
