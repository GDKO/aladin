#!/usr/bin/env perl

# To do
# Multi-fasta
# Parallel
###

use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;
use List::Util 'sum';
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($samfile,$fastafile,$prefix) = ("","","corrected");

GetOptions (
 "s|samfile:s" => \$samfile,
 "f|fastafile:s" => \$fastafile,
 "p|prefix:s" => \$prefix,
);

if (not $fastafile or not $samfile) {
  print STDERR "Usage:
  -s|samfile <samfile>
  -f|fastafile <fastafile>
  -p|prefix <string>
";
exit;
}

my @ref_bases;
my $header;
my $seqio = Bio::SeqIO->new(-file => $fastafile, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
    my $string = $seq->seq;
    $header = $seq->id;
    @ref_bases = split("",$string);
}

open IN, "$samfile";
my %ref_positions;
my %ins_positions;
my $i = 0;
my $k = 0;

while (my $line=<IN>) {
  next if $line=~/^\@/;
  my @p =split(/\t/,$line);
  my $read_name = $p[0];
  my $sam_flag = $p[1];
  my $read_start = $p[3];
  my $cigar = $p[5];
  my $read = $p[9];
  if ($sam_flag != 0 && $sam_flag != 16){
    next;
  } #avoid multimapping

  # Parse CIGAR string
  my @type_counts = split(/\D+/,$cigar);
  my @types = split(/\d+/,$cigar);
  shift @types;
  $i = 0;
  my @cs;
  while ($i<scalar(@types)) {
    $k = 0;
    while ($k<$type_counts[$i]) {
      push @cs, $types[$i];
      $k++;
    }
    $i++;
  }

  # Create hash with M/D and I
  my @rb = split("",$read);
  $i = 0;
  $k = 0;
  while ($i<scalar(@cs)) {
    if ($cs[$i] eq "M") {
      $ref_positions{$read_start+$i}{$rb[$k]}++;
    }
    elsif ($cs[$i] eq "D") {
      $ref_positions{$read_start+$i}{"-"}++;
      $k--;
    }
    elsif ($cs[$i] eq "I") {
      $ins_positions{$read_start+$i}{$read_name}.=$rb[$k];
      $read_start--;
    }
    else {
      $read_start--;
    }
    $i++;
    $k++;
  }
}

my %insertions;
my %ins_base;
# Count the insertions
foreach my $key (keys %ins_positions) {
  foreach my $skey (keys %{$ins_positions{$key}}){
    $ins_base{$key}{$ins_positions{$key}{$skey}}++;
  }
  $insertions{$key} = sum values %{$ins_base{$key}};
}


open CORRFASTA, ">$prefix.fasta";
open CORRLST, ">$prefix.lst";

my $len = scalar(@ref_bases)+1;

my $final_ref = "";
$i = 1;
while ($i<$len) {

  my $nocov = 0;
  unless (defined($ref_positions{$i})) {$ref_positions{$i}{$ref_bases[$i-1]}=1;print CORRLST "[N] at position $i\n";$nocov=1;} # No LR cov

  my @bases = (sort {$ref_positions{$i}{$b} <=> $ref_positions{$i}{$a}} keys %{$ref_positions{$i}})[0..1];
  my $base = $bases[0];
  my $sbase = $bases[1];
  my $total_match = sum values %{$ref_positions{$i}};

  unless (defined($insertions{$i})) {$insertions{$i}=0;} # No insertions at the position
  my $perc = int(($insertions{$i}/$total_match)*100);
  #if ($insertions{$i}>0) {print CORRLST "pos $i $total_match:" . Dumper($ins_base{$i});}
  if ($perc > 50) {
    my $max_value = (sort {$ins_base{$i}{$b} <=> $ins_base{$i}{$a}} keys %{$ins_base{$i}})[0];
    $final_ref .= uc($max_value);
    print CORRLST "[I] at position $i: $max_value\n";
  }

  # Check that deletion is the most probable
  if ($base eq "-") {
    if (int(($ref_positions{$i}{$base}/$total_match)*100) > 50) {
      $base = "";
      print CORRLST "[D] at position $i\n";
    }
    else {
      $base = uc($sbase);
    }
  }

  # Write the substitution output
  if ($base ne $ref_bases[$i-1] && $base ne "") {
    print CORRLST "[S] at position $i: $ref_bases[$i-1] to $base\n";
  }
  # Lowercase if no coverage
  if ($nocov) {$final_ref .= lc($base);}
  else {$final_ref .= uc($base);}
  $i++;


}

print CORRFASTA ">$header\n$final_ref\n";
