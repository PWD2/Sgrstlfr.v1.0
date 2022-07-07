#!/usr/bin/perl -w

use strict;
use Data::Printer;
use Cwd qw/getcwd/;

my $usage=<<USAGE;

    Example: perl $0 <nodebarcode> <name> <outdir>

USAGE

my ($barcode_stat,$name,$outdir) = @ARGV;

die $usage if (!$ARGV[0]);
$name ||='test';
$outdir ||=getcwd;

$/="\@";
my %set;
open INA,$barcode_stat or die $!;
<INA>;
while(<INA>){
    chomp;
    my ($id,$dis,$bili,$node) = split/\n/,$_;
    my ($barcode,$edge_nod) = split/\s+/,$id;
    my $edge_num = split/\,/,$edge_nod;
    my $nod_num = split/\,/,$node;
    push @{$set{$edge_num}{$nod_num}{$barcode}},($edge_num,$edge_nod,$nod_num,$node,$bili);
}
$/="\n";
close INA;

open OUTA,">$outdir/$name.sort.barcode.stat" or die $!;
foreach my $keys (sort {$a<=>$b} keys %set){
    my %new_hash = %{$set{$keys}};
    foreach my $keys2 (sort {$a<=>$b} keys %new_hash){
        my %new_hash_two = %{$new_hash{$keys2}};
        next if ($keys2 < $keys - 1);
        foreach my $keys3 (keys %new_hash_two){
            my @this_data_set = @{$new_hash_two{$keys3}};
            print OUTA "$keys3\t".join("\t",@this_data_set)."\n";
        }
    }
}

close OUTA;


