#!/usr/bin/perl -w

use strict;
use Data::Printer;
use List::MoreUtils qw/uniq/;
use Array::Utils qw(:all);
use Getopt::Long;
use Cwd qw/getcwd/;

my $usage=<<USAGE;

    This script is order to repair the LFR inforamtion !!! 
    
    Example: perl $0 <nodebarcode> [-option]

        -a          :[node|none] allow to plus the node and edge to complete the LFR [default:node]
        -n          : the file name
        -o          : the directory 

USAGE

my $barcode_stat = $ARGV[0];
my ($allow,$name,$outdir);

GetOptions(
    "a:s" => \$allow,
    "n:s" => \$name,
    "o:s" => \$outdir
);

die $usage if (!$ARGV[0]);

$allow ||='node';
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
    my @nod_set = split/\,/,$node;
    my $nod_num = scalar(@nod_set);
    next if ($nod_num == 1);
    my @edge_set = split/\,/,$edge_nod;
    my $edge_num = scalar(@edge_set);
    if ($allow eq 'none'){
        if ($edge_num - $nod_num <= 1){
            push @{$set{$edge_num}{$nod_num}{$barcode}},($edge_num,$edge_nod,$nod_num,$node,$bili);
        }
    }elsif($allow eq 'node'){
        if ($edge_num - $nod_num == 1){
            push @{$set{$edge_num}{$nod_num}{$barcode}},($edge_num,$edge_nod,$nod_num,$node,$bili);
        }else{
            my @n_edge =uniq(ntoedge(\@nod_set));
            my $n_edge_num = scalar(@n_edge);
            if ($n_edge_num - $nod_num == 1){
                my @new_bili_line;
                my @new_edge_line;
                foreach my $aedge (sort {(split/_/,$a)[0] <=> (split/_/,$b)[0]} @n_edge){
                    push @new_bili_line,($aedge.':1');
                    push @new_edge_line,$aedge; 
                }
                my $edge_lin = join(",",@new_edge_line);
                my $bili_lin = join(";",@new_bili_line);
                push @{$set{$n_edge_num}{$nod_num}{$barcode}},($n_edge_num,$edge_lin,$nod_num,$node,$bili_lin);
            }
        }
    }else{
        print 'please input the corret parmaters !!!'."\n";
        exit;
    }
}
$/="\n";
close INA;

open OUTA,">$outdir/$name.sort.barcode.stat" or die $!;
foreach my $keys (sort {$a<=>$b} keys %set){
    my %new_hash = %{$set{$keys}};
    foreach my $keys2 (sort {$a<=>$b} keys %new_hash){
        my %new_hash_two = %{$new_hash{$keys2}};
        next if ($keys2 == 1);   #repair this !
        foreach my $keys3 (keys %new_hash_two){
            my @this_data_set = @{$new_hash_two{$keys3}};
            print OUTA "$keys3\t".join("\t",@this_data_set)."\n";
        }
    }
}

close OUTA;

#####################################
sub ntoedge{
    my $input = $_[0];
    my @a_nod_set;
    foreach my $tnod (@$input){
        my ($l,$r) = split/\//,$tnod;
        my $left_sight = judge_odd_even($l,'l');
        my $right_sight = judge_odd_even($r,'r');
        push @a_nod_set,($left_sight,$right_sight);
    }
    my @sort_edge = re_sort_edge(\@a_nod_set);
    return @sort_edge;
}

###################################
sub judge_odd_even{
    my ($input,$direct) = @_;
    my $yu = $input % 2;
    my $new_link;
    if ($yu == 1 ){
        my $asight = $input - 1;
        if ($direct eq 'l'){
            $new_link = $asight.'_'.$input;
        }else{
            $new_link = $input.'_'.$asight;
        }
    }else{
        my $asight = $input + 1;
        if ($direct eq 'r'){
            $new_link = $input.'_'.$asight;
        }else{
            $new_link = $asight.'_'.$input;
        }
    }
    return $new_link;
}

############################################
sub re_sort_edge{
    my $input = $_[0];
    my @new_set;
    foreach my $edge (@$input){
        my ($l,$r)=(split/_/,$edge)[0,1];
        if ($l > $r){
            my $new_edge = $r.'_'.$l;
            push @new_set,$new_edge;
        }else{
            push @new_set,$edge;
        }
    }
    return @new_set;
}