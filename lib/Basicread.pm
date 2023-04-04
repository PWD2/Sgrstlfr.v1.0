package Basicread;

###This module is purpose to read the syn-yeast-chr file.

use strict;
use warnings;
use Exporter;
use List::Util qw/min max/;

use vars qw($VERSION @ISA @EXPORT_OK);

our @EXPORT_OK = qw(read_region read_edge read_coverage read_fasta mersh_edge_region read_sub_list read_pathlist);

our $VERSION  = 1.00;
our @ISA =qw(Exporter);

##### This sub is input the `region` and output the '1'=> ('IXR',1,700) ######
sub read_region{
    my $input = $_[0];
    my %loxp_region;
    open INA,$input or die $!;
    while(<INA>){
        chomp;
        next if  (/^\#/);
        my ($chrid,$start,$end,$in_dex)=(split)[0,3,4,8];
        @{$loxp_region{$in_dex}}= ($chrid,$start,$end);
    }
    close INA;
    return %loxp_region;
}

##### This sub is input the `edge` and output the '0_1'=>(1,1) #######
sub read_edge{
    my $input = $_[0];
    my %edge;
    open INB,$input or die $!;
    while(<INB>){
        chomp;
        next if (/^\#/);
        my @set = split/\s+/,$_;
        if (@set == 4){
            @{$edge{$set[0]}} = ($set[1],$set[3]);
        }elsif (@set == 5){
            @{$edge{$set[0]}} = ($set[1],$set[4]);
        }else{
            print "Wrong!!!The edge file line not is 4 or 5 !!!\n";
            exit;
        }
    }
    close INB;
    return %edge;
}

##### This sub is input the `[samtools depth -a *.sort.bam > depth.txt]` and output the 'IXR'=>`(1,2,1,1,3)` ######
sub read_coverage{
    my ($input,$chrid) = @_;
    open INC,$input or die $!;
    my @coverage;
    while(<INC>){
        chomp;
        my ($chr,$point,$depth)=split/\s+/,$_;
        if ($chr eq $chrid){
            push @coverage,$depth;
        }
    }
    close INC;
    return @coverage;
}

#### This sub is read the fasta file and return the %fasta ################
sub read_fasta{
    my ($input,$chr)=@_;
    open IND,$input or die $!;
    $/=">";
    <IND>;
    my $syn_seq;
    while(<IND>){
        chomp;
        my ($i_d,$s_eq)=(split/\n/,$_,2)[0,1];
        if ($i_d =~ /$chr/){
            $s_eq =~ s/\r/\n/g;
            $s_eq =~ s/\n//g;
            $syn_seq = $s_eq;
        }
    }
    $/="\n";
    unless($syn_seq){
        print '! Not find the appoint chrmosome name seqence. Please corret the fasta file and chrmosome name!'."\n";
        exit;
    }
    close IND;
    return $syn_seq;
}

######## mersh the edge and region information. it can get the new loxp region for other use #####
sub mersh_edge_region{
    my ($edge,$region)=@_;
    my %edge_hash = read_edge($edge);
    my %region_hash = read_region($region);
    my %result;
    foreach my $keys (keys %edge_hash){
        my @aset = @{$edge_hash{$keys}};
        my @index_set =split/\,/,$aset[1]; 
        if (@index_set == 1){
            my @region_range = @{$region_hash{$index_set[0]}};
            my $length = $region_range[2] - $region_range[1] +1;
            @{$result{$keys}} =($aset[0],$length,$region_range[1],$region_range[2],$aset[1]); 
        }else{
            my $max_index = max(@index_set);
            my $min_index = min(@index_set);
            my $max_point = ${$region_hash{$max_index}}[2];
            my $min_point = ${$region_hash{$min_index}}[1];
            my $length = $max_point - $min_point +1 ;
            @{$result{$keys}} =($aset[0],$length,$min_point,$max_point,$aset[1]);
        }
    }
    return %result;
}

####### read the true genome path ###########
sub read_pathlist{
    my $pathlist = $_[0];
    open INF,$pathlist or die $!;
    my $corret_path;
    while(<INF>){
        chomp;
        next if (/^\#/);
        if (/edge path/){
            $corret_path = (split/\s+/,$_)[2];
        }
    }
    unless ($corret_path){
        print "The pathlist file is emtpy? Not corrt path!\n";
        exit;
    }else{
        return $corret_path;
    }
}

1;


