#!/usr/bin/perl -w

use strict;
use List::MoreUtils qw/uniq/;
use List::Util qw/sum/; 

my ($name,$outdir,$confict)=@ARGV;

my $usage=<<USAGE;

Description: Link the edge and the nod ! Though produce the dot graph!
Usage:perl $0 <name|list> <outdir> <confict>

USAGE

die $usage if (@ARGV < 2);

$confict ||='0';

print '======start:'.localtime.'======='."\n";

my @list;
if ($name =~ /\.txt/){
    open INA,$name or die $!;
    while(<INA>){
        chomp;
        push @list,$_;
    }
    close INA;
}else{
    push @list,$name;
}

foreach my $filename (@list){
    open INB,"$outdir/$filename.nod" or die $!;
    my %nod_all;
    my %nod_hash_set;
    while(<INB>){
        chomp;
        next if (/^\#/);
        my ($a_n,$num) = split/\s+/,$_;
        $nod_hash_set{$a_n} = '1';
        my ($le,$ri) = split/\//,$a_n;
        push @{$nod_all{$le}},($a_n);
        my $ri_an = $ri.'/'.$le;
        push @{$nod_all{$ri}},($ri_an);
    }
    close INB;
    my %edge_cn_hash;
    my %edge_index_hash;
    my @uniq_edge_set;
    my $j=1;
    open INC,"$outdir/$filename.edg" or die $!;
    while(<INC>){
        chomp;
        next if (/^\#/);
        my ($edge,$cn) = (split,/\s+/,$_)[0,1];
        next if ($cn == 0);
        $edge_cn_hash{$edge} = $cn;
        # push @origin_edge,$edge if ($j==1);
        push @uniq_edge_set,$edge if ($cn == 1);
        my ($l,$r) = split/_/,$edge;
        $edge_index_hash{$l} = ($edge);
        my $b_edge = $r.'_'.$l;
        $edge_index_hash{$r} = ($b_edge);
        $j++;
    }
    close INC;
    my $genome_length = sum(values %edge_cn_hash);
    open OUTA,">$outdir/$filename.info" or die $!;
<<<<<<< HEAD
=======
    open OUTB,">$outdir/$filename.dot" or die $!;
    print OUTB 'digraph dot_graph{'."\n";
    print OUTB '  rankdir=LR'."\n";
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
    my %dot_set;
    my $stat_number=0;
    my @origin_edge = ($uniq_edge_set[0]);
    while(1){
        my $edge1 = shift @origin_edge;
        my @this_solution_length = split/\//,$edge1;
        if (@this_solution_length >= int($genome_length/2)){
            my $edg_c = $genome_length - scalar(@this_solution_length);
            my $nod_c = control_confict($edge1,\%nod_hash_set);
            if ($nod_c + $edg_c <= $confict){
                $stat_number++;
                print OUTA "$edge1\t$edg_c\t$nod_c\n";
            }
        }
        my $edge_rindex = (split/_/,$edge1)[-1];
        if (exists $nod_all{$edge_rindex}){
            my @this_nod = @{$nod_all{$edge_rindex}};
            foreach my $vert (@this_nod){
                my $new_nod = (split/\//,$vert)[1];
                if (exists $edge_index_hash{$new_nod}){
                    my $new_path = $edge1.'/'.$edge_index_hash{$new_nod};
                    if (control_the_cnv($new_path,\%edge_cn_hash) eq 'T'){
<<<<<<< HEAD
=======
                            # my $line = ("  \"$this_solution_length[-1]\"->\"$edge_index_hash{$new_nod}\"\[label=\"$vert\"\]");
                            my $line = ("  \"$this_solution_length[-1]\"->\"$edge_index_hash{$new_nod}\"");
                            $dot_set{$line} = '1';
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
                            push @origin_edge,$new_path;
                    }
                }
            }
        }
        if (@origin_edge == 0){
            last;
        }
    }
<<<<<<< HEAD
    print 'We get the solution number: '.$stat_number."\n";
    close OUTA;
=======
    print OUTB join("\n",keys %dot_set)."\n";
    print OUTB '}'."\n";
    print 'We get the solution number: '.$stat_number."\n";
    close OUTA;
    close OUTB;
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
}

print '======end:'.localtime.'======='."\n";

##########sub control_the_cnv($new_path,\%edge) ####################
sub control_the_cnv{
    my ($this_path,$this_hash)=@_;
    my @set = re_sort_edge($this_path);
    my @uniq =uniq(@set);
    my @stat;
    my $check = 'T';
    foreach my $re_ed (@uniq){
        my $num = grep {$_ eq $re_ed} @set;
        if ($num <= ${$this_hash}{$re_ed}){
            next; 
        }else{
            $check = 'F';
            last;
        }
    }
    return $check;
}

########re-sort the sub_path edge set likes 11_10→10_11 ##########
sub re_sort_edge{
    my $input = $_[0];
    my @array = split/\//,$input;
    my @new_set;
    foreach my $edge (@array){
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

#################sub control_confict #######################
sub control_confict{
    my ($edg,$nod_set) = @_;
    my @nodset = split/_/,$edg;
    my $first_poin = shift @nodset; 
    my $last_poin = pop @nodset;
    my $local;
    if ($first_poin > $last_poin){
        $local = $last_poin.'/'.$first_poin;
    }else{
        $local = $first_poin.'/'.$last_poin;
    }
    my %double_nod = %$nod_set;
    delete $double_nod{$local} if (exists $double_nod{$local});
    foreach my $nd (@nodset){
        my($ln,$rn)=split/\//,$nd;
        if ($ln > $rn){
            $nd = $rn."/".$ln;
        }
        if ($double_nod{$nd}){
            delete $double_nod{$nd} if (exists $double_nod{$nd});
        }
    }
    my $nod_confict= values %double_nod;
    return $nod_confict;
}

