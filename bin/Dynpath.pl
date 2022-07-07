#!/usr/bin/perl -w

use strict;
use List::MoreUtils qw/uniq/;
use List::Util qw/sum/; 

my ($nod,$edge_list,$outdir,$name)=@ARGV;

my $usage=<<USAGE;

Usage:perl $0 <nod> <edge> <outdir> <name>

USAGE

die $usage if (@ARGV < 4);

print '======start:'.localtime.'======='."\n";
my %nod_set;
open INB,$nod or die $!;
while(<INB>){
    chomp;
    my $set = (split/\s+/,$_)[0];
    $nod_set{$set} = 1;
    
}
close INB;

my %edge;    
my @edge_set;
open INC,$edge_list or die $!;
while(<INC>){
    chomp;
    next if (/^\#/);
    my ($edg,$num)=(split/\s+/,$_)[0,1];
    if ($num > 0){
        $edge{$edg}=$num;
        push @edge_set,$edg;
    }
}
close INC;

open OUTA,">$outdir/$name.info" or die $!;
my $origin = $edge_set[0];
my %path;
$path{$origin} ='1';
while(1){
    my @list = ();
    @list = keys %path;
    foreach my $edge1 (@list){
        my ($edge_confict,$nod_confict)=control_confict($edge1,\%edge,\%nod_set);
        my $confict = $edge_confict+$nod_confict;
        ###perfect set ##########################################
        if ($confict == 0 ){
            print OUTA "0\t$edge1\n";
            delete $path{$edge1};       
        }
        my @index = split/_/,$edge1;
        my $edge_lindex = $index[0];
        my $edge_rindex = $index[-1];
        my $new_path;
        foreach my $vert (keys %nod_set){
            my ($nod_lindex,$nod_rindex)=split/\//,$vert;
            my $new_nod;
            if ($edge_rindex == $nod_lindex or $edge_rindex == $nod_rindex){
                $new_nod = $nod_rindex if ($edge_rindex == $nod_lindex);
                $new_nod = $nod_lindex if ($edge_rindex == $nod_rindex);
                foreach my $edge2 (@edge_set){
                    my @index2=split/_/,$edge2;
                    my $lindex =$index2[0];
                    my $rindex = $index2[-1];
                    my $new_edge;                      
                    if ($lindex == $new_nod or $rindex == $new_nod){
                        $new_edge = $edge2 if ($lindex == $new_nod);
                        $new_edge = $rindex.'_'.$lindex if ($rindex == $new_nod);
                        $new_path = $edge1.'/'.$new_edge;
                        if (control_the_cnv($new_path,\%edge) eq 'T'){
                                $path{$new_path} = '1';
                        }
                    }
                }
            }
        }
        delete $path{$edge1};        
    }
    last if (keys %path == 0);
}

print '======end:'.localtime.'======='."\n";
##########sub control_the_cnv($new_path,\%edge) ####################
sub control_the_cnv{
    my ($new_path,$hash)=@_;
    my @stat;
    my @new = reverse_edge($new_path);
    my @uniq =uniq(@new);
    foreach my $re_ed (@uniq){
        my $num = grep {$_ eq $re_ed} @new;
        if ($num <= ${$hash}{$re_ed}){
            push @stat,'1'; 
        }else{
            push @stat,'0';
        }
    }
    if (grep {$_ eq '0'} @stat){
        return 'F';
    }else{
        return 'T';
    }
}

#################sub control_confict #######################
sub control_confict{
    my ($edg,$edge,$nod_set) = @_;
    my @nodset = split/_/,$edg;
    my $first_poin = $nodset[0]; 
    shift @nodset;
    my $last_poin = $nodset[-1];
    pop @nodset;
    my @new_edge_set = reverse_edge($edg);
    my @erro_edge;
    ##stat the edge confitct
    my %double_edge = %$edge;
    my @uniq_arr = uniq(@new_edge_set);
    foreach my $de2 (@uniq_arr){
        my $num = grep {$_ eq $de2} @new_edge_set;
        if ($num == ${$edge}{$de2}){
            push @erro_edge,'0';
        }elsif($num < ${$edge}{$de2}){
            my $err_num = ${$edge}{$de2} - $num;
            push @erro_edge,$err_num;
        }else{
            my $err_num = $num - ${$edge}{$de2};
            push @erro_edge,$err_num;
        }
        delete $double_edge{$de2} if (exists $double_edge{$de2});
    }
    #print Dumper(%double_edge);
    my @edgeleft = values %double_edge;
    push @erro_edge,@edgeleft;
    my $edge_confict = sum(@erro_edge);
    ##stat the nod confict ##
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
    #print Dumper(%double_nod);
    my $nod_confict= values %double_nod;
    return ($edge_confict,$nod_confict);
}

########re-sort the sub_path edge set likes 11_10→10_11 ##########
sub reverse_edge{
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



