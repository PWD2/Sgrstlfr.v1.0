package Classbarcode;

use strict;
use warnings;
use Exporter;
use Array::Utils qw(:all);

use vars qw($VERSION @ISA @EXPORT_OK);

our @EXPORT_OK = qw(class_barcode same_class_barcode);

our $VERSION  = 1.10;
our @ISA =qw(Exporter);

###########This sub is class the barcode by the edge and nod number ####
sub class_barcode{
    my %filter_hash;
    my ($barcode,$sort_edge,$sort_nod,$stand) = @_;
    my $edge_num = scalar(@$sort_edge);
    my $nod_num = scalar(@$sort_nod);
    open INA,$barcode or die $!;
    while(<INA>){
        chomp;
        my ($id,$this_edge_num,$ed,$this_nod_num,$no,$bili)=split/\s+/,$_;
        my @this_edge = split/\,/,$ed;
        my @this_no = split/\,/,$no;
        if ($this_edge_num == $edge_num + $stand){
            if (intersect(@this_edge,@$sort_edge) == $edge_num){
                if (intersect(@this_no,@$sort_nod) == $nod_num){
                    @{$filter_hash{$id}} = ($ed,$no,$bili);
                }
            }
        }
    }
    close INA;
    return %filter_hash;
}

############################################################################
sub same_class_barcode{
    my %filter_hash;
    my ($barcode,$sort_edge,$sort_nod) = @_;
    open INA,$barcode or die $!;        
    my $edge_number = scalar(@$sort_edge);
    my $nod_number = scalar(@$sort_nod);
    while(<INA>){
        chomp;
        my ($id,$this_edge_number,$ed,$this_nod_number,$no,$bili)=split/\s+/,$_;
        my @this_edge = split/\,/,$ed;
        my @this_no = split/\,/,$no;
        if ($this_edge_number == $edge_number + 1){
            if ($this_edge_number - $this_nod_number <= 0){
                if (intersect(@this_edge,@$sort_edge) == $edge_number){
                    if (intersect(@this_no,@$sort_nod) == $nod_number){
                        @{$filter_hash{$id}} = ($ed,$no,$bili);
                    }
                }
            }
        }
    }
    close INA;
    return %filter_hash;
}








