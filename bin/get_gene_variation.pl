#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#use Data::Printer;
use List::Util qw/max min/;

my ( $loxpregfile,$encodfile, $prefix, $outdir) = @ARGV;

die "\nUsage: perl $0 <loxpregion> <encod> <name> <outdir>\n\n" if (!$loxpregfile||!$encodfile);

open LXPR, "$loxpregfile" || die $!; # loxpreg
open ENC,  "$encodfile"   || die $!;

open OUTV,"> $outdir/$prefix.variation.xls" || die $!;
open OUTVCP,"> $outdir/$prefix.segmcp.vari.xls" || die $!;
#open OUTVGE, "> $prefix.gene.vari.xls" || die $!;

my %loxpreg;
my $lastregid=0;

while ( <LXPR> ) {
	chomp;
	my ($chrid, $star, $END, $indx) = (split /\s+/,$_)[0,3,4,8];
	$loxpreg {$indx} = [$chrid, $star, $END];

    $lastregid = $indx;
}
close LXPR;
#print Dumper %loxpreg;
my %loxp_halfsite;
my %halfsite2segm;
my $_initial_hs = 0;
my $_last_hs ;

foreach my $_segmid ( 1 .. $lastregid ) {
    my $_half_site_L = $_segmid * 2 - 2;
    my $_half_site_R = $_half_site_L + 1;
    $loxp_halfsite{$_segmid} = [ $_half_site_L, $_half_site_R ];

    $halfsite2segm {$_half_site_L} = [0, $_segmid];
    $halfsite2segm {$_half_site_R} = [1, $_segmid];
    $_last_hs = $_half_site_R;
}

my %samencods;
my @sample_ids;

while (<ENC>) {
    chomp;
    next if $_ =~ /^#/;
	my ($samid,$encod);
	my @this_set =split;
	if (@this_set == 2){
    	$samid = $this_set[0];
		$encod = $this_set[1];
	}else{
		$samid = $this_set[0];
		$encod = $this_set[2];
	}
    push @sample_ids, $samid;

    push @{$samencods{$samid}}, $encod;
}
close ENC;
#print Dumper %samencods;

my %sgemcp;

foreach my $_sample_id (@sample_ids) {
    my @_hs_paths;

    my %_cp_tmp;
    foreach my $indx ( 1..$lastregid ){
	   $_cp_tmp{$indx} = 0;
    }

    foreach my $_encod (@{$samencods{$_sample_id}}) {
        my @idx = split /,/, $_encod;
        my @_tmp_hs_path;
        foreach my $_segm_idx (@idx) {
            if ($_segm_idx > 0) {
                push @_tmp_hs_path, $loxp_halfsite{$_segm_idx};
                $_cp_tmp{$_segm_idx} ++;
            }else {
                my $_segm_idx_tmp = abs($_segm_idx);
                my @_hs_tmp = reverse ( @{$loxp_halfsite{$_segm_idx_tmp}});

                #print Dumper @_hs_tmp;
                push @_tmp_hs_path, \@_hs_tmp;

                $_cp_tmp{$_segm_idx_tmp} ++;
            }
        }
        push @_hs_paths, \@_tmp_hs_path;
    }
    $sgemcp{$_sample_id} = \%_cp_tmp;

    #print Dumper %_cp_tmp;
    #print Dumper @_hs_paths;
    # identify breakpoint hs...
    my @_hs_scbsite;
    my %novel_junction;

    my %_hs_path_direction;  ## check direction only on the condition of one copy ,
    foreach my $i (0..@_hs_paths-1) {
        my $_hs_path = shift @{$_hs_paths[$i]};
        my ($_initial_start, $_up_end) =  @{$_hs_path};
        #print "_initial_start, _up_end: $_initial_start, $_up_end\n";
        my $_up_direction = "+" ;
        $_up_direction = "-" if  $_initial_start > $_up_end;

        $_hs_path_direction{$_initial_start} = $_up_direction;
        $_hs_path_direction{$_up_end} = $_up_direction;

        foreach $_hs_path (@{$_hs_paths[$i]}) {
            my ($_hs_start,  $_hs_end) = @{$_hs_path};
            #print "_hs_start, _hs_end: $_hs_start, $_hs_end\n";
            my $_tmp_;
            if ( $_up_direction eq "+" ) {
                $_tmp_ = $_up_end + 1;
            }else {
                $_tmp_ = $_up_end - 1;
            }
            #print "$_tmp_==$_hs_start\n";

            # #######
            if ($_hs_start < $_hs_end) {
                $_up_direction = "+";
            }else { $_up_direction = "-"; }

            $_hs_path_direction{$_hs_start} = $_up_direction;
            $_hs_path_direction{$_hs_end}   = $_up_direction;

            if ($_hs_start == $_tmp_) {
                $_up_end = $_hs_end;
                next;
            }else {
                push @_hs_scbsite, ($_up_end, $_hs_start);

                my @_sortpoint = sort {$a <=> $b} ($_up_end, $_hs_start);
                push @{$novel_junction{$_up_end}},   \@_sortpoint;
                push @{$novel_junction{$_hs_start}}, \@_sortpoint;
                $_up_end = $_hs_end;
            }
        }

        unless ( $_up_end == $_last_hs || $_up_end == $_initial_hs || $_initial_start == $_initial_hs){
            push @_hs_scbsite, ($_up_end, $_initial_start);

            my @_sortpoint = sort {$a <=> $b} ($_up_end, $_initial_start);
            push @{$novel_junction{$_up_end}},        \@_sortpoint;
            push @{$novel_junction{$_initial_start}}, \@_sortpoint;
        }
    }

    #print Dumper %_hs_path_direction;
    #print Dumper  @_hs_scbsite;
    #print Dumper  %novel_junction;

    if (@_hs_scbsite == 0){
        print "$_sample_id does NOT have any variation.\n";
        next;
    }

    ## fragment _hs_path ...
    my $_up_start = 0;
    my @new_hs_segm;

    my %_uniq;
    foreach my $_hs ( sort {$a <=> $b} @_hs_scbsite ){
         if ( $halfsite2segm {$_hs} -> [0] == 1 ){

             next if $_up_start >= $_hs ;
             push @new_hs_segm, [$_up_start, $_hs];
             $_uniq{$_up_start} = 1;
             $_up_start = $_hs + 1;
         }else{
             my $_tmp_hs = $_hs - 1;
             next if $_up_start >= $_tmp_hs;

             push @new_hs_segm, [$_up_start, $_tmp_hs];
             $_up_start = $_hs;
         }
    }

    push @new_hs_segm, [$_up_start,$_last_hs];

    #print Dumper @new_hs_segm;
    ## justify the variation ...
	#p(@new_hs_segm);exit;
    foreach my $_hs_segm (@new_hs_segm) {
       
	   	my ($_half_site_L, $_half_site_R) = @{$_hs_segm};
		#p($_half_site_L);p($_half_site_R);exit;
        ## verify copy_number ..
		#p(@{$halfsite2segm{$_half_site_L}});
		#p(@{$halfsite2segm{$_half_site_R}});exit;
        my ($_ix_L, $_segm_start) = @{$halfsite2segm {$_half_site_L}};
        my ($_ix_R, $_segm_end) = @{$halfsite2segm {$_half_site_R}};

        my $loxpsegm_encod = join (",", $_segm_start..$_segm_end);
        my $_hs_segm_cp = $sgemcp{$_sample_id}{$_segm_start};
        my $_ck_check = 1;
        foreach my $_ix ( $_segm_start .. $_segm_end ){
            if ($_hs_segm_cp ==  $sgemcp{$_sample_id}{$_ix}){
                next;
            }else{
                $_ck_check = 0;
                print "The fragmentation of breakpoints make an error !!\n";
            }
        }

        ## extract variation ...
        my ($varitype, $varijunct, $complex );
        my @junction_infos;
        if ($_ck_check == 1) {  # copy number check passed ...

            if (exists $novel_junction{$_half_site_L}){
                push @junction_infos, @{$novel_junction{$_half_site_L}};
            }

            if (exists $novel_junction{$_half_site_R}){
                push @junction_infos, @{$novel_junction{$_half_site_R}};
            }

            if ( $_hs_segm_cp == 0 ){
                $varitype = "DEL";

                if ( @junction_infos == 0 ) {
                    my ($_tmp_L, $_tmp_R);
                    if ($_half_site_L == $_initial_hs) {
                        $_tmp_L = $_half_site_R + 1;
                        push @junction_infos, @{$novel_junction{$_tmp_L}};
                        $complex = 1;

                    } elsif ($_half_site_R == $_last_hs ){
                        $_tmp_R = $_half_site_L - 1;
                        push @junction_infos, @{$novel_junction{$_tmp_R}};
                        $complex = 1;

                    } else {
                        $_tmp_L = $_half_site_R + 1;
                        $_tmp_R = $_half_site_L - 1;
                        push @junction_infos, @{$novel_junction{$_tmp_L}};
                        push @junction_infos, @{$novel_junction{$_tmp_R}};

                        if ($novel_junction{$_tmp_L}->[0]->[0] == $_tmp_R && $novel_junction{$_tmp_L}->[0]->[1] == $_tmp_L){
                            $complex = 0;
                        }else{
                            $complex = 1;
                        }
                    }
                }
            }elsif( $_hs_segm_cp > 1 ) {
                $varitype = "DUP";
                # $complex ? details Inverted-Tandem/Tandem DUP or others ..
                $complex = 1;

            }elsif( $_hs_segm_cp == 1 ){
                # need direction to check the variations ...
                #print "_half_site_L: $_half_site_L\n" unless exists $_hs_path_direction{$_half_site_L};

                if ( $_hs_path_direction{$_half_site_L} eq "+" ){
                    $varitype = "NA";

                }else{
                    $varitype = "INV";
                    if (@junction_infos == 2) {
                        #print Dumper @junction_infos;

                        my ( $_link_L1, $_link_R1 ) = @{$junction_infos[0]};
                        my ( $_link_L2, $_link_R2 ) = @{$junction_infos[1]};

                        if ( abs($_link_L1 - $_link_L2) == 1 && abs($_link_R1 - $_link_R2) == 1 ){
                            $complex = 0;
                        }else{
                            $complex = 1;
                        }
                    }elsif (@junction_infos == 0){
                        $varitype = "NA";

                    } else {
                        $complex = 1;

                    }
                }
                #only junction related to gene.
            }
            ###########
            #output variation ..
        }else{ next; }

        if  ($varitype eq "NA"){

            foreach my $_segm_id ( $_segm_start .. $_segm_end ) {
                print OUTVCP "$_sample_id\t$loxpreg{$_segm_start}->[0]\t$loxpreg{$_segm_start}->[1]\t$loxpreg{$_segm_end}->[2]\t$varitype\t" . "-\t";
                print OUTVCP "-" . "\t" . "-" . "\t$_segm_id\t$sgemcp{$_sample_id}{$_segm_id}\n";
            }
            next;
        }

        my %_uniq_junction;
        my @_uniq_varijunct;
        my @segm_junct;
        #print Dumper @junction_infos;

        foreach my $_junct (@junction_infos) {
            my ($_link_L, $_link_R) = @{$_junct};

            my $_junction = "$_link_L/$_link_R";

            unless ($_uniq_junction{$_junction}) {
                push @_uniq_varijunct, $_junction;
                $_uniq_junction{$_junction} = 1;
            }
        }
        my @_sort_uniq_varijunct = sortjunct ( @_uniq_varijunct );

        foreach my $_varijunct ( @_sort_uniq_varijunct ) {
            my ($_link_L, $_link_R) =  split /\//, $_varijunct;
            my ($_side_L, $_segmid_L) = @{$halfsite2segm{$_link_L}};
            my ($_side_R, $_segmid_R) = @{$halfsite2segm{$_link_R}};

            my $_segm_junct;
            if ($_side_L == 0) {$_segm_junct = "-$_segmid_L"; } else { $_segm_junct = "$_segmid_L";   }
            if ($_side_R == 0) {$_segm_junct .= "/$_segmid_R";} else { $_segm_junct .= "/-$_segmid_R";}

            push @segm_junct, $_segm_junct;
        }

        print OUTV "$_sample_id\t$loxpreg{$_segm_start}->[0]\t$loxpreg{$_segm_start}->[1]\t$loxpreg{$_segm_end}->[2]\t$varitype\t"."$complex\t";
        print OUTV join (",", @_uniq_varijunct) . "\t" . join (",",@segm_junct) ."\t$loxpsegm_encod\n";

        foreach my $_segm_id ( $_segm_start .. $_segm_end ) {
            print OUTVCP "$_sample_id\t$loxpreg{$_segm_start}->[0]\t$loxpreg{$_segm_start}->[1]\t$loxpreg{$_segm_end}->[2]\t$varitype\t" . "$complex\t";
            print OUTVCP join (",", @_uniq_varijunct) . "\t" . join (",",@segm_junct) . "\t$_segm_id\t$sgemcp{$_sample_id}{$_segm_id}\n";
        }
    }
}

sub sortjunct{
    my @_in_junct  = @_;
    my @_new_juct;
    my %_junct;
    foreach my $junct_info (@_in_junct){
      	my @_index = split /\//,$junct_info;

        foreach  my $i ( 0..@_index-1 ) {
           $_index[$i] = -1 * $_index[$i] if  $_index[$i] < 0;
        }

        my $_junct_min_index = min (@_index);
        my $_junct_sum_index = 0;
        foreach (@_index){ $_junct_sum_index += $_; }

        $_junct{$_junct_min_index}{$_junct_sum_index} = $junct_info;
    }

  foreach my $_index (sort{$a<=>$b} keys %_junct){
       foreach my $_sum_index (sort{$a<=>$b} keys %{$_junct{$_index}}){
       push  @_new_juct,$_junct{$_index}{$_sum_index};
   }
 }
  return @_new_juct;
}

=head
my %genidx;
my @genlst;

while(<IN2>) { # GFF ...
	chomp;
	my ($star,$END,$ori,$genescript) = ( split /\s+/,$_ )[3,4,6,8];

    my @_info =  split /;/, $genescript;
    my $gename = (split /Parent=/, $_info[2])[1];
    $gename =~ s/_mRNA//g;

    push @genlst,$gename;
    my $ck = 0;
    foreach (keys %loxpreg) {
        if ($loxpreg{$_}[0] <= $star && $loxpreg{$_}[1] >= $END) {
            $genidx{$gename} = [ $_, 0 ];
            $ck = 1;
            last;
        }
    }

    if ($ck ==0) {print "Not Found $gename in any loxp region\n";}
}



close IN1;
close IN2;
close IN3;
close OUT;
close OUTS;
=cut

