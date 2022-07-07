#!/usr/bin/perl -w

use strict;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../perlib";
use Getopt::Long;
use Pod::Usage;
use Data::Printer;

my $usage=<<USAGE;

    Description: This script is used to restructure the synthethic yeast SCRaMbLE genome by stLFR technology !
                 v1.0 2022.06.12 wendingpang
    Usage: perl $0 [-option]

            --fa            [file]  : <fasta> the unSCRaMbLE genome FASTA [bwa index must!]
            --fq1           [file]  : <fastq1> the stLFR FASTQ1
            --fq2           [file]  : <fastq2> the stLFR FASTQ2 
            --lxpreg        [file]  : <loxp_region> the loxp region file        
            --blist         [file]  : <barcode_list> the barcode_list from stLFR-seq [default]
            --rfcvg         [file]  : <refcoverage> the seqence of unSCRaMbLE genome's depth

            --n             : the output median file and directory name  [defacult: test]
            --o             : the output file directory and the median-file directory [defacult: ./]
            --step          : [all|assembly|others] which steps you want to run Sgrstlfr.
                                    all:run the all [defacult]
                                    assembly: only run the module7 to repair the scaffold and genome assembly
                                    others: you appoint which step to run (likes 4,5,6,7) 
            <step7-2>
            --min_edge      : :the minnum of the barcode include edge [defacult:3]
            --max_edge      : the maxnum of the barcode include edge [defacult:7] 
            --p_nod         : the left and right nod number from diff nod start to stat bs [defacult:1]  
            --solu_way      : [max | plus ] the solution use barcode support way (max: get the max bs path or plus:get the bs more than 0 path) [defacult:max]
            
            <step7-1>
            --start_nod     : the scaffold assembly start type.(defacult auto find uniq edge as start or appoint mutipe nod) [default:undef]
            --winfo         : [F | T ]in the scaffold assembly if produce the all complete edge-nod link set [defacult:F]
            --s_len         : the contig length of start to scaffold assembly [defacult:5000]
            --s_way         : [only | all] only:only produce all scaffold set not include filter genome path all:produce all result [defacult:all]
            --fix_dup       : if from the info file find the duplication type is erro in the scaffold assembly. 
                              you can direct repair in this paramter likes (0_1,T,T;2_3,F,T) [defacult:undef]
            --bar_edge      : the process of assembly scaffold you use barcode max edge number (N+this paramter) [defaclut:1] 
            --s_bs          : the barcode_support of sub_path (more than Zero) [defacult:3]
            
            <step6>
            --overlap       : the alon-reads coverage one edge precentage [deafacult:0.8]
            --isz           : the PE seqence average insert_size [defacult:500]
            --f_bar         : the number of every edge include reads. [defacult:4].
            
            <step5>
            --cutlen        : cut each region into subregion for Mahalanobis Distance(MD) analysis [default: 500] #500 model only used
            --mcycle        : the maximal cycle for estimate copy number [default: 10] 
            --mdep          : the maximal average sequencing deth of copy number = 0 [default: 10] 
            --minread       : the minimal split reads for supporting a breakpoint [default: 2]

            <step2>
            --s_l           : the low quality threshold [defacult:5]
            --s_q           : the low quality rate [defacult:0.5]
            --s_n           : the N rate threshold [defacult:0.05]
            --Q           : the quality system 1:64 2:33 [defacult:2]
            
            <basic inforamtion>
            --chrtype      : [liner | cycle ]the chrmosome type (liner or cycle) [defacult: cycle] 
            --chrid        : the chrid of want to restructure [defacult:IXR_BACseq]
            --t             : the threads about this script run [defacult: 4]
            
            <tools>
            --st            : the tools path of samtools [defacult]
            -soke           : the tools path of soapnuke [defacult]
            -bwa            : the tools path of bwa [defacult]

            --h|-help       : display this help 

USAGE


# pre-defined
my ($fasta,$fastq1,$fastq2,$region,$barcode_list,$refcoverage);
my ($filename,$outdir,$step,$threads,$chrid,$chr_type,$help);
my ($minread,$mdep,$mcycle,$cutlen);
my ($insert_size,$overlap,$filter_barcode);
my ($start_nod,$scaffold_way,$winfo,$start_len,$fix_dup,$bar_edge_num,$f_barcode_num);
my ($the_min_barcode_edge,$the_max_barcode_edge,$bs_edge_num,$solution_way);
my ($s_l,$s_q,$s_n,$s_Q);
my ($samtools,$bwa,$soapnuke);

if ( !GetOptions(
    "fa:s" =>\$fasta,
    "fq1:s" =>\$fastq1,
    "fq2:s" =>\$fastq2,
    "lxpreg:s" =>\$region,
    "blist:s" =>\$barcode_list,
    "rfcvg:s" =>\$refcoverage,
    "n:s" =>\$filename,
    "o:s" =>\$outdir,
    "step:s" =>\$step,
    "t:s" =>\$threads,
    "chrid:s" =>\$chrid,
    "chrtype:s" =>\$chr_type,
    "minread:s" => \$minread,
    "mdep:s" => \$mdep,
    "mcycle:s" => \$mcycle,
    "cutlen:s" => \$cutlen,
    "isz:s" =>\$insert_size,
    "overlap:s" =>\$overlap,
    "f_bar:s" =>\$filter_barcode,
    "start_nod:s" =>\$start_nod,
    "s_way:s" =>\$scaffold_way,
    "winfo:s" =>\$winfo,
    "s_len:s" =>\$start_len,
    "fix_dup:s" =>\$fix_dup,
    "bar_edge:s" =>\$bar_edge_num,
    "s_bs:s" =>\$f_barcode_num,
    "min_edge:s" =>\$the_min_barcode_edge,
    "max_edge:s" =>\$the_max_barcode_edge,
    "p_nod:s" =>\$bs_edge_num,
    "solu_way:s" =>\$solution_way,
    "s_l:s" => \$s_l,
    "s_q:s" => \$s_q,
    "s_n:s" => \$s_n,
    "Q:s" => \$s_Q,
    "st:s" =>\$samtools,
    "soke:s" =>\$soapnuke,
    "bwa:s" => \$bwa,
    "help|h!" => \$help
)){
    pod2usage ({
        -message => "!Failed to parse commmand line.\n$usage\n",
        -verbose => 1,
        -exitval =>1
    });
}

if ($help){
    pod2usage({
        -verbose =>0,
        -exitval =>0,
        -message => "$usage\n"
    });
}

# basic information
$filename ||='test';
$outdir ||=getcwd;
$threads ||=4;
$chrid ||='IXR_BACseq';
$chr_type ||='cycle';
$barcode_list ||="$Bin/stLFR_barcode_split/barcode_list.txt";
# step 5
$minread ||=2;
$mdep ||=10;
$mcycle ||=10;
$cutlen ||=500;
# step 6
$insert_size ||=500;
$overlap ||=0.8;
$filter_barcode ||=4;
# step 7-1
$start_nod ||='undef';
$scaffold_way ||='all';
$winfo ||='F';
$start_len ||=5000;
$fix_dup ||='';
$bar_edge_num ||=1;
$f_barcode_num ||=3;
# step7-2
$the_min_barcode_edge ||=3;
$the_max_barcode_edge ||=7;
$bs_edge_num ||=1;
$solution_way ||='max';
# step2
$s_l ||=5;
$s_q ||=0.5;
$s_n ||=0.05;
$s_Q ||=2; 
#tools path 
$samtools ||= "/hwfssz1/ST_BIOCHEM/P18Z10200N0255/PMO/F14ZQSYJSY1726/pangwending/some_tools/samtools-1.13/bin/samtools";
$bwa ||="/share/app/bwa/0.7.17/bwa";
$soapnuke ||= "/share/app/soapnuke/2.1.0/SOAPnuke";
# step information
my ($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8);
$step ||='all';

print '***** Start run the [Sgrstlfr] *****'."\n";
print localtime()."\n";
my $start_time = time();
my @a_step_set = qw/0 0 0 0 0 0 0 0/;
($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8) = @a_step_set;
if ($step eq 'all'){
    ($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8) = (1,1,1,1,1,1,1,1);
    print '##### Start all the step!,loading.........'."\n";
}elsif ($step eq 'assembly'){
    $step7 = 1;
    print 'Only start the step7!, loading......'."\n";
}else{
    my @step_set = split/\,/,$step;
    print '##### Start the step:'.join(",",@step_set)."\t".'loading.....'."\n";
    foreach my $anum (@step_set){
        if ($anum =~ /[1-8]/){
            splice(@a_step_set,$anum-1,1,1);
        }else{
            print '!!!The parameter -step extists not a number in (1-8) ! '.$anum."\n";
            exit;
        }

    }
    ($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8) = @a_step_set;
    # p(@a_step_set);exit;
}


# script path
my $split_barcode = "$Bin/bin/split_barcode.pl";
my $class_mapping_result = "$Bin/bin/class_mapping_result.pl";
my $get_unmap_reads_split = "$Bin/bin/get_unmap_reads_split.pl";
my $getvariareads = "$Bin/bin/getvariareads.pl";
my $extnodedge = "$Bin/bin/extnodedge.pl";
my $get_barcode_split_node = "$Bin/bin/get_barcode_split_node.pl";
my $get_barcode_pe_node = "$Bin/bin/get_barcode_pe_node.pl";
my $mersh_basic_file = "$Bin/bin/mersh_basic_file.pl";
my $mersh_check_node = "$Bin/bin/mersh_check_node.pl";
my $stat_barcode_information = "$Bin/bin/stat_barcode_information.pl";
my $arrange_barcode_information = "$Bin/bin/arrange_barcode_information.pl";
my $exact_link = "$Bin/bin/exact_link.pl";
my $global_filter = "$Bin/bin/global_filter.pl";
my $translate_path = "$Bin/bin/translate_path.pl";

# check file and dependencies
die "The script split_barcode.pl is not found in $split_barcode!\n" unless -s $split_barcode;
die "The script class_mapping_result.pl is not found in $class_mapping_result!\n" unless -s $class_mapping_result;
die "The script get_unmap_reads_split.pl is not found in $get_unmap_reads_split!\n" unless -s $get_unmap_reads_split;
die "The script getvariareads.pl is not found in $getvariareads!\n" unless -s $getvariareads;
die "The script extnodedge.pl is not found in $extnodedge!\n" unless -s $extnodedge;
die "The script get_barcode_split_node.pl is not found in $get_barcode_split_node!\n" unless -s $get_barcode_split_node;
die "The script get_barcode_pe_node.pl is not found in $get_barcode_pe_node!\n" unless -s $get_barcode_pe_node;
die "The script mersh_basic_file.pl is not found in $mersh_basic_file!\n" unless -s $mersh_basic_file;
die "The script mersh_check_node.pl is not found in $mersh_check_node!\n" unless -s $mersh_check_node;
die "The script stat_barcode_information.pl is not found in $stat_barcode_information!\n" unless -s $stat_barcode_information;
die "The script arrange_barcode_information.pl is not found in $arrange_barcode_information!\n" unless -s $arrange_barcode_information;
die "The script exact_link.pl is not found in $exact_link!\n" unless -s $exact_link;
die "The script global_filter.pl is not found in $global_filter!\n" unless -s $global_filter;
die "The script translate_path.pl is not found in $translate_path!\n" unless -s $translate_path;

# mkdir father-directory
system ("mkdir -p $outdir/$filename");
my $file_path = "$outdir/$filename";
#### module 1 split barcode 
if ($step1 == 1){
    die "ERROR: The stLFR-seq file fastq1 not found in $fastq1!\n" unless -s $fastq1;
    die "ERROR: The stLFR-seq file fastq2 not found in $fastq2!\n" unless -s $fastq2;
    print '!Start to run the step1,loading......'."\n";
    `perl $split_barcode $barcode_list $fastq1 $fastq2 split_reads $file_path > $file_path/split_barcode.log 2>$file_path/split_barcode.err`;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step1 split_barcode is not finished! Pleas check the running commond!"."\n";
    }
    print '#The step1 done!'."\n";
}

#### module 2 clean data 
if ($step2 == 1){
    die "ERROR: The split_barcode results split_reads.1.fq.gz not found in $file_path/split_reads.1.fq.gz!\n " unless -s "$file_path/split_reads.1.fq.gz";
    die "ERROR: The split_barcode results split_reads.2.fq.gz not found in $file_path/split_reads.2.fq.gz!\n" unless -s "$file_path/split_reads.2.fq.gz";
    die "The tools of soapnuke path is not found in $soapnuke!\n" unless -s $soapnuke;
    print '!Start to run the step2,loading......'."\n";
    `$soapnuke filter -l 5 -q 0.1 -n 0.01 -1 $file_path/split_reads.1.fq.gz -2 $file_path/split_reads.2.fq.gz -C split_clean_reads.1.fq.gz -D split_clean_reads.2.fq.gz -o $file_path `;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step2 clean_data is not finished! Pleas check the following run commond!"."\n";
        print "$soapnuke filter -l $s_l -q $s_q -n $s_n -Q $s_Q -1 $file_path/split_reads.1.fq.gz -2 $file_path/split_reads.2.fq.gz -C split_clean_reads.1.fq.gz -D split_clean_reads.2.fq.gz -o $file_path\n";
        exit;
    }
    print "#The step2 done!\n";
}

#### module 3 alignment to unSCRaMbLE genome
if ($step3 == 1){
    die "ERRO: The split and clean results split_clean_reads.1.fq.gz not found in $file_path/split_clean_reads.1.fq.gz!\n " unless -s "$file_path/split_clean_reads.1.fq.gz";
    die "ERRO: The split and clean results split_clean_reads.2.fq.gz not found in $file_path/split_clean_reads.2.fq.gz!\n" unless -s "$file_path/split_clean_reads.2.fq.gz";
    die "ERRO: The reference of fasta file is not found in $fasta!\n" unless -s $fasta;
    die "The tools of bwa path is not found in $bwa!\n" unless -s $bwa;
    die "The tools of samtools path is not found in $samtools!\n" unless -s $samtools;
    print '!Start to run the step3,loading......'."\n";
    `$bwa mem -t $threads $fasta $file_path/split_clean_reads.1.fq.gz $file_path/split_clean_reads.2.fq.gz |$samtools sort -> $file_path/$filename.sort.bam `;
    `$samtools index $file_path/$filename.sort.bam `;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step3 alignment is not finished! Pleas check the following run commond!"."\n";
        print "$bwa mem -t $threads $fasta $file_path/split_clean_reads.1.fq.gz $file_path/split_clean_reads.2.fq.gz |$samtools sort -> $file_path/$filename.sort.bam\n";
        exit;
    }
    print "#The step3 done!\n";
    
}

#### module 4 stat the mapping result and split the unmapping reads include loxp
if ($step4 == 1){
    die "ERROR: The alignment results $filename.sort.bam not found in $file_path/$filename.sort.bam!\n" unless -s "$file_path/$filename.sort.bam";
    die "ERROR: The unSCRaMbLE genome fasta file not found in $fasta!\n" unless -s $fasta;
    die "The tools of bwa path is not found in $bwa!\n" unless -s $bwa;
    die "The tools of samtools path is not found in $samtools!\n" unless -s $samtools;
    print '!Start to run the step4,loading......'."\n";
    `$samtools depth -a $file_path/$filename.sort.bam | grep "$chrid" > $file_path/$filename.depth `;
    `perl $class_mapping_result $file_path/$filename.sort.bam -n $filename -chr $chrid -s $samtools -o $file_path >$file_path/stat.bam.log `;
    `$bwa mem $fasta $file_path/$filename.split.fq > $file_path/$filename.split.sam `;
    `perl $getvariareads $file_path/$filename.split.sam -n $filename -s $samtools -o $file_path `;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the alignment is not finished! Pleas check the median file adn running commond !"."\n";
        exit;
    }
    print "#Step 4 done!"."\n";
}

#### module 5 get the edge and node 
if ($step5 == 1){
    die "ERROR: The step4 produce  file $filename.gvr not found in $file_path/$filename.gvr!\n" unless -s "$file_path/$filename.gvr";
    die "ERROR: The step4 produce coverage file $filename.depth not found in $file_path/$filename.depth!\n" unless -s "$file_path/$filename.depth";
    die "ERROR: The alignment results $filename.sort.bam not found in $file_path/$filename.sort.bam!\n" unless -s "$file_path/$filename.sort.bam";
    die "ERROR: The reference of unSCRaMbLE genome coverage file not found in $refcoverage!\n" unless -s $refcoverage;
    die "ERROR: The loxpsym region file not found in $region!\n" unless -s $region; 
    print '!Start to run the step5,loading......'."\n";
    `perl $extnodedge $file_path/$filename.depth $file_path/$filename.sort.bam $region $file_path/$filename.gvr -n $filename -chrid $chrid -s $samtools -refcover $refcoverage -minread $minread -mdep $mdep -mcycle $mcycle -cutlen $cutlen -o $file_path >extnodedge.log`;
    # print "perl $extnodedge $file_path/$filename.depth $file_path/$filename.sort.bam $region $file_path/$filename.gvr -n $filename -chrid $chrid -s $samtools -refcover $refcoverage -minread $minread -mdep $mdep -mcycle $mcycle -cutlen $cutlen -o $file_path\n";
    print "The erro number:$?\n";  #6400 edge canclute erro repair the mdep  
      if ($? != 0){
        print "ERROR: the step5:get the edge and node is not finished! Pleas check the extnodedge.log and running commmond !!!"."\n";
        exit;
    }
    print "#Step 5 done!"."\n";
}

#### module 6 prepare the basic co-barcode information and filter the barcode
if ($step6 == 1){
    die "ERROR: The step5 produce $filename.edg file not found in $file_path/$filename.edg!\n" unless -s "$file_path/$filename.edg";
    die "ERROR: The step5 produce $filename.nod file not found in $file_path/$filename.nod!\n" unless -s "$file_path/$filename.nod";
    die "ERROR: The step4 produce $filename.depth file not found in $file_path/$filename.depth!\n" unless -s "$file_path/$filename.depth";
    die "ERROR: The step4 produce $filename.synmap.gz not found in $file_path/$filename.synmap.gz!\n" unless -s "$file_path/$filename.synmap.gz";
    die "ERROR: The step4 produce $filename.split.sam not found in $file_path/$filename.split.sam!\n" unless -s "$file_path/$filename.split.sam";
    die "ERROR: The loxpsym region file not found in $region!\n" unless -s $region;
    print '!Start to run the step6,loading......'."\n";
    `perl $mersh_basic_file $file_path/$filename.edg $region $file_path/$filename.depth $chrid $filename $file_path `;
    die "ERROR: The step6 produce $filename.total file not found in $file_path/$filename.total!\n" unless -s "$file_path/$filename.total";
    `perl $get_barcode_pe_node $file_path/$filename.edg $region $file_path/$filename.synmap.gz -n $filename -isz $insert_size -ct $chr_type -overlap $overlap -o $file_path `;
    die "ERROR: The step6 produce $filename.pe.nod not found in $file_path/$filename.pe.nod!\n" unless -s "$file_path/$filename.pe.nod";
    `perl $get_barcode_split_node $file_path/$filename.edg $region $file_path/$filename.split.sam -s $samtools -t $chr_type -n $filename -o $file_path `;
    die "ERROR: The step6 produce $filename.split.nod not found in $file_path/$filename.split.nod!\n" unless -s "$file_path/$filename.split.nod";
    `perl $mersh_check_node $file_path/$filename.nod $file_path/$filename.pe.nod $file_path/$filename.split.nod $filename $file_path `;
    die "ERROR: The step6 produce $filename.mersh.nod not found in $file_path/$filename.mersh.nod!\n" unless -s "$file_path/$filename.mersh.nod";
    `perl $stat_barcode_information $file_path/$filename.edg $region $file_path/$filename.synmap.gz $file_path/$filename.mersh.nod -n $filename -overlap $overlap -f $filter_barcode -o $file_path `;
    `perl $arrange_barcode_information $file_path/$filename.node.barcode $filename $file_path `;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step6:get the pre-file is not finished! Pleas check the median file and run commond!!!"."\n";
        exit;
    }
    print "#Step 6 done!"."\n";
}

#### module 7 get the scaffold result and find final exact genemo result! 
if ($step7 == 1){
    die "ERROR: The step5 produce $filename.edg file not found in $file_path/$filename.edg!\n" unless -s "$file_path/$filename.edg";
    die "ERROR: The step5 produce $filename.nod file not found in $file_path/$filename.nod!\n" unless -s "$file_path/$filename.nod";
    die "ERROR: The step6 produce barcode file $filename.sort.barcode.stat not found in $file_path/$filename.sort.barcode.stat!\n" unless -s "$file_path/$filename.sort.barcode.stat";
    die "ERROR: The step6 produce basic file total file $filename.total not found in $file_path/$filename.total!\n" unless -s "$file_path/$filename.total";
    print '!Start to run the step7,loading......'."\n";
    `perl $exact_link $file_path/$filename.sort.barcode.stat $file_path/$filename.nod $file_path/$filename.edg $file_path/$filename.total -n $filename -o $file_path -chr_type $chr_type -b $bar_edge_num -fix $fix_dup -fp $f_barcode_num -i $winfo -len $start_len -start_nod $start_nod -way $scaffold_way > $file_path/$filename.scaf.log `;
    `perl $global_filter $file_path/$filename.result $file_path/$filename.sort.barcode.stat -n $filename -o $file_path -min $the_min_barcode_edge -max $the_max_barcode_edge -wtnod $bs_edge_num -way $solution_way >$file_path/$filename.solu.log`;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step7:get the scaffold and exact genome is not finished! Pleas check the scaffold log!!!"."\n";
        exit;
    }
    print "#Step 7 done!"."\n";
}

#### module 8 transfer the edge-nod path genome to fasta!
if ($step8 == 1){
    die "ERROR: The step7 produce $filename.final not found in $file_path/$filename.final!\n" unless -s "$file_path/$filename.final";
    die "ERROR: The loxpsym region file not found in $region!\n" unless -s $region;
    die "ERROR: The step5 produce $filename.edg file not found in $file_path/$filename.edg!\n" unless -s "$file_path/$filename.edg";
    die "ERROR: The unSCRaMbLE fasta file is not found in $fasta!\n" unless -s $fasta;
    print '!Start to run the step8,loading......'."\n";
    `perl $translate_path $file_path/$filename.final $region $file_path/$filename.edg $fasta -c $chrid -n $filename -o $outdir `;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step8:tranfer the genome to fasta is not finished! Pleas check the following run commond!"."\n";
        print "perl $translate_path $file_path/$filename.final $region $file_path/$filename.edg $fasta -c $chrid -n $filename -o $outdir\n";
        exit;
    }
    print "#Step 8 done!"."\n";
}
my $end_time = time();
my $run_time = $end_time - $start_time;
print $run_time."\n";
print '############### All work is finished!!! ##################'."\n";
print localtime()."\n";






