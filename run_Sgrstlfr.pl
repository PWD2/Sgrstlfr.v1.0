#!/usr/bin/perl -w

use strict;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../perlib";
use List::Util qw/pairs/;
use Getopt::Long;
use Pod::Usage;

my $usage=<<USAGE;

    Description: This script is used to reconstruct the synthethic yeast SCRaMbLE genome by stLFR technology !
                 v1.0 2023.03.24 pangwending\@genomics.cn wangyun\@genomics.cn
    Usage: perl $0 [-option] <h|help>

            --fa            [file]  : <fasta> the unSCRaMbLE genome FASTA [alignment index must!]
            --fq1           [file]  : <fastq1> the stLFR FASTQ1 [inculde:raw_data、split_data、split_clean_data]
            --fq2           [file]  : <fastq2> the stLFR FASTQ2 [inculde:raw_data、split_data、split_clean_data]       
            --rfcvg         [file]  : <refcoverage> the seqence of unSCRaMbLE genome's single site depth
            <basic information [must]>
            --chrtype       : [liner | cycle ] the chrmosome type [default: cycle]
            --chrid         : the chrid what want to restructure [default: IXR_BACseq]
            --t             : the threads about this script run [default: 4]
            --n             : the output median file and directory name  [default: test]
            --o             : the output file directory and the median-file directory [default: ./]
            --step          : [all|assembly|others] which steps you want to run Sgrstlfr.
                                    all:run the all [default]
                                    assembly: only run the module7 to repair the scaffold and genome assembly
                                    others: you can appoint which step will run (likes 4,5,6,7)
          ----------------------------------------------------------------------------------------------------------------------------------- 
            <step7>
            --min_edge      : the minnum of the barcode include edge [default:3]
            --max_edge      : the maxnum of the barcode include edge [default:12] 
            --sp_minus      : [number | max] the barcode support decline precent of two solution [default:max]  
            --start_nod     : the scaffold assembly start type.(default auto find uniq edge as start or appoint mutipe nod) [default:undef]
            --bias          : allow the erro number in the scaffold set [default:1]
            --s_len         : the contig length of start to scaffold assembly [default:5000]
            --s_way         : [only | all] only:only produce all scaffold set not include filter genome path all:produce all result [default:all]
            --fix_dup       : if from the info file find the duplication type is erro in the scaffold assembly. 
                              you can direct repair in this paramter likes (0_1,T,T;2_3,F,T) [default:undef]
            --bar_edge      : the process of assembly scaffold you use barcode max edge number (N+this paramter) [defaclut:1] 
            --s_bs          : [auto | others ] the barcode_support of sub_path (more than Zero) [default:auto]
            --r_type        : [all | first ]the allow of the produce fasta file number [default:first]
            --c             : the solution confict allow [default:0] 
            <step6>
            --overlap       : the alon-reads coverage one edge precentage [deafacult:0.8]
            --isz           : the PE seqence average insert_size [default:500]
            --f_bar         : the number of every edge include reads. [default:4].
            --m_q           : the mapping quality of all alignment [deafacult:50]
            <step5>
            --cutlen        : cut each region into subregion for Mahalanobis Distance(MD) analysis [default: 500] #500 model only used
            --mcycle        : the maximal cycle for estimate copy number [default: 10] 
            --mdep          : the maximal average sequencing deth of copy number = 0 [default: auto] 
            --minread       : the minimal split reads for supporting a breakpoint [default: 5]
            <step4 svg graph>
            -p_svg          : if produce the svg graph of the split-mapping result [default:T]
            -score          : score pattern calculating cumulative score,match score:mismath penalty [1:2]
            -mlen           : the min length of match and unmatch [default:20]
            -d              : the min mapping distance or recongnized as recombination [default:30]
            -t_rpe          : threshold to identify as a breakpoint of rpe [default: 2]
            -sd_rpe         : maximal breakpoint coordinate difference allowed of rpe [default:5]
            -t_rse          : threshold to identify as a breakpoint of rse [default:1]
            -sd_rse         : maximal breakpoint coordinate difference allowed rse [defaclut:20]
            -wind           : depth graph length of a window [default: 10]
            -rdsn		    : minimal read number to surport a breakpoint [default 2]
            <step2>
            --s_l           : the low quality threshold [default:5]
            --s_q           : the low quality rate [default:0.1]
            --s_n           : the N rate threshold [default:0.01]
            --f_other       : the SOAPnuke other option [likes: Q,2,f,NNN,r,NNN ..]
            <tools>
            --st            : the tools path of samtools [default]
            -soke           : the tools path of soapnuke [default]
            -bwa            : the tools path of bwa [default]
            -botie2         : the tools path of bowtie2 [default]
            --h|-help       : display this help 

USAGE


# pre-defined
my ($fasta,$fastq1,$fastq2,$refcoverage);
my ($filename,$outdir,$step,$threads,$chrid,$chr_type,$m_q,$help);
my ($minread,$mdep,$mcycle,$cutlen);
my ($insert_size,$overlap,$filter_barcode);
my ($start_nod,$scaffold_way,$bias,$start_len,$fix_dup,$bar_edge_num,$f_barcode_num);
my ($the_min_barcode_edge,$the_max_barcode_edge,$sp_minus,$rtype,$confict);
my ($score,$alen,$dt,$tp,$sdp,$ts,$sds,$wind,$rdsn,$p_svg);
my ($s_l,$s_q,$s_n,$f_other);
my ($samtools,$bwa,$soapnuke,$bamdeal,$bowtie2);

if ( !GetOptions(
    "fa:s" =>\$fasta,
    "fq1:s" =>\$fastq1,
    "fq2:s" =>\$fastq2,
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
    "m_q:s" =>\$m_q,
    "start_nod:s" =>\$start_nod,
    "s_way:s" =>\$scaffold_way,
    "bias:s" =>\$bias,
    "s_len:s" =>\$start_len,
    "fix_dup:s" =>\$fix_dup,
    "bar_edge:s" =>\$bar_edge_num,
    "s_bs:s" =>\$f_barcode_num,
    "min_edge:s" =>\$the_min_barcode_edge,
    "max_edge:s" =>\$the_max_barcode_edge,
    "sp_minus:s" =>\$sp_minus,
    "r_type:s" => \$rtype,
    "c:s" => \$confict,
    "p_svg:s" => \$p_svg,
    "score:s" => \$score,
    "mlen:s" => \$alen,
    "d:s" => \$dt,
    "t_rpe:s" => \$tp,
    "sd_rpe:s" => \$sdp,
    "t_rse:s" => \$ts,
    "sd_rse:s" => \$sds,
    "wind:s" => \$wind,
    "rdsn:s" => \$rdsn,
    "s_l:s" => \$s_l,
    "s_q:s" => \$s_q,
    "s_n:s" => \$s_n,
    "f_other:s" => \$f_other,
    "st:s" =>\$samtools,
    "soke:s" =>\$soapnuke,
    "bwa:s" => \$bwa,
    "botie2:s" => \$bowtie2,
    "bdl:s" => \$bamdeal,
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
my $barcode_list ="$Bin/stLFR_barcode_split/barcode_list.txt";
# step 5
$minread ||=5;
$mdep ||='auto';
$mcycle ||=10;
$cutlen ||=500;
# step 6
$insert_size ||=500;
$overlap ||=0.8;
$filter_barcode ||=4;
$m_q ||=50;
# step 7
$start_nod ||='undef';
$scaffold_way ||='all';
$bias ||=1;
$start_len ||=5000;
$fix_dup ||='';
$bar_edge_num ||=1;
$f_barcode_num ||='auto';
$the_min_barcode_edge ||=3;
$the_max_barcode_edge ||=12;
$sp_minus ||='max';
$rtype ||='first';
$confict ||='0';
#step4 
$p_svg ||='T';
$score ||='1:2';
$alen ||=20;
$dt ||=30;
$tp ||=2;
$sdp ||=5;
$ts ||=1;
$sds ||=20;
$wind ||=10;
$rdsn ||=2;
# step2
$s_l ||=5;
$s_q ||=0.1;
$s_n ||=0.01;
my $soap_paramter = '';
if ($f_other){
    my @option_set = split/\,/,$f_other;
    foreach my $pairs (pairs @option_set){
        $soap_paramter .='-'.${$pairs}[0].' '.${$pairs}[1].' ';
    }
}
#tools path 
$samtools ||= "/hwfssz1/ST_BIOCHEM/P18Z10200N0255/PMO/F14ZQSYJSY1726/pangwending/some_tools/samtools-1.13/bin/samtools";
$bwa ||="/share/app/bwa/0.7.17/bwa";
$bowtie2 ||="/share/app/bowtie2/2.4.1/bowtie2";
$soapnuke ||= "/share/app/soapnuke/2.1.0/SOAPnuke";
$bamdeal = "$Bin/BamDeal-master/bin/BamDeal_Linux";
`chmod a+x $bamdeal`;
# step information
my ($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8);
$step ||='all';

if ($step eq 'all'){
    die $usage unless -s $fastq1;
}

print '***** Start run the [Sgrstlfr] *****'."\n";
print 'Start: '.localtime()."\n";
my $start_time = time();
my @a_step_set = qw/0 0 0 0 0 0 0 0/;
($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8) = @a_step_set;
if ($step eq 'all'){
    ($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8) = (1,1,1,1,1,1,1,1);
    print '##### Start all the step!,loading.........'."\n";
}elsif ($step eq 'assembly'){
    $step7 = 1;$step8=1;
    print 'Only start the step7-8!, loading......'."\n";
}else{
    my @step_set = split/\,/,$step;
    print '##### Start the step:'.join(",",@step_set)."\t".'loading.....'."\n";
    foreach my $anum (@step_set){
        if ($anum =~ /[1-8]/){
            splice(@a_step_set,$anum-1,1,1);
        }else{
            print '!!!The parameter -step not a number in (1-8) ! '.$anum."\n";
            exit;
        }
    }
    ($step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8) = @a_step_set;
}


# script path
my $split_barcode = "$Bin/bin/split_barcode.pl";
my $class_mapping_result = "$Bin/bin/class_mapping_result.pl";
my $get_loxp_region = "$Bin/bin/get_loxp_region.pl";
my $classify ="$Bin/bin/classify.pl";
my $bkptlog2gvr = "$Bin/bin/bkptlog2gvr.pl";
my $bpktdepth = "$Bin/bin/bpktdepth.pl";
my $drawdepth_splitmap = "$Bin/bin/drawdepth_splitmap.pl";
my $targetcoverage = "$Bin/bin/targetcoverage.pl";
my $getvariareads = "$Bin/bin/getvariareads.pl";
my $extnodedge = "$Bin/bin/extnodedge.pl";
my $fix_edge = "$Bin/bin/fix_edge.pl";
my $draw_dot = "$Bin/bin/draw_dot.pl";
my $get_barcode_node = "$Bin/bin/get_barcode_node.pl";
my $mersh_basic_file = "$Bin/bin/mersh_basic_file.pl";
my $repair_barcode ="$Bin/bin/repair_barcode.pl";
my $exact_link = "$Bin/bin/exact_link.pl";
my $Dynpath = "$Bin/bin/Dynpath.pl";
my $stat_length_sp = "$Bin/bin/stat_length_sp.pl";
my $translate_path = "$Bin/bin/translate_path.pl";
my $get_gene_variation = "$Bin/bin/get_gene_variation.pl";
my $SCBINdotplot = "$Bin/bin/SCBINdotplot.pl";


# check file and dependencies
die "The script split_barcode.pl is not found in $split_barcode!\n" unless -s $split_barcode;
die "The script class_mapping_result.pl is not found in $class_mapping_result!\n" unless -s $class_mapping_result;
die "The script get_loxp_region.pl is not found in $get_loxp_region!\n" unless -s $get_loxp_region;
die "The script classify.pl is not found in $classify!\n" unless -s $classify;
die "The script bkptlog2gvr.pl is not found in $bkptlog2gvr!\n" unless -s $bkptlog2gvr;
die "The script bpktdepth.pl is not found in $bpktdepth!\n" unless -s $bpktdepth;
die "The script drawdepth_splitmap.pl is not found in $drawdepth_splitmap!\n" unless -s $drawdepth_splitmap;
die "The script targetcoverage.pl is not found in $targetcoverage!\n" unless -s $targetcoverage;
die "The script extnodedge2.pl is not found in $extnodedge!\n" unless -s $extnodedge;
die "The script draw_dot.pl is not found in $draw_dot!\n" unless -s $draw_dot;
die "The script get_barcode_node.pl is not found in $get_barcode_node!\n" unless -s $get_barcode_node;
die "The script mersh_basic_file.pl is not found in $mersh_basic_file!\n" unless -s $mersh_basic_file;
die "The script repair_barcode.pl is not found in $repair_barcode!\n" unless -s $repair_barcode;
die "The script exact_link.pl is not found in $exact_link!\n" unless -s $exact_link;
die "Thes script Dynpath.pl is not found in $Dynpath!\n" unless -s $Dynpath;
die "The script stat_length_sp.pl is not found in $stat_length_sp!\n" unless -s $stat_length_sp;
die "The script translate_path.pl is not found in $translate_path!\n" unless -s $translate_path;
die "The script get_gene_variation.pl is not found in $get_gene_variation!\n" unless -s $get_gene_variation;
die "The script SCBINdotplot.pl is not found in $SCBINdotplot!\n" unless -s $SCBINdotplot;
# mkdir father-directory

system ("mkdir -p $outdir/$filename");
my $file_path = "$outdir/$filename";
#### module 1 split barcode 
if ($step1 == 1){
    die "ERROR: The stLFR-seq file fastq1 not found in $fastq1!\n" unless -s $fastq1;
    die "ERROR: The stLFR-seq file fastq2 not found in $fastq2!\n" unless -s $fastq2;
    die "ERROR: The barcode list file not found in $barcode_list !\n" unless -s $barcode_list;
    print '!Start to run the step1,loading......'."\n";
    `perl $split_barcode $barcode_list $fastq1 $fastq2 split_reads $file_path > $file_path/split_barcode.log 2>$file_path/split_barcode.err`;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step1 split_barcode is not finished! Pleas check the running commond!"."\n";
        exit;
    }
    print '#The step1 done!'."\n";
}

#### module 2 clean data 
if ($step2 == 1){
    die "The tools of soapnuke path is not found in $soapnuke!\n" unless -s $soapnuke;
    if ($step1 != 0){
        die "ERROR: The split_barcode results split_reads.1.fq.gz not found in $file_path/split_reads.1.fq.gz!\n " unless -s "$file_path/split_reads.1.fq.gz";
        die "ERROR: The split_barcode results split_reads.2.fq.gz not found in $file_path/split_reads.2.fq.gz!\n" unless -s "$file_path/split_reads.2.fq.gz";
        print 'Detect the split_reads.1.fq.gz split_reads.2.fq.gz !!!'."\n";
        print '!Start to run the step2,loading......'."\n";
        `$soapnuke filter -l $s_l -q $s_q -n $s_n -T $threads $soap_paramter -1 $file_path/split_reads.1.fq.gz -2 $file_path/split_reads.2.fq.gz -C split_clean_reads.1.fq.gz -D split_clean_reads.2.fq.gz -o $file_path `;
        print "The erro number:$?\n";
        if ($? != 0){
            print "ERROR: the step2 clean_data is not finished! Pleas check the following run commond!"."\n";
            print "$soapnuke filter -l $s_l -q $s_q -n $s_n -T $threads $soap_paramter -1 $file_path/split_reads.1.fq.gz -2 $file_path/split_reads.2.fq.gz -C split_clean_reads.1.fq.gz -D split_clean_reads.2.fq.gz -o $file_path\n";
            exit;
        }
    }else{
        die "ERROR: The stLFR-seq split_data.1 not found in $fastq1!\n" unless -s $fastq1;
        die "ERROR: The stLFR-seq split_data.2 not found in $fastq2!\n" unless -s $fastq2;
        print "Use the input $fastq1 $fastq2 to filter !!!\n";
        print '!Start to run the step2,loading......'."\n";
        `$soapnuke filter -l $s_l -q $s_q -n $s_n -T $threads $soap_paramter -1 $fastq1 -2 $fastq2 -C split_clean_reads.1.fq.gz -D split_clean_reads.2.fq.gz -o $file_path `;
        print "The erro number:$?\n";
        if ($? != 0){
            print "ERROR: the step2 clean_data is not finished! Pleas check the following run commond!"."\n";
            print "$soapnuke filter -l $s_l -q $s_q -n $s_n -T $threads $soap_paramter -1 $fastq1 -2 $fastq2 -C split_clean_reads.1.fq.gz -D split_clean_reads.2.fq.gz -o $file_path\n";
            exit;
        }
    }
    print "#The step2 done!\n";
}

#### module 3 alignment to unSCRaMbLE genome
if ($step3 == 1){
    die "ERROR: The reference of fasta file is not found in $fasta!\n" unless -s $fasta;
    die "The tools of bwa path is not found in $bwa!\n" unless -s $bwa;
    die "The tools of samtools path is not found in $samtools!\n" unless -s $samtools;
    if ($step2 != 0){
        die "ERROR: The split and clean results split_clean_reads.1.fq.gz not found in $file_path/split_clean_reads.1.fq.gz!\n " unless -s "$file_path/split_clean_reads.1.fq.gz";
        die "ERROR: The split and clean results split_clean_reads.2.fq.gz not found in $file_path/split_clean_reads.2.fq.gz!\n" unless -s "$file_path/split_clean_reads.2.fq.gz";
        print 'Detect the split_clean_reads.1.fq.gz and split_clean_reads.2.fq.gz !!!'."\n";
        print '!Start to run the step3,loading......'."\n";
        `$bwa mem -t $threads $fasta $file_path/split_clean_reads.1.fq.gz $file_path/split_clean_reads.2.fq.gz |$samtools sort -@ $threads -> $file_path/$filename.sort.bam `;
        `$samtools index $file_path/$filename.sort.bam `;
        print "The erro number:$?\n";
        if ($? != 0){
            print "ERROR: the step3 alignment is not finished! Pleas check the following run commond!"."\n";
            print "$bwa mem -t $threads $fasta $file_path/split_clean_reads.1.fq.gz $file_path/split_clean_reads.2.fq.gz |$samtools sort -@ $threads -> $file_path/$filename.sort.bam\n";
            exit;
        }
    }else{
        die "ERROR: The stLFR-seq Clean_data.1 not found in $fastq1!\n" unless -s $fastq1;
        die "ERROR: The stLFR-seq Clean_data.2 not found in $fastq2!\n" unless -s $fastq2;
        print "Use the input fastq file $fastq1 $fastq2 to alignment!!!\n";
        print '!Start to run the step3,loading......'."\n";
        `$bwa mem -t $threads $fasta $fastq1 $fastq2 |$samtools sort -@ $threads -> $file_path/$filename.sort.bam `;
        `$samtools index $file_path/$filename.sort.bam `;
        print "The erro number:$?\n";
        if ($? != 0){
            print "ERROR: the step3 alignment is not finished! Pleas check the following run commond!"."\n";
            print "$bwa mem -t $threads $fasta $fastq1 $fastq2 | $samtools sort -@ $threads -> $file_path/$filename.sort.bam\n";
            exit;
        }
        print "#The step3 done!\n";
    }
}

#### module 4 stat the mapping result and split the unmapping reads include loxp
if ($step4 == 1){
    die "ERROR: The alignment results $filename.sort.bam not found in $file_path/$filename.sort.bam!\n" unless -s "$file_path/$filename.sort.bam";
    die "ERROR: The unSCRaMbLE genome fasta file not found in $fasta!\n" unless -s $fasta;
    die "The tools of bamdeal path is not found in $bamdeal!\n" unless -s $bamdeal;
    die "The tools of samtools path is not found in $samtools!\n" unless -s $samtools;
    print '!Start to run the step4,loading......'."\n";
    `perl $get_loxp_region $fasta $chrid $filename $file_path`;
    `$bamdeal statistics Coverage -i $file_path/$filename.sort.bam -r $fasta -q 0 -o $file_path/$filename -w 5000 `;
    `perl $class_mapping_result $file_path/$filename.sort.bam -n $filename -chr $chrid -q $m_q -s $samtools -o $file_path`;
    `$bwa mem $fasta $file_path/$filename.split.fq > $file_path/$filename.split.sam`;
    if ($p_svg eq 'T'){
        print 'Start to produce the split-mapping svg graph! loading........'."\n";
        my $chr_len = `cat $file_path/$filename.list.loxpreg |tail -1 |awk '{print \$5}'`;
        chomp($chr_len);
        `$bowtie2 -x $fasta -q $file_path/$filename.split.fq -k 100 > $file_path/$filename.bowsplit.sam `;
        `perl $targetcoverage $file_path/$filename.depthsite.fa.gz -spid $filename -wind $wind -chrid $chrid -outdir $file_path`;
        `perl $classify $file_path/$filename.bowsplit.sam -n $filename -chrid $chrid -l $chr_len -score $score -mlen $alen -d $dt -o $file_path `;
        `perl $bpktdepth $file_path/$filename.multi.rpe.gz -n $filename -sp $filename -rdstype rpe -uniq -t $tp -sd $sdp -o $file_path `;
        `perl $bpktdepth $file_path/$filename.multi.rse.gz -n $filename -sp $filename -rdstype rse -uniq -t $ts -sd $sds -o $file_path `;
        `perl $bkptlog2gvr $file_path/$filename.rse.bkpt.log -n $filename -o $file_path -type rse -rdsn $rdsn`;
        `perl $bkptlog2gvr $file_path/$filename.rpe.bkpt.log -n $filename -o $file_path -type rpe -rdsn $rdsn`;
        `perl $drawdepth_splitmap -n $filename -isz $insert_size -refid $chrid -reflen $chr_len -loxpftr $file_path/$filename.lxpftr.lst -cvg $file_path/$filename.tar.depth -rpe $file_path/$filename.all.rpe.bkpt.gvr -rse $file_path/$filename.all.rse.bkpt.gvr -o $file_path`;
        print 'Get svg graph is ok !'."\n";
        `rm $file_path/$filename.rse.bkpt.log $file_path/$filename.rpe.bkpt.log $file_path/$filename.multi.rpe.gz $file_path/$filename.multi.rse.gz $file_path/$filename.tar.depth $file_path/$filename.bowsplit.sam`;
        `rm $file_path/$filename.bpk.rpe.xls $file_path/$filename.bpk.rse.xls`;
    }
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step4 is not finished! Pleas check the median file and the running commond !"."\n";
        exit;
    }
    
    print "#Step 4 done!"."\n";
}

#### module 5 get the edge and node 
if ($step5 == 1){
    die "ERROR: The step4 produce  file $filename.split.sam not found in $file_path/$filename.split.sam!\n" unless -s "$file_path/$filename.split.sam";
    die "ERROR: The step4 produce coverage file $filename.depthsite.fa.gz not found in $file_path/$filename.depthsite.fa.gz!\n" unless -s "$file_path/$filename.depthsite.fa.gz";
    die "ERROR: The step4 produce $filename.synmap.gz not found in $file_path/$filename.synmap.gz!\n" unless -s "$file_path/$filename.synmap.gz";
    die "ERROR: The step4 produce $filename.list.loxpreg not found in $file_path/$filename.list.loxpreg!\n" unless -s "$file_path/$filename.list.loxpreg";
    die "ERROR: The alignment results $filename.sort.bam not found in $file_path/$filename.sort.bam!\n" unless -s "$file_path/$filename.sort.bam";
    die "ERROR: The reference of unSCRaMbLE genome coverage file not found in $refcoverage!\n" unless -s $refcoverage;
    print '!Start to run the step5,loading......'."\n";
    `perl $getvariareads $file_path/$filename.split.sam -n $filename -s $samtools -o $file_path`;
    if ($mdep eq 'auto'){
        die "ERROR: The seqence depth stat file not found in $file_path/$filename.stat! You could input the int number of mdep paramaters or produce\n" unless -s "$file_path/$filename.stat";
        my $mean_depth = `cat $file_path/$filename.stat |tail -1 |awk '{print \$6}'`;
        chomp ($mean_depth);
        if ($chr_type eq 'cycle'){
            $mean_depth = int($mean_depth/2);
        }
        print '#Use the average seqence depth '.$mean_depth.' to estimation the cnv!!!'."\n";
        `perl $extnodedge $file_path/$filename.depthsite.fa.gz $file_path/$filename.sort.bam $file_path/$filename.list.loxpreg $file_path/$filename.gvr -refcover $refcoverage -n $filename -chrid $chrid -chrtype $chr_type -s $samtools -minread $minread -mdep $mean_depth -mcycle $mcycle -cutlen $cutlen -o $file_path > $file_path/extnodedge.log`;
    }else{
        `perl $extnodedge $file_path/$filename.depthsite.fa.gz $file_path/$filename.sort.bam $file_path/$filename.list.loxpreg $file_path/$filename.gvr -refcover $refcoverage -n $filename -chrid $chrid -chrtype $chr_type -s $samtools -minread $minread -mdep $mdep -mcycle $mcycle -cutlen $cutlen -o $file_path > $file_path/extnodedge.log`;
    }
    system("perl $draw_dot $file_path/$filename.edg $file_path/$filename.nod $filename $file_path");
    print "The erro number:$?\n";  #6400 edge calculte erro repair the mdep  
    if ($? != 0){
        print "ERROR: the step5 get the edge and node is not finished! Pleas check the extnodedge.log and running commmond !!!"."\n";
        exit;
    }
    
    print "#Step 5 done!"."\n";
}

#### module 6 prepare the basic co-barcode information and filter the barcode
if ($step6 == 1){
    die "ERROR: The step5 produce $filename.edg file not found in $file_path/$filename.edg!\n" unless -s "$file_path/$filename.edg";
    die "ERROR: The step5 produce $filename.nod file not found in $file_path/$filename.nod!\n" unless -s "$file_path/$filename.nod";
    die "ERROR: The step4 produce $filename.depthsite.fa.gz file not found in $file_path/$filename.depth!\n" unless -s "$file_path/$filename.depthsite.fa.gz";
    die "ERROR: The step4 produce $filename.synmap.gz not found in $file_path/$filename.synmap.gz!\n" unless -s "$file_path/$filename.synmap.gz";
    die "ERROR: The step4 produce $filename.split.sam not found in $file_path/$filename.split.sam!\n" unless -s "$file_path/$filename.split.sam";
    print '!Start to run the step6,loading......'."\n";
    `perl $mersh_basic_file $file_path/$filename.edg $file_path/$filename.list.loxpreg $file_path/$filename.depthsite.fa.gz $chrid $filename $file_path `;
    `perl $get_barcode_node $file_path/$filename.edg $file_path/$filename.nod $file_path/$filename.list.loxpreg $file_path/$filename.split.sam $file_path/$filename.synmap.gz -n $filename -isz $insert_size -t $chr_type -overlap $overlap -o $file_path `;
    `perl $repair_barcode $file_path/$filename.node.barcode $file_path/$filename.notnode.barcode $file_path/$filename.edg $file_path/$filename.nod $filename $file_path`;
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
    die "ERROR: Thes step6 produce $filename.merge.barcode not found in $file_path/$filename.merge.barcode!\n" unless -s "$file_path/$filename.merge.barcode";
    print '!Start to run the step7,loading......'."\n";
    `perl $exact_link $file_path/$filename.sort.barcode.stat $file_path/$filename.nod $file_path/$filename.edg $file_path/$filename.total -n $filename -o $file_path -chr_type $chr_type -b $bar_edge_num -cf $confict -fix $fix_dup -fp $f_barcode_num -i $bias -len $start_len -start_nod $start_nod -way $scaffold_way > $file_path/$filename.scaf.log `;
    unless (-s "$file_path/$filename.info"){
        `perl $Dynpath $filename $file_path $confict`;
    }
    `perl $stat_length_sp $file_path/$filename.result $file_path/$filename.info $file_path/$filename.total $file_path/$filename.merge.barcode -n $filename -chr_type $chr_type -min $the_min_barcode_edge -max $the_max_barcode_edge -t $threads -sp_minus $sp_minus -o $file_path`;
    if (-s "$file_path/$filename.2.dot"){
        `dot -Tpdf $file_path/$filename.2.dot -o $file_path/$filename.2.dot_graph.pdf`;
    }
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step7:get the scaffold and exact genome is not finished! Pleas check the scaffold log!!!"."\n";
        exit;
    }
    print "#Step 7 done!"."\n";
}

#### module 8 transfer the edge-nod path genome to fasta!
if ($step8 == 1){
    die "ERROR: The step7 produce $filename.path not found in $file_path/$filename.path!\n" unless -s "$file_path/$filename.path";
    die "ERROR: The step5 produce $filename.edg file not found in $file_path/$filename.edg!\n" unless -s "$file_path/$filename.edg";
    die "ERROR: The unSCRaMbLE fasta file is not found in $fasta!\n" unless -s $fasta;
    print '!Start to run the step8,loading......'."\n";
    `perl $translate_path $file_path/$filename.path $file_path/$filename.list.loxpreg $file_path/$filename.edg $fasta -c $chrid -n $filename -o $file_path -a $rtype `;
    `perl $get_gene_variation $file_path/$filename.list.loxpreg $file_path/$filename.path.index $filename $file_path`;
    `perl $SCBINdotplot -prefix $filename -Xencod $file_path/$filename.ref.encod -Yencod $file_path/$filename.path.index -outdir $file_path`;
    print "The erro number:$?\n";
    if ($? != 0){
        print "ERROR: the step8:transfer the genome to fasta is not finished! Pleas check the following run commond!"."\n";
        print "perl $translate_path $file_path/$filename.path $file_path/$filename.list.loxpreg $file_path/$filename.edg $fasta -c $chrid -t $chr_type -n $filename -a $rtype -o $file_path\n";
        exit;
    }
    print "#Step 8 done!"."\n";
}
print 'End: '.localtime()."\n";
print '############### All work is finished!!! ##################'."\n";





