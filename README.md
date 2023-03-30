Descrption: Sgrstlfr(Synthetic genome reconstrution by single tube long fragment reads)

##This software is order to assembly the synthetic yeast chromosome by the stLFR seqence !!!

##The SCRaMbLE[1] genome has a complex structure, which brings great difficulty to conventional assembly methods. Moreover, due to the limitations of the short read length of Syngenor, which is currently developed based on next-generation sequencing technology, 60% of the SCRaMbLE genome cannot be accurately analyzed [2]. However, fully reconstructing the SCRaMbLE genome is of great significance for the downstream application research of synthetic yeast genome. Therefore, in order to solve the above technical problems, we use the virtual long molecular interval sequencing technology developed by BGI, stLFR[3] (single-tube long-read sequencing technology), this technology uses the method of next-generation sequencing to sequence the same barcode ( barcode) long molecules are sequenced, so that the short reads from the same long molecule after sequencing contain the same barcode, that is, the short reads of the same barcode have an adjacent relationship (adjacent relationship), so to some extent, It is equivalent to "long molecule". This software is based on the principle of the proximity relationship of long molecules, which solves the problem of path branching caused by variation in the genome assembly process, and achieves the purpose of accurately assembling the SCRaMbLE genome.


**The simple useage:
The program only needs to contain the following basic parameters to run

--fa ; The fasta file of the unSCRaMbLE genome that has been indexed.
--fq1 : 1-end fastq file of stLFR paired-end sequencing.
--fq2 : 2-terminal fastq file of stLFR paired-end sequencing.
--rfcvg : Sequencing depth file for the unSCRaMbLE genome
--n : The output process file and the file name of the result. [Note: The pre-file name of the process file search]
--o : The output process file and the final result path. [Note: The path to find the process file]
--step : The steps of the program running, divided into three types [all | assembly |others].
all: run all steps [default]; assembly: only run module 7, the process file containing modules 4-6 is required! ! ! ;others: specify the corresponding run steps (likes 4,5,6,7)
--chrtype : type of synthetic chromosome, linear or circular [default: cycle]
--chrid : synthetic chromosome name [default:IXR_BACseq]
--t : Number of threads for software comparison to run [default: 4]
--tools : Enter the relevant software path: including bwa, bowtie2, bamdeal, soapnuke, etc.

**For example: 
-> perl run_Sgrstlfr.pl -fa BY4741chr9RD_SynIXR.fa -fq1 test1.fq.gz -fq2 test2.fq.gz -rfcvg JS94.depthsite.fa.gz -n test -o outdir -chrid IXR_BACseq -chrtype cycle -t 8 -step all
-> perl run_Sgrstlfr.pl -fa BY4741chr9RD_SynIXR.fa -fq1 test1.fq.gz -fq2 test2.fq.gz -rfcvg JS94.depthsite.fa.gz -n test -o outdir -chrid IXR_BACseq -chrtype cycle -t 8 -step 6,7 -sp_minus max 


Description: This script is used to reconstruct the synthethic yeast SCRaMbLE genome by stLFR technology !
                 v1.0 2023.03.24 
    Usage: perl run_Sgrstlfr.pl [-option] <h|help>

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

**Reference:

[1]	DYMOND J, BOEKE J. The Saccharomyces cerevisiae SCRaMbLE system and genome minimization [J]. Bioeng Bugs, 2012, 3(3): 168-71.

[2]	SHEN Y, STRACQUADANIO G, WANG Y, et al. SCRaMbLE generates designed combinatorial stochastic diversity in synthetic chromosomes [J]. Genome Res, 2016, 26(1): 36-49.

[3]	WANG O, CHIN R, CHENG X, et al. Efficient and unique cobarcoding of second-generation sequencing reads from long DNA molecules enabling cost-effective and accurate sequencing, haplotyping, and de novo assembly [J]. Genome Res, 2019, 29(5): 798-808.

**and some reference software and Thans for the software

SOAPnuke:https://github.com/BGI-flexlab/SOAPnuke。

BWA:https://github.com/lh3/bwa

Bamdeal:https://github.com/BGI-shenzhen/BamDeal

stLFR_barcode_split:https://github.com/BGI-Qingdao/stLFR_barcode_split

****if have some bug please contact:
pangwending\@genomics.cn    wangyun\@genomics.cn
