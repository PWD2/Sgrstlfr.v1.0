#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;

Describe: identify breakpoint in a reads; 
          classify 3 type of read with bisection by split mapping...
          read type: syn-syn syn-wt wt-wt 
         
Author: Yun Wang, wangyun\@genomic.cn
Version: 1.0, Date: 2014-03-27

Usage:perl $0 <split.map> [-option] 

	-chrid        		: synthetic chromosome id [defaclut:IXR_BACseq] 
	-l         	 	: synthetic sequence length (start..end) [defaclut:1..100371] 
	-score        		: score pattern calculating cumulative score,match score:mismath penalty [1:2] 
	-mlen           	: minimal length of match or mismatch allowed  [default: 10]   
	-d             		: minimal mapping distance recongnized as recombination events [default: 100]  
	-cinfo          	: circle chromosome information [default:chrmt:85779,2-micron:6318,IXR_BACseq:100371]
	-n     		  	: Sample ID/prefix of output file   [defaclut:test]
	-o            	 	: Directory of output [default: ./]   
	-h            	 	: show this message

Example: perl $0 JS599.split.map -n JS599 -score 1:2 -l 1..100371 -mlen 10 -d 100 -o ./ 
   
Note: # mapping with gaps is illigle in this programe,we wil ignore it...         
     indel was recognized as Mismatch...
     Priority:--
      1> no recombination events;  x100
      2> chromosome internal recombination; x50 
      3> chromosome external recombination with wild-type chromosome; x25
      4> chromosome external recombination with synthetic chromosome; x5

USAGE

my ($prefix,$synchr,$synscale,$score,$minlen,$redis,$circhr,$outdir,$help);
GetOptions(
   "n:s" => \$prefix,
   "chrid:s" => \$synchr,
   "l:s" => \$synscale,
   "score:s" => \$score,
   "mlen:s" => \$minlen,
   "d:s" => \$redis,
   "cinfo:s" => \$circhr,
   "o:s"  => \$outdir,
   "h"  => \$help,
);

die $usage if (!$ARGV[0] ||$help);

$score ||= "1:2";
$prefix ||= "test";
$minlen ||= 5;
$redis ||= 100;
$synchr ||= "IXR_BACseq";
$synscale ||= "1..100371";
$outdir ||= getcwd;
$circhr ||= "chrmt:85779,2-micron:6318,IXR_BACseq:100371";

my %circhrinfo; 
foreach (split /,/,$circhr){
  my ($circhr_id,$circhr_len)= split /:/,$_;
  $circhrinfo{$circhr_id} = $circhr_len;
}

my ($match_score,$mismatch_score) = split /:/,$score;
my $gap_score = 3;
my ($synstar,$synend) = split /\.\./,$synscale;
#print "##$synstar,$synend\n";

if ($ARGV[0] =~ /.gz$/){
  	open IN,"gzip -dc $ARGV[0] |" || die $!;
}else{
	open IN, $ARGV[0] ||die $!; #open SAM alignment file
}
	
open OUT1, "| gzip >$outdir/$prefix.multi.rpe.gz" ||die $!; #JSnum.rpe
open OUT2, "| gzip >$outdir/$prefix.multi.rse.gz" ||die $!; #JSnum.rse
#open OUT2, ">$outdir/$prefix.rse" ||die $!; #JSnum.rse

my (%rds1,%rds2);
my %mapnum;
my ($prdsid,$pid);
my ($bestrds1,$bestrds2);
my ($bestscore1,$bestscore2) = (0,0);
my $optinfo ="";
my $optscore = 0;
my ($readsid,$strand,$chrid,$site,$seq,$macht,$md);

while(1){
	  my $Filine=<IN>;
	  next if ($Filine =~ /^\@/);
    chomp $Filine;
    my @line = split /\s+/,$Filine;	

    if ($line[1]!=4){		
		  if ($line[1] == 16 || $line[1] == 272) {
			  $strand = "-";
		  }elsif($line[1] == 0 || $line[1] == 256){
			  $strand = "+";
		  }	  
	    ($readsid,$chrid,$site,$macht,$seq,$md) = @line[0,2,3,5,9,18];
	    
	    if(exists $mapnum{$readsid}){
	    	$mapnum{$readsid} ++;
	    }else{$mapnum{$readsid} =1;}
	    
	    $md = $line[17] unless($md=~ /^MD/);   	    
	    $prdsid = (split /\#/,$readsid)[0];
      #$pid = (split /_/,$readsid)[0]; ##for pair-end sequencing reads...
      $pid= $prdsid;
      
      my $maplen = length ($seq);  
      next if $maplen < $minlen;
      my $mismatch_num  = 0;
      $mismatch_num  ++ while ($md =~/[ATCGN]/gi);
      my $gaps_num = 0;
      $gaps_num ++ while($macht =~ /[ID]/gi);
      
      my $rdsscore = $match_score * ($maplen  - $mismatch_num) - $mismatch_score * $mismatch_num -  $gaps_num*$gap_score;
      my $mapinfo = "$readsid\t$strand\t$chrid\t$maplen\t$site\t$macht";
        
	      if ($readsid =~/\/1$/){
	         push @{$rds1{$readsid}{"info"}},$mapinfo; 
	         push	@{$rds1{$readsid}{"score"}},$rdsscore; 
	       	 $bestrds1 = "$mapinfo"."\t-\t5\t-";
	       	 $bestscore1 = $rdsscore;
	      }elsif($readsid =~/\/2$/){
	         push @{$rds2{$readsid}{"info"}}, $mapinfo; 
	       	 push @{$rds2{$readsid}{"score"}} , $rdsscore; 
	       	 $bestrds2 = "$mapinfo"."\t-\t5\t-";;
	       	 $bestscore2 = $rdsscore; 	       	 	       	       
	      }else {print "reads ID error:  ectopic ID -- $readsid\n";}	    
	    last; 
	 }
}

while (<IN>){
    chomp;
    my @line = split /\s+/,$_;	
    if ($line[1]!=4){		
		  if ($line[1] == 16 || $line[1] == 272) {
			  $strand = "-";
		  }elsif($line[1] == 0 || $line[1] == 256){
			  $strand = "+";
		  }
		    
	    ($readsid,$chrid,$site,$macht,$seq,$md) = @line[0,2,3,5,9,18];
	    
	    if(exists $mapnum{$readsid}){
	    	$mapnum{$readsid} ++;
	    }else{$mapnum{$readsid} =1;}    
	    
	    $md = $line[17] unless($md=~ /^MD/);    
	     
	    my $prdsidnew = (split /\#/,$readsid)[0];
      my $pidnew = (split /_/,$readsid)[0]; 
     
      
      my $maplen = length ($seq);  
      next if $maplen < $minlen;
      my $mismatch_num = 0;
      $mismatch_num  ++ while ($md =~/[ATCGN]/gi);
      
      my $gaps_num = 0;
      $gaps_num ++ while($macht =~ /[ID]/gi);      
      my $rdsscore = $match_score * ($maplen  - $mismatch_num) - $mismatch_score * $mismatch_num -  $gaps_num*$gap_score;
      my $mapinfo = "$readsid\t$strand\t$chrid\t$maplen\t$site\t$macht";            
	 
      
      if ($prdsidnew ne $prdsid){
      ##OUTPUT
         
         ($optinfo,$optscore) = identify ($optinfo,$optscore,\%rds1,\%rds2,\%mapnum);   
         print OUT1 "$optinfo\n" if defined $optinfo; 
         
         unless(defined $optinfo){
            if (defined $bestrds1 ||$bestrds2){
               if($bestscore1 >= $bestscore2){
            	    my $optial_id = (split /\s+/,$bestrds1)[0];
            	    if (exists $mapnum{$optial_id}){
         	  	       print OUT2 "$bestrds1\t$mapnum{$optial_id}\n";
         	        }else{
         	  	       print "warning: Not found the number of mapping location...\n";
         	        }  
            	    #print OUT2 "$bestrds1\n";
               }else{
               	  my $optial_id = (split /\s+/,$bestrds2)[0];
            	    if (exists $mapnum{$optial_id}){
         	  	       print OUT2 "$bestrds2\t$mapnum{$optial_id}\n";
         	        }else{
         	  	       print "warning: Not found the number of mapping location...\n";
         	        }  
               	#print OUT2 "$bestrds2\n";
               }
            }
         }
     ###initialization        
         undef %rds1;
         undef %rds2;
         undef $optinfo;
         $optscore = 0; 
         
         undef %mapnum;
         $mapnum{$readsid} =1;
         
         ($bestrds1,$bestrds2) =("",""); 
         ($bestscore1,$bestscore2) = (0,0);      
         
         $prdsid = $prdsidnew; 
	       $pid = $pidnew;
	       
	       if ($readsid =~/\/1$/){
	         push @{$rds1{$readsid}{"info"}}, $mapinfo; 
	       	 push @{$rds1{$readsid}{"score"}}, $rdsscore; 
	       	 if ($bestscore1 < $rdsscore){
	       	   $bestrds1 = "$mapinfo"."\t-\t5\t-";
	       	   $bestscore1 = $rdsscore;
	       	 }
	       }elsif($readsid =~/\/2$/){
	        push  @{$rds2{$readsid}{"info"}}, $mapinfo; 
	       	push  @{$rds2{$readsid}{"score"}}, $rdsscore; 
	       	 if ($bestscore2 < $rdsscore){
	       	   $bestrds2 = "$mapinfo"."\t-\t5\t-";
	       	   $bestscore2 = $rdsscore;
	       	   $bestrds2 = "$mapinfo"."\t-\t5\t-";
	       	   $bestscore2 = $rdsscore;
	       	 }	      
	       }else {print "reads ID error:  ectopic ID -- $readsid\n";}	          
      
      }else{
        if ($pidnew ne $pid){
           ($optinfo,$optscore) = identify ($optinfo,$optscore,\%rds1,\%rds2,\%mapnum);
           $pid = $pidnew; 
           undef %rds1;
           undef %rds2;
           
           
           if ($readsid =~/\/1$/){
	           push @{$rds1{$readsid}{"info"}}, $mapinfo; 
	         	 push @{$rds1{$readsid}{"score"}}, $rdsscore; 
	       	   if ($bestscore1 < $rdsscore){
	       	     $bestrds1 = "$mapinfo"."\t-\t5\t-";
	       	     $bestscore1 = $rdsscore;
	       	   }
	         }elsif($readsid =~/\/2$/){
	           push @{$rds2{$readsid}{"info"}}, $mapinfo; 
	         	 push @{$rds2{$readsid}{"score"}},$rdsscore; 
	       	   if ($bestscore2 < $rdsscore){
	       	      $bestrds2 = "$mapinfo"."\t-\t5\t-";
	       	      $bestscore2 = $rdsscore;
	       	   }	      
	         }else{print "reads ID error:  ectopic ID -- $readsid\n";}	          
        
        }else{
	         if ($readsid =~/\/1$/){
	         push @{$rds1{$readsid}{"info"}}, $mapinfo; 
	         push @{$rds1{$readsid}{"score"}}, $rdsscore; 
	       	   if ($bestscore1 < $rdsscore){
	       	     $bestrds1 = "$mapinfo"."\t-\t5\t-";
	       	     $bestscore1 = $rdsscore;
	       	   }
	       	 }elsif($readsid =~/\/2$/){
	         push @{$rds2{$readsid}{"info"}}, $mapinfo; 
	         push @{$rds2{$readsid}{"score"}}, $rdsscore; 
	       	   if ($bestscore2 < $rdsscore){
	       	      $bestrds2 = "$mapinfo"."\t-\t5\t-";
	       	      $bestscore2 = $rdsscore;
	       	     }       	 	       	       
	         }else {print "reads ID error:  ectopic ID -- $readsid\n";}	          
              	
       }
        
    }
            
	  
  }
}

($optinfo,$optscore) = identify ($optinfo,$optscore,\%rds1,\%rds2,\%mapnum);            
## LAST OUTPUT...
print OUT1 "$optinfo" if defined $optinfo; 

unless(defined $optinfo){
  if (defined $bestrds1 ||$bestrds2){
      if($bestscore1 >= $bestscore2){
          my $optial_id = (split /\s+/,$bestrds1)[0];
          if (exists $mapnum{$optial_id}){
         	 	print OUT2 "$bestrds1\t$mapnum{$optial_id}\n";
         	}else{
         	 	print "warning: Not found the number of mapping location...\n";
         	}  
             #print OUT2 "$bestrds1\n";
       }else{
           my $optial_id = (split /\s+/,$bestrds2)[0];
           if (exists $mapnum{$optial_id}){
         	   print OUT2 "$bestrds2\t$mapnum{$optial_id}\n";
           }else{
           	 print "warning: Not found the number of mapping location...\n";
         	 }
       }  
  }
}

close IN;
close OUT1;
close OUT2;
###########==========================================#########
###########==========================================#########

sub maprec #recongnize mapped reads without recombination event...
{
	my ($mapchc,$mapscore);
	my $type;
	my $original;
	my $recbinfo;
	my $prior;
	my ($map1,$map2,$scor1,$scor2,$map_info) = @_; 
	my ($strand1,$chr1,$maplen1,$site1) = (split /\s+/,$map1)[1,2,3,4];
  my ($strand2,$chr2,$maplen2,$site2) = (split /\s+/,$map2)[1,2,3,4];   
  $mapscore = $scor1 + $scor2;
  #print "$map1\n$map2\n##$scor1,$scor2\n"; 
  #print "$strand1,$chr1,$maplen1,$site1\n$strand2,$chr2,$maplen2,$site2\n";

  if ($chr1 eq $chr2 ){
  	if ($chr1 =~ /$synchr/){
  		if(($site1 >= $synstar&&$site1<=$synend)&&($site2 >= $synstar&&$site2<=$synend)){
  		  $type = "syn_syn";
  	  }elsif(($site1 >= $synstar&&$site1<=$synend)||($site2 >= $synstar&&$site2<=$synend)){
  	    $type = "syn_wt";
  	  }else{$type = "wt_wt";}
  	
  	}else{$type = "wt_wt";}
    
    if (($strand1 eq $strand2)&& ( ($strand1 eq "-" && $site1>$site2 && abs($site1 -$site2 - $maplen2) <= $redis) || ($strand1 eq "+" && $site1<=$site2&& abs($site2 -$site1-$maplen1) <=  $redis) )  ){
        $original = 1;
      	$recbinfo = "NO";
        $prior = 100;        
    }else{
    	 if(exists $circhrinfo{$chr1}){
    	    if($strand1 eq $strand2){
    	    	 if(($site1>$site2 && $circhrinfo{$chr1} - $site1 +$site2 <= $maplen1 + $redis) || ($site1<=$site2&&$circhrinfo{$chr1} - $site2 +$site1 <= $maplen2 + $redis)){
    	          $original = 1;
      	        $recbinfo = "NO";
                $prior = 100;   
    	      }
    	    }
        }
       unless (defined $recbinfo){    	  
        $original = 2;
      	$recbinfo = "REC";
        $prior = 50;         
       }
    }
          # $mapchc .= "$map1\t$type\t1\tNO\n$map2\t$type\t1\tNO\n" ;                   	 

  }else{
  	if (($chr1 =~ /$synchr/ && ($site1 >= $synstar&&$site1<=$synend)) ||($chr2 =~ /$synchr/ && ($site2 >= $synstar&&$site2<=$synend))){
  		 $type = "syn_wt";
       $original = 4;
       $recbinfo = "REC";
       $prior = 5;          		   	
  	}else{
 		   $type = "wt_wt";
       $original = 3;
       $recbinfo = "REC";
       $prior = 25;      		 
  		}  	
  } 
  
  $mapscore = $mapscore*$prior;
  
  my $id1 = (split /\s+/,$map1)[0];
  my $id2 = (split /\s+/,$map2)[0];
  
  $mapchc = "$map1\t$type\t$original\t$recbinfo\t${$map_info}{$id1}\n$map2\t$type\t$original\t$recbinfo\t${$map_info}{$id2}" ;  
  
  return  ($mapchc,$mapscore);
}



sub identify
{
	my ($opt_info,$opt_score,$read1,$read2,$map_num) = @_;
	
	
	if ( %{$read1} && %{$read2}){
	   foreach (keys %{$read1}){
	       my $rds1id  = $_;
	       my $rds2id  = $_;
         $rds2id =~s/\/1$/\/2/;   
         
         next unless (exists ${$read2}{$rds2id}{info}); ##debug item...
         my $num1 = @{$read1->{$_}->{info}};
	       my $num2 = @{$read2->{$rds2id}->{info}};
	       foreach my $i(0..$num1-1){
	          foreach my $j (0..$num2-1){
	            my ($opt_info_tmp,$opt_score_tmp) = maprec (${$read1->{$rds1id}->{info}}[$i],${$read2->{$rds2id}->{info}}[$j],${$read1->{$rds1id}->{score}}[$i],${$read2->{$rds2id}->{score}}[$j],\%{$map_num});
	            if ($opt_score_tmp > $opt_score){
	            	 $opt_info = $opt_info_tmp ;
	               $opt_score = $opt_score_tmp ;
	            }
	          }
	          
	       } 
	   }
  
  }

  return ($opt_info,$opt_score);
}



##===end===##


