#! /usr/bin/perl -w

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;
Describe: This script is used to draw a standard reads mapping svg figure for 
          validation of the reconstructed reference. .
   
  Author: Yun Wang, wangyun\@genomic.cn
 Version: 1.0, Date: 2014-04-01

Usage :perl $0 [-option]

        -loxpftr    [file]      :  loxp feature of start and end point  
        -cvg        [file]      :  average sequence coverage for appoint region 
        -rpe        [file]      :  abnormal mapped split reads (*.rpe.bkpt.sort.gvr) in reference chromosome
        -rse        [file]      :  single split END mapping (*.rse.bkpt.sort.gvr) 
        -spr                    :  pair_end mapping info from normal and abnormal SOAP results --mapping on synthetic reference
        -n                      :  prefix of output file [defacult:test]
        -isz                    :  insert size [default:500] 
        -refid                  :  id of synthetic reference [defacult:IXR_BACseq]
        -reflen                 :  length of synthetic chromosome [defacult:100371]
        -o                      :  output directory [default: ./]
        -help                   :  show this help message


USAGE

my ($samplename,$insertsize,$refid,$reflen,$loxpftr,$cvg,$gvr,$spr,$rse,$outdir,$help);
GetOptions(
    "n:s"     => \$samplename,
    "isz:s"     => \$insertsize,
    "refid:s"      => \$refid,
    "reflen:s"     => \$reflen,
    "loxpftr:s"    => \$loxpftr,
    "cvg:s"        => \$cvg,  
    "rpe:s"        => \$gvr,
    "spr:s"        => \$spr,
    "rse:s"        => \$rse,
    "o:s"     => \$outdir,
    "help"       => \$help,
);

die $usage if (!$refid||!$cvg||!$gvr||$help);

$samplename ||='test';
$outdir ||=getcwd;
$refid ||= "IXR_BACseq";
$insertsize ||= 500;
my $version = "v1.0";
my $software = basename ($0); 
my $scal_fig = 10;

print "#$samplename\n";
print "#drawing SVG\n";

open CVG, $cvg ||die $!;
if (defined $spr){
  open SPR, $spr ||die $! 
}

my $height = estimateY($gvr,60);

$height = $height + 50;
print "###$height\n";

my $width = int($reflen/10)+200;

open GVR, $gvr ||die $!;
$outdir ||= getcwd;
open OUT, ">$outdir/$samplename.svg" ||die $!;
print OUT "<svg height=\"$height\" width=\"$width\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

Drawscript(5,35,"Coverage");
while (<CVG>){
	chomp;
	my ($site,$depth) = (split /\s+/,$_)[1,3];
  Drawdpline(50,$site,$depth);
}
my $axipos = 58;
Drawscript(5,54,"$refid");

Drawaxis ($axipos,$reflen,$height,1);

my $loxpY = 54;
if (defined $loxpftr){

#   Drawscript(5,60,"Loxp site");
 	open FTR, $loxpftr ||die $!;
  while (<FTR>){
	   next if $_ =~ /^#/;
	   my ($star,$END) = (split /\s+/,$_)[3,4];
	 
	   Drawloxp($star,$loxpY);
	}
}

my $scripty = $loxpY + 50;
Drawscript(5,$scripty,"Split reads");
my @recpos;
my $posy = $scripty + 30;
my $suposy = 0;
my $max;
my ($direct,$lchrid,$lstar,$lend,$lori,$rchrid,$rstar,$rend,$rori);
while(1){
  my $lingvr = <GVR>;
  last unless defined $lingvr;
  chomp $lingvr;
  ($direct,$lchrid,$lstar,$lend,$lori,$rchrid,$rstar,$rend,$rori) = (split /\s+/,$lingvr)[0,2,3,4,5,6,7,8,9];
   if ($lchrid eq $refid && $rchrid eq $refid){
     Drawsplitrds($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy);
     $max = max($lstar,$lend,$rstar,$rend,);
     push @recpos, [$max,$posy];
     last;
   }elsif($lchrid eq $refid){
     my ($tfrstar,$tfrend,$tfscript);
   
     $tfrstar = $lend + 100 * $scal_fig;
     $tfrend = $rend - $rstar+1 + $tfrstar;
     $tfscript = "$rchrid"."_$rstar"." $rstar $rend $rori";
     Drawsplitrds($direct,$lstar,$lend,$lori,$tfrstar,$tfrend,$rori,$posy);
     
     print "chc A..$tfscript\n"; 
     my $maxtf =  max($tfrstar,$tfrend); 
     $max = max($lstar,$lend,$tfrstar,$tfrend);
     push @recpos, [$max,$posy];
     
     $maxtf =  $maxtf/$scal_fig + 100;  	    
     Drawscript($maxtf+14,$posy+8,"$tfscript","star");
     last;
   }elsif($rchrid eq $refid){
    	my ($tflstar,$tflend,$tfscript);
    	#if ($rori eq "+"){    	 	  	
    	$tflend =  $rstar - 100* $scal_fig;
    	$tflstar = $tflend - ($lend-$lstar+1) ;     	      
    	#}
    	$tfscript = "$lchrid"."_$lend"." $lstar $lend $rori";  
    	 	 
    	Drawsplitrds($direct,$tflstar,$tflend,$lori,$rstar,$rend,$rori,$posy); 
    	print "chc B..$tfscript\n"; 
    	my $mintf = min($tflstar,$tflend);
    	$max = max($tflstar,$tflend,$rstar,$rend);
      push @recpos, [$max,$posy];
      
    	$mintf =  $mintf/$scal_fig + 100;      	    
    	Drawscript($mintf-8,$posy+8,"$tfscript","end");
      last;
   }
 }

while (<GVR>){
	chomp;
	my $line = $_;
	($direct,$lchrid,$lstar,$lend,$lori,$rchrid,$rstar,$rend,$rori)= (split /\s+/,$line)[0,2,3,4,5,6,7,8,9];
  if ($lchrid eq $refid && $rchrid eq $refid){
     $posy += 16;
     my $minpos = min ($lstar,$lend,$rstar,$rend);
     my $maxpos = max ($lstar,$lend,$rstar,$rend);
     if($minpos <$max){
        Drawsplitrds($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy);
        push @recpos,[$maxpos,$posy];
        $max = $maxpos;	
     }else{
        my @recposR = reverse @recpos;
        foreach my $record (@recposR){
          if(${$record}[0] < $minpos){
            $posy= ${$record}[1];
          }else{
            last;
          }
        }
       #@recpos = ();
       Drawsplitrds($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy);
       push @recpos,[$maxpos,$posy];
       $max = $maxpos;
      }
     $suposy = $posy if $posy > $suposy; 
  }elsif($lchrid eq $refid){
     my ($tfrstar,$tfrend,$tfscript);
     $posy += 16;
     $tfrstar = $lend + 100 * $scal_fig;
     $tfrend = $rend - $rstar+1 + $tfrstar;
     
     my $minpos = min ($lstar,$lend,$tfrstar,$tfrend);
     my $maxpos = max ($lstar,$lend,$tfrstar,$tfrend);
     
      $tfscript = "$rchrid"."_$rstar"." $rstar $rend $rori";
      if($minpos <$max){
        print "chc C..$tfscript\n";
        Drawsplitrds($direct,$lstar,$lend,$lori,$tfrstar,$tfrend,$rori,$posy); 
        my $maxtf =  max($tfrstar,$tfrend);
        $maxtf =  $maxtf/$scal_fig + 100;
        Drawscript($maxtf+14,$posy+8,"$tfscript","star");
        push @recpos,[$maxpos,$posy];
        
        $max = $maxpos;	
     }else{
  	    my @recposR = reverse @recpos;
  	    foreach my $record (@recposR){
  	 	    ${$record}[0] < $minpos?$posy= ${$record}[1]:last; 
  	    }
  	    #@recpos = ();	  
        print "chc D..$tfscript\n";
        Drawsplitrds($direct,$lstar,$lend,$lori,$tfrstar,$tfrend,$rori,$posy); 
        push @recpos,[$maxpos,$posy];
        my $maxtf =  max($tfrstar,$tfrend);
        $maxtf =  $maxtf/$scal_fig + 100;
        Drawscript($maxtf+8,$posy+8,"$tfscript","star");
        $max = $maxpos;
      }
     $suposy = $posy if $posy > $suposy;
     
   }elsif($rchrid eq $refid){
    	my ($tflstar,$tflend,$tfscript);
    	#if ($rori eq "+"){
    	$tflend =  $rstar - 100* $scal_fig;
    	$tflstar = $tflend - ($lend-$lstar+1);
    	#}
    	  
      my $minpos = min ($tflend,$tflstar,$rstar,$rend);
      my $maxpos = max ($tflend,$tflstar,$rstar,$rend);
    	
    	$tfscript = "$lchrid"."_$lend"." $lstar $lend $rori"; 
    	if($minpos <$max){
    	  
    	  Drawsplitrds($direct,$tflstar,$tflend,$lori,$rstar,$rend,$rori,$posy); 
      	my $mintf = min($tflstar,$tflend); 
      	$mintf =  $mintf/$scal_fig + 100; 
    	  Drawscript($mintf-8,$posy+8,"$tfscript","end");
        print "chc E..$tfscript\n";
        $max = $maxpos;
        push @recpos, [$max,$posy];
      }else{
        my @recposR = reverse @recpos;
  	    foreach my $record (@recposR){
  	 	    ${$record}[0] < $minpos?$posy= ${$record}[1]:last; 
  	    }
  	    #@recpos = ();	 
        print "chc F..$tfscript\n";
        Drawsplitrds($direct,$tflstar,$tflend,$lori,$rstar,$rend,$rori,$posy); 
        push @recpos,[$maxpos,$posy];
        my $mintf = min($tflstar,$tflend); 
        $mintf =  $mintf/$scal_fig +100;
    	  Drawscript($mintf-8,$posy+8,"$tfscript","end");
        $max = $maxpos;
      }
      $suposy = $posy if $posy > $suposy;
   }
}


print "The maximal y-axis coordinate of splitted reads:$suposy\n";

if($suposy <200){
   $suposy = 100;
}

$posy = $suposy;
$axipos = $posy  + 30;

#$axipos = $height - 50;

Drawaxis ($axipos,$reflen,$height,0);
$max = 0;
my @recpose;

$posy =  $axipos;

$suposy = $posy + 20;
if (defined $rse){
   open RSE,"$rse" || die $!;
   while (<RSE>){
     chomp;
     my ($rdsid,$chrid,$star,$ends,$ori) =(split /\s+/,$_)[1,2,3,4,5];
     if ($chrid eq $refid){
       my $minsite = min($star,$ends);
       my $maxsite = max($star,$ends);
       my ($tfstar,$tfends,$rds_info);
       
       if ($minsite > $max){
         $posy = $suposy;
       }else{
         $posy += 16;
       }
        $max = $maxsite if $maxsite > $max;
        
       if ($rdsid =~ /\/1$/){
          $rds_info = 1;
       }elsif($rdsid =~ /\/2$/){
          $rds_info = 2;
       }
       Drawrse($direct,$star,$ends,$ori,$rds_info,$posy);
     }
  }
}
print "The maximal y-axis coordinate of splitted reads: $posy\n";
=head
#$axipos = ($height-450)+80;
$axipos = $suposy + 100;
Drawaxis ($axipos,$reflen,$height,0);

$posy = $axipos - 8;
$scripty = $axipos - 70;
Drawscript(5,$scripty," ");
=cut

if (defined $spr){
  while(<SPR>){
    chomp;
    my ($site1,$ori1,$site2,$ori2) = (split /\s+/,$_)[2,3,6,7];
    my $ry = 50;
    Drawsprds($site1,$ori1,$site2,$ori2,$posy,$ry);
  }
}

####################################

$scripty = $height - 20;
my $loctime = scalar (localtime); 
Drawscript(25,$scripty,"## Created by $software $version at $loctime.");
print OUT "</svg>";

###======sub function=========##### 

### draw description #####
sub Drawscript{
	 my ($posx,$posy,$script,$anchor) = @_;
	 if (defined $anchor){
	   print OUT "<text fill=\"rgb(47,79,79)\" font=\"Arial\" font-size=\"14\" font-weight=\"regular\" text-anchor=\"$anchor\" x=\"$posx\" y=\"$posy\">$script </text>\n";
	 }else{
	   print OUT "<text fill=\"rgb(47,79,79)\" font=\"Arial\" font-size=\"14\" font-weight=\"bold\" x=\"$posx\" y=\"$posy\">$script </text>\n";    
   }
}
 
### draw sequence coverage #### 
sub Drawdpline{
	my ($starY,$starX,$depth) = @_;
	$starX = ($starX/10)+100;
	my $endY =$starY*(1-$depth);  
	print OUT "<line style=\"fill: rgb(0,0,255); fill-opacity: 1.0; stroke: rgb(0,0,255); stroke-opacity: 1.0; stroke-width: 0.5; stroke-linecap: square\" x1=\"$starX\" x2=\"$starX\" y1=\"$starY\" y2=\"$endY\" />\n";
}

##### draw coordinate axis and background#####
sub Drawaxis{
	my ($pos,$seqlen,$height,$oper)= @_;
	my $starym = $pos - 4;
	my $endym = $pos + 4;
	my $starys = $pos - 2;
	my $endys = $pos + 2;
	my $codin = 0;
	my $mx = int($seqlen/10) + 120;
	my $mxa = $mx - 4;
	
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"$mx\" y1=\"$pos\" y2=\"$pos\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$starym\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$endym\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$mx\" x2=\"$mxa\" y1=\"$pos\" y2=\"$starym\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$mx\" x2=\"$mxa\" y1=\"$pos\" y2=\"$endym\" />\n";
			
 
	my $n = int ($seqlen/1000); 

	foreach my $i(0..$n){
		 my $starxm = $i*100+100;
		 my $TextX = $starxm-5;
		 my $TextY = $pos +16;
		 my $ik ="$i"."k";
		   
		 print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$starym\" y2=\"$endym\" />\n"; 
		 print OUT "<line style=\"fill: rgb(211,211,211); fill-opacity: 1.0; stroke: rgb(211,211,211); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$pos\" y2=\"$height\" />\n" if $oper == 1;
	   print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"11\" font-weight=\"normal\" x=\"$TextX\" y=\"$TextY\">$ik </text>\n";
	       last if ($i== $n);
	       foreach my $j(1..9){
	       	  my $starxs = (($i*1000)+ $j*100)/10 + 100;
	       	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
		        print OUT "<line style=\"fill: rgb(224,255,255); fill-opacity: 1.0; stroke: rgb(224,255,255); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$pos\" y2=\"$height\" />\n"if $oper == 1;
	         }
	   }
  my $m = int (($seqlen%1000)/100) + 1;
  foreach my $k (1..$m){
	   my $starxs = (($n*1000)+ $k*100)/10 + 100;
	   print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
		 print OUT "<line style=\"fill: rgb(224,255,255); fill-opacity: 1.0; stroke: rgb(224,255,255); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$pos\" y2=\"$height\" />\n"if $oper == 1;
	 }       	
}	

#### Draw Loxp feature ####
sub Drawloxp {
	 my ($site,$posy) = @_;
	 my $posx = $site/10+100;
	 print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
 }

### Draw structure variation ###  
sub Drawvariate {
	 my @para = @_;
	 my ($star,$END,$tstar,$tend) = (0,0,0,0);
	 my ($type,$direc); 
	 if (@para == 4){
	 	  ($posy,$type,$star,$END) = @para;
	 	  my $posx = $star/10+100;
	 	  my $width = ($END - $star + 1)/10;	 	  
	 	  if ($type eq "DEL"){	 	  	
	 	  	print OUT "<rect style=\"fill: rgb(105,105,105); fill-opacity: 1.0; stroke: rgb(105,105,105); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"$width\" x=\"$posx\" y=\"$posy\" />\n";
	 	  }elsif($type eq "INV"){
	 	  	print OUT "<rect style=\"fill: rgb(165,42,42); fill-opacity: 1.0; stroke: rgb(165,42,42); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"$width\" x=\"$posx\" y=\"$posy\" />\n";
      }else{print "lsv format error\n";}
   }else{
   	  ($posy,$type,$star,$END,$tstar,$tend,$direc) = @para;
   	  
   	  my $posx = $star/10+100;
	 	  my $width = ($END - $star + 1)/10;
	 	  my ($pathsx,$pathex) = (0,0);
	 	  
      if($END <= $tstar){
	 	   	 $pathsx = $END/10+100;
	 	     $pathex = $tstar/10+100;
	 	  }else{	 	     
	 	   	 $pathsx = ($star/10)+100;
	 	     $pathex = ($tstar/10)+100;
	 	   	}
	 	   	
	 	  
	 	  my $script = "$direc";
	 	  my $arrowl = $pathex-3;
	 	  my $arrowr = $pathex+3;
	 	  my $posym = $posy + 4;
	 	  my $posymd = $posy + 2;
	 	  my $posyd = $posy + 8;
	 	  my $posyu = $posy - 2;
	 	  my $rgb;
	 	  
	 	  if($type eq "TDUP"){
	 	   $rgb = "25,25,112";#midnightblue
	 	   
	 	  }elsif($type eq "NTDUP"){
	 	   $rgb = "255,165,0";# orange
	 	  }elsif($type eq "D&I"){
	 	   $rgb = "148,0,211";#purple
	 			   		 	      	
      }else{
      print "Weird structure variationn !!\n"	
     }
      print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"$width\" x=\"$posx\" y=\"$posy\" />\n";
      print OUT "<line style=\"stroke-dasharray: 2; fill: none; fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathsx\" x2=\"$pathex\" y1=\"$posym\" y2=\"$posym\" />\n";
      print OUT "<line style=\"fill: rgb(25,25,112); fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathex\" x2=\"$pathex\" y1=\"$posyu\" y2=\"$posyd\" />\n";
     	print OUT "<line style=\"fill: rgb(25,25,112); fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathex\" x2=\"$arrowl\" y1=\"$posyd\" y2=\"$posymd\" />\n";
      print OUT "<line style=\"fill: rgb(25,25,112); fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathex\" x2=\"$arrowr\" y1=\"$posyd\" y2=\"$posymd\" />\n";
      print OUT "<text fill=\"rgb(119,136,153)\" font=\"Arial\" font-size=\"15\" font-weight=\"bold\" x=\"$arrowl\" y=\"$posyu\">$script </text>\n";    
     
  }
}

##### Draw loxp Reads ######
 
sub Drawloxprds{
	  my ($site,$posy,$fowrds,$revrds) = @_;
	  my $posx = $site/10+100;
	  
	  if($fowrds > 0){
	  my $posxs = $posx - 6;
	  my $posxe = $posx + 8;
	  my $posxa = $posxe - 3;
	  my $posyu = $posy + 1;
	  my $posyl = $posy + 4;
		my $posyd = $posy + 7;
	  my $posxt = $posx + 9;
	  my $posyt = $posy + 8;
	  
	  print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxe\" y1=\"$posyl\" y2=\"$posyl\" />\n";
	  print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxe\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyu\" />\n";
    print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxe\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyd\" />\n";
    print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
    print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"10\" font-weight=\"normal\" x=\"$posxt\" y=\"$posyt\">$fowrds </text>\n";    
  }
###    
    if ($revrds>0){
    $posy = $posy + 12;	
    my $posxs = $posx - 8;
    my $posxa = $posxs + 3;
	  my $posxe = $posx + 6;
	  my $posyu = $posy + 1;
	  my $posyl = $posy + 4;
	  my $posyd = $posy + 7;
	  my $posxt = $posx + 9;
	  my $posyt = $posy + 8;
    print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxe\" y1=\"$posyl\" y2=\"$posyl\" />\n";
	  print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyu\" />\n";
    print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyd\" />\n";
    print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
    print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"10\" font-weight=\"normal\" x=\"$posxt\" y=\"$posyt\">$revrds </text>\n";    
   }
}

##### Draw splitted loxp reads ###
sub Drawrds {
	my ($star,$loxpend,$ori,$posy) = @_;
  my ($starx,$endx,$arxo,$arx,$aruy,$ardy); 
  
  $aruy = $posy - 3;
  $ardy = $posy + 3;	
	if ($star>$loxpend){
		  $starx = ($loxpend/10) + 100;
		  $endx = ($star/10 +100) + 10;
	}else{       
		  $starx = ($star/10+100) - 10;
		  $endx = ($loxpend/10) + 100;
	}
	
	if ($ori eq "+"){
		  $arxo = $endx;
		  $arx = $arxo - 3;
		}else{
		  $arxo = $starx;
	    $arx = $arxo + 3;
	 }

#Only uesd for regulate picture!!
	 if ($star>$loxpend && $ori eq "-"){
	 	   $arxo = $arxo + 3;
	 	   $arx = $arxo + 3;
	 	   }  
##END
	 
	print OUT "<line style=\"fill: rgb(75,0,130); fill-opacity: 1.0; stroke: rgb(75,0,130); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starx\" x2=\"$endx\" y1=\"$posy\" y2=\"$posy\" />\n";
  print OUT "<line style=\"fill: rgb(75,0,130); fill-opacity: 1.0; stroke: rgb(75,0,130); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arxo\" x2=\"$arx\" y1=\"$posy\" y2=\"$aruy\" />\n";
  print OUT "<line style=\"fill: rgb(75,0,130); fill-opacity: 1.0; stroke: rgb(75,0,130); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arxo\" x2=\"$arx\" y1=\"$posy\" y2=\"$ardy\" />\n";
}

	  
sub Drawsplitrds{
	 my ($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy) = @_;
	 my $rgb  = "0,0,0";
	 if($direct eq "R"){
	 	  $rgb = "165,42,42";
	 	}else {$rgb = "106,90,205";}
	 
	 my $loxpl = $lend/10 +100;
	 my $loxpr = $rstar/10 +100; 
	 
	 print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"1.5\" x=\"$loxpl\" y=\"$posy\" />\n";
	 print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"1.5\" x=\"$loxpr\" y=\"$posy\" />\n";
	 $posy = $posy + 4;
	 if ($direct eq "F"){
	 print OUT "<line style=\"stroke-dasharray: 2; fill: none; fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpl\" x2=\"$loxpr\" y1=\"$posy\" y2=\"$posy\" />\n";
	}else{
	 print OUT "<line style=\"stroke-dasharray: 2; fill: none; fill-opacity: 1.0; stroke: rgb(255,192,203); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpl\" x2=\"$loxpr\" y1=\"$posy\" y2=\"$posy\" />\n";	
		}
	 Drawrds ($lstar,$lend,$lori,$posy);
	 Drawrds ($rend,$rstar,$rori,$posy);
	} 
	
### estimate the maximal y-axis coordinate ####
sub estimateY{
  my ($gvr,$posy) = @_;
  #print "$posy\n";  
 
  my $posysv = 0;
  open GVR, $gvr ||die $!;
  my @recpos;
  my $max = 0;
   my $suposy = $posy; 
  while (<GVR>){
	chomp;
	my $line = $_;
	($direct,$lchrid,$lstar,$lend,$lori,$rchrid,$rstar,$rend,$rori)= (split /\s+/,$line)[0,2,3,4,5,6,7,8,9];
  if ($lchrid eq $refid && $rchrid eq $refid){
     $posy += 16;
     my $minpos = min ($lstar,$lend,$rstar,$rend);
     my $maxpos = max ($lstar,$lend,$rstar,$rend);
     if($minpos <$max){
        push @recpos,[$maxpos,$posy];
        $max = $maxpos;	
     }else{
  	    my @recposR = reverse @recpos;
  	    foreach my $record (@recposR){
  	 	    if (${$record}[0] < $minpos){
  	 	       $posy = ${$record}[1];
  	 	    }else{
  	 	      last;
  	 	    }
  	    }
  	    #@recpos = ();
  	    push @recpos,[$maxpos,$posy];
  	    $max = $maxpos;
      }
     $suposy = $posy if $posy > $suposy; 
  }elsif($lchrid eq $refid){
     my ($tfrstar,$tfrend,$tfscript);
     $posy += 16; 
     $tfrstar = $lend + 100 * $scal_fig;
     $tfrend = $rend - $rstar+1 + $tfrstar; 
     
     my $minpos = min ($lstar,$lend,$tfrstar,$tfrend);
     my $maxpos = max ($lstar,$lend,$tfrstar,$tfrend);
     
      if($minpos <$max){          
        push @recpos,[$maxpos,$posy];
        $max = $maxpos;	
     }else{
  	    my @recposR = reverse @recpos;
  	    foreach my $record (@recposR){
  	 	    if(${$record}[0] < $minpos){
  	 	      $posy= ${$record}[1];
  	 	    }else{last;} 
  	    }
  	    #@recpos = ();	          
        push @recpos,[$maxpos,$posy]; 
        $max = $maxpos;
      }
     $suposy = $posy if $posy > $suposy; 
     
   }elsif($rchrid eq $refid){
    	my ($tflstar,$tflend,$tfscript);
    	#if ($rori eq "+"){    	 	  	
    	$tflend =  $rstar - 100* $scal_fig;
    	$tflstar = $tflend - ($lend-$lstar+1);
    	#}
    	     
      my $minpos = min ($tflend,$tflstar,$rstar,$rend);
      my $maxpos = max ($tflend,$tflstar,$rstar,$rend);
    	    
    	if($minpos <$max){
        $max = $maxpos;
        push @recpos, [$max,$posy];
      }else{
        my @recposR = reverse @recpos;
  	    foreach my $record (@recposR){
  	 	    if(${$record}[0] < $minpos){
  	 	      $posy= ${$record}[1];
  	 	    }else {last;} 
  	    }
  	    #@recpos = ();
        push @recpos,[$maxpos,$posy]; 
        $max = $maxpos;
      }
      $suposy = $posy if $posy > $suposy;
   }
  }
 
  
 if (defined $rse){
   open RSE,"$rse" || die $!;
   my $max = 0;
   my $maxposy = $suposy;
   while (<RSE>){
     chomp;
     my ($rdsid,$chrid,$star,$ends,$ori) =(split /\s+/,$_)[1,2,3,4,5];
     if ($chrid eq $refid){
       my $minsite = min($star,$ends);
       my $maxsite = max($star,$ends);
       my ($tfstar,$tfends,$rds_info);
       
       if ($minsite > $max){
         $posy = $suposy;
       }else{
         $posy += 16;
       }
        $max = $maxsite if $maxsite > $max;
        
       if ($rdsid =~ /\/1$/){
          $rds_info = 1; 
       }elsif($rdsid =~ /\/2$/){
          $rds_info = 2;
       }
       $maxposy = $posy if  $posy >  $maxposy;
     }	
  }
  
  $suposy = $maxposy;
  
}

close RSE;
 return $suposy;
}
    	
sub Drawsprds{
	my ($site1,$ori1,$site2,$ori2,$posy,$ry,$rgb) = @_;
	my $min = min($site1,$site2);
	my $max = max($site1,$site2);
	my $rgba;
	if($ori1 eq $ori2){
	 	$rgba = "165,42,42";
	}elsif(($max-$min)>=500){$rgba = "106,90,205";
	}else{$rgba = "128,128,128";}
   
  $rgb ||=	$rgba;
	
	my ($starx,$endx,$rx); 
	 $starx = ($min/10) + 100;
	 $endx = ($max/10) + 100;
	 $rx = ($endx-$starx)/2;
	 print OUT "<path d=\"M$starx,$posy,A$rx,$ry,0,0,1,$endx,$posy\" style=\"fill: none; fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 0.8; stroke-linecap: square\"/>\n" if ($max-$min)>=500;
}


sub Drawrse{
	 my ($direct,$star,$ends,$ori,$rdsinfo,$posy) = @_;
	 my $rgb  = "0,0,0";
   my $bpk_site;
   my $tfend = $ends/10 +100;
   my $tfstar = $star/10 +100;
   if($rdsinfo == 1){	 
	   $bpk_site = $tfend;
	 }elsif($rdsinfo == 2){
	   $bpk_site = $tfstar;
	 }
	 
	 print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"1.5\" x=\"$bpk_site\" y=\"$posy\" />\n";

	 $posy = $posy + 4;

   if($rdsinfo == 1){	 
	    Drawrds ($star,$ends,$ori,$posy);
	 }elsif($rdsinfo == 2){
	   Drawrds ($ends,$star,$ori,$posy);  
	 }
}
