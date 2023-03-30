#! /usr/bin/perl -w
use strict;
use FindBin qw ($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;
#use Bio::SeqIO;

my $usage=<<USAGE;
Describe: SCB grams of SCRaMbLEd synthetic chrosome....
   
  Author: Yun Wang, wangyun\@genomic.cn
 Version: 2.0, Date: 2023-3-17

Usage :perl $0 [option]
        -prefix    <str> :  prefix of output file [default:test]
        -Xencod    <str> :  encod order of x-axis 
        -Yencod    <str> :  encod order of Y-axis
        -outdir          :  output directory [default: ./]
        -help            :  show this help message

USAGE

my ($prefix,$Yencod,$Xencod,$outdir,$help);
GetOptions(
    "prefix:s"      => \$prefix,
    "Yencod:s"      => \$Yencod,
    "Xencod:s"      => \$Xencod,
    "outdir:s"     => \$outdir,
    "help"       => \$help,
);

$prefix ||= "test";
$outdir ||= getcwd;

die $usage if (!$Yencod || !$Xencod);

my $encod = $Xencod;
open ENC, "$encod" || die $!;

my @Xordcods;
my $Xid;

while (<ENC>){
   chomp;

   my ($id,$encod);
   my @set = split;
   
   if (@set == 3){
	    ($id,$encod) = ($set[0],$set[2]);
   }else{
	    ($id,$encod) = ($set[0],$set[1]);
	 }
    
    $Xid = $id;
    @Xordcods = split /,/,$encod;
}

my $max_segm = @Xordcods;
close ENC;

#################
$encod = $Yencod;

open ENC, "$encod" || die $!;
####Start to draw svg figure.....
my $intseg = 8;
my $x_border = 30;
my $y_border = 30;

while (<ENC>){
   chomp;
   my $posY = $y_border;
   my ($id,$encod);
   my @set = split /\s+/;
   
   if (@set == 3){
	    ($id,$encod) = ($set[0],$set[2]);
   }else{
	    ($id,$encod) = ($set[0],$set[1]);
	 }
   
   my @ordcods = split /,/,$encod;
   
   open OUT, ">$outdir/$id.$prefix.svg" ||die $!;
   my $width = $max_segm  * $intseg +  $x_border + $x_border;
   my $height =  @ordcods  * $intseg +  $y_border + $y_border;
   print OUT "<svg height=\"$height\" width=\"$width\" fill = \"#e0ffff\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

   
   my $posX =  ($max_segm * $intseg + $x_border + $y_border)/2;
   Drawscript($posX, $posY, "$id", "left");
   
   ###Drawaxis
   my $x1 = $x_border  +  20;
   my $x2 = $max_segm * $intseg + $x1;
   $posY = @ordcods * $intseg + $y_border;
   
   print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x1\" y1=\"$posY\" x2=\"$x2\" y2=\"$posY\" />\n";
   print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x1\" y1=\"$y_border\" x2=\"$x1\" y2=\"$posY\" />\n";

   my $textY = $posY + 24;
   Drawscript($posX,$textY, "$Xid", "left");
  
   $textY = $posY + 12;
   foreach my $segm (1..$max_segm){
     my $y1 = $posY;
     my $y2 = $y1 + 4;
     my $x = $segm * $intseg - $intseg/2 + $x1;
     my $X_segm = $Xordcods[$segm-1];
     
     print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x\" y1=\"$y1\" x2=\"$x\" y2=\"$y2\" />\n";
     print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"7\" font-weight=\"normal\" text-anchor=\"middle\" x=\"$x\" y=\"$textY\">$X_segm</text>\n";
   }
   
   my $n = 1;
   my $textX = $x1 - 16;
   my ($tempx,$tempy);
   foreach my $segm (@ordcods){
     my $y = $posY - $n * $intseg + $intseg/2;
     my $x2 = $x1;
     my $x1 = $x1 - 4;
     
     my $textY = $y + 2;
     print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x1\" y1=\"$y\" x2=\"$x2\" y2=\"$y\" />\n";
     print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"7\" font-weight=\"normal\" text-anchor=\"left\"  x=\"$textX\" y=\"$textY\">$segm</text>\n";
     
     $tempy = $y;
     
   
     my $_check_point = 0;
     foreach my $k (0..@Xordcods-1){         
           if (abs ($Xordcods[$k]) == abs ($segm) ){
               $_check_point = 1;
               my $x = abs($k+1) * $intseg - $intseg/2 +  $x_border  +  20;
               if ($Xordcods[$k] * $segm > 0){
                   print OUT "<circle style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 0; stroke-linecap: square\"  cx=\"$x\" cy=\"$y\" r=\"2.5\" />\n";
           
               }else{
                  print OUT "<circle style=\"fill: rgb(30,144,225); fill-opacity: 1.0; stroke: rgb(30,144,225); stroke-opacity: 1.0; stroke-width: 0; stroke-linecap: square\"  cx=\"$x\" cy=\"$y\" r=\"2.5\" />\n";
               }
        
              $tempx = $x;
           }
     }
        
     unless ($_check_point == 1){
   
        my $x = $tempx + $intseg/2;
        my $y = $y;
        my $x1 = $x-2;
        my $x2 = $x+2;
        my $y1 = $y-2;
        my $y2 = $y+2;
        print OUT "<line style=\"fill: rgb(220,20,60); fill-opacity: 1.0; stroke: rgb(220,20,60); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x1\" y1=\"$y\" x2=\"$x2\" y2=\"$y\" />\n";
        print OUT "<line style=\"fill: rgb(220,20,60); fill-opacity: 1.0; stroke: rgb(220,20,60); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x\" y1=\"$y1\" x2=\"$x\" y2=\"$y2\" />\n";
     }
     
     $n ++;
   }
   
print OUT "</svg>\n";

}

### draw description #####
sub Drawscript{
   my ($posx,$posy,$script,$anchor) = @_;
   if (defined $anchor){
     print OUT "<text fill=\"rgb(0,0,0)\" font=\"Calibri\" font-size=\"9\" font-weight=\"normal\" text-anchor=\"$anchor\" x=\"$posx\" y=\"$posy\">$script </text>\n";
   }else{
     print OUT "<text fill=\"rgb(0,0,0)\" font=\"Calibri\" font-size=\"9\" font-weight=\"normal\" x=\"$posx\" y=\"$posy\">$script </text>\n";    
   }
}

=head
##### draw coordinate axis and background#####
sub Drawaxis{
  my ($pos,$starx,$endx)= @_;
  my $starym = $pos - 4;
  my $endym = $pos + 4;
  my $starys = $pos - 2;
  my $endys = $pos + 2;
  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starx\" x2=\"$endx\" y1=\"$pos\" y2=\"$pos\" />\n";
  #print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starx\" x2=\"$starx\" y1=\"$starym\" y2=\"$endym\" />\n";
  #print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$endx\" x2=\"$endx\" y1=\"$starym\" y2=\"$endym\" />\n";

  
  my $chrpos = 0;
  my $i;
  while ($chrpos < $synreflen){
     my $starxm = $chrpos/$pixelbase + $border_distance;
     my $TextX = $starxm;
     my $TextY = $pos - 8;
     $i = $chrpos/1000; 
     my $ik ="$i"."k";
      
     print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$starym\" y2=\"$endym\" />\n"; 
     print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"7\" font-weight=\"normal\" text-anchor=\"middle\"  x=\"$TextX\" y=\"$TextY\">$ik </text>\n";
     if ($chrpos == 0){
       $chrpos += 50000 ;
       next; 
     }
      
     foreach my $j(1..4){
        my $starxs = (($i*1000 - 50000)+ $j*10000)/$pixelbase + $border_distance;
        print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
        }
     
     $chrpos += 50000;
     }
     
  my $m = int (($synreflen%50000)/10000) ;
  
  foreach my $k (1..$m){
     my $starxs = ($i*1000+ $k*10000)/$pixelbase + $border_distance;
     print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
   }
}	

=cut
