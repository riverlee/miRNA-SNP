#!/usr/bin/perl
############################################################
#Center for Quantitative Science, Vanderbilt University
#System Bioinformatics Engineer
#Jiang Li
#riverlee2008@gmail.com
#############################################################
use strict;
use warnings;

#define variants
my $isdownload=0; #have already downloaded the neccessary file
my $miRNAdatFile="miRNA.dat";
my $premiRNAgffFile="hsa.gff";
my $premiRNA2miRNAFile="human_pre-miRNA_to_miRNA.txt";
my $dbSNPFile="/2TB/reference/dbsnp/vcf/00-All.vcf"; #dbsnp135

my $outfile="miRNA_SNP_details.txt";

my %premiRNA; #Store premiRNA's information
my %position2premiRNA; #key are genome positions, chr->position->premiRNA = 1


#If you don't have downloaded the mirbase data and dbsnp vcf file
if($isdownload){
    downloadAndTar();
}



#1) parsing premiRNA gff file
print info(),"Parsing pre-miRNA gff file \n";
parsingPremiRNA($premiRNAgffFile,\%premiRNA,\%position2premiRNA);

#2) parsing miRNA.dat file to get the pre-miRNA mature miRNA relationship and premiRNA's sequence
#Now the %premiRNA has more information with two more keys 'mature' and 'seq'
print info(),"Parsing miRNA.dat file \n";
parsingmiRNAdat($miRNAdatFile,$premiRNA2miRNAFile,\%premiRNA);

#3) Read dbSNP vcf file and determine whether this SNPs located in the pre-miRNA
open(IN,$dbSNPFile) or die $!;
open(OUT,">$outfile") or die $!;
print OUT join "\t",("#pre-miRNA ID","pre-miRNA ACC","chr","start","end","strand","SNP ID","SNP Loc","SNP ref","SNP alt","type","details","#miRNAs","mature miRNAs","","premiRNA Seq\n");
while(<IN>){
    next if (/^#/);
    s/\r|\n//g;
    my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info) = split "\t";
    next if (! exists($position2premiRNA{$chr}->{$pos}));  #only those SNPs whose position are in the pre-miRNAs
    foreach my $premirna (sort keys %{$position2premiRNA{$chr}->{$pos}}){
        my $premiRNARef = $premiRNA{$premirna};
        my ($type,$details) = getSNPmiRNAType($pos,$premiRNARef);
        my @matures = sort keys %{$premiRNA{$premirna}->{'mature'}};
        print OUT join "\t",($premirna,$premiRNA{$premirna}->{'ACC'},
                            $premiRNA{$premirna}->{'chr'},
                            $premiRNA{$premirna}->{'start'},
                            $premiRNA{$premirna}->{'end'},
                            $premiRNA{$premirna}->{'strand'},
                            $id,$pos,$ref,$alt,$type,$details,
                            scalar(@matures));
       
       my $seq = $premiRNA{$premirna}->{'seq'};
       foreach my $m (@matures){
           my $tstart=$premiRNA{$premirna}->{'mature'}->{$m}->{'start'};
           my $tend=$premiRNA{$premirna}->{'mature'}->{$m}->{'end'};
           my $chrStart=$premiRNA{$premirna}->{'mature'}->{$m}->{'chrStart'};
           my $chrEnd=$premiRNA{$premirna}->{'mature'}->{$m}->{'chrEnd'};
           my $tacc=$premiRNA{$premirna}->{'mature'}->{$m}->{'ACC'};
           my $tid=$m;
           my $tstr=$tstart."-".$tend.":".$chrStart."-".$chrEnd.":".$tacc.":".$tid.":".substr($seq,$tstart-1,$tend-$tstart+1);
	       print OUT "\t",$tstr;
	   }
	   if(scalar(@matures)==2){
		    print OUT "\t$seq\n";
		}else{
		    print OUT "\t\t$seq\n";
		}

    }
}


#Download the input data and tar them
sub downloadAndTar{
    my $dbsnpVcfURL="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/v4.0/00-All.vcf.gz";
    my $miRNADatURL="ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz";
    my $premiRNAgffURL="ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff";
    
    `wget $dbsnpVcfURL`;
    `wget $miRNADatURL`;
    
    `gunzip 00-All.vcf.gz`;
    `gunzip miRNA.dat.gz`;
    $dbSNPFile="00-All.vcf";  #rename it;
    
}


#since the SNP is in the premiRNA,
#This is used to determine whether the SNP is in the loop, upstream of mature miRNA, in the seed of mature miRNA etc.
sub getSNPmiRNAType{
    my($pos,$ref) = @_;
    my $return="";
    if(scalar(keys %{$ref->{'mature'}})==2){
        #could be in the loop
        #sort it's mature miRNA by chromosome location small to large
        my @mature_mirnas = sort {$ref->{'mature'}->{$a}->{'chrStart'} <=> $ref->{'mature'}->{$b}->{'chrStart'}} keys %{$ref->{'mature'}};
        #To see whether in the loop
        if($ref->{'mature'}->{$mature_mirnas[0]}->{'chrEnd'}< $pos && $ref->{'mature'}->{$mature_mirnas[1]}->{'chrStart'}>$pos){
            #it is in the loop
            #return Type: means the distance to miRNA1's 3' prime is x while the distance to miRNA2's 5 prime is y
            #loop|miRNA1:three:x|miRNA2:five:y
            $return="loop|";
            if($ref->{'strand'} eq "+"){
                my $threeDistance=$pos - $ref->{'mature'}->{$mature_mirnas[0]}->{'chrEnd'};
                my $fiveDistance=$ref->{'mature'}->{$mature_mirnas[1]}->{'chrStart'}-$pos;
                $return.=$mature_mirnas[0].":three:".$threeDistance."|".$mature_mirnas[1].":five:".$fiveDistance;
                return ("loop",$return);
            }else{
                my $fiveDistance=$pos- $ref->{'mature'}->{$mature_mirnas[0]}->{'chrEnd'};
                my $threeDistance=$ref->{'mature'}->{$mature_mirnas[1]}->{'chrStart'}-$pos;
                $return.=$mature_mirnas[1].":three:".$threeDistance."|".$mature_mirnas[0].":five:".$fiveDistance;
                return ("loop",$return);
            }    
        }
     }
     
     #If not in the loop then could be one of upstream, in the seed, in mature miRNA but not in the seed, downstream
     #sort it's mature miRNA by the distance to SNPs from small to larger, since it could not be located in the two miRNAs
     my @mature_mirnas = sort { abs($pos-$ref->{'mature'}->{$a}->{'chrStart'}) <=> abs($pos-$ref->{'mature'}->{$b}->{'chrStart'})} keys %{$ref->{'mature'}};
     foreach my $mirna (@mature_mirnas){
        #to see whether it is located in the miRNA/seed
        if($ref->{'mature'}->{$mirna}->{'chrStart'}<=$pos && $ref->{'mature'}->{$mirna}->{'chrEnd'}>=$pos){
            if($ref->{'strand'} eq "+"){ #deal with + strand
                my $fiveDistance=$pos-$ref->{'mature'}->{$mirna}->{'chrStart'};
                if($ref->{'mature'}->{$mirna}->{'chrStart'}<$pos && ($ref->{'mature'}->{$mirna}->{'chrStart'}+7)>=$pos){
                    $return="seed|$mirna:five:".$fiveDistance;
                    return ("seed",$return);
                }else{
                    $return="inButnotSeed|$mirna:five:".$fiveDistance;
                    return ("inButnotSeed",$return);
                }
            }else{                      #deal with  - strand
                my $fiveDistance=$ref->{'mature'}->{$mirna}->{'chrEnd'}-$pos;
                if($ref->{'mature'}->{$mirna}->{'chrEnd'}>$pos && ($ref->{'mature'}->{$mirna}->{'chrEnd'}-7)<=$pos){
                    $return="seed|$mirna:five:".$fiveDistance;
                    return ("seed",$return);
                }else{
                    $return="inButnotSeed|$mirna:five:".$fiveDistance;
                  return ("inButnotSeed",$return);
                }
            }
        }
        
        #if not located in the miRNA itself, determine upstream or downstream
        if($pos<$ref->{'mature'}->{$mirna}->{'chrStart'}){
            if($ref->{'strand'} eq "+"){ #deal with + strand
                my $fiveDistance=$ref->{'mature'}->{$mirna}->{'chrStart'}-$pos;
                $return="upstream|$mirna:five:".$fiveDistance;
                return ("upstream",$return);
            }else{
                my $threeDistance=$ref->{'mature'}->{$mirna}->{'chrStart'}-$pos;
                $return="downstream|$mirna:three:$threeDistance";
                return ("downstream",$return);
            }
        }else{
            if($ref->{'strand'} eq "+"){
                my $threeDistance = $pos-$ref->{'mature'}->{$mirna}->{'chrEnd'};
                $return="downstream|$mirna:three:$threeDistance";
                return ("downstream",$return);
            }else{
                my $fiveDistance = $pos-$ref->{'mature'}->{$mirna}->{'chrEnd'};
                $return="upstream|$mirna:five:$fiveDistance";
                return ("upstream",$return);
            }
        }
     }
}






#parse the pre-miRNA gff file, use pre-miRNA ID as key to store the dat
#items include ACC, chr,start,end,strand
sub parsingPremiRNA{
	my($in,$ref,$positionRef) = @_;
	open(IN,$in) or die $!;
	while(<IN>){
		next if (/^#|^$/);
		s/\r|\n//g;
		my($chr,undef,undef,$start,$end,undef,$strand,undef,$ids) = split "\t";
		my $acc="";
		my $id="";
		if($ids=~/ACC="(\w+)"; ID="([\w-]+)";/){
			$acc=$1;$id=$2;
		}
		$ref->{$id}->{'ACC'}=$acc;
		$ref->{$id}->{'chr'}=$chr;
		$ref->{$id}->{'start'}=$start;
		$ref->{$id}->{'end'}=$end;
		$ref->{$id}->{'strand'}=$strand;
		
		for(my $a=$start;$a<=$end;$a++){
		    $positionRef->{$chr}->{$a}->{$id}=1;
		}
	}
	close IN;
}




#this code is used to pasring the miRNA.dat file
#@parameters
#$infile: miRNA.dat 
#$outfile: store the relationship of pre-miRNA and miRNA, columns are #pre-miRNA ID","pre-miRNA ACC","chr","start","end","strand","#miRNAs","mature miRNAs\n"
#$ref:  premiRNA hash reference
sub parsingmiRNAdat{
	my ($infile,$outfile,$ref) = @_;
	open(IN,$infile) or die $!;
	open(OUT,">$outfile") or die $!;
	print OUT join "\t",("#pre-miRNA ID","pre-miRNA ACC","chr","start","end","strand","#miRNAs","mature miRNAs\n");
	my $string="";
	while(<IN>){
		if(/^\/\//){
			my ($acc,$id,$seq,@matures)=doit($string);
			if(exists($ref->{$id})){ #only parse those human pre-miRNA to mature miRNAs
				#first print OUT 
				print OUT join "\t",($id,$acc,$ref->{$id}->{'chr'},$ref->{$id}->{'start'},$ref->{$id}->{'end'},$ref->{$id}->{'strand'},scalar(@matures));
				$seq=uc($seq);
				$ref->{$id}->{'seq'}=uc($seq);
				foreach my $tmpRef (@matures){
					my ($tstart,$tend,$tacc,$tid) = @{$tmpRef};
					$ref->{$id}->{'mature'}->{$tid}->{'ACC'}=$tacc;
					$ref->{$id}->{'mature'}->{$tid}->{'start'}=$tstart; #this start is based on pre-miRNA, from 5'-> 3'
					$ref->{$id}->{'mature'}->{$tid}->{'end'}=$tend;
					my ($chrStart,$chrEnd) = determineRegion($ref->{$id}->{'start'},$ref->{$id}->{'end'},$ref->{$id}->{'strand'},$tstart,$tend);
					$ref->{$id}->{'mature'}->{$tid}->{'chrStart'}=$chrStart;
					$ref->{$id}->{'mature'}->{$tid}->{'chrEnd'}=$chrEnd;
										
					my $tstr=$tstart."-".$tend.":".$chrStart."-".$chrEnd.":".$tacc.":".$tid.":".substr($seq,$tstart-1,$tend-$tstart+1);
					print OUT "\t",$tstr;
				}
				if(scalar(@matures)==2){
				    print OUT "\t$seq\n";
				}else{
				    print OUT "\t\t$seq\n";
				}
			}
			$string="";
		}else{
			$string.=$_;
		}
	}
}


#To determine the miRNA's chromosome region based on it relative position on pre-miRNA and pre-miRNA's chromosome location
#Strand information is important here cause if strand=+, it is straitforward while if strand=-, Since miRNA's relative position to pre-miRNA is based on 
#pre-miRNA's 5'->3', that is 5's chromosome location is the "start" while the 3's is the "end"
#for strand -, reverse->comp
sub determineRegion{
    my($chrStart,$chrEnd,$strand,$start,$end) = @_;
    my $returnStart=0;
    my $returnEnd=0;
    if($strand eq "+"){
        #straightforward
        $returnStart = $chrStart+$start-1;
        $returnEnd = $chrStart+$end-1;
    }else{
        $returnStart = $chrEnd-$end+1;
        $returnEnd = $chrEnd-$start+1;
    }
    return ($returnStart,$returnEnd);
}


#actuall code to parse the miRNA.dat file
#FT   miRNA           16..36
#FT                   /accession="MIMAT0000002"
#FT                   /product="cel-lin-4"
#FT                   /evidence=experimental
#FT                   /experiment="cloned [1,3-5]"
sub doit{
	my $string=$_[0];
	$string=~/^ID\s+([\w-]+)\s/m;
	my $hairpin=$1;  #which is the pre-miRNA ID
	$string=~/^AC\s+(\w+)/m;
	my $ac=$1;       #which is the pre-miRNA ACC
	my @mature_mirna;
	
	#parsing miRNA info
	while($string=~/FT\s+?miRNA\s+(\d+)\.\.(\d+)\nFT\s+\/accession="(\w+)"\nFT\s+\/product="(.*?)"/mg){
		#1 start ; 2 end; 3 acc; 4 id
		push @mature_mirna,[$1,$2,$3,$4];
	}
	
	#pasring sequence
	my $seq="";
	pos($string)=0;
	if($string=~/SQ\s+?Sequence[\w\W]+other;([\S\s\r\n]+)/m){
	    $seq=$1;
	}
	$seq=~s/\s|\d//g;
	
	return ($ac,$hairpin,$seq,@mature_mirna);
}





sub info{
	return "[",scalar(localtime),"]";
}

