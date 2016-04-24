#!uer/bin/perl 
use warnings;
use strict;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Getopt::Long;
sub usage{
	print <<USAGE;
usage:
	perl $0 -listMode1 listfile -report -project -prefix -QCRaw -Resequencelist -outdir -QN -exists
		1. -listMode1 listfile
			Column 1:sampleName
			Column 2:sample library
			Column 3:the time of beginning sequence
			Column 4:Flowcell number
			Column 5:lane
			Column 6:Rawdir(eg:/share/fqdata032/solexa)
		2. -report the path of report output
		3. -prefix 1
		4. -project projectname(eg:HK11165_HUMqvwR)
		5. -help 
		6. -QCRaw
		7. -Resequencelist
		8. -NRate 3-10
		9. -maxErrorRate 2
		10.-contGC 37-49
		11. -maxproGC 0.1
		12. -minQ20   83
		13. -minQ30 75
		14. -librarylist
		15. -outdir
		16. -Nnumber 5
		17. -Lowquality 0.5
		18. -QN 
		19. -listMode2 listfile
			Column 1:lanelib
			Column 2:Rawdir
		20. -finsih

USAGE
}
#===========================parameters===========================#
my($help,$listMode1,$listMode2,$prefix,$project,$report,$QCRaw,$Resequencelist,$NRate,$maxErrorRate,$contGC,$maxproGC,$minQ20,$minQ30,$librarylist,$QN,$finish);
my($FilterRawdir,$Nnumber,$Lowquality,$outdir);
my%hashexists;
GetOptions(
	"help"=>\$help,
	"listMode1:s"=>\$listMode1,
	"listMode2:s"=>\$listMode2,
	"librarylist:s"=>\$librarylist,
	"report:s"=>\$report,
	"project:s"=>\$project,
	"prefix:s"=>\$prefix,
	"QCRaw"=>\$QCRaw,
	"Resequencelist"=>\$Resequencelist,
	"NRate:s"=>\$NRate,
	"maxErrorRate:s"=>\$maxErrorRate,
	"contGC:s"=>\$contGC,
	"maxproGC:s"=>\$maxproGC,
	"minQ20:s"=>\$minQ20,
	"minQ30:s"=>\$minQ30,
	"FilterRawdir:s"=>\$FilterRawdir,
	"outdir:s"=>\$outdir,
	"Nnumber:s"=>\$Nnumber,
	"Lowquality:s"=>\$Lowquality,
	"QN"=>\$QN,
	"finish:s"=>\$finish,
	);
if($help){ &usage();exit 0;}
unless(defined $project){&usage();exit 0;}
unless(defined $prefix){$prefix="01"};
if($QCRaw){
	unless (defined $report){$report="./"};
#	unless (defined $FilterRawdir){$FilterRawdir="./FilterRawdir";mkdir($FilterRawdir)||die $!};
	}
if($Resequencelist){
	unless (defined $outdir){$outdir="./"};
	unless (defined $Nnumber){$Nnumber="5"};
	unless (defined $Lowquality ){$Lowquality="0.5"};
}
#===============================================default==============================================#
if($Resequencelist){
	my$Resequencelistfile="$outdir/$prefix.$project.RE.list";
	open RESE,">","$Resequencelistfile" or die $!;
}
if($finish){
	open EX,"$finish" or die $!;
	while(<EX>){
		chomp;
		my$existslane=(split /\s+/,$_)[0];
		$hashexists{$existslane}=1;
	}
}
#===========================================do QC Rawdata============================================#
my($QCRawreport,$mincontGC,$maxcontGC,$firstNRate,$secondNRate);
if($QCRaw){
	my$QCRawreport="$report/$prefix.$project.QC_raw.report";
	my$disqualification="$report/$prefix.$project.disqualification";
	open DIS,">","$disqualification" or die $!;
	open QCRAWRE,">","$QCRawreport" or die $!;	
	unless (defined $NRate){$NRate="3-10"};
	unless (defined $contGC){$contGC="37-49"};
	unless (defined $maxErrorRate){$maxErrorRate="2"};
	unless (defined $maxproGC){$maxproGC="0.1"};
	unless (defined $minQ20){$minQ20="83"};
	unless (defined $minQ30){$minQ30="75"};
	($mincontGC,$maxcontGC)=(split /\-/,$contGC)[0,1];
	($firstNRate,$secondNRate)=(split /\-/,$NRate)[0,1];
	print QCRAWRE"Samplename\tPath\tN_Rate($firstNRate%-$secondNRate%)\tError($maxErrorRate%)\tcontentGC($mincontGC%-$maxcontGC%)\tproportionGC($maxproGC)\tQ20($minQ20%)\tQ30($minQ30%)\tparticulars\n";
	}
#==============================================Main process===========================================#
my($sampleName,$samplelibrary,$Flowcell,$lane,$Rawdir,$time,$RawPath);
my%hashlib;
if($listMode1){
	open LIST,"$listMode1" or die $!;
}elsif($listMode2){
	if($librarylist){
	open LIB,"$librarylist" or die $!;
	open LIST,"$listMode2" or die $!;
	while(<LIB>){
		chomp;
		my($l_sample,$l_library)=(split /\t/,$_)[0,1];
		unless(exists $hashlib{$l_library}){
			$hashlib{$l_library}=$l_sample;
			}
		}
	}else{
	 print "WARNING: Need the file of library information list ,you can get from BMS!\n";
	&usage();exit 0;
	}
}
#=========================================================================================
	while(<LIST>){
	chomp;
	my@list=split;
	if($listMode1){
	$sampleName=$list[0];
	$samplelibrary=$list[1];
	$time=$list[2];
	$Flowcell=$list[3];
	$lane=$list[4];
	$Rawdir=$list[5];
	$RawPath="$Rawdir/$project/*/*\_$Flowcell\_L$lane\_$samplelibrary/*\_$Flowcell\_L$lane\_$samplelibrary\_1.fq.gz";
	}elsif($listMode2){
	my$lanelibT=$list[0];
	($time,$Flowcell,$lane,$samplelibrary)=(split /\_/,$lanelibT)[0,2,3,4];
	$lane=~s/L//g;
	$Rawdir=$list[1];
	if(exists $hashlib{$samplelibrary}){
		$sampleName=$hashlib{$samplelibrary};
		}else{
		print "WARNING:the library  is wrong\n";
		}
	$RawPath="$Rawdir/$project/*/$lanelibT/$lanelibT\_1.fq.gz";
}
	my@fqlist = glob($RawPath);
	if(@fqlist==0){
		print "WARNING:No FASTQ file for $sampleName,$samplelibrary,$Flowcell,$Rawdir/$project\n";
		next;
		}
	foreach my$path(@fqlist){
		my$fq_gz=basename($path);my@temp=split(/\_/,$fq_gz);my$lanelib=join"_",@temp[0..4];
	if($finish){
		if(exists $hashexists{$lanelib}){
			print "the $lanelib have been exists\n";
			next;
		}
		}
		my$rawLaneDir=dirname($path);
		my$pooling=basename(dirname($rawLaneDir));
		my$projectT=basename(dirname(dirname($rawLaneDir)));
		my$raw_pe1="$rawLaneDir/$lanelib\_1.fq.gz";
		my$raw_pe2="$rawLaneDir/$lanelib\_2.fq.gz";
		my$adapterlist1="$rawLaneDir/1.adapter.list.gz";
		my$adapterlist2="$rawLaneDir/2.adapter.list.gz";
	if($Resequencelist){
		unless($QCRaw){
		print RESE"$sampleName\t$raw_pe1\t$raw_pe2\t$adapterlist1\t$adapterlist2\t$lanelib\t$Nnumber\t$Lowquality\n";}
	}
#==================================================QC Rawfile==========================================#
	if($QCRaw){
	my($NRateflag,$error_rateflag,$contGCflag,$maxproGCflag,$minQ20flag,$minQ30flag)=(0,0,0,0,0,0);
	my($NRatevalue,$errorvalue,$contGCvalue,$maxproGCvalue,$minQ20value,$minQ30value)=(0,0,0,0,0,0);
	my($GC_number,$position_number)=(0,0);
	my@Nlist1=();
	my@Nlist2=();
	my@Ntemp=();
	my$Nflag;
	my$check_pe1="$rawLaneDir/1.fqcheck";
	my$check_pe2="$rawLaneDir/2.fqcheck";
	my@particulars=("congratulations");
	open CH1,"$check_pe1" or die $!;
	open CH2,"$check_pe2" or die $!;
		while(<CH1>){
			chomp;
			if($_=~/^base/){
			$position_number++;
			my@check1=split;
			my$pos1=$check1[1]-1;
			my$GC_D_value1=abs($check1[3]-$check1[4]);
			my$AT_D_value1=abs($check1[2]-$check1[5]);
			if($GC_D_value1>2 or $AT_D_value1>2){
				$GC_number++;
			}
			if($check1[6]>=$firstNRate and $check1[6]<$secondNRate){
				push @Nlist1,$pos1;
				if($NRateflag==0){
					$NRatevalue=$check1[6];
					$NRateflag=1;
					}else{
					my($maxNRatevalue)=sort{$b <=> $a}($NRatevalue,$check1[6]);
					$NRatevalue=$maxNRatevalue;
					}
			}elsif($check1[6]>$secondNRate){
				push @Nlist1,$pos1;	
					$NRateflag=2;
					my($maxNRatevalue)=sort{$b <=> $a}($NRatevalue,$check1[6]);
					$NRatevalue=$maxNRatevalue;
			}
			}elsif($_=~/^Error/){
				my$QC=<CH1>;
				my($error_rate,$GC,$Q20,$Q30)=(split /\s+/,$QC)[0,1,2,3];
				if($error_rate>$maxErrorRate){
					$errorvalue=$error_rate;
					$error_rateflag=1;
					}
				if($GC>$maxcontGC or $GC<$mincontGC ){
					$contGCvalue=$GC;
					$contGCflag=1;
					}
				if($Q20<$minQ20 ){
					$minQ20value=$Q20;
					$minQ20flag=1;
					}
				if($Q30<$minQ30 ){
					$minQ30value=$Q30;
					$minQ30flag=1;
					}
			}
		}
		while(<CH2>){
                        chomp;
                        if($_=~/^base/){
                        $position_number++;
			my@check2=split;
			my$pos2=$check2[1]-1;
                        my$GC_D_value2=abs($check2[3]-$check2[4]);
			my$AT_D_value2=abs($check2[2]-$check2[5]);
			if($GC_D_value2>2 or $AT_D_value2>2){
				$GC_number++;
			}
		   if($check2[6]>=$firstNRate and $check2[6]<$secondNRate){
                               	push @Nlist2,$pos2;
				 if($NRateflag==0){
                                $NRatevalue=$check2[6];	
				$NRateflag=1;
                        	}else{	
				my($maxNRatevalue)=sort{$b <=>$a}($NRatevalue,$check2[6]);
				$NRatevalue=$maxNRatevalue;
				}
                        }elsif($check2[6]>$secondNRate){
				push @Nlist2,$pos2;
				$NRateflag=2; 
                                my($maxNRatevalue)=sort{$b <=>$a}($NRatevalue,$check2[6]);
				$NRatevalue=$maxNRatevalue;
                        }
                	}elsif($_=~/^Error/){
				my$QC=<CH2>;
                              	 my($error_rate,$GC,$Q20,$Q30)=(split /\s+/,$QC)[0,1,2,3];
                                if($error_rate>$maxErrorRate){
                                            my($maxerrorvalue)=sort{$b <=> $a}($error_rate,$errorvalue);
					    $errorvalue=$maxerrorvalue;
					    $error_rateflag=1;		
                                        }
                                if($GC>$maxcontGC or $GC<$mincontGC){
                                             if($contGCflag==1){
						my$temp=$contGCvalue;
						$contGCvalue=join":",$temp,$GC;
						}
                                        }
                                if($Q20<$minQ20){
						if($minQ20flag==1){
						my($minQ20)=sort{$a<=>$b}($Q20,$minQ20value);
                                                $minQ20value=$minQ20;
                                                }else{
                                                $minQ20value=$Q20;
						$minQ20flag=1;
						}
                                        }
                                if($Q30<$minQ30){
						if($minQ30flag){
					       my($minQ30)=sort{$a<=>$b}{$Q30,$minQ30value};
                                               $minQ30value=$minQ30;
						}else{
					       $minQ30value=$Q30;
					       $minQ30flag=1;
						}
                                        }
					
				}
			}
		my$Nlisttemp1=join ",",@Nlist1;my$Nlisttemp2=join ",",@Nlist2;
		if(@Nlist1){my$Nlist_fq1=join":",1,$Nlisttemp1;push @Ntemp,$Nlist_fq1;}
		if(@Nlist2){my$Nlist_fq2=join":",2,$Nlisttemp2;push @Ntemp,$Nlist_fq2;}
		if(@Ntemp){
			$Nflag=join ";",@Ntemp;
		}
		my$proGC=$GC_number/$position_number;
		if($proGC>$maxproGC){
					$maxproGCvalue=$proGC;
					$maxproGCflag=1;
			}
	if($NRateflag or $error_rateflag or $contGCflag or $minQ20flag or $minQ30flag or $maxproGCflag){
		      if($NRateflag){
			print DIS"$sampleName\t$lanelib\t$Rawdir\t$Nflag\n";
			}else{
			 print DIS"$sampleName\t$lanelib\t$Rawdir\ttemp\n";
			if($QN){
			print RESE"$sampleName\t$raw_pe1\t$raw_pe2\t$adapterlist1\t$adapterlist2\t$lanelib\t$Nnumber\t$Lowquality\n";
				}
			}
			my@particulars;
		if($NRateflag==1){
			push @particulars,"you just need remove all of the N base that in the fastq";	
			}elsif($NRateflag==2){
			push @particulars,"you need the tile !!";
			}else{
				#push @particulars,"the rate of N base is ok .";
				$NRatevalue="OK";
			}
		if($error_rateflag){
				push @particulars,"Unfortunately,the rate of error base is too high";
			}else{
				#push @particulars,"the Error rate is too ok.";
				$errorvalue="OK";
			}
		if($maxproGCflag){
				push @particulars,"the GC is unequal";
			}else{
				$maxproGCvalue="OK";
			}
		if($contGCflag){
				push @particulars,"the content of GC is unreasonable.Maybe the data is unavailable";
			}else{
				#push @particulars,"the content of GC is reasonable .";
				$contGCvalue="OK";
			}
		if($minQ20flag){
				push @particulars,"Q20 is disqualification";
			}else{
				#push @particulars,"Q20 is ok .";
				$minQ20value="OK";
			}
		if($minQ30flag){
				push @particulars,"Q30 is  disqualification";
			}else{
				#push @particulars,"Q30 is ok .";
				$minQ30value="OK";
			}
		if($maxproGCflag){
				push @particulars,"GC is unequal";
			}else{
				#push @particulars,"the proportion of GC is ok .";
				$maxproGCvalue="OK";
			}
		my$particulars_list=join ";",@particulars;
		print QCRAWRE"$sampleName\t$lanelib\t$NRatevalue\t$errorvalue\t$contGCvalue\t$maxproGCvalue\t$minQ20value\t$minQ30value\tUnfortunately\n";	
		print "$particulars_list\n";
		}else{
		print QCRAWRE"$sampleName\t$lanelib\tOK\tOK\tOK\tOK\tOK\tOK\tcongratulations\n";
		if($Resequencelist ){
		print RESE"$sampleName\t$raw_pe1\t$raw_pe2\t$adapterlist1\t$adapterlist2\t$lanelib\t$Nnumber\t$Lowquality\n";
		}
		}
}
#==========================================================================================================#
		}
	}

