#!/usr/local/bin/perl -w

use Getopt::Long;
use Bio::DB::Sam;
use Data::Dumper;
use strict;

#use lib '/package/cbrg/lib/';
#use CapSequm_utils;

&GetOptions
(
	#"n=s"=>\my $data_name,
	#"t=s"=> \my $data_tag,
	"r=s"=> \my $ref,
	"b=s"=> \my $build,
	"chr=s"=> \my $chr_lengths_file,
	"bed=s"=> \my $coord_data,
	"o=i"=>\my $oligo_size,
	"c=s"=>\my $cutsite
);

$| = 1;

#my $tmp_data_dir = '/t1-data/wwwtmp/CapSequm/' . $data_tag . '/';
#my $coord_data =  $tmp_data_dir . 'input.bed';	
my $cutsite_length = length($cutsite);

#my $output_file_root = $tmp_data_dir . $data_name;

#open (OUTPUT,  ">$output_file_root\_Walk\_loops.mfa") or die "cannot open mfa file: $! ";
#open (OUTPUT1, ">$output_file_root\_Walk\_loops.bed") or die "cannot open bed file: $! ";	
#open (OUTPUT2, ">$output_file_root\_Not_Designed.txt") or die "cannot open txt file: $! ";
#open (OUTPUT3, ">$output_file_root\_Fragments.bed") or die "cannot open bed file: $! ";
#open (OUTPUT4, ">$output_file_root\_Fragments.gff") or die "cannot open gff file: $! ";
#open (OUTPUT5, ">$output_file_root\_tmp.txt") or die "open temp file: $! ";

open (OUTPUT,  ">oligos.mfa") or die "cannot open mfa file: $! ";
open (OUTPUT1, ">oligos.bed") or die "cannot open bed file: $! ";	
open (OUTPUT2, ">Not_Designed.txt") or die "cannot open txt file: $! ";
open (OUTPUT3, ">fragments.bed") or die "cannot open bed file: $! ";
open (OUTPUT4, ">fragments.gff") or die "cannot open gff file: $! ";
open (OUTPUT5, ">tmp.txt") or die "open temp file: $! ";


# get build info:
#my ($genome_file, $chr_lengths_file, $organism, $scientific_name) = &CapSequm_utils::_get_build_data_files($build);

# Store chromosome length of build to prevent cutsite serach falling of end of chromosome
#my $chr_lengths_file = "/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt";
my %chrlength;
open (SIZES, $chr_lengths_file) or die "couldn't open chr_lengths_file: $chr_lengths_file: $!\n";	
while (<SIZES>)
{
	chomp;
	my ($gchr, $gsize)= split(/\t+/);
	unless ($gchr =~ /M/gi)
	{   
		$chrlength{$gchr} = $gsize;
	}		
}	
close SIZES;


# Open design coordinate file and store:
my %StoredCoor;
my %designed_Probes;
my %output1_data;
my %Redisigners;

my $TempGeneID = 1;
my %lookup;


#my $coord_data = '/t1-data1/WTSA_Dev/jkerry/CaptureC/DDownes/CapsequmInput_2.txt';
open(INFO, $coord_data) or die "couldn't open coord_data file: $coord_data: $!\n";
while (<INFO>)
{
	chomp;
	if(/Gene/gi){next;}
	
	my ($chr, $TSS, $TSS2, $Gene) = split(/\s+/);

	my $Newline = "$chr" . "\t" . "$TSS" . "\t" . "$TSS2" . "\t" . "$TempGeneID";
	$StoredCoor{$chr}{$TSS}=$Newline;	

	$lookup{$TempGeneID} = $Gene;
	print OUTPUT5 "$TempGeneID\t$Gene\n";
	$TempGeneID++;
}
close INFO;
close OUTPUT5;

# Load Genome:
#my $fai = Bio::DB::Sam::Fai->load("$genome_file");
my $fai = Bio::DB::Sam::Fai->load("$ref");

# Run through design Hash for analysis:
foreach my $Storedchr (sort keys %StoredCoor)
{
	foreach my $StoredTss (sort by_number keys %{$StoredCoor{$Storedchr}})
	{
		my ($Tchr, $Tstart, $Tstop, $Target) = split(/\s+/,$StoredCoor{$Storedchr}{$StoredTss});
		
		my $OriginalId = $lookup{$Target};
		
		my $GenomicLeftStart;
		my $GenomicLeftStop;
		
		my $InitialStart = $StoredTss;
		my $InitialCheckStart = ($StoredTss - ($cutsite_length-1));
		my $InitialStop = ($StoredTss + ($cutsite_length-1));
		
		my $SToredChrLength = $chrlength{$Storedchr};
		
		my $fai_location_Check = $Storedchr . ':' . $InitialCheckStart . '-' . $InitialStop;
		my $subseqcheck = $fai->fetch($fai_location_Check);
		
		
		# Check coordinate is not in a cutsite if so skip and report to Redesigners
		if ($subseqcheck =~ /$cutsite/gi)
		{
			$Redisigners{$Storedchr}{$StoredTss}= "$Tchr\t$Tstart\t$Tstop\t$OriginalId\t is a $cutsite cutsite $subseqcheck!";
			next;
		}
		
		
		# Walk left and find cutsite, unless fall of the start of a chromosome if so skip and report to Redesigners		
		my $LeftMoveStart = $InitialStart;
		my $subseqLeft = 'NNNN';
		my $LFlag = 'unmatched';
		
		until (($LFlag eq 'match') || ($LFlag eq 'start'))
		{
			$LeftMoveStart --;
			my $LeftMoveEnd = ($LeftMoveStart + ($cutsite_length-1));
			
			my $fai_locationL = $Storedchr . ':' . $LeftMoveStart . '-' . $LeftMoveEnd;
			$subseqLeft = $fai->fetch($fai_locationL);
		
			if ($subseqLeft =~ /$cutsite/gi)
			{
				$GenomicLeftStart = 	$LeftMoveStart;		
				$GenomicLeftStop = 	$LeftMoveStart + $oligo_size;
				$LFlag = 'match';
			}
			elsif ($LeftMoveStart == 1)
			{
				$Redisigners{$Storedchr}{$StoredTss}= "$Tchr\t$Tstart\t$Tstop\t$Target\t is the start of a Chromosome!";
				$LFlag = 'start';
			}
		}

		if ($LFlag eq 'start')
		{
			#print"skipping chr start\n";
			next;
		}

		
		# Walk Right and find cutsite, unless fall of the end of a chromosome if so skip and report to Redesigners
		my $GenomicRightStart;
		my $GenomicRightStop;
		my $RightMoveStart = $InitialStart;
		my $RightMoveEnd;
		my $subseqRight = 'NNNN';
		my $RFlag = 'unmatched';
		
		until (($RFlag eq 'match') || ($RFlag eq 'end'))
		{
			$RightMoveStart++;
			$RightMoveEnd = ($RightMoveStart + ($cutsite_length-1));
			
			my $fai_locationR = $Storedchr . ':' . $RightMoveStart . '-' . $RightMoveEnd;
			$subseqRight= $fai->fetch($fai_locationR);
			
			if ($subseqRight =~ /$cutsite/gi)
			{
				$GenomicRightStart = $RightMoveEnd - $oligo_size;
				$GenomicRightStop = $RightMoveEnd;
				$RFlag = 'match';
			}
			elsif ($RightMoveEnd == $SToredChrLength)
			{
				$Redisigners{$Storedchr}{$StoredTss}= "$Tchr\t$Tstart\t$Tstop\t$OriginalId\t is the end of a Chromosome";
				$RFlag = 'end';
			}
		}
		
		if ($RFlag eq 'end')
		{
			#print"skipping chr end\n";
			next;
		}
		
		
		# Reformat chr name and ID and generate unique parseble IDS		
		my $targetchr = $Storedchr;
		$targetchr =~ s/chr//gi;		
		my $LeftID  = $targetchr . ':' . $GenomicLeftStart . '-' . $GenomicLeftStop . '_' . $Target;
		my $RightID = $targetchr . ':' . $GenomicRightStart . '-' . $GenomicRightStop . '_' . $Target;

		
		# Check if design is not redundant (two specified points in a single fragment).  If redundant keep one and report the redundancy of the other to Redesigners		
		my $RedFlag = 'unique';
		
		if (exists $designed_Probes{$targetchr}{$GenomicLeftStart}{$GenomicLeftStop}{$GenomicRightStart}{$GenomicRightStop})
		{
			my $storedTarget = $designed_Probes{$targetchr}{$GenomicLeftStart}{$GenomicLeftStop}{$GenomicRightStart}{$GenomicRightStop};
			$RedFlag ='redundant';
			$Redisigners{$targetchr}{$StoredTss} = "$Tchr\t$Tstart\t$Tstop\t$OriginalId\t is redundant to target $storedTarget";
		}
		else
		{
			$designed_Probes{$targetchr}{$GenomicLeftStart}{$GenomicLeftStop}{$GenomicRightStart}{$GenomicRightStop}=$OriginalId;
		}
		
		if ($RedFlag eq 'redundant')
		{
			#print"skipping redundant design\n";
			next;
		}


		# Extract the sequence of the complete targeted fragment and check its size.  If it is too small skip and report to Redesigners
		my $FragCoor = $Storedchr . ':' . $LeftMoveStart . '-' . $RightMoveEnd;
		my $nseq = $fai->fetch($FragCoor);
		
		$nseq =~ tr/gatc/GATC/;
		my $SizeFlag = 'OK';
		
		my $size = length($nseq);
		
		if ($size<$oligo_size)
		{
			$Redisigners{$Storedchr}{$StoredTss} = "$Tchr\t$Tstart\t$Tstop\t$OriginalId\t fragment is too small, $size bp!";
			$SizeFlag = 'TooSmall';
		}
		
		if ($SizeFlag eq 'TooSmall')
		{
			#print"skipping too small fragment\n";
			next;
		}	
		
		my $EndStart = $size - $oligo_size;
		my $fragment1 = substr($nseq, 0, $oligo_size);
		my $fragment2 = substr($nseq, $EndStart, $oligo_size);
		
		unless(exists($chrlength{$targetchr}))
		{
			$targetchr = 'chr' . $targetchr;
			$LeftID = 'chr' . $LeftID;
			$RightID = 'chr' . $RightID;
		}
		
		print OUTPUT ">$LeftID\n$fragment1\n>$RightID\n$fragment2\n";
		print OUTPUT3 "$targetchr\t$LeftMoveStart\t$RightMoveEnd\t$OriginalId\n";
		print OUTPUT4 "$targetchr\tFragExtract.pl\t$build\t$LeftMoveStart\t$RightMoveEnd\t\.\t\.\t\.\tName=$OriginalId\;LeftID=$LeftID;RightID=$RightID\n";
		
		$output1_data{$targetchr}{$GenomicLeftStart} = "$targetchr\t$GenomicLeftStart\t$GenomicLeftStop\t$LeftID\n$targetchr\t$GenomicRightStart\t$GenomicRightStop\t$RightID\n";
	}
}
close OUTPUT;
close OUTPUT3;
close OUTPUT4;


foreach my $SortChr (sort keys %output1_data)
{
	foreach my $SortPosition (sort by_number keys %{$output1_data{$SortChr}})
	{		
		print OUTPUT1 $output1_data{$SortChr}{$SortPosition};
	}
}
close OUTPUT1;


foreach my $Errorchr (sort keys %Redisigners)
{
	foreach my $ErrorTss (sort by_number keys %{$Redisigners{$Errorchr}} )
	{		
		my $StoredErrorCode = $Redisigners{$Errorchr}{$ErrorTss};
		
		print OUTPUT2 "$StoredErrorCode\n";
	}
}
close OUTPUT2;


exit(0);

#--------------------------		
		
sub by_number
{
	$a <=> $b;
}
