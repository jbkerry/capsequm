#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;
use File::Copy;
use Data::Dumper;
use IO::Handle;

use Env::Modulecmd;
Env::Modulecmd::load (qw(ucsctools));

use lib '/package/cbrg/lib/';
use CapSequm_utils;

&GetOptions (
		"blat=s"=> \my $file,
		"maf=s"=> \my $MafFile,
		"t=s"=> \my $data_tag,
		"n=s"=> \my $data_name,
		"b=s"=> \my $build
	     );

$| = 1;

my ($genome_file, $chr_lengths_file, $organism, $scientific_name) = &CapSequm_utils::_get_build_data_files($build);	


my $tmp_data_dir = '/t1-data/wwwtmp/CapSequm/' . $data_tag . '/';
my $www_tmp_dir = '/public/CaptureC/';
my $www_tmp_data_dir = $www_tmp_dir . $data_tag . '/';  

# full path to repeatmasker file:
my $repeatmasker_file = $MafFile . '.out';

# full paths to output files:
my $wig_file = $tmp_data_dir . $data_name . '.wig';
my $MIG_gff_file = $tmp_data_dir . $data_name . '_MIG.gff'; 
my $Oligos_gff_file = $tmp_data_dir . $data_name . '_Oligos.gff'; 
my $bigwig_file = $tmp_data_dir . $data_name . '.bw';

my $lookup = $tmp_data_dir . $data_name . '_tmp.txt';

open (OUTPUT, ">$wig_file") or die "cannot open file $wig_file: $!\n";
open (OUTPUT2, ">$MIG_gff_file") or die "cannot open file $MIG_gff_file: $!\n";
open (OUTPUT3, ">$Oligos_gff_file") or die "cannot open file $Oligos_gff_file: $!\n";

my %depthgauge;
my %Maffed;
my %used;
my %oligoCoverage;
my %OligoValue;
my %OligoSize;
my %StoredSeq;
my %dup_counts;
my %ProblemSequences;     # repeatmasker



# loop to load look up table generated by Fragextract
my %LoopUp;
open (LOOKUP, $lookup);
while (<LOOKUP>)
{
	chomp;	
	my ($LookUpId, $RealID)=split(/\t+/);
	$LoopUp{$LookUpId} = $RealID;
}
close LOOKUP;


my $MAfLine = 0;
my $MafHeader;
my $MafSeq;

open (MAF, $MafFile);	
while (<MAF>)
{
	chomp;	    
	$MAfLine++;
	
	if ($MAfLine == 1)
	{
		$MafHeader = $_;
		$MafHeader =~ s/>//g;

		my ($Mchr, $Mstart, $Mstop, @Mrest) = split(/[\:\-\_]/, $MafHeader);
		my $Mcounter = $Mstart;
		$Maffed{$MafHeader}=0;
		until ($Mcounter == $Mstop + 1)
		{
			$depthgauge{$Mchr}{$Mcounter}=0;
			$oligoCoverage{$Mchr}{$Mcounter}++;
			$Mcounter++;
		}
	}
	
	if ($MAfLine == 2)
	{
		$MafSeq = $_;
		$StoredSeq{$MafHeader}=$MafSeq;
		$MAfLine=0;
	}
}	
close MAF;


my $linenumber = 0;
open (BLAT, $file);	
while (<BLAT>)
{
	chomp;
	$linenumber++;
	if ($linenumber <= 5) { next; }
	
	my ($match, $mismatch, $repmatch, $ns, $qgap, $qgapbases, $tgap, $tgapbases, $strand, $query, $qsize, $qstart, $qend, $tname, $tsize, $tstart, $tend, $blockcount, $blocksize, $qtarts, $tstarts)=split(/\s+/);
	if ($tname =~ /_/g){next;}
	
	my $percent = ($qsize / 100);	
	my $percent_match = ($match / $percent);	
	if ($percent_match >= 70) { $dup_counts{$query}++; }

	my ($chr, $start, $stop, @rest)=split(/[\:\-\_]/, $query);		
	
	unless (exists $used{$query})
	{
		my $counter = $start;			
		until ($counter == ($start + $qsize + 1))
		{
			$oligoCoverage{$chr}{$counter}++;
			$counter++; 
		} 
	}
	
	$Maffed{$query}++;
	$used{$query}++;
	$OligoValue{$query} = 0;		
	
	my $loopstart = $start + $qstart;
	my $counter = $qstart;
	
	until ($counter == ($qend + 1))
	{
		$depthgauge{$chr}{$loopstart}++;
		$counter++;
		$loopstart++;		
	}
}
close BLAT;



my $LineCounter=0;
open (INFO3, $repeatmasker_file);
while (<INFO3>)
{
	chomp;    
	$LineCounter++;    
	unless ($LineCounter > 3){next;}
	
	my ($stuffer, $mscore, $percDiv, $percDel, $percIns, $Query, $QPositionSt, $QPositionStp, $left, $strand, $RepType, $RepClass, $RPositionST, $RpPositionStp, $RepLeft, $RepID)=split(/\s+/);
	
	my $match_length = $QPositionStp - $QPositionSt;
	
	my $stored_sequence = $StoredSeq{$Query};
	my $nasty_seq = substr($stored_sequence, $QPositionSt, $QPositionStp);
	
	$ProblemSequences{$Query} = $match_length . '_' . $nasty_seq . '_' . $RepType . '_' .$RepClass;
}
close INFO3;


foreach my $StoredID (keys %Maffed)
{	
	my ($storedchr, $storedstr, $storedstp, @rest2)=split(/[\:\-\_]/, $StoredID); 
	my $posCounter = $storedstr;

	my $size = $storedstp - $storedstr;
	$OligoSize{$StoredID}=$size;

	until ($posCounter == $storedstp + 1)
	{
		if (exists $depthgauge{$storedchr}{$posCounter})
		{
			$OligoValue{$StoredID}+=$depthgauge{$storedchr}{$posCounter};
			$posCounter++;
		}
		else
		{
			print "Warning: entry missing from blat file\n";   ### ...and I quote, "Fuckity Fuck!"
		}
	}
}



foreach my $Did (sort keys %OligoValue)
{
	my $density = $OligoValue{$Did}/$OligoSize{$Did};
	my $GC=0;
	my $CG=0;
	my $DupFlag = 'null';
	my $BlatFlag = 'null';
	my $dupCount = 0;
	
	###Captures info from MAF header Oligo, coor, Target name, Fragment Coor, Fragment Type (Core or left right flank), loop mumber
	my ($Dchr, $Dstr, $Dstp, $TargetNumber)=split(/[\:\-\_]/, $Did);

	$Dchr =~ s/chr//gi;
	my $StoredOligoSeq = $StoredSeq{$Did};

	$GC++ while ($StoredOligoSeq =~ /G|C/ig);
	$CG++ while ($StoredOligoSeq =~ /CG/ig);
	my $SeqLength = length($StoredOligoSeq);
	
	my $GCPercent = $GC ? ($GC/($SeqLength/100)) : 0;
	my $CGPercent = $CG ? ($CG/($SeqLength/100)) : 0;
	
	if (exists $dup_counts{$Did})
	{	
		$dupCount = $dup_counts{$Did};
		if ($dupCount > 1) { $DupFlag = 'TRUE'; }
		else { $DupFlag = 'FALSE'; }
	}	
	if (exists $used{$Did})  { $BlatFlag = 'TRUE'; }
	else { $BlatFlag = 'FALSE'; }

	#   Reform new id with real gene name for output;
	my $StoredRealID = $LoopUp{$TargetNumber};
	my $NewID = "$Dchr" . ':' . "$Dstr" . '-' . "$Dstp" . '_' . "$StoredRealID";

	print OUTPUT2 "chr$Dchr\tDepthGauge_MIG\t$data_name\t$Dstr\t$Dstp\t.\t.\t.\tNOFILTER_FragID=$NewID\; Oligo_Size=$OligoSize{$Did}\; Density=$density\; Percent_GC=$GCPercent\; Percent_CG=$CGPercent\; NOFILTER_Sequence=$StoredOligoSeq\; Duplication=$DupFlag\; Duplicates=$dupCount\; BLATed=$BlatFlag\; ";
	
	if (exists $ProblemSequences{$Did})
	{
		my ($StoredNLength, $StoredNsequence, $StoredType, $StoredClass) = split(/\_/, $ProblemSequences{$Did});
        	print OUTPUT2 "SimpleRepeat=TRUE\; SRepeatLength=$StoredNLength\; NOFILTER_RepeatSequence=$StoredNsequence\; RepeatType=$StoredType\; RepeatClass=$StoredClass\;\n";
	}
	else
	{
		print OUTPUT2 "SimpleRepeat=FALSE\; SRepeatLength=0\; NOFILTER_RepeatSequence=NA\; RepeatType=NA\; RepeatClass=NA\;\n";
	}
	
	print OUTPUT3 "chr$Dchr\tDepthGauge_Oligos\t$data_name\t$Dstr\t$Dstp\t.\t.\t.\tName=$NewID\;\n";	
}
close OUTPUT2;
close OUTPUT3;

foreach my $dChr(sort keys %depthgauge)
{
	print OUTPUT "variableStep  chrom=$dChr\n";
	
	foreach my $dcoor(sort numerically keys %{$depthgauge{$dChr}})
	{		    
		my $depth = $depthgauge{$dChr}{$dcoor};
		
		unless (exists $oligoCoverage{$dChr}{$dcoor})
		{
			print "warning $dChr $dcoor has no coverage value\n";
		}
	   
		my $coverage = $oligoCoverage{$dChr}{$dcoor};

		my $AdjustedDepth = $depth/$coverage;
		my $fullAdjust = sprintf("%.10f", $AdjustedDepth);
	    
		print OUTPUT "$dcoor\t$fullAdjust\n";
	}
}
close OUTPUT;




# convert wig -> bigwig
my $wigToBigWig_cmd = "wigToBigWig -clip $wig_file $chr_lengths_file $bigwig_file";
print "wigToBigWig cmd: $wigToBigWig_cmd\n";
system($wigToBigWig_cmd) == 0 or die "couldn't bigwig file: $wig_file: $!\n";

my $www_bigwig_file = $bigwig_file;
$www_bigwig_file =~ s/$tmp_data_dir/$www_tmp_data_dir/;
print "move bigwig cmd: move($bigwig_file, $www_bigwig_file)\n";
move($bigwig_file, $www_bigwig_file) or die "cannot move files: $bigwig_file, $www_bigwig_file: $!";

my $www_MIG_gff_file = $MIG_gff_file;
$www_MIG_gff_file =~ s/$tmp_data_dir/$www_tmp_data_dir/;
print "move MIG gff cmd: move($MIG_gff_file, $www_MIG_gff_file)\n";
move($MIG_gff_file, $www_MIG_gff_file) or die "cannot move files: $MIG_gff_file, $www_MIG_gff_file: $!";



exit(0);


#--------------------------


sub numerically
{
	$a <=> $b;
}