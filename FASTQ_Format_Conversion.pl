#!/usr/bin/perl
# Created by Sumeet Gupta (sgupta@wi.mit.edu) - 21st June 2016

# Converts FASTQ read header/indentifier format

# Automatically identifies the format of the read id and changes it between the 2 following 2 formats
# read id format 1: @HWUSI-EAS100R:6:73:941:1973#0/1
# read id format 2: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

# FASTQ files from WI Genome core end with ";1" or ";2" in format 1 specified above - This script can handle this variation

use warnings;
use strict;
use Carp;

my $randomrunid = int(rand(10)*100);

my $USAGE = "Usage: perl FASTQ_Format_Conversion.pl <FASTQ file>\n\n";
$USAGE .= "Example1: perl FASTQ_Format_Conversion.pl s_1_sequence.txt.tar.gz > Output.filename\n\n";
$USAGE .= "Example2: perl FASTQ_Format_Conversion.pl s_1_sequence.txt.gz > Output.filename\n\n";
$USAGE .= "Example3: perl FASTQ_Format_Conversion.pl s_1_sequence.txt > Output.filename\n\n";
$USAGE .= "Example3: perl FASTQ_Format_Conversion.pl s_1_sequence.fq > Output.filename\n\n";

if( (!defined($ARGV[0])) || (scalar(@ARGV) != 1)) { 
	print $USAGE."\n"; 
	exit; 
}

$ARGV[0] =~ s/[\n\r]//g;
my $pipecmd = compressioncheck($ARGV[0]);

if($pipecmd !~ m/cat/){
	print $pipecmd;
	exit;
}

open(my $PIPEIN, '-|', $pipecmd) or die "Opening pipe [$pipecmd]: $!\n";

while(my $newline = <$PIPEIN>) {

	$newline =~ s/[\n\r]//g;
	my $id = $newline;
	
	my $formatconv;
	if($id =~ m/\#/){
		$formatconv = 1;
	}
	elsif($id =~ m/\s/) {
		$formatconv = 2;
	}
	
	my $newid = changeidformat($id,$formatconv,$randomrunid);
	print $newid."\n";
	
	$newline = <$PIPEIN>;
	print $newline;
	$newline = <$PIPEIN>;
	$newline =~ s/[\n\r]//g;
	$newid = changeidformat($newline,$formatconv,$randomrunid);
	print $newid."\n";
	$newline = <$PIPEIN>;
	print $newline;

}

close ($PIPEIN);

sub changeidformat {
	
	my ($id,$convtype,$runid) = @_;
	
	my $formatedid;
	my @seqarray = ();
	my @idstats = ();
	my @barcode = ();
	
	if($convtype == 1){
		
		@seqarray = split(/\#/, $id);
		@barcode = split(/\//, $seqarray[1]);
		@idstats = split(/:/, $seqarray[0]);
		my $readnumber;
		my $passfiler;
		if($barcode[1] =~ m/;/){
			my @readandpf = split(/;/, $barcode[1]);
			$readnumber = $readandpf[0];
			$passfiler ="Y";
			if($readandpf[1] != 1){
				$passfiler = "N";
			}
		}
		else {
			$readnumber = $barcode[1];
			$passfiler ="Y";
		}
		
		$formatedid = $idstats[0].":".$runid.":Unknown:".$idstats[1].":".$idstats[2].":".$idstats[3].":".$idstats[4]." ".$readnumber.":".$passfiler.":0:".$barcode[0];
		
	}
	elsif($convtype == 2){
		
		@seqarray = split(/\s/, $id);
		@idstats = split(/:/, $seqarray[0]);
		@barcode = split(/:/, $seqarray[1]);
		
		$formatedid = $idstats[0].":".$idstats[3].":".$idstats[4].":".$idstats[5].":".$idstats[6]."#".$barcode[3]."/".$barcode[0];
		
	}
	
	return $formatedid;
	
}

sub compressioncheck {
	
	my ($file) = @_;
	my $pipecmd;
		
	if($file =~ m/.tar.gz$/){
		$pipecmd = "zcat $file | tar  -O -xf -";
	}
	elsif($file =~ m/.gz$/) {
		$pipecmd = "zcat $file"; 
	}
	elsif($file =~ m/.txt$/ || $file =~ m/.fq$/) {
		$pipecmd = "cat $file"; 
	}
	else {
		$pipecmd = "================\n\nFilename does not end with txt or .gz or tar.gz or fq\n\n================\n\nExiting Attempting Format Change!!!!!!!\n\n================\n\n\n";
	}
	
	return $pipecmd;
	
}

