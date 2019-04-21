#!/usr/bin/perl
use strict;
{
    opendir(my $DIR,"Perihelions") || die;
    my %angles=();
    while (my $angle=readdir($DIR)){
	if ($angle=~/^[0-9,.E]+$/ && $angle=~ /[0-9]/)
	{	$angles{$angle}=1;}
    }
    for my $angle (sort {$a<=>$b} keys %angles){
	opendir(my $D,"Perihelions/$angle/") || die;
	print "$angle\n";
	open(my $IN,"<Perihelions/$angle/schw_perihelion.dat") || die;
	while(my $line= <$IN>){
	    print $line;
	}
	close $IN;
    }
}
