#!/usr/bin/env perl
#
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my $de_info = get_de_info( $$options{ 'file' } );
print_info( $de_info );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'file|f=s@', 'help|h' );
	unless( $$options{ 'file' } ) {
		my $usage = "$0 <--file|-f> [--file|-f ...]";
		print $usage, "\n";
	}
	return $options;
}

sub get_de_info {
	my( $file_list ) = @_;
	my $info = {};
	foreach my $file( @$file_list ){
		my( $comp ) = ( basename( $file ) =~ /(.+)\.deseq\.csv/ );
		open( FH, "<$file" ) or die "Error opening the file, $file, $!\n";
		my $header = <FH>;
		while( my $line = <FH> ) {
			chomp $line;
			my( $id, $mean, $logfc, $lfsce, $stat, $p, $p_adj ) = split( ",", $line );
			next if( $p_adj eq 'NA' );
			if( $p_adj <= 0.1 ) {
				if( $logfc <= -2 ) {
					$$info{ $comp }{ 'down_p1_log2' } += 1;
				} elsif( $logfc >= 2 ) {
					$$info{ $comp }{ 'up_p1_log2' } += 1;
				}
				if( $logfc <= -1 ) {
					$$info{ $comp }{ 'down_p1_log1' } += 1;	
				} elsif( $logfc >= 1 ) {
					$$info{ $comp }{ 'up_p1_log1' } += 1;
				}
			}
			if( $p_adj <= 0.05 ) {
				if( $logfc <= -2 ) {   
                                        $$info{ $comp }{ 'down_p05_log2' } += 1;
                                } elsif( $logfc >= 2 ) {
                                        $$info{ $comp }{ 'up_p05_log2' } += 1;
                                }
				if( $logfc <= -1 ) {
                                        $$info{ $comp }{ 'down_p05_log1' } += 1;
                                } elsif( $logfc >= 1 ) {
                                        $$info{ $comp }{ 'up_p05_log1' } += 1;
                                }
			}
			if( $p_adj <= 0.01 ) {
				if( $logfc <= -2 ) {   
                                        $$info{ $comp }{ 'down_p01_log2' } += 1;
                                } elsif( $logfc >= 2 ) {
                                        $$info{ $comp }{ 'up_p01_log2' } += 1;
                                }
				if( $logfc <= -1 ) {
                                        $$info{ $comp }{ 'down_p01_log1' } += 1;
                                } elsif( $logfc >= 1 ) {
                                        $$info{ $comp }{ 'up_p01_log1' } += 1;
                                }
			}
			$$info{ $comp }{ 'down_p1_log2' } += 0;
			$$info{ $comp }{ 'up_p1_log2' } += 0;
			$$info{ $comp }{ 'down_p05_log2' } += 0;
			$$info{ $comp }{ 'up_p05_log2' } += 0;
			$$info{ $comp }{ 'down_p01_log2' } += 0;
			$$info{ $comp }{ 'up_p01_log2' } += 0;
			$$info{ $comp }{ 'down_p1_log1' } += 0;
                        $$info{ $comp }{ 'up_p1_log1' } += 0;
                        $$info{ $comp }{ 'down_p05_log1' } += 0;
                        $$info{ $comp }{ 'up_p05_log1' } += 0;
                        $$info{ $comp }{ 'down_p01_log1' } += 0;
                        $$info{ $comp }{ 'up_p01_log1' } += 0;
		}
		close FH or die "Error closign the file, $file, $!\n"; 
	}
	return $info;
}

sub print_info {
	my( $info ) = @_;
	my @pvals = qw(up_p1_log2 down_p1_log2 up_p05_log2 down_p05_log2 up_p01_log2 down_p01_log2 up_p1_log1 down_p1_log1 up_p05_log1 down_p05_log1 up_p01_log1 down_p01_log1);
	my $header = join(",", ("Comparison", @pvals));
	print STDOUT $header,"\n";
	foreach my $comp( keys %$info ) {
		my @vals = ();
		foreach my $p( @pvals ) {
			push @vals, $$info{$comp}{$p} || 0;
		}
		print STDOUT join(",",($comp,@vals)),"\n";
	}
}
