#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my ( $matrix, $files ) = get_matrix( $options );
print_matrix( $matrix, $options, $files );
exit $?;


sub parse_options {
	my $options = {};
	GetOptions( $options, 'file|f=s@', 'column|l', 'htseq|t', 'cufflinks|c', 'remove_ERCC_ids|e', 'help|h' );
	unless( $$options{ 'file' } ) {
		my $usage = "$0 <file|-f> [--file|-f] [--column|-l <1 or 2 or 3; default=3>] [--htseq|-t] [--cufflinks|-c] [--remove_ERCC_ids|-e]";
		print STDERR $usage, "\n";
		exit 1;
	}
	unless( $$options{ 'column' } ) {
		$$options{ 'column' } = 3;
	}
	return $options;
}

sub get_matrix {
	my( $options ) = @_;
	my $basenames = [];
	my $matrix = {};
	my @files = @{ $$options{ 'file' } };
	foreach my $file( @files ) {
		open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
		my $file_base  = basename( $file );
		if( $$options{ 'cufflinks' } ) {
			$file_base = basename( dirname( $file ) );
		} elsif( $$options{ 'htseq' } ) {
			$file_base =~ s/\.htseq\.counts//;
		} else {
			$file_base =~ s/\.counts\.tab//;
		}
		push @$basenames, $file_base;
		while( my $line = <FH> ) {
			chomp $line;
			if( $$options{ 'cufflinks' } ) {
				my @array_of_vals = split( "\t", $line );
				#my( $gene_id, $fpkm ) = ( $array_of_vals[ 0 ], $array_of_vals[ 9 ] );
				my( $gene_id, $fpkm ) = ( $array_of_vals[ 4 ], $array_of_vals[ 9 ] );
				if( not substr( $gene_id, 0, 2 ) eq '__' ) {
					if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
						# do nothing
					} else {
						$$matrix{ $gene_id }{ $file_base } = $fpkm;
					}
				}
			} elsif( $$options{ 'htseq' } ) {
				my( $gene_id, $count ) = split( "\t", $line );
				if( not substr( $gene_id, 0, 2 ) eq '__' ) {
                                        if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
                                                # do nothing
                                        } else {
						$$matrix{ $gene_id }{ $file_base } = $count;
					}
				}
			} else {
				my @array_of_vals = split( "\t", $line );
				my( $gene_id, $count ) = @array_of_vals[ 0, $$options{ 'column' } ];
                                if( not substr( $gene_id, 0, 2 ) eq 'N_' ) {
                                        if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
                                                # do nothing
                                        } else {
                                                $$matrix{ $gene_id }{ $file_base } = $count;
                                        }
                                }

			}
		}
		close FH or die "Error in closing the file, $file, $!\n";
	}
	return $matrix, $basenames;
}

sub print_matrix {
	my( $matrix, $options, $files ) = @_;
	print STDOUT "Gene_ID", ",", join( ",", @$files ), "\n";
	foreach my $gene_id( keys %$matrix ) {
		my @counts = ();
		foreach my $file_base( @$files ) { 
			if( exists $$matrix{ $gene_id }{ $file_base } ) {
				push @counts, $$matrix{ $gene_id }{ $file_base };
			}
			else {
				push @counts, 0;
			}
		}
		print STDOUT join( ",", ( $gene_id, @counts ) ), "\n";
	}
}
