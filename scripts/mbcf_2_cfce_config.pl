#!/usr/bin/env perl
#
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use YAML::Tiny;

my $options = parse_options();
my $yaml = YAML::Tiny -> read( $$options{ 'config' } );
my $meta_info = undef;
($yaml, $meta_info ) = get_meta_info( $yaml, $$options{'metasheet'} );
print_meta_info( "metasheet.csv", $meta_info );
$yaml -> write( "config.yaml" );
exit $?;

sub print_meta_info {
	my( $meta_file, $meta_info ) = @_;
	open( OFH, ">$meta_file" ) or die "Error in writing to the file, $meta_file, $!\n";
	print OFH join( "\n", @$meta_info );
	close OFH or die "Error in closingt he file, $meta_file, $!\n";
}

sub parse_options {
	my $options = {};
	GetOptions( $options, 'config|c=s', 'metasheet|m=s', 'help|h' );
	unless( $$options{ 'config' } or $$options{ 'metasheet' } ) {
		print STDERR "Usage: $0 <--config|-c> <--metasheet|-m>\n";
		exit 1;
	}
	return $options;
}

sub get_meta_info {
	my( $yaml, $meta_file ) = @_;
	my $meta_info = [];
	open( FH, "<$meta_file" ) or die "Error in opening the file, $meta_file, $!\n";
	my $header = <FH>;
	chomp $header;
	my @head = split( ",", $header );
	push @$meta_info, join(",", @head[1..scalar @head - 1] );
	my $sample_info = {};
	while( my $line = <FH> ) {
		chomp $line;
		if( $line =~ /^\s+$/ ) {
			next;
		} else {
			my( $left_mate, $sample, @rest ) = split( ",", $line );
			push @$meta_info, join( ",", ( $sample, @rest ) );
			my $right_mate = $left_mate;
			$right_mate =~ s/_R1_/_R2_/;
			my $cur_dir = getcwd();
			if( -f './data/' . $right_mate ) {
				$$sample_info{ 'samples' }{  $sample } = [ './data/' . $left_mate, './data/' . $right_mate ];
			} else {
				$$sample_info{ 'samples' }{ $sample } = [ './data/' . $left_mate ];
			}
		}
	}
	close FH or die "Error in closing the file, $meta_file, $!\n";
	$$yaml[0]->{'samples'} = $$sample_info{'samples'};
	return ( $yaml, $meta_info );
}
