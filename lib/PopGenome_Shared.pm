package PopGenome_Shared;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use YAML::Tiny;
use Bio::SeqIO;
#----------------------------------- other subroutines ---------------------------------------------
sub LOADREF{
	my ($reference) = @_;
	my %genome;
	open REF, "$reference";
	my $id;
	while (<REF>){
		chomp;
		if (/\>(\S+)/){
			$id = $1;
			$genome{seq}{$id}="";
		}
		elsif(/(\w+)/){
			$genome{seq}{$id}.=$1;
		}
	}
	close REF;
	$genome{sumlen}=0;
	foreach my $id(keys %{$genome{seq}}){
		$genome{len}{$id}=length($genome{seq}{$id});
		$genome{sumlen} += $genome{len}{$id};
	}
	return \%genome;
}

1;