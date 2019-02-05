#!/root/miniconda3/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

$| = 1;    # Flush output
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].cds.fasta" );
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].pep.fasta" );

my $file_fasta = $ARGV[0];
my $db = Bio::DB::Fasta->new($file_fasta);

my @mRNA;
my $mRNA_name;
my $frame;
my $flag = 0;

my %site;
my (%genetic_code) = get_genetic_code();

my ($codon2syn_site, $codon2non_site) = get_codon2syn_site_AND_codon2non_site();
my %codon2syn_site = %$codon2syn_site;
my %codon2non_site = %$codon2non_site;

my @atgc=('A', 'T', 'G', 'C');
foreach my $codon1(keys %genetic_code){
    my @a = (0,0,0);
    my @codon = split //, $codon1;
    for (my $i=0; $i<3;$i++){
        my @temp = @codon;
        foreach my $base (@atgc){
            next if ($base eq $codon[$i]);
            $temp[$i] = $base;
            my $codon2 = join '', @temp;
            if ($genetic_code{$codon1} eq $genetic_code{$codon2}){
                $a[$i]=$a[$i]+1/3;
            }
        }
    }
    $site{$codon1}=\@a;
}

my $n=0;
foreach (keys %site){
    $n++;
    if ($genetic_code{$_} eq '_'){
        #print "$n\t$_:$genetic_code{$_} $site{$_}->[0], $site{$_}->[1], $site{$_}->[2]\n";
    }
    else{
 #       print "$n\t$_:$genetic_code{$_} $site{$_}->[0], $site{$_}->[1], $site{$_}->[2]\tsyn:$codon2syn_site{$_}\tnonsyn:$codon2non_site{$_}\n" if (($site{$_}->[0] + $site{$_}->[1] + $site{$_}->[2]) != $codon2syn_site{$_});
    }
}


open GFF, "<$ARGV[1]" or die $!;
while ( my $line = <GFF> ) {
	next unless ($line =~ /CDS|mRNA/);
	chomp $line;
	my @array = split( "\t", $line );
	my $type = $array[2];
	next if ( $type eq 'exon' or $type eq 'UTR' );
	if ( (( $type eq 'mRNA' ) and ( $flag > 0 ))or((eof) and ( $flag > 0) ) ) {

		# Collect CDSs and extract sequence of the previous mRNA
		my $mRNA_seq;
        my $chr = "";
		foreach my $coord (@mRNA) {
			my @cds_coord = split( " ", $coord );
			my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
			$mRNA_seq .= $cds_seq;
            $chr = $cds_coord[0];
		}	
		
		if(defined $mRNA_seq){
			my $output_nucleotide = Bio::Seq->new(
				-seq        => $mRNA_seq,
				-id         => $mRNA_name,
				-display_id => $mRNA_name,
				-alphabet   => 'dna',
			);

			if ($frame eq '-') {$output_nucleotide = $output_nucleotide->revcom();} 
			my $output_protein = $output_nucleotide->translate();
		
			########## output cds and protein sequences ###################
			$outfile_cds->write_seq($output_nucleotide);
			$outfile_pep->write_seq($output_protein);		
			###############################################################	
            my @temp;
            my ($s,$end,$i);
            $s=0; $end =0;$i=0;
            foreach my $coord (@mRNA){
                my @cds_coord = split( " ", $coord );
                $end = $s + $cds_coord[2]-$cds_coord[1];
                for($i=$s;$i <= $end;$i++){
                    $temp[$i]=$i-$s+$cds_coord[1];
                }
                $s=$i;
            }
            my @i2pos;
            if ($frame eq '-'){
                @i2pos = reverse @temp;
            }else{
                @i2pos = @temp;
            }

			###############################################################
			my $cds = $output_nucleotide->seq();
			my $len_cds = length $cds;
			my @base;

#			print ">$mRNA_name\n";
			foreach my $index ( 0 .. $len_cds - 1){
        		$base[$index] = substr($cds,$index,1);
    		}
    		my $triple_switch = 1;
    		my $codon = $base[0];

			foreach my $index ( 1 .. $len_cds - 1){
        		my $b = $base[$index];

				if ($b ne '-'){
					if ($triple_switch < 3){
						#print $codon1, '.=', $b1, "\n"; sleep 1;
						$codon .= $b;
						$triple_switch ++;
					}
					if ($triple_switch == 3){
                        $codon = uc($codon);
						print $mRNA_name,"\t","|$codon|","\t",$chr,"\t",$i2pos[$index-2],"\t",$site{$codon}->[0],"\n" if ($codon !~ /N/);
                        print $mRNA_name,"\t","|$codon|","\t",$chr,"\t",$i2pos[$index-1],"\t",$site{$codon}->[1],"\n" if ($codon !~ /N/);
                        print $mRNA_name,"\t","|$codon|","\t",$chr,"\t",$i2pos[$index],"\t",$site{$codon}->[2],"\n" if ($codon !~ /N/);
						$codon = '';
						$triple_switch = 0;

					}
        		}
    		}
    		#print "\n";

			###############################################################
		}

		# Now initialize the next mRNA
		my @attrs = split( ";", $array[8] );
		$attrs[0] =~ s/ID=//;
		$mRNA_name = $attrs[0];
		$frame=$array[6];
		@mRNA = (); # Empty the mRNA
		$flag ++;
	}

	elsif ( $type eq 'mRNA' ) {    # First mRNA
		my @attrs = split( ";", $array[8] );
		$attrs[0] =~ s/ID=//;
		$mRNA_name = $attrs[0];
		$frame=$array[6];
		$flag = 1;
	}
	elsif ( $type eq 'CDS' ) {
		my $cds_coord = $array[0] . " " . $array[3] . " " . $array[4];
		push( @mRNA, $cds_coord );
	}
}

close GFF;

sub codon2aa {
    my ($codon) = @_;
    $codon = uc $codon;
    my (%genetic_code) = get_genetic_code();
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }elsif($codon =~ /N/){
    	$genetic_code{$codon}="X";
    	return $genetic_code{$codon};
    }else{
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}

# a hash connecting codon and amino acid
sub get_genetic_code {
    my (%genetic_code) = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
    return %genetic_code;
}

# get synonymous and non-synonymous site for codons
sub get_codon2syn_site_AND_codon2non_site {
    my %codon2syn_site = ();
    my %codon2non_site = ();
    %codon2syn_site = (
        'AAA' => 1/3,
        'AAC' => 1/3,
        'AAG' => 1/3,
        'AAT' => 1/3,
        'ACA' => 1,
        'ACC' => 1,
        'ACG' => 1,
        'ACT' => 1,
        'AGA' => 5/6,
        'AGC' => 1/3,
        'AGG' => 2/3,
        'AGT' => 1/3,
        'ATA' => 2/3,
        'ATC' => 2/3,
        'ATG' => 0,
        'ATT' => 2/3,
        'CAA' => 1/3,
        'CAC' => 1/3,
        'CAG' => 1/3,
        'CAT' => 1/3,
        'CCA' => 1,
        'CCC' => 1,
        'CCG' => 1,
        'CCT' => 1,
        'CGA' => 3/2,
        'CGC' => 1,
        'CGG' => 4/3,
        'CGT' => 1,
        'CTA' => 4/3,
        'CTC' => 1,
        'CTG' => 4/3,
        'CTT' => 1,
        'GAA' => 1/3,
        'GAC' => 1/3,
        'GAG' => 1/3,
        'GAT' => 1/3,
        'GCA' => 1,
        'GCC' => 1,
        'GCG' => 1,
        'GCT' => 1,
        'GGA' => 1,
        'GGC' => 1,
        'GGG' => 1,
        'GGT' => 1,
        'GTA' => 1,
        'GTC' => 1,
        'GTG' => 1,
        'GTT' => 1,
        'TAC' => 1,
        'TAT' => 1,
        'TCA' => 1,
        'TCC' => 1,
        'TCG' => 1,
        'TCT' => 1,
        'TGC' => 1/2,
        'TGG' => 0,
        'TGT' => 1/2,
        'TTA' => 2/3,
        'TTC' => 1/3,
        'TTG' => 2/3,
        'TTT' => 1/3
    );
    foreach my $codon (keys %codon2syn_site){
        $codon2non_site{$codon} = 3 - $codon2syn_site{$codon};
    }
    return (\%codon2syn_site, \%codon2non_site);
}
