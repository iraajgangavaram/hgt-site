#!/usr/bin/env perl
# Run LGTSeek's own read classification (LGT::LGTSeek->bwaPostProcess) on a
# donor BAM and a host BAM made from the same reads in the same order.
# Usage: run_lgtseek.pl <lgtseek lib dir> <donor.bam> <host.bam> <output dir> <samtools>
use strict;
use warnings;

my ($lib, $donor_bam, $host_bam, $out_dir, $samtools) = @ARGV;
die "usage: $0 lib donor.bam host.bam out_dir samtools\n" unless $samtools;
unshift @INC, $lib;
require LGT::LGTSeek;

my $lgtseek = LGT::LGTSeek->new({ samtools_bin => $samtools, output_dir => $out_dir, verbose => 0 });
my $result = $lgtseek->bwaPostProcess({
    donor_bam     => $donor_bam,
    host_bam      => $host_bam,
    output_dir    => $out_dir,
    output_prefix => 'lgtseek',
});
for my $class (sort keys %{ $result->{counts} || {} }) {
    my $n = $result->{counts}{$class} // 0;
    print "$class\t$n\n";
}
