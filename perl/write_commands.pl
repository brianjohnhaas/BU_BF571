#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

my $matlab_runner = "$FindBin::Bin/matlab_runner.pl";

foreach my $file (<*.m>) {
    print "$matlab_runner $file\n";
}

exit(0);

