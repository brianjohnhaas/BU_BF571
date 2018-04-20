#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Process_cmd;

my $usage = "usage: $0 matlab.m\n\n";

my $matlab_script = $ARGV[0] or die $usage;


# matlab -nodesktop -nodisplay -nosplash -r "run('/seq/RNASEQ/bf5seg/small_test/test.m');quit"


main: {

    $matlab_script = &ensure_full_path($matlab_script);

    my $cmd = "matlab -nodesktop -nodisplay -nosplash -r \"run('$matlab_script');quit\" ";

    &process_cmd($cmd);

    exit(0);
}

