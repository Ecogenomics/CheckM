#!/usr/bin/perl -w
# Compares bootstrap and length for splits that match
use strict;
use lib "$ENV{HOME}/Genomics/Perl/modules";
use MOTree;
use Getopt::Long;

{
    my $usage = "CompareToBootstrap.pl -tree treeFile1 -boot treeFile2 [-debug]\n"
	. " Compares bootstrap replicate trees in file2 to the tree in file1\n"
	. " Outputs a new version of treeFile1 with bootstrap proportions\n";

    my $file1 = undef;
    my $file2 = undef;
    my $debug = 0;
    if (@ARGV==2 && -e $ARGV[0] && -e $ARGV[1]) {
	($file1,$file2) = @ARGV;
    } else {
	die $usage unless GetOptions('tree=s' => \$file1,
				     'boot=s' => \$file2,
				     'debug' => \$debug)
	    && defined $file1 
	    && defined $file2
	    && @ARGV == 0;
    }
    open(FILE1,"<",$file1) || die "Cannot read $file1";
    open(FILE2,"<",$file2) || die "Cannot read $file2";
    my $tree1 = MOTree::new(fh => \*FILE1);
    die "Cannot read a tree from $file1" unless $tree1;

    # remove leading and trailing underscores (an artefact of using rose)
    foreach my $node (@{ $tree1->depthfirst() }) {
	if ($tree1->is_Leaf($node)) {
	    my $id = $tree1->id($node);
	    $id =~ s/^_+//;
	    $id =~ s/_+$//;
	    $tree1->set_id($node,$id);
	}
    }

    my $nDesc1 = $tree1->countAllLeaves(); # node => number of leaves beneath it
    my $totLeaf1 = $nDesc1->{$tree1->get_root_node()};
    print STDERR "$totLeaf1 leaves in input\n" if $debug;

    my $nTree = 0;
    my %support = (); # internal node in tree1 => count of times seen
    while (my $tree2 = MOTree::new(fh => \*FILE2)) {
	$nTree++;
	print STDERR "Read bootstrap $nTree\n" if $debug;

	foreach my $node (@{ $tree2->depthfirst() }) {
	    if ($tree2->is_Leaf($node)) {
		my $id = $tree2->id($node);
		if (defined $id) {
		    $id =~ s/^_+//;
		    $id =~ s/_+$//;
		    $tree2->set_id($node,$id);
		}
	    }
	}
	print STDERR "Tweaked ids for bootstrap $nTree\n" if $debug;
	my $match = $tree1->Compare($tree2);
	print STDERR "Matched tree $nTree\n" if $debug;
	foreach my $node1 (keys %$match) {
	    $support{$node1}++;
	}
    }
    die "No trees in $file2" if $nTree == 0;
    print STDERR "Read $nTree trees in $file2\n";
    foreach my $node1 (@{ $tree1->depthfirst() }) {
	next if $tree1->is_Leaf($node1);
	if ($nDesc1->{$node1} > 1 && $nDesc1->{$node1} < $totLeaf1 - 1) {
	    my $nsupport = $support{$node1} || 0;
	    $tree1->set_id($node1,sprintf("%.3f",$nsupport/$nTree));
	} else {
	    $tree1->set_id($node1,"");
	}
    }
    print $tree1->toNewick() . "\n";
}
