#!/usr/bin/env perl
###############################################################################
#
#    bin-contamination
#    
#    Check for multiple copies of known single-copy marker genes
#
#    Copyright (C) 2012 Connor Skennerton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;

#CPAN modules

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
my $global_options = checkParams();
unless($global_options->{'quiet'}) {
    printAtStart();
}

if($global_options->{'verbose'}) {
    $global_log_commands = 1;
}

# list of phylosift markers
my %markers;
$markers{'PMPROK00003'}="50S Ribosomal Protein L14P (rpl14P)";
$markers{'PMPROK00014'}="O-sialoglycoprotein endopeptidase (gcp)";
$markers{'PMPROK00015'}="30S Ribosomal Protein S3P (rps3P)";
$markers{'PMPROK00019'}="50S Ribosomal Protein S5 (rps5P)";
$markers{'PMPROK00020'}="50S Ribosomal Protein L6P (rpl6P)";
$markers{'PMPROK00022'}="50S Ribosomal Protein L1P (rpl1P)";
$markers{'PMPROK00024'}="50S Ribosomal Protein L22P (rpl22P)";
$markers{'PMPROK00025'}="30S Ribosomal Protein S19P (rps19P)";
$markers{'PMPROK00028'}="30S Ribosmal Protein S2P (rps2P)";
$markers{'PMPROK00029'}="50S Ribosomal Protein L11P (rpl11P)";
$markers{'PMPROK00031'}="GTP-binding signal recognition particle G-domain protein (srp54)";
$markers{'PMPROK00034'}="50S Ribosomal Protein L4P (rpl4P)";
$markers{'PMPROK00041'}="30S Ribsomal Protein S17P (rps17P)";
$markers{'PMPROK00047'}="CTP synthase, UTP-ammonia lyase (pyrG)";
$markers{'PMPROK00048'}="50S Ribosomal Protein L2P (rpl2P)";
$markers{'PMPROK00050'}="30S Ribosomal Protein S15P/S13e (rps15P)";
$markers{'PMPROK00051'}="30S Ribosomal Protein S9P/S13 (rps9P)";
$markers{'PMPROK00052'}="50S Ribsomal Protein L18P (rpl18P)";
$markers{'PMPROK00053'}="50S Ribosomal Protein L5P (rpl5P)";
$markers{'PMPROK00054'}="30S Ribosomal Protein S7P (rps7P)";
$markers{'PMPROK00060'}="Translation elongation factor aEF-2";
$markers{'PMPROK00064'}="Translation initiation factor aIF-2";
$markers{'PMPROK00067'}="50S Ribosomal Protein L24P (rpl24P)";
$markers{'PMPROK00068'}="30S Ribosomal Protein S11P (rps11P)";
$markers{'PMPROK00069'}="50S Ribosomal Protein L10E (rp10E)";
$markers{'PMPROK00071'}="50S Ribosomal Protein L16/L10E";
$markers{'PMPROK00074'}="30S Ribosomal Protein S8P (rps8E)";
$markers{'PMPROK00075'}="50S Ribosomal Protein L3P (rpl3P)";
$markers{'PMPROK00081'}="30S Ribosomal Protein S13P (rps13P)";
$markers{'PMPROK00086'}="Phenylalanyl-tRNA synthetase, beta subunit (pheT)";
$markers{'PMPROK00087'}="Phenylalanyl-tRNA synthetase, alpha subunit (pheS)";
$markers{'PMPROK00092'}="50S Ribosomal Protein L15P (rpl15P)";
$markers{'PMPROK00093'}="30S Ribosomal Protein S10P (rps10P)";
$markers{'PMPROK00094'}="30S Ribosomal Protein S12P (rps12P)";
$markers{'PMPROK00097'}="Triosephosphate isomerase (tpiA)";
$markers{'PMPROK00106'}="tRNA-ribosyltransferase - tRNA-guanine transglycosylase";
$markers{'PMPROK00123'}="rdgB/HAM1 protein";
$markers{'PMPROK00126'}="Ribonuclease HII";
my %colours = (
    FG_BLACK=>"\033[0;30m",
    FG_RED=>"\033[0;31m",
    FG_GREEN=>"\033[0;32m",
    FG_YELLOW=>"\033[0;33m",
    FG_BLUE=>"\033[0;34m",
    FG_MAGENTA=>"\033[0;35m",
    FG_CYAN=>"\033[0;36m",
    FG_WHITE=>"\033[0;37m",
    FG_LTBLACK=>"\033[1;30m",
    FG_LTRED=>"\033[1;31m",
    FG_LTGREEN=>"\033[1;32m",
    FG_LTYELLOW=>"\033[1;33m",
    FG_LTBLUE=>"\033[1;34m",
    FG_LTMAGENTA=>"\033[1;35m",
    FG_LTCYAN=>"\033[1;36m",
    FG_LTWHITE=>"\033[1;37m",
    BG_BLACK=>"\033[0;40m",
    BG_RED=>"\033[0;41m",
    BG_GREEN=>"\033[0;42m",
    BG_YELLOW=>"\033[0;43m",
    BG_BLUE=>"\033[0;44m",
    BG_MAGENTA=>"\033[0;45m",
    BG_CYAN=>"\033[0;46m",
    BG_WHITE=>"\033[0;47m",
    FG_BD=>"\033[1m",
    FG_UL=>"\033[4m",
    NOCOLOR=>"\033[0m",
);

my $prefix = overrideDefault('PS_temp','prefix');
if ($global_options->{'runPS'}) {
    foreach my $binfile (@ARGV) {
        # run phylosift search and phylosift align
        checkAndRunCommand('phylosift',[['search', "--output=$prefix/$binfile", $binfile, ">/dev/null"]], DIE_ON_FAILURE);
        checkAndRunCommand('phylosift',[['align', "--output=$prefix/$binfile", $binfile, ">/dev/null"]], DIE_ON_FAILURE);
    }
}

my $outdir;
opendir($outdir, $prefix) || die $!;
my @resultdirs = grep { $_ !~ /^\./ } readdir($outdir);
for my $resultdir (@resultdirs) { 
# check for multiple marker sequences
# search through the output of each of the marker files
    unless ($global_options->{'quiet'} ) {
        print "$resultdir\n";
    }
    my %marker_counts;

    foreach my $marker (keys %markers) {
        # check that the marker has been found
        if (-e "$prefix/$resultdir/alignDir/$marker.fasta") {
            $marker_counts{$marker} = `grep -c '>' $prefix/$resultdir/alignDir/$marker.fasta`;
            chomp $marker_counts{$marker};
        } else {
            $marker_counts{$marker} = undef;
        }
    }
    my $single_count = 0;
    my $zero_count = 0; 
    my$multi_count = 0;
    while (my($marker,$count) = each %marker_counts) {
        if(! defined $count || $count == 0) {
            $zero_count++;
            unless($global_options->{'quiet'}){ 
                print colourize("$marker $markers{$marker} 0",$global_options->{colour},'FG_YELLOW'),"\n";
            }
        } elsif ($count == 1) {
            $single_count++;
            unless($global_options->{'quiet'}){
                print colourize("$marker $markers{$marker} 1",$global_options->{colour},'FG_GREEN'),"\n";
            }
        } else {
            $multi_count++;
            unless($global_options->{'quiet'}){
                print colourize("$marker $markers{$marker} $count",$global_options->{colour},'FG_RED'),"\n";
            }
        }
    }
    if($global_options->{'quiet'}){
        print colourize($resultdir ,$global_options->{colour},'FG_MAGENTA')," ";
    }
    print colourize($single_count,$global_options->{colour},'FG_GREEN'), ' ';
    print colourize($zero_count, $global_options->{colour}, 'FG_YELLOW'), ' ';
    print colourize($multi_count,$global_options->{colour}, 'FG_RED'), "\n";

}
exit;
######################################################################
# CUSTOM SUBS

######################################################################

sub colourize {
    my ($exp, $is_colouring, $color) = @_;
    if($is_colouring) {
        return $colours{$color}.$exp.$colours{NOCOLOR};
    } else {
        return $exp;
    }
}


######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "runPS+","quiet","prefix:s", "colour+", "verbose+");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to 
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME


   bin-contamination

=head1 COPYRIGHT

   copyright (C) Connor Skennerton

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   bin-contamination is a wrapper around phylosift for checking metagenomic
   contig sets for multiple copies of any of the marders provided by pyhlosift.
   The output can be a simple count of the number of single, duplicate or 
   unfound genes, or a more detailed list of those markers.

=head1 SYNOPSIS

    bin-contamination  [-help|h] [-runPS] [-prefix STRING] [-quiet]

      [-help -h]                   Displays basic usage information
      [-runPS]                     Run Phylosift on the input files to generate results
      [-prefix STRING]             Output folder for the results
      [-quiet]                     Output only the aggregate counts for each bin
      [-verbose]                   Print alot of stuff to screen
      [-noColour]                  Don't output anything colourful
         
=cut
