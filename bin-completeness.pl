#!/usr/bin/env perl
###############################################################################
#
#    bin-completeness
#    
#    <one line to give the program's name and a brief idea of what it does.>
#
#    Copyright (C) Connor Skennerton
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
use File::Path qw(make_path);
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
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################
my $prefix = overrideDefault('BC_out', 'prefix');
# call prodigal
if ($global_options->{callORFs}) {
    for my $infile (@ARGV) {
        unless (-e "$prefix/$infile") {
            make_path("$prefix/$infile");
        }
        checkAndRunCommand("prodigal", [{
                                        '-a' => "$prefix/$infile/prodigal_out.faa",
                                        '-i' => "$infile"
                                        },
                                        ["-c","-q", ">/dev/null"]
                                    ],
                           DIE_ON_FAILURE);
    }
}

# call Hmmer3

if ($global_options->{runHmmer}) {
    for my $infile (@ARGV) {
        unless (-e "$prefix/$infile") {
            make_path("$prefix/$infile");
        }
        checkAndRunCommand('hmmsearch', 
                [[
                '--tblout',
                "$prefix/$infile/hmmer_out.txt", 
                $global_options->{"hmmerDB"}, 
                "$prefix/$infile/prodigal_out.faa", 
                ">$prefix/$infile/hmmer_out.hmmer3" 
                ]], 
            DIE_ON_FAILURE);
    }
}

my $outdir;
opendir($outdir, $prefix) || die $!;
my @resultsdir = grep { $_ !~ /^\./ } readdir($outdir);

# parse the tabular output from hmmer
for my $infile (@resultsdir) {
    my %found_genes;
    print "$infile\n";
    my $hmmerout = openRead("$prefix/$infile/hmmer_out.txt");
    while (<$hmmerout>) {
        next if (/^#/);
        chomp;
        my @columns = split /\s+/;
        if ($columns[4] <= 1e-10) {
            $found_genes{$columns[2]}++
        }
    }
    while(my($name,$count) = each %found_genes) {
        print "$name $count\n";
    }
    print scalar keys %found_genes, "\n";
}

######################################################################
# CUSTOM SUBS
######################################################################

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "runHmmer+", "callORFs+", "prefix:s", "hmmerDB:s");
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

    bin-completeness

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

   Call orfs, run Hmmer3 to find 111 "essential" genes.  The default behaviour
   is to look in the directory pointed to by -prefix and parse the output results
   from Hmmer.  However the user can also specify to input files of either contigs
   or ORFs, in which case the ORFs will be called automatically and run through Hmmer

=head1 SYNOPSIS

    bin-completeness  [-help|h] [-callORFs] [-runHmmer] [-hmmerDB FILE] -prefix STRING [<input files>]

      [-help -h]                   Displays basic usage information
      [-callORFs]                  Run callORFs script to call the orfs, implies -runHmmer
      [-runHmmer]                  Run Hmmer3 on the called ORFs
      [-hmmerDB FILE]              File of HMMs for hmmer
      -prefix STRING               Output directory for results
      [<input files>]              If either -callORFs -runHmmer are used the input files must be given
         
=cut

