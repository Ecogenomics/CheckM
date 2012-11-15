.                                                                             
.                888b     d888  .d8888b.   .d8888b.  888    d8P                 
.                8888b   d8888 d88P  Y88b d88P  Y88b 888   d8P                  
.                88888b.d88888 888    888        888 888  d8P                 
.                888Y88888P888 888             .d88P 888d88K                  
.                888 Y888P 888 888         .od888P"  8888888b                 
.                888  Y8P  888 888    888 d88P"      888  Y88b                
.                888   "   888 Y88b  d88P 888"       888   Y88b               
.                888       888  "Y8888P"  888888888  888    Y88b              
.                                                                            

Overview
=========

MetaChecka is a metagenomic QA toolset. It uses sets of reference
genes to test for bin completeness and contamination

Installation
=========

Should be as simple as

    pip install MetaChecka2000

Data preparation and running MetaChecka2000
=========

MetaChecka2000 expects to see a folder containing one or more multiple fasta 
files containg assembled and binned contigs, one file for each bin.

You need to:

  1. Extract genes from the contigs
  2. Use hmmer to match these genes to an existing database
  3. Do some counting

Use: metachecka2000 -h for more detailed help

Licence and referencing
=========

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/MetaChecka2000

This software is currently unpublished. Please contact me at m_dot_imelfort_at_uq_dot_edu_dot_au for more information about referencing this software.

Copyright © 2012 Michael Imelfort, Connor Skennerton. See LICENSE.txt for further details.
