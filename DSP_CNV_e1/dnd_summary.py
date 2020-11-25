#!/usr/bin/python3
"""
dnd_summary.py 

DnD: Dsp ngs Dpp
Digital Spatial Profiling Next Generation Sequencing Data Processing Pipeline Summary Generator
Version: 0.1

DM: Zach Norgaard
Last Updated: 07January2019
"""
# General Packages
import logging
import argparse
import glob
import os

# Custom Packages 
import dndFxns

# Logging Pattern
logging.basicConfig(level=logging.DEBUG, format=('%(asctime)s - %(levelname)s '
    '- %(message)s'))

# Arguments
parser = argparse.ArgumentParser()

# Not a list of files so the command line character limit can be avoided.
parser.add_argument('-dcd', '--DCC_Directory', default='.', 
                    help=('Directory containing DCC Files to be summarized.'
                          ' Default is the current directory.'))

parser.add_argument('-o', '--output', default='',
                    help=('Root name of output files.'))

args = parser.parse_args()

# Get List of DCC files
dccfiles = glob.glob(os.path.join(args.DCC_Directory, '*.dcc'))

# Read Files
logging.debug('Reading %d dcc files.' % len(dccfiles))
dccobjs = [ dndFxns.dccfile(dcc) for dcc in dccfiles ]

# Create Deduplicated Counts Table and save
logging.debug('Creating deduplicated counts table.')
df = dndFxns.generateCountTable(dccobjs)
dedupname = 'dedup.txt' if args.output=='' else '%s_dedup.txt' % args.output
df.to_csv(dedupname, sep='\t', index=False)

# Create Read Summary Table and save
logging.debug('Creating processing summary table.')
sumdf = dndFxns.pipelineSummaryTable(dccobjs)
summaryname = 'summary.txt' if args.output=='' else '%s_summary.txt' % args.output
sumdf.to_csv(summaryname, sep='\t', index=False)

logging.debug('DONE')
