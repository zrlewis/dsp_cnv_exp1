"""
lockedFxns

This is the set of resources, functions, and classes for 
the final customer facing DSP-NGS Data Processing Pipeline.

znorgaard
Last Updated: 07January2019

Functions:
safteryfirst = verify string contains no 'dangerous' symbols'
file_removal = remove list of files without raising critical errors
fileLines = count lines in a file
fastqReads = count sequence reads in fastq file
samSummary = provide unaligned, aligned, and specific target alignment counts
samToDCC = convert SAM to DCC-like file
splitSam = split sam file based on dsp tag id
mergeSAM = merge SAM files
dccToTSV = convert a dcc file to a tsv file (no summary information, just counts)
generateCountTable = take a list of DCC files and generate a table of counts
pipelineSummaryTable = take a list of DCC files and generate a read summary table

Classes:
summaryinfo = object containing sample summary information
dccfile = object containing information of dcc / dcc-like file

"""

import csv
import datetime
import gzip
import collections
import re
import os
import logging
import pandas as pd

# Logging Stuff
logging.basicConfig(level=logging.DEBUG, format=('%(asctime)s - %(levelname)s '
    '- %(message)s'))

file_version = '0.01'
soft_version = '0.01'

# safteryfirst
#    Function to verify a string is safe for execution
#    INPUT:
#        somstring = str
#    OUTPUT:
#        None or Exception if ['&', '|', ';'] detected in string
def safetyfirst(somestring):
    dangerlist = ['&', '|', ';']
    assert not any(x in somestring for x in dangerlist), (
        'Dangerous symbol detected. Execution halted.')

# file_removal
#    Function to remove a file
#    INPUT:
#        filenames = list of file names to remove as strings
#    OUTPUT:
#        NA
def file_removal(filenames):
    for fh in filenames:
        try:
            os.remove(fh)
        except:
            logging.warning('Could not remove %s.' % fh)

# fileLines
#     Function to count lines in a file.
#     INPUT:
#        fName = file name of gzipped or text file
#     OUTPUT:
#        Number of lines in file as int
def fileLines(fName):
    if fName.endswith('.gz'):
        # Open gzipped file
        with gzip.open(fName, 'rb') as fh:
            for i, _line in enumerate(fh):
                pass
    else:
        with open(fName) as fh:
            for i, _line in enumerate(fh):
                pass
    return i + 1
    
# fastqReads
#     Function to count reads in fastq file
#     INPUT:
#        fName = fastq file name (gzipped supported)
#     OUTPUT:
#        Number of sequencing reads as int
def fastqReads(fName):
    try:
        return fileLines(fName) / 4
    except:
        return 0

# samSummary
#    Function to count number of aligned reads in sam file
#        and number of reads aligned to each target
#    INPUT:
#        sam = sam file name
#    OUTPUT:
#        unalinged = number of reads which did not align to
#            a target sequence as an int
#        aligned = number of sequences which did align to
#            a target sequence as an int
#        targets = dictionary with targets (str) as keys and
#            aligned read counts (int) as values
def samSummary(sam):
    unaligned = 0; aligned = 0; targets = {}
    with open(sam) as fh:
        for line in csv.reader(fh, delimiter='\t'):
            if line[0].startswith('@'): 
                pass
            elif line[2] == '*':
                unaligned += 1
            else:
                aligned += 1
                targ = line[2]
                targets[targ] = targets.setdefault(targ, 0) + 1
    return unaligned, aligned, targets

# samToDCC
#    Function to convert a sam file to a nanostring DCC-like file
#    INPUT:
#        samFile = list of sam files to convert
#        dccfname = name of dcc file to output
#        sampID = name of sample for dcc sample attributes section
#        sumObj = instance of summaryInfo object 
#        overcounts = comma separated string of analytes that were
#            deduplicated at hamming distance 0 because they
#            were detected at a frequency greater than the max
#            umi count threshold for standard deduplication
#        owner = file owner
#        comments = any extra information
#        seqKit = sequencing kit used
#        dateTag = date of dcc generation
#    OUTPUT:
#        writes dcc file
def samToDCC(samFile, dccfname, sumObj,
             seqSetID, trimOpts, flash2Opts, 
             umiExtractOpts, bowtie2Opts, umiDedupOpts, umiLimit,
             overcounts = '', plateid='1012207777777', wellid='A01', 
             dateTag=datetime.datetime.today().strftime('%Y-%m-%d')):
    tags = {}
    with open(samFile) as fh:
        for line in csv.reader(fh, delimiter='\t'):
            if line[0].startswith('@'): continue
            dspTag = line[2]
            if dspTag == '*': continue
            tags[dspTag] = tags.setdefault(dspTag, 0) + 1
    out = open(dccfname, 'w')
    out.write(("<Header>\n"
              "FileVersion,%s\n"
              "SoftwareVersion,%s\n"
              "Date,%s\n"
              "</Header>\n\n"
              
              "<Scan_Attributes>\n"
              "ID,%s\n"
              "Plate_ID,%s\n"
              "Well,%s\n"
              "</Scan_Attributes>\n\n"
              
              "<NGS_Processing_Attributes>\n"
              "seqSetID,%s\n"
              "tamperedIni,No\n"
              "trimGaloreOpts,%s\n"
              "flash2Opts,%s\n"
              "umiExtractOpts,%s\n"
              "bowtie2Opts,%s\n"
              "umiDedupOpts,%s\n"
              "umiLimit,%s\n"
              "Raw,%d\n"
              "Trimmed,%d\n"
              "Stitched,%d\n"
              "Aligned,%d\n"
              "Overcounts,%s\n"
              "</NGS_Processing_Attributes>\n\n"
              
              "<Code_Summary>\n"
              "%s\n"
              "</Code_Summary>\n\n"
              
              "SOMEHASH100000000000" % ( file_version, soft_version, dateTag,
                  sumObj.fRoot, plateid, wellid, seqSetID, trimOpts, flash2Opts, 
                  umiExtractOpts, bowtie2Opts, umiDedupOpts, umiLimit,
                  sumObj.vals['Raw'], sumObj.vals['Trimmed'], 
                  sumObj.vals['Stitched'], sumObj.vals['Aligned'], 
                  overcounts,
                  '\n'.join(['%s,%d' % (tag, tCount) for tag, tCount in 
                             tags.items()])
                  )))
    
# splitSam
#    splits a sam file in to seperate sam files for each tag
#    INPUT:
#        samFile = name of sam file as string
#    OUTPUT:
#        outFiles = dictionary with sam file names as strings
#            for keys and max count for individual UMI assigned
#            to the file as an int (used to determine if deduplication
#            at a hamming distance > 0 is feasible).
#        alnCount = count of total aligned reads as int
def splitSAM(samFile):
    # Storage Dictionary
    linesByTag = {}
    # Loop through file
    with open(samFile) as fh:
        # Header Storage 
        header = ''
        alnCount = 0
        umiCounts = {}
        for line in fh:
            # Store then skip header lines
            if line.startswith('@'):
                header += line 
                continue
            # Extract dsptagID and umi sequence
            qId, _flg, dsptagID = line.rstrip().split('\t')[0:3]
            if dsptagID == '*':
                continue
            # Update summary information
            else:
                alnCount += 1
                umi = qId.split('_')[-1]
                umiCounts.setdefault(dsptagID, {})
                umiCounts[dsptagID][umi] = umiCounts[dsptagID].setdefault(umi, 0) + 1
                # Add line to dictionary
                linesByTag.setdefault(dsptagID, [])
                linesByTag[dsptagID].append(line.rstrip())
    
    # Write outputs
    outFiles = {}
    for dsptagID, tagLines in linesByTag.items():
        outName = samFile.replace('.sam', '_%s.sam' % dsptagID)
        outFiles[outName] = len(umiCounts[dsptagID])
        out = open(outName, 'w') 
        out.write('%s%s' % (header, '\n'.join(tagLines)))
        out.close()
    return outFiles, alnCount

# mergeSAM
#    Merges SAM files keeping header from first file in list
#    INPUT:
#        samList = list of full filepaths to SAM files to merge
#    OUTPUT:
#        outSam = full filepath to desired output location
def mergeSAM(samList, outSam):
    # First grep grabs header from first file
    # Second grep grabs sequence alignments from all files
    os.system('grep ^@ %s > %s' % (samList[0], outSam))
    for fh in samList:
        os.system('grep -v ^@ %s >> %s' % (fh, outSam))
        
# dccToTSV
#    Convert DCC File to TSV File
#    INPUT:
#        dccfilename = name and full path to dcc file as string
#        tsvfilename = name and full path to desired tsv file
#            output
#    OUTPUT:
#        tsv file at specified location with two columns
#            one: analyte/target/tile name as specified in dcc file
#            two: deduplicated count contained in dcc file
def dccToTSV(dccfilename, tsvfilename):
    dccObj = dccfile(dccfilename)
    dccObj.tsvwrite(tsvfilename)

# generateCountTable
#    Generate a table of counts for a list of dcc files
#    INPUT:
#        dccobjs = list of dcc file objects
#    OUTPUT: 
#        df = pandas data frame with first column being Sample_ID
#            and remaining columns corresponding to analyte/target/tile
def generateCountTable(dccobjs):
    rstnames = list(set(
        [ analyte for dcc in dccobjs for analyte in dcc.codesum.keys() ]
        ))
    rstnames.sort()
    df = pd.DataFrame()
    for dcc in dccobjs:
        orderedCounts = [dcc.codesum.get(rst, [0]*99)[dcc.countcol] for rst in rstnames]
        df[dcc.scanattribs['ID']] = orderedCounts
    
    df = df.transpose()
    df.columns = rstnames
    df.insert(0, 'Sample_ID', df.index.values)
    df.index = range(len(df.index))
    
    return df

# pipelineSummaryTable
#    Generate a table summarizing the number of reads
#        present at each step in the analysis
#    INPUT:
#        dccobjs = list of dcc file objects
#    OUTPUT:
#        pandas data frame
def pipelineSummaryTable(dccObjs):
    df = pd.DataFrame()
    for dcc in dccObjs:
        dedupCount = sum(dcc.analytecounts().values())
        dsum = dcc.getSummary()
        df[dsum.fRoot] = dsum.orderedCounts() + [dedupCount]
        
    df = df.transpose()
    df.columns = dsum.steps + ['Unique']
    df.insert(0, 'Sample_ID', df.index.values)
    df.index = range(len(df.index))
    
    return df
    

# summaryInfo
#    Object containing summary information
#    ATTRIBUTES:
#        self.fRoot = file name identifier as string
#        self.steps = list of steps to summarize as strings
#        self.vals = dictionary with keys as steps and values
#            as the summary value
#    METHODS:
#        __init__ = instance initialization
#        update = update step summary value
#        __str__ = convert to single summary string
class summaryInfo(object):
    # __init__
    #    Object Initialization. Defines sample name and 
    #    sets initial summary values to 0. 
    #    INPUT:
    #        self = duh
    #        fRoot = file name, just used for reporting
    #    OUTPUT:
    #        NA
    steps = ['Raw', 'Trimmed', 'Stitched', 'Aligned']
    def __init__(self, fRoot):
        self.fRoot = fRoot
        self.vals = { x: 0  for x in self.steps }
        
    # update
    #    Update specific summary value
    #    INPUT:
    #        self = duh
    #        step = string matching one of self.steps
    #        stepFile = file for summary stat
    #        sampLog = logger to report to
    #        forceVal = known summary value to override file analysis
    #    OUTPUT:
    #        boolean. True if summary value > 0. False if summary value == 0.
    def update(self, step, stepFile, forceVal=0):
        if step not in self.steps:
            raise Exception('%s not recognized as summary step.' % step)
        # Check if a forced value is desired
        if forceVal != 0:
            v = forceVal
        # Process fastq file for sequence count
        elif stepFile.endswith(('fq', 'fastq', 'fastq.gz')):
            v = fastqReads(stepFile)
        # Other methods unsupported
        elif stepFile.endswith('sam'):
            aln, _unaln, _targs = samSummary(stepFile)
            v = aln
        else:
            v = 0
            raise Exception('File type is not supported for summary info. '
                            '%s %s' % (step, stepFile))
        self.vals[step] = v
        # Check step value and return boolean.
        if v == 0:
            logging.debug('No reads found in %s.' % stepFile)
            return False
        else:
            return True
    
    # orderedCounts
    #    Return list of counts ordered by step (1st to last)
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        list of counts at each step in order
    def orderedCounts(self):
        return [ self.vals[stp] for stp in self.steps ]
    
    # __str__
    #    String generation definition.
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        String of summary information.
    def __str__(self):
        return 'File: %s, %s' % (
            self.fRoot, ', '.join([ "%s: %d" % (step, self.vals[step]) 
                                  for step in self.steps ]))
        
# dccfile
#    Object representing a DCC
#    ATTRIBUTES:
#        self.header contains DCC file header information as ordered 
#            dictionary with labels as keys and values as values
#        self.scanattribs contains DCC file scan attribute information as 
#            ordered dictionary with labels as keys and values as values
#        self.ngsattribs contains DCC file ngs processing attribute information  
#            as an ordered dictionary with labels as keys and values as values
#        self.codelabs is a list of column names for the code summary section
#        self.codesum is an ordered dictionary of code summary rows with
#            keys as analyte names and values are the corresponding row as a 
#            list
#        self.namecol
#        self.countcol
#    METHODS:
#        __init__ = instance initialization
#        importDCC = read DCC file and populate attributes
#        listattributes = return string of column labels for code summary 
#            section
#        listanalytes = return string of analytes (line separted)
#        getcodeval = return value for a specific analyte in a specific column
#            in the code summary section
#        analytecounts = return dictionary of analytes and counts
#        addcounts = combine counts from two dccfile objects.
#        __str__ = convert to single DCC/DCC-like string
class dccfile(object):
    # __init__
    #    Object Initialization. Initiates attributes and sets values if
    #        DCC file name is supplied.
    #    INPUT:
    #        self = duh
    #        dccfileName = full path to DCC/DCC-like file if you want to
    #            populate attributes on instance creation
    #    OUTPUT:
    #        NA
    def __init__(self, readfromfile=''):
        # Initialize empty attributes
        self.header = collections.OrderedDict()
        self.scanattribs = collections.OrderedDict()
        self.ngsattribs = collections.OrderedDict()
        self.codelabs = []
        self.codesum = collections.OrderedDict()
        self.countcol = 0
        self.namecol = 0
        # If a file was supplied populate attributes
        if readfromfile.lower().endswith('dcc'):
            self.importDCC(readfromfile)
        elif readfromfile.endswith('sam'):
            self.importSAM(readfromfile)
        else:
            pass
            
    # importDCC
    #    Read an DCC/DCC-like file and populate attributes
    #    INPUT:
    #        self = duh
    #        dccfileName = full path to DCC file
    #    OUTPUT:
    #        None. Populates attributes.
    def importDCC(self, dccfileName):
        # Open file
        with open(dccfileName) as fh:
            # Flag variables to identify location in file
            section = ''
            codeSummaryStarted = False
            # Loop through file
            for line in fh:
                line = line.rstrip()
                # Skip blank lines
                if line == '':
                    continue
                # Check if section end and reset section
                if line.startswith('</'):
                    section = ''
                    continue
                # Check if section start and get section name
                elif line.startswith('<'):
                    section = re.sub('<|>', '', line)
                    continue
                # Split line by ','
                line = line.split(',')
                # Store line as header info if in header section
                if section == 'Header':
                    self.header[line[0]] = line[1]
                # Store line as sample info if in sample attributes section
                elif section == 'Scan_Attributes':
                    self.scanattribs[line[0]] = line[1]
                # Store line as lane info if in late attributes section
                elif section == 'NGS_Processing_Attributes':
                    if len(line)<2:
                        print('Skipping line '+str(line[0])+' to avoid version-specific formatting differences.')
                    else:
                        self.ngsattribs[line[0]] = line[1]
                # Need to do special stuff in code summary section
                elif section == 'Code_Summary':
                    # If this is the first line of the code summary section
                    # store as column labels
                    try:
                        if not codeSummaryStarted:
                            self.codelabs = line
                            codeSummaryStarted = True
                            self.namecol = line.index('Name')
                            self.countcol = line.index('Count')
                            logging.debug("Using old DCC File Format")
                        # If this is not the first line of the code summary section
                        # store the line info in the ordered dictionary
                        else:
                            nm = line[self.namecol]
                            if nm in self.codesum:
                                raise Exception('Repeat analyte in Code_Summary')
                            self.codesum[nm] = line

                    except ValueError: 
                        logging.debug("Using new DCC File Format")
                        codeSummaryStarted = True
                        self.codelabs = ['Name', 'Count']
                        self.namecol = 0
                        self.countcol = 1 
                        nm = line[self.namecol]
                        if nm in self.codesum:
                            raise Exception('Repeat analyte in Code_Summary')
                        self.codesum[nm] = line
                # Do you even know where we are?
                else:
                    logging.warning('Unrecognized Section.')
        
    # listattributes
    #    Return string of column labels for summary section
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        comma separted string of column labels for summary section
    def listattributes(self):
        return ', '.join(self.codelabs)
    
    # listanalytes
    #    Return string of analytes (one analyte per line)
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        line separted string of analytes in code summary section
    def listanalytes(self):
        return '\n'.join(self.codesum.keys())
    
    # getcodeval
    #    Return value for specific analyte and column in code summary section
    #    INPUT:
    #        self = duh
    #        analyte = analyte name as string
    #        attribute = code summary column label as string
    #    OUTPUT:
    #        value of specific analyte and column in code summary section
    #            as string
    def getcodeval(self, analyte, attribute):
        # get desired column as index
        idx = self.codelabs.index(attribute)
        return self.codesum[analyte][idx]
    
    # analytecounts
    #    Return dictionary of analytes and counts
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        Dictionary with keys as analyte names and values as counts
    def analytecounts(self):
        return {analyte:int(self.getcodeval(analyte, 'Count')) for analyte in 
                self.codesum.keys()}
        
    # addcounts
    #    Combine counts from two dccfile objects
    #    INPUT:
    #        self = duh
    #        dccObj = a dccfile object
    #    OUTPUT:
    #        dccfile objects with all attributes of self
    #        except with counts of dccObj added to counts
    #        from self
    def addcounts(self, dccObj):
        countIdx = self.codelabs.index('Count')
        nameIdx = self.codelabs.index('Name')
        altCounts = dccObj.analytecounts()
        for analyte, c in altCounts.items():
                if analyte not in self.codesum.keys():
                    self.codesum[analyte] = [''] * len(self.codelabs)
                    self.codesum[analyte][nameIdx] = analyte
                    self.codesum[analyte][countIdx] = c
                else:
                    self.codesum[analyte][countIdx] = int(self.getcodeval(
                        analyte, 'Count')) + c
        
    # windowswrite
    #    Write official windows formatted DCC
    #    INPUT:
    #        self = duh
    #        outFileName = full path to desired file output as string
    #    OUTPUT:
    #        written windows formatted DCC
    def windowswrite(self, outFileName):
        out = open(outFileName, 'w')
        out.write(('%s' % self).replace('\n', '\r\n'))
        out.close()
    
    # unixwrite
    #    Write official unix formatted DCC
    #    INPUT:
    #        self = duh
    #        outFileName = full path to desired file output as string
    #    OUTPUT:
    #        written unix formatted DCC
    def unixwrite(self, outFileName):
        out = open(outFileName, 'w')
        out.write(('%s' % self))
        out.close()
        
    # tsvwrite
    #    Write tsv file for codesummary section
    #    INPUT:
    #        self = duh
    #        outFileName = full path to desired file output as string
    #    OUTPUT:
    #        written tsv formatted codesummary section    
    def tsvwrite(self, outFileName, sep='\t'):
        out = open(outFileName, 'w')
        out.write('Analyte%sCount\n%s' % (sep,
            '\n'.join(['%s%s%s' % (analyte, sep, self.codesum[analyte][
                self.countcol]) for analyte in self.codesum ])))
        out.close()
        
    # getSummary
    #    Generate Summary Object for dcc object
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        summaryInfo object populated with values from
    #            ngsattribs if available.
    def getSummary(self):
        samp = self.scanattribs['ID']
        sumObj = summaryInfo(samp)
        for stp in sumObj.steps:
            try:
                sumObj.vals[stp] = int(self.ngsattribs[stp])
            except:
                logging.warning('Could not find %s in NGS_Processing_'
                                'Attributes section of dcc object for %s.' 
                                % (stp, samp))
        return sumObj
    
    # countlessdcc
    #    Fill out a DCC object without any counts
    #    INPUT:
    #        self = duh
    #        sampID = name of sample for dcc sample attributes section
    #        sumObj = instance of summaryInfo object 
    #        plateid = unique DSP plate ID
    #        wellid = well name on DSP plate ID for specified isolate 
    #        owner = file owner
    #        comments = any extra information
    #        seqKit = sequencing kit used
    #        dateTag = date of dcc generation
    
    def countlessdcc(self, sumObj, comments='', 
                     plateid='P1012207777777', wellid='A01', owner='Nanostring',
                     seqKit='', 
                     dateTag=datetime.datetime.today().strftime('%Y-%m-%d')):
        # Fill Information
        self.header = collections.OrderedDict([
            ("FileVersion", file_version),
            ("SoftwareVersion", soft_version)
            ])
        
        self.scanattribs = collections.OrderedDict([
            ("ID", sumObj.fRoot), ("Plate_ID", plateid),
            ("Well", wellid), ("Owner", owner),
            ("Comments", comments), ("Date", dateTag)
            ])
        
        self.ngsattribs = collections.OrderedDict(
            [("Sequencing_Kit", seqKit)] +
            [(nm, sumObj.vals[nm]) for nm in sumObj.steps] +
            [('Overcounts', '')]
            )
    
    # __str__
    #    String generation definition.
    #    INPUT:
    #        self = duh
    #    OUTPUT:
    #        Full string of DCC / DCC-like file (unix new line characters)
    def __str__(self):
        return (
            '<Header>\n' +
            '\n'.join(['%s,%s' % (k,v) for k, v in self.header.items()]) +
            '\n</Header>\n\n<Scan_Attributes>\n' +
            '\n'.join(['%s,%s' % (k,v) for k, v in self.scanattribs.items()]) +
            '\n</Scan_Attributes>\n\n<NGS_Processing_Attributes>\n' +
            '\n'.join(['%s,%s' % (k,v) for k, v in self.ngsattribs.items()]) +
            '\n</NGS_Processing_Attributes>\n\n<Code_Summary>\n' +
            ','.join(self.codelabs) + '\n' +
            '\n'.join([ ','.join(v) for v in self.codesum.values() ]) +
            '\n</Code_Summary>\n'
            )
