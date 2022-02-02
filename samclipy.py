#!/usr/bin/env python

import sys
import re
import argparse

def mainArgs():
    parser = argparse.ArgumentParser(description='Filter SAM file for clipped alignments.',prog='samclipy')
    # Input options
    parser.add_argument('samfile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    # parser.add_argument('--fai',type=str,required=True,help='Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`')
    parser.add_argument('--minClip',type=float,default=1,help='Required clip length (left + right). If >= 1 total number; if < 1 relative to read length.')
    parser.add_argument('--invert',default=False,action='store_true',help='Output only soft-clipped SAM records and ignore the good ones.')
    args = parser.parse_args()
    return args

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    sys.exit(error_type)

def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')

def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (len,operator)
    """
    CIGARlist = list()
    for x in re.findall('[0-9]*[A-Z|=]',SAM_CIGAR):
        CIGARlist.append((int(re.findall('[0-9]*', x)[0]),re.findall('[A-Z]|=', x)[0]))
    #174M76S --> [(174,M),(76,S)]
    #96S154M --> [(96,S),(154,M)]
    return CIGARlist

def checkClips(SAM_CIGAR):
    """
    Get lengths of soft-clipped blocks from either end of an alignment given a CIGAR string.
    """
    leftClipLen = None
    rightClipLen = None
    CIGARlist = splitCIGAR(SAM_CIGAR)
    # Check if first segment is soft-clipped
    if CIGARlist[0][1] == "S" :
        leftClipLen = int(CIGARlist[0][0])
    # Check if last segment is soft-clipped
    if CIGARlist[-1][1] == "S" :
        rightClipLen = int(CIGARlist[-1][0])
    return (leftClipLen,rightClipLen)

def lenCIGAR(SAM_CIGAR):
    """
    Calculate length of alignment in reference sequence as sum of 
    match, read-deletion, splice, mismatch, and read-match block values.
    Ignore read-insertions, padding, hard and soft clip blocks.
    """
    alnLen = 0
    CIGARlist = splitCIGAR(SAM_CIGAR)
    for x in CIGARlist: # i.e. = [(174,M),(76,S)]
        if x[1] in set(['D','M','N','X','=']):
            alnLen += x[0]
    #Ignore operators in set('P','H','S','I')
    return alnLen

def fai_to_dict(fai_file):
    ref_len_dict = {}
    with open(fai_file, 'r') as fai:
        for line in fai:
            contig = line.split('\t')[0]
            len = line.split('\t')[1]
            ref_len_dict[contig] = len
    return ref_len_dict

def main():
    # Get cmd line args
    args = mainArgs()
    # SAM line index keys
    SAM_RNAME = 2
    SAM_POS   = 3
    SAM_CIGAR = 5
    SAM_TLEN  = 8
    SAM_SEQ   = 9
    # log counters
    log_total = 0
    log_removed = 0
    log_passed = 0
    log_kept = 0
    log_header = 0
    # check input
    if args.minClip <= 0:
        error(f"minClip must be greater than 0")

    # len_dict = fai_to_dict(fai)
    # log(f"Found {len(len_dict)} sequences in {fai}.")

    # iterate over sam file
    for line in args.samfile:
        line = line.strip() 
        if line.startswith('@'):
            log_header += 1
            print(line)
        else:
            log_total += 1
            sam_rec = line.split('\t')
            if "S" in sam_rec[SAM_CIGAR] and not "H" in sam_rec[SAM_CIGAR]:
                leftClipLen,rightClipLen = checkClips(sam_rec[SAM_CIGAR])
                clipSum = sum(filter(None, [leftClipLen, rightClipLen]))
                if args.minClip >= 1: # check absolute number of soft-clipped positions 
                    if clipSum < args.minClip:
                        # S but not to much in sum
                        log_passed += 1
                        if not args.invert:
                            print(line)
                    else:
                        # to much soft-clipping
                        log_removed += 1
                        if args.invert:
                            print(line)
                else: # check soft-clipped positions relative to alignment length
                    readLen = len(sam_rec[SAM_SEQ])
                    relativeClipLen = clipSum/readLen
                    if relativeClipLen < args.minClip:
                        log_passed += 1
                        if not args.invert:
                            print(line)
                    else:
                        log_removed += 1
                        if args.invert:
                            print(line)
            else:
                log_kept += 1
                if not args.invert:
                    print(line)
    log(f"Total SAM records {log_total}, clipped and removed {log_removed}, clipped and passed {log_passed}, not clipped {log_kept}")
    log(f"Header contained {log_header} lines")
    if args.invert:
        log(f"Printed {log_removed} clipped SAM records")
    else:
        log(f"Printed {log_passed + log_kept} SAM records")
    log('Done.')

if __name__ == '__main__':
    main()