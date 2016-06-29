#!/usr/bin/python

import os
import sys
import argparse

if __name__ == '__main__':

#metavar='' is the text shown after then option argument
    parser = argparse.ArgumentParser(
            description='Template python code with argument parsing',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''\
Created 2011-09-12, updated 2011-09-12, Nanjiang Shu

Examples:
    argParser.py -h
''')
    parser.add_argument('filenamelist', metavar='FILE', nargs='*', 
            help='supply one or more filenames')
    parser.add_argument('-l', metavar='FILE', dest='listfile', 
            help='provide a file with a list of filenames')
    parser.add_argument('-o' , metavar='OUTFILE', dest='outfile', 
            help='output the result to outfile')
    parser.add_argument('-outpath', metavar='DIR', dest='outpath', 
            help='output the result to the dir outpath')
    parser.add_argument('-v', dest='verbose', nargs='?', type=int, default=0, const=1, 
            help='show verbose information, (default: 0)')

    args = parser.parse_args()
    
    listfile= args.listfile
    outpath=args.outpath
    verbose=args.verbose
    filenamelist=args.filenamelist

    if listfile != None:
        fpin=open(listfile,"r")
        filenamelist += fpin.read().split()
        fpin.close(); 

    for file in filenamelist:
        print file

