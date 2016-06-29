#!/usr/bin/python
# Description:
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -l LISTFILE   Set the listfile
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2009-06-08, updated 2013-03-14, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    outfile = ""
    fileListFile = ""
    fileList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o", "-outfile"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

#     if infile == "":
#         print >> sys.stderr, "infile not set"
#         print >> sys.stderr, usage_short
#         return 1
#     elif not os.path.exists(infile):
#         print >> sys.stderr, "infile %s does not exist"%(infile)
#         return 1

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    print "outfile = %s"%(outfile)
    print "Number of ID =", len(fileList)
    print "They are:"
    for i in xrange(len(fileList)):
        print "%d --> %s" %(i, fileList[i])
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
