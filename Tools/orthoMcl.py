#!/usr/bin/env python

'''
Oct 10, 2017: Pasi Korhonen, The University of Melbourne

Simplifies running orthoMCL with a wrapper and pre-checks the known
formatting issues with FASTA headers to avoid failure in later stages of the run.


'''

import os, sys, optparse, getpass
from multiprocessing import Process, Pipe
from Utils import Base, Fasta

#################################################
def options():
    parser = optparse.OptionParser('example: %prog -i residues1.fa,residues2.fa,...,residuesN.fa -l LB1,LB2,...,LBN -p 1,1,...,1 -e 1e-5 -s 0.5')
    parser.add_option('-a', '--ip', dest='ip', help='MySQL server IP address', metavar='IPADDRESS', default='127.0.0.1')
    parser.add_option('-d', '--dir', dest='wd', help='The directory, in which the FASTA files for the analysis are copied to', metavar='DIR', default='TmpOrthoMcl')
    parser.add_option('-i', '--filenames', dest='filenames', help='Proteome, CDS or transcript FASTA files of species (separated by commas)', metavar='FILES', default='')
    parser.add_option('-l', '--labels', dest='labs', help="Respective labels for proteomes following the order of FASTA files (separated by commas)", metavar='LABELS', default='')
    parser.add_option('-p', '--positions', dest='positions', help="Positions of a unique identifier in FASTA header separated by |. Default position is 1 (separated by commas).", metavar='POSITIONS', default='')
    parser.add_option('-T', '--threads', dest='pCnt', help='Number of parallel threads (default is half of the capacity but >= 1)', metavar='THREADS', default='0')
    parser.add_option('-e', '--evalue', dest='evalue', help="E-value for BLAST run. Default is 1e-5. Use always E-value <= 1e-5 and 1e-X format only!", metavar='EVALUE', default='1e-5')
    parser.add_option('-s', '--similarity', dest='sim', help="Required similarity (0 .. 1) in mcl algorithm. Default if 0.5", metavar='SIM', default='0.5')
    parser.add_option('-m', '--minlen', dest='minlen', help="Required minimum lenght of a sequence. Default is 20.", metavar='MINLEN', default='20')
    parser.add_option('-b', '--noblast', dest='skipBlast', action='store_true', help="Skip BLAST (used to rerun mcl using different E-value and similarity settings)", default=False)
    parser.add_option('-n', '--nucl', dest='nucl', action='store_true', help="The residues in sequences represent nucleotides instead of proteins", default=False)
    options, args = parser.parse_args()
    if options.filenames == '' or options.labs == '':
        parser.print_help()
        print '\nE.g.: orthoMcl -i proteome1.fa,proteome2.fa -l Tax,Tvi -p 4,4 -e 1e-5'
        print "Results will be collected to 'Results' directory in groups.txt file."
        print "Note! Labels must be exactly 3 characters long."
        sys.exit(-1)
    return options

#################################################
def checkMySqlVersion(wd, mySqlIpAddress):
    '''
    '''
    base = Base()
    with open("%s/version.sql" %wd, 'w') as handle:
        handle.write('SHOW VARIABLES LIKE "%version%";\n')
    mySqlVer = base.shell("mysql -h %s -P3306 --protocol tcp --user=root --password=password < %s/version.sql" %(mySqlIpAddress, wd), myStdout = True)
    verStr = ""
    found = False
    for line in mySqlVer.stdout:
        items = line.strip().split()
        for i in xrange(len(items)):
            if "inno" in items[i]:
                try:
                    verStr = items[i+1].strip()
                    found = True
                except IndexError:
                    break
        if found == True: break
    try:
        verList = [int(item) for item in verStr.split('.')]
    except ValueError:
        print >> sys.stderr, "### Fatal error: Could not read mysql version (%s). Exiting..." %verStr
        sys.exit(-1)
    if len(verList) != 3:
        print >> sys.stderr, "### Fatal error: Could not read mysql version (%s). Exiting..." %verStr
        sys.exit(-1)
    return verList


#################################################
def checkResidue(fastaFile):
    '''
    '''
    retVal = "nucleotides"
    try:
        limit = 100
        fasta = Fasta(fastaFile)
        for i in xrange(fasta.cnt()):
            if i > limit: break
            seq = fasta.seq(i).upper()
            for item in seq:
                if item not in ['A', 'T', 'C', 'G', 'N']:
                    retVal = "amino acids"
                    break
    except IOError:
        print >> sys.stderr, "### Fatal error: file %s not found. Exiting..." %fastaFile
        sys.exit(-1)
    return retVal


#################################################
def checkUniqueIds(fastaFile):
    '''
    '''
    fasta = Fasta(fastaFile)
    if fasta.cnt() != len(set(fasta.headers)):
        print >> sys.stderr, "### Fatal error: FASTA sequence identifiers are not unique in %s. Exiting..." %fastaFile
        print >> sys.stderr, "### Probably position for this file is given wrong..."
        sys.exit(-1)


#################################################
def checkIdLen(fastaFile):
    '''
    '''
    fasta = Fasta(fastaFile)
    for i in xrange(fasta.cnt()):
        seqId = fasta.headers[i].split()[0]
        if len(seqId) > 56:
            print >> sys.stderr, "### Fatal error: FASTA sequence identifier %s is too long in %s. Exiting..." %(seqId, fastaFile)
            sys.exit(-1)


#################################################
def createOrthoMclConfigFile(wd, userName, eValue, similarity, mySqlIpAddress):
    '''
    '''
    eValue = eValue.split('e')[1]
    similarity = int(float(similarity) * 100.0)
    with open("%s/orthomcl.config" %wd, 'w') as handle:
        handle.write("# this config assumes a mysql database named 'orthomcl'.  adjust according\n")
        handle.write("# to your situation.\n")
        handle.write("dbVendor=mysql\n")
        #handle.write("dbConnectString=dbi:mysql:database=ortho%s;host=%s;port=3306\n" %(userName, os.environ["MYSQLHOST"]))
        handle.write("dbConnectString=dbi:mysql:database=ortho%s;host=%s;port=3306\n" %(userName, mySqlIpAddress))
        handle.write("dbLogin=ortho%s\n" %userName)
        handle.write("dbPassword=password\n")
        handle.write("similarSequencesTable=SimilarSequences\n")
        handle.write("orthologTable=Ortholog\n")
        handle.write("inParalogTable=InParalog\n")
        handle.write("coOrthologTable=CoOrtholog\n")
        handle.write("interTaxonMatchView=InterTaxonMatch\n")
        handle.write("percentMatchCutoff=%d\n" %similarity)
        handle.write("evalueExponentCutoff=%s\n" %eValue)
        handle.write("oracleIndexTblSpc=NONE\n")


#################################################
def createMySqlScripts(wd, userName, ver):
    '''
    '''
    with open("%s/createDb.sql" %wd, 'w') as handle:
        if ver[0] > 5 or (ver[0] == 5 and ver[1] > 7) or (ver[0] == 5 and ver[1] == 7 and ver[2] > 5):
            handle.write("CREATE USER IF NOT EXISTS 'ortho%s'@'%%' IDENTIFIED BY 'password';\n" %userName)
        handle.write("CREATE DATABASE ortho%s;\n" %userName)
        #handle.write("GRANT SELECT,INSERT,UPDATE,DELETE,CREATE VIEW,CREATE,INDEX,DROP on ortho%s.* TO 'ortho%s'@'%%';\n" %(userName, userName))
        handle.write("GRANT ALL PRIVILEGES ON ortho%s.* TO 'ortho%s'@'%%';\n" %(userName, userName))
        handle.close()
        handle = open("%s/dropDb.sql" %wd, 'w')
        handle.write("drop database if exists ortho%s;\n" %userName)


#################################################
def callShell(base, cmdStr, dummy = None):
    '''
    '''
    base.shell(cmdStr)


#################################################
def main():
    '''
    '''
    opts = options() # files contains exactly two PE files

    pCnt = int(opts.pCnt)
    if pCnt == 0:
        pCnt = int(float(multiprocessing.cpu_count()) / 2.0 + 0.5)
    eValue = opts.evalue
    similarity = opts.sim
    minlen = opts.minlen
    files = ["%s/%s" %(opts.wd.strip('"').strip("'").rstrip('/'), myFile) for myFile in opts.filenames.strip().split(',')]
    labels = opts.labs.strip().split(',')
    if len(labels) != len(set(labels)):
        print >> sys.stderr, "### Fatal error: duplicate labels found. Exiting..."
        sys.exit(-1)
    if len(files) != len(set(files)):
        print >> sys.stderr, "### Fatal error: duplicate fasta file names found. Exiting..."
        sys.exit(-1)
    positions = None
    if opts.positions != "":
        positions = opts.positions.strip().split(',')
    if positions == None:
        positions = []
        for i in xrange(len(files)):
            positions.append("1")
    if len(files) != len(labels):
        print >> sys.stderr, "### Fatal error: number of files does not match with the number of labels. Exiting..."
        sys.exit(-1)
    if len(positions) != len(labels):
        print >> sys.stderr, "### Fatal error: number of labels does not match with the number of positions of the ids. Exiting..."
        sys.exit(-1)
    for lab in labels:
        if len(lab) != 3:
            print >> sys.stderr, "### Fatal error: labels have to be exactly three characters long. Exiting..."
            sys.exit(-1)

    base = Base()
    wd = "Results"
    wdFasta = "%s/Fasta" %wd
    base.shell("rm -rf Results")
    base.createDir(wd)
    logHandle = open("%s/log.txt" %wd, 'w')
    base.setLogHandle(logHandle)
    base.createDir(wdFasta)
    #userName = getpass.getuser()
    #createOrthoMclConfigFile(wd, userName, eValue, similarity)
    #createMySqlScripts(wd, userName)
    verList = checkMySqlVersion(wd, opts.ip)
    createOrthoMclConfigFile(wd, "root", eValue, similarity, opts.ip)
    createMySqlScripts(wd, "root", verList)

    requiredMolType = "amino acids"
    if opts.nucl == True:
        requiredMolType = "nucleotides"
    for myFile in files:
        molType = checkResidue(myFile)
        if requiredMolType != molType:
            print >> sys.stderr, "### Fatal error: files have to all be %s. Exiting..." %requiredMolType
            print >> sys.stderr, "### File %s failed and was %s." %(myFile, molType)
            sys.exit(-1)

    base.shell("rm -f %s/*.fasta" %wd)
    base.shell("rm -f %s/*.fasta" %wdFasta)

    for i in xrange(len(files)):
        myLab, myFile, myPos = labels[i], files[i], positions[i]
        if myFile == "%s.fasta" %myLab:
            print >> sys.stderr, "### Fatal error: orthoMCL produces same filenames that you already have. Please rename your fasta files e.g. to .fa instead of .fasta. Exiting..."
            sys.exit(-1)
        base.shell("orthomclAdjustFasta %s %s %s" %(myLab, myFile, myPos))
        checkUniqueIds("%s.fasta" %myLab)
        checkIdLen("%s.fasta" %myLab)
        base.shell("mv -f %s.fasta %s" %(myLab, wdFasta))

    if opts.skipBlast == False:
        base.shell("orthomclFilterFasta %s %s 20" %(wdFasta, minlen))
        base.shell("mv -f poorProteins.* %s" %wd)

    # Blast all against all
    if opts.skipBlast == False:
        if opts.nucl == False:
            base.shell("makeblastdb -in goodProteins.fasta -dbtype prot")
        else:
            base.shell("makeblastdb -in goodProteins.fasta -dbtype nucl")
        base.shell("cp goodProteins.fasta %s/" %wd)
    blastEvalue = eValue
    if float(blastEvalue) < 1e-5: blastEvalue = "1e-5"
    if opts.skipBlast == False:
        if opts.nucl == False:
            base.shell("blastp -db goodProteins.fasta -query goodProteins.fasta -outfmt 6 -evalue %s -num_threads %d > %s/goodProteins.blast" %(blastEvalue, pCnt, wd))
        else:
            base.shell("blastn -db goodProteins.fasta -query goodProteins.fasta -outfmt 6 -evalue %s -num_threads %d > %s/goodProteins.blast" %(blastEvalue, pCnt, wd))
    base.shell("""awk '{if ($11<=%s) print $0}' %s/goodProteins.blast | grep -v "^#" > %s/filtered.blast""" %(eValue, wd, wd))
    base.shell("mv -f goodProteins.* %s/" %wd)

    base.shell("orthomclBlastParser %s/filtered.blast %s > %s/similarSequences.txt" %(wd, wdFasta, wd))
    # Prepare database
    base.shell("mysql -h %s -P 3306 --protocol tcp --user=root --password=password < %s/dropDb.sql" %(opts.ip, wd))
    base.shell("mysql -h %s -P 3306 --protocol tcp --user=root --password=password < %s/createDb.sql" %(opts.ip, wd))
    base.shell("orthomclInstallSchema %s/orthomcl.config" %wd)
    base.shell("orthomclLoadBlast %s/orthomcl.config %s/similarSequences.txt" %(wd, wd))
    # Identify potential orthologs
    base.shell("orthomclPairs %s/orthomcl.config %s/orthomclPairs.log cleanup=no" %(wd, wd))
    base.shell("rm -rf pairs")
    base.shell("rm -rf %s/pairs" %wd)
    base.shell("orthomclDumpPairsFiles %s/orthomcl.config" %wd)
    base.shell("mv -f pairs %s" %wd)
    # Group the orthologs
    base.shell("mcl mclInput --abc -I 2.0 -o mclOutput")
    base.shell("orthomclMclToGroups OWN_ 1 < mclOutput > %s/groups.txt" %wd)
    base.shell("mv -f mclInput %s" %wd)
    base.shell("mv -f mclOutput %s" %wd) 
    logHandle.close()
       

#################################################
if __name__ == "__main__":
    main()
