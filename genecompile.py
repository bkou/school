#!/usr/bin/python
"""genecompiler.py takes many genomes and compiles them into one list of genes while trying to avoid repeats."""

import os
import argparse
from string import strip, split, join

# Blast result Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

def main(startfiles, files, outfile, evalue, minlength, minid):
    """Main is main.
-Checks to make sure output file does not exist to avoid overwriting data.
-Takes the seed lists and adds all genes to output file.
-Adds the genes in the other files iteratively check to see if they already exist in the master list and don't add to master list if it does exist."""
    try:
        #this is here to avoid overwriting someone else's file. I know it's annoying, sorry.
        if(os.path.exists(outfile)):
            raise Exception("Output file already exists. Delete it or choose a different name.")

        #add genes from seed files
        masterf = open(outfile,'a')
        for f in startfiles:
            pnames,pproteins = getProteins(f)
            for name,seq in zip(pnames,pproteins):
                blastprint(name,seq,masterf)
        masterf.close()

        print "Listing removed genes:"
        #add genes from other files
        for f in files:
            pnames,pproteins = getProteins(f)
            for name,seq in zip(pnames,pproteins):
                if os.path.getsize(outfile) == 0:
                    masterf = open(outfile,'a')
                    blastprint(name,seq,masterf)
                    masterf.close()
                else:
                    blastnappend(name, seq, outfile, minlength, minid, evalue)

    except:
        cleantmpfiles()
        raise
    cleantmpfiles()

def cleantmpfiles():
    """Deletes temporary files created by this program"""
    if os.path.isfile("tmpquery.fa"):
        os.remove("tmpquery.fa")

def blastnappend(pname, pseq, outfile, minlength, minid, evalue):
    """blastnappend blast a gene against the master list to see if it's already in there, and if it isn't, adds gene to master list
-creates a temp file with a gene and blasts against master list.
-checks to see if length is greater than minlength, identity is higher than min identity, and if it is, doesn't add to list. Else adds."""
    queryf = open('tmpquery.fa','w')
    queryf.write(">" + pname + "\n")
    queryf.write(pseq)
    queryf.close()

    blastresults = os.popen("blast2 -p blastp -i tmpquery.fa -j " + outfile + " -m 8 -e " + str(evalue))

    masterf = open(outfile,'a')

    for line in blastresults:
        line = split(line[:-1])
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
        if float(line[-10]) > minid and (int(line[-3]) - int(line[-4]) >  minlength or int(line[-5]) - int(line[-6]) > minlength):
            print line[0], "query:", line[-6], "-", line[-5], "subject:", line[-4], "-", line[-3], "identity:", line[-10] + "% evalue:", line[-2]
            #print line
            masterf.close()
            return

    blastprint(pname,pseq,masterf)
    masterf.close()

def blastprint(name,seq,masterf):
    """Writes to file"""
    masterf.write(">" + name + "\n")
    masterf.write(seq)
    masterf.write("\n\n")

def getProteins(subjectf):
    """getProteins returns the name of all the proteins in the comment file, as well as the protein sequences themselves in a dictionary with the name as they key."""
    f = open(subjectf,'rU')
    names = []
    genomes = []
    gen = ""
    first = True
    for line in f:
        if line.startswith(">"):
            if first != True:
                genomes.append(gen)
                gen = ""
            else:
                first = False
            names.append(strip(line[1:-1]))
        else:
            gen += strip(line)
    genomes.append(gen)
    #genomes = {name:genome for name, genome in zip(names,genomes)}
    return names, genomes

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='gene compiler')
    parser.add_argument('-s', action="append",dest="startfiles",# default=[],
                        help='files to be added whole to master list')
    parser.add_argument('-f', action="append",dest="files",# default=[],
                        help='fasta files of sequences')
    parser.add_argument('-o', action="store",dest="outfile",default="master.faa",
                        help='output file, default is master.fa')
    parser.add_argument('-e', action="store",dest="evalue",default=1e-5,
                        help='evalue for blast statement')
    parser.add_argument('-l', action="store",dest="minlength",default=60,
                        help='minimum length for two peptides to be considered the same')
    parser.add_argument('-i', action="store",dest="minid",default=90,
                        help='minimum %identity for two peptides to be considered the same')
    args = parser.parse_args()

    main(args.startfiles,args.files,args.outfile,args.evalue,args.minlength,args.minid)
