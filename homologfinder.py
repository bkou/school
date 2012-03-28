#!/usr/bin/python

"""homologfinder.py is a program that takes in pre-computed BLAST results and attempts to compile multiple BLAST hits into one big hit as if they were multiple exons for the same mRNA transcript."""

import sys #for taking in argv
import os #for opening, reading, writing files
import argparse
from string import strip, split, join

# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

def main(queryf, subjectf,blastf):

    #get names and sequences of genomes from contigs
    gnames,genomes = getNames(subjectf)
    #get names and sequences of proteins as dictionary data structure
    proteins = getProteins(queryf,gnames)
    #read blast results from
    loadHits(blastf,queryf,proteins,gnames)

    #find exon combos and return them
    findcombos(proteins, gnames, genomes,queryf)

class Seq():
    """Array structure for each protein to store blast hits and percent coverage"""
    def __init__(self,name,length,gnames):
        #name of sequence in fasta file comment line
        self.name = name
        #length of fasta sequence
        self.length = length
        #creating empty array for each gene to be filled with hits between protein and gene
        self.posframe = {name:[] for name in gnames}
        self.negframe = {name:[] for name in gnames}
        #percent coverage of protein by each gene
        self.posscore = {genename:0 for genename in gnames}
        self.negscore = {genename:0 for genename in gnames}
    def addpos(self,line,gname):
        """adds hit to list of positively framed hits"""
        self.posframe[gname].append(line)
    def addneg(self,line,gname):
        """adds hit to list of negatively framed hits"""
        self.negframe[gname].append(line)
    def getpos(self,gname):
        """gets list of positively framed hits"""
        return self.posframe[gname]
    def getneg(self,gname):
        """gets list of negatively framed hits"""
        return self.negframe[gname]
    def getlen(self):
        """returns length of protein"""
        return self.length
    def getposscore(self,gname):
        """returns score in positive direction."""
        return self.posscore[gname]
    def getnegscore(self,gname):
        """returns score in negative direction."""
        return self.negscore[gname]

def findcombos(proteins, gnames,genomes,queryf):
"""This method finds sequences of exons in sequential order on the protein that are also in sequential order on the gene, and then prints the results using the method blastnprint. This was written before I had learned the dynamic programming algorithm for longest exon sequence length, so it uses the greedy algorithm that continuously takes the next non-overlapping exon."""
    for key in sorted(proteins):
        length = proteins[key].getlen()
        for name in gnames:
            hits = []
            for line in proteins[key].getneg(name):
                hits.append(line)
            if len(hits) != 0:
                hits = sorted(hits, key=lambda line: line[8],reverse=True)
                combolist = []
                combo = [hits[0]]
                for line1,line2 in zip(hits[:-1], hits[1:]):
                    if line1[7] <= line2[6] and int(line1[8]) - int(line2[9]) < 1500:
                        combo.append(line2)
                    else:
                        if len(combo) != 0 and scorer(combo, length) == 2:
                            combolist.append(combo)
                        combo = [line2]
                if len(combo) != 0 and scorer(combo,length) == 2:
                    combolist.append(combo)
                if len(combolist) > 0:
                    print "Protein: " + key
                    print ">length:", length
                    print "Contig: " + name
                    print ">length:", len(genomes[name])
                    for combo in combolist:
                        blastnprint(combo,genomes[name],queryf,key)

        for name in gnames:
            hits = []
            for line in proteins[key].getpos(name):
                hits.append(line)
            if len(hits) != 0:
                hits = sorted(hits, key=lambda line: line[8])
                combolist = []
                combo = [hits[0]]
                for line1,line2 in zip(hits[:-1], hits[1:]):
                    if line1[7] <= line2[6] and int(line2[9]) - int(line1[8]) < 1500:
                        combo.append(line2)
                    else:
                        if len(combo) != 0 and scorer(combo,length) == 2:
                            combolist.append(combo)
                        combo = [line2]
                if len(combo) != 0 and scorer(combo,length) == 2:
                    combolist.append(combo)
                if len(combolist) > 0:
                    print "Protein: " + key
                    print ">length:", length
                    print "Contig: " + name
                    print ">length:", len(genomes[name])
                    for combo in combolist:
                        blastnprint(combo,genomes[name],queryf,key)

def blastnprint(combo,genome,queryf,key):
    """This takes the gene and protein with 50%+ identity and prints information on it."""
    gname = ""
    pname = ""
    for line in combo:
        print "exon  :" + '\t' + '\t'.join([line[0],line[1],line[6],line[7],line[8],line[9]])
    gname = line[1] 
    pname = line[0]
    seq = ""
    for line in sorted(combo, key=lambda c: c[8]):
        if int(line[8]) < int(line[9]):
            start, end = line[8], line[9]
            start, end = int(start)-1, int(end)
            seq += (genome[start:end])
            #print seq
        else:
            start, end = line[9], line[8]
            start, end = int(start)-1, int(end)
            seq += (genome[start:end])[::-1]
            tempseq = ""
            for i in range(len(seq)):
                if seq[i] == 'A':
                    tempseq += 'T'
                if seq[i] == 'T':
                    tempseq += 'A'
                if seq[i] == 'G':
                    tempseq += 'C'
                if seq[i] == 'C':
                    tempseq += 'G'
            seq = tempseq
            #print seq
    f = open('tmpsubject.fa','w')
    f.write(seq)
    f.close()
    #f1 = open('tmpsubject.fa','w')
    f=open(queryf,'rU')
    record = False
    aaseq = ""
    for line in f:
        if line.startswith('>'):
            if strip(line[1:-1]) == key:
                record = True
            else:
                record = False
        elif record == True:
            aaseq += strip(line)
    f.close()
    f = open('tmpquery.fa','w')
    f.write(">" + pname + "\n")
    f.write(aaseq)
    f.close()
    f = os.popen("blast2 -p tblastn -i tmpquery.fa -j tmpsubject.fa -m 8 -e 1e-5")
    for line in f:
        line = split(line)
        line[0],line[1] = pname,gname
        print "RESULT:" + '\t' + '\t'.join([line[0],line[1],line[6],line[7],line[8],line[9]]) + '\n'

    f = open('tmpsubject.fa','w')
    f.write(">" + gname + "\n")
    f.write(translate(seq))
    f.close()
    good = False
    f = os.popen("blast2 -p blastp -i tmpquery.fa -j tmpsubject.fa -e 1e-5")
    for line in f:
        if line.startswith(" Score"):
            good = True
        if line.startswith("Lambda"):
            break
        if good:
            print line

def translate(nseq):
    """Takes nucletide sequence and returns translated amino acid sequence"""
    aseq = ""
    for n1,n2,n3 in zip(nseq[::3],nseq[1::3],nseq[2::3]):
        codon = n1+n2+n3
        if   codon in ["ATT","ATC","ATA"]:
            aseq += "I"
        elif codon in ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"]:
            aseq += "L"
        elif codon in ["GTT", "GTC", "GTA", "GTG"]:
            aseq += "V"
        elif codon in ["TTT", "TTC"]:
            aseq += "F"
        elif codon in ["ATG"]:
            aseq += "M"
        elif codon in ["TGT", "TGC"]:
            aseq += "C"
        elif codon in ["GCT", "GCC", "GCA", "GCG"]:
            aseq += "A"
        elif codon in ["GGT", "GGC", "GGA", "GGG"]:
            aseq += "G"
        elif codon in ["CCT", "CCC", "CCA", "CCG"]:
            aseq += "P"
        elif codon in ["ACT", "ACC", "ACA", "ACG"]:
            aseq += "T"
        elif codon in ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]:
            aseq += "S"
        elif codon in ["TAT", "TAC"]:
            aseq += "Y"
        elif codon in ["TGG"]:
            aseq += "W"
        elif codon in ["CAA", "CAG"]:
            aseq += "Q"
        elif codon in ["AAT", "AAC"]:
            aseq += "N"
        elif codon in ["CAT", "CAC"]:
            aseq += "H"
        elif codon in ["GAA", "GAG"]:
            aseq += "E"
        elif codon in ["GAT", "GAC"]:
            aseq += "D"
        elif codon in ["AAA", "AAG"]:
            aseq += "K"
        elif codon in ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]:
            aseq += "R"
        elif codon in ["TAA", "TAG", "TGA"]:
            aseq += "?"
    return aseq
 
def score(percent):
    if percent == 0:
        return 0
    if percent > 0.5:
        return 2
    else:
        return 1

def scorer(combo,length):
    unhit = [1, length]
    for line in combo:
        insert(unhit,int(line[6]),int(line[7]))
    num = length
    for start,end in zip(unhit[::2],unhit[1::2]):
        num -= end - start + 1
    return score(float(num)/float(length))

def insert(unhit, start,end):
    """insert takes a list representing unhit ranges, and a start/stop for a hit area.
    It then detects what hit areas have been hit and removes them from the unhit list."""
    end += 1
    spos = 0 #len(unhit)
    epos = 0 #len(unhit)
    for i in range(len(unhit)):
        if start <= unhit[i]:
            break
        spos += 1
    for i in range(len(unhit)):
        if end <= unhit[i]:
            break
        epos += 1
    start -= 1
    if spos%2 == 0 and epos%2 == 0:
        unhit[spos:epos] = []
    elif spos%2 == 0 and epos%2 == 1:
        unhit[spos:epos] = [end]
    elif spos%2 == 1 and epos%2 == 0:
        unhit[spos:epos] = [start]
    elif spos%2 == 1 and epos%2 == 1:
        unhit[spos:epos] = [start,end]

def loadHits(blastf,queryf,proteins,gnames):
    """loadHits will distribute hits with appropriate score to their respective instances
in the protein dictionary. It will also detect frame.

Note that the reverse frame is not accurate.
It will not give the correct frame, but it will be able to distinguish between frames."""
    #print "Getting blast results from NCBI blast 2.0 process"
    #f = os.popen("blast2 -p tblastn -m 8 -i " + queryf + " -j " + subjectf)
    f = open(blastf,'rU')
    for line in f:
        line = split(line[:-1])
        frame = int((line)[9])-int((line)[8])
        if frame < 0:
            frame = -1 - int((line)[8])%3
        else:
            frame = int((line)[8])%3 + 1
        
        line.append(str(frame))

        line[1] = gnames[ int((line[1])[:-2])-1 ]

        if frame > 0:
            proteins[line[0]].addpos(line,line[1])
        else:
            proteins[line[0]].addneg(line,line[1])

def getProteins(queryf,gnames):
    """getProteins initializes the dictionary that will hold information for each protein.
Initially it will contain just the name and length."""
    #print "Reading nucleotide sequence"
    f = open(queryf,'rU')
    proteins = {}
    first = True
    aaseq = ""
    for line in f:
        if line.startswith('>'):
            if first is True:
                first = False
            else:
                proteins[name] = (Seq(name,len(aaseq),gnames))
                aaseq=""
            name = line[1:-1]
        aaseq += strip(line)
    proteins[name] = (Seq(name,len(aaseq),gnames))
    return proteins

def getNames(subjectf):
    """getGenome returns an uninterrupted string representing the subject genome"""
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
    genomes = {name:genome for name, genome in zip(names,genomes)}
    return names, genomes

if __name__ == "__main__":
    #argument parser
    parser = argparse.ArgumentParser(description='whatever ujj said')
    parser.add_argument('-q', action="store",dest="queryf",
                        help='fasta file containing PROTEIN sequences')
    parser.add_argument('-s', action="store",dest="subjectf",
                        help='fasta file containing GENE sequences')
    parser.add_argument('-b', action="store",dest="blastf",
                        help='tabulated blast output')
    args = parser.parse_args()

    #checking to see if the necessary arguments were passed through command line
    argsok = True
    if not args.queryf:
        print "Query file not specified"
        argsok = False
    if not args.subjectf:
        print "Subject file not specified"
        argsok = False
    if not args.queryf:
        print "BLAST result file not specified"
        argsok = False
    if not argsok:
        print "Usage: ./homologfinder.py -q <queryfile> -s <subjectfile> -b <blastresultfile>"
        print "queryfile      : PROTEIN fasta file"
        print "subjectfile    : GENE fasta file"
        print "blastresultfile: file with recorded output from tblastn between protein query and genome subect file"
        print "Exiting program"
        exit()
        
    main(args.queryf, args.subjectf, args.blastf)
