#!/usr/bin/python
"""pfinder is a program that takes a protein fragment, blasts it against a proteome of a similar species, take the top hit from the results, and blast that protein back against the proteome of the species that the protein fragment belongs to. The purpose of this is to find a similar protein in the same genome that may be more complete.

This program simply outputs to stdout. To save the results to a file, use the > operator. Check the tcdb wiki article "Redirecting Standard Input/Output" for more info on saving the results to a file."""
import os
import argparse
from string import strip, split

# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

def main(queryf,otherfs,proteinf,evalue):
    """Main is main.
-Calls getProteins to record protein sequences.
-Calls otherlister to get list of similar proteins in other species.
-Calls blastiterator to blast back to original species and print results."""
    try:
        pnames,pproteins = getProteins(proteinf)
        qnames,qproteins = getProteins(queryf)
        qproteins["No Hits"] = ""

        othernames = []
        otherproteins = []
        for otherf in otherfs:
            onames,oproteins = getProteins(otherf)
            oproteins["No Hits"] = ""
            othernames.append(onames)
            otherproteins.append(oproteins)
        for name in pnames:
            printheader(name, pproteins[name])

            for otherf in otherfs:
                printotherheader(otherf)
                oname = otherfinder(name,pproteins[name],otherf,onames,evalue)
                blaster(oname,oproteins[oname],name,pproteins[name],queryf,evalue, qnames, qproteins)

    except:
        cleantmpfiles()
        raise

    cleantmpfiles()

def cleantmpfiles():
    """Deletes temporary files created by this program"""
    if os.path.isfile("tmpquery.fa"):
        os.remove("tmpquery.fa")
    if os.path.isfile("tmpsubject.fa"):
        os.remove("tmpsubject.fa")

def blaster(oname, oseq, pname, pseq, queryf,evalue,qnames,qproteins):
    """blastiterator takes the list of similar proteins and blasts them back to the original genome one by one.
-Calls getProteins to get the sequences of the original species.
-Iterates through the list of similar proteins in other species, gets the top hit, and blasts fragment against it. And prints results"""
   
    hitnames = samefinder(oname, oseq, qnames, queryf, evalue)
    for hitname in hitnames:
        blastnprint(hitname, qproteins[hitname], oname, oseq, pname, pseq, evalue)

def printotherheader(otherf):
    print "****************************************"
    print "HITS USING FILE " + otherf + ":"

def printheader(pname, pseq):
    print "============================================================"
    print "============================================================"
    print "============================================================"
    print "============================================================"
    print "PROTEIN FRAGMENT:"
    print ">" + pname
    print ">Length:", len(pseq)
    print pseq
    print "============================================================"
    print 

def blastnprint(hitname, hitseq, oname, oseq, pname, pseq, evalue):
    """blastother creates a file containing only the protein found similar to a protein on the fragment list. It then blasts that file against the database. The results are printed to stdout.
-makes temporary query and subject files to blast single protein, and blasts them.
-Outputs protein and simliar protein names, lengths, sequences
-outputs blast results"""
    f = open('tmpquery.fa','w')
    f.write(">" + hitname + "\n")
    f.write(hitseq)
    f.close()

    f = open('tmpsubject.fa','w')
    f.write(">" + oname + "\n")
    f.write(oseq)
    f.close()

    print "TOP MATCHED PROTEIN FROM SIMILAR GENOME:"
    if oname == "No Hits":
        print "No Hits"
    else:
        print ">" + oname
        print ">Length:", len(oseq)
        print oseq
        print

        print "TOP MATCHED PROTEIN FROM SAME GENOME:"
        if hitname == "No Hits":
            print "No Hits"
        else:
            print ">" + hitname
            print ">Length:", len(hitseq)
            print hitseq
            print

            print "BLAST RESULTS BETWEEN SIMILAR GENOME (subject) AND SAME GENOME (query):"
            print

            f = os.popen("blast2 -p blastp -i tmpquery.fa -j tmpsubject.fa")#-e " + str(evalue))
            blastprint(f)
            f.close()

def blastprint(f):
    good = False
    for line in f:
        if line.startswith("Lambda") or line.startswith(">"):
            if good == False:
                good = True
            else:
                break
        if good:
            print line,

def samefinder(oname, oseq, qnames, queryf, evalue):
    """samefinder takes the homolog on the other genome and finds the top blast hit against the same genome.
-makes new file and inserts single query """

    if oname == "No Hits":
        return ["No Hits"]
    f = open('tmpquery.fa','w')
    f.write(">" + oname + "\n")
    f.write(oseq)
    f.close()

    f = os.popen("blast2 -p blastp -m 8 -i tmpquery.fa -j " + queryf + " -e " + str(evalue))

    results = []
    for line in f:
        results.append(proteinfromnumber(line, qnames))
    f.close()

    results = list(set(results))

    return results

def otherfinder(pname,pseq,otherf,othernames,evalue):
    """otherlister creates a list of best blast hits in the similar genome corresponding to the list of protein fragments."""
    
    f = open('tmpquery.fa',"w")
    f.write(">" + pname + "\n")
    f.write(pseq + "\n")
    f.close()

    f = os.popen("blast2 -p blastp -m 8 -e " + str(evalue) + " -i tmpquery.fa -j " + otherf)
    
    try:
        return proteinfromnumber(f.next(),othernames)
    except StopIteration:
        return ("No Hits")

def proteinfromnumber(line, nameslist):
"""blast2 sometimes decides to output a number instead of the protein name when outputting tabulated results. I'm not sure why this is, but I've noticed that the number corresponds to the number in sequence of the protein in the file. This method goes back and gets the name from the number"""
    line = split(line[:-1])
    try:
        if len(line[1]) > 2:
            result = (nameslist[int(line[1][:-2])-1])
        else:
            result = (nameslist[int(line[2][:-2])-1])
    except ValueError:
        result = (nameslist[int(line[2][:-2])-1])
    return result


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
    genomes = {name:genome for name, genome in zip(names,genomes)}
    return names, genomes

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='whatever ujj said')
    parser.add_argument('-p', action="store",dest="proteinf",
                        help='fasta file containing one sequence that is being looked for')
    parser.add_argument('-q', action="store",dest="queryf",
                        help='fasta file of sequences, same species as query protein')
    parser.add_argument('-s', action="append",dest="otherfs",
                        help='fasta file of sequences, similar species')
    parser.add_argument('-e', action="store",dest="evalue", default=1e-2,
                        help='fasta file of sequences, similar species')
    args = parser.parse_args()

    main(args.queryf, args.otherfs,args.proteinf, args.evalue)
