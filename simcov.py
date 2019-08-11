#!/usr/bin/env python

import sys
import math
import os.path
import random
import numpy as np
from BIutils.BImisc import Output, genOpen
from BIutils.BIscript import Script

def reverseComplement(seq):
    bases = {'A': 'T',
             'C': 'G',
             'G': 'C',
             'T': 'A',
             'a': 't',
             'c': 'g',
             'g': 'c',
             't': 'a'}
    rc = [ bases[b] for b in seq ]
    rc.reverse()
    return rc

class CovStats():
    maxCov = 0
    numCovBases = 0
    counts = []

    def __init__(self):
        self.counts = [ [0, 0],   # 25th percentile
                        [0, 0],   # 50th percentile
                        [0, 0],    # 75th percentile
                        [1, 0],
                        [5, 0],
                        [10, 0],
                        [20, 0],
                        [30, 0],
                        [50, 0],
                        [100, 0] ]

    def pass1(self, vector, size):
        self.vectorSize = size
        for i in range(size):
            x = vector[i]
            if x > 0:
                self.numCovBases += 1
                if x > self.maxCov:
                    self.maxCov = x
        self.counts[0][0] = int(np.around(self.maxCov * 0.25))
        self.counts[1][0] = int(np.around(self.maxCov * 0.50))
        self.counts[2][0] = int(np.around(self.maxCov * 0.75))

    def pass2(self, vector, size):
        for i in range(size):
            x = vector[i]
            for c in self.counts:
                if x >= c[0]:
                    c[1] += 1

    def report(self, scale, genomesize):
        for i in range(3, len(self.counts)):
            bc = self.counts[i][1] / scale
            pct = bc / genomesize
            sys.stdout.write("Bases covered at {:3}X: {:,}bp ({:.1f}%)\n".format(self.counts[i][0], int(bc), 100.0 * pct))
            
        bc = self.counts[0][1] / scale
        pct = bc / genomesize
        sys.stdout.write("25th percentile - bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[0][0], int(bc), 100.0 * pct))
        bc = self.counts[1][1] / scale
        pct = bc / genomesize
        sys.stdout.write("50th percentile - bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[1][0], int(bc), 100.0 * pct))
        bc = self.counts[2][1] / scale
        pct = bc / genomesize
        sys.stdout.write("75th percentile - bases covered at {}X: {:,}bp ({:.1f}%)\n".format(self.counts[2][0], int(bc), 100.0 * pct))

class SNPsite():
    truepos = 0
    ref = ""
    alt = ""
    all1freq = 0.0
    coverage = 0
    all1count = 0
    all2count = 0

    def __init__(self):
        self.all1freq = random.random()

    def setRefAllele(self, a):
        self.ref = a
        while True:
            b = random.choice('ACGT')
            if b != a:
                self.alt = b
                break

    def genAllele(self):
        if random.random() <= self.all1freq:
            return self.ref
        else:
            return self.alt

    def seen(self):
        self.coverage += 1
        if random.random() <= self.all1freq:
            self.all1count += 1
        else:
            self.all2count += 1

    def sqerror(self):
        x = self.all1freq - 1.0 * self.all1count / self.coverage
        return x * x
        
class SNPstats():
    positions = []
    snps = {}
    counts = []
    ndetected = 0
    
    def __init__(self, size, sitefreq):
        nsites = int(sitefreq * size)
        sys.stderr.write("Simulating {} SNP positions.\n".format(nsites))
        self.positions = np.random.choice(size, nsites, replace=False)
        self.positions.sort()
        for p in self.positions:
            self.snps[p] = SNPsite()
        self.counts = [ [1, 0],
                        [5, 0],
                        [10, 0],
                        [20, 0],
                        [30, 0],
                        [50, 0],
                        [100, 0] ]

    def addread(self, start, end):
        for i in range(start, end):
            if i in self.snps:
                self.snps[i].seen()

    def snpCoverage(self):
        for p in self.positions:
            snp = self.snps[p]
            x = snp.coverage
            if x:
                self.ndetected += 1
                for c in self.counts:
                    if x >= c[0]:
                        c[1] += 1

    def avgError(self):
        errs = []
        for p in self.positions:
            snp = self.snps[p]
            if snp.coverage > 0:
                errs.append(snp.sqerror())
        if errs:
            return np.mean(errs)
        else:
            return None
                        
    def report(self):
        self.snpCoverage()
        sys.stdout.write("\n=== SNPs ===\n")
        nsnps = len(self.positions)
        sys.stdout.write("Simulated: {}\n".format(nsnps))
        sys.stdout.write("Detected: {} ({:.1f}%)\n".format(self.ndetected, 100.0 * self.ndetected / nsnps))
        for c in self.counts:
            sys.stdout.write("  at {}X: {} ({:.1f}%)\n".format(c[0], c[1], 100.0 * c[1] / nsnps))
        sys.stdout.write("Average squared error: {}\n".format(self.avgError()))
        sys.stdout.write("\n")

def usage0():
    sys.stdout.write("""simcov.py - Short reads simulator.

This program can be used in three different modes. The mode can be selected
using one of the command-line options -C, -R, and -S, or by invoking this 
program through a symlink with the appropriate name. The following table 
describes each mode and how to select it:

Opt | Symlink  | Description
----+----------+---------------------------------------
 -C | simcov   | Simulate short-read coverage (default)
 -S | simseq   | Generate a random reference sequence
 -R | simreads | Generate simulated short reads

Note that only the basename of the link is used, so for example both `simseq' 
and  `simseq.py' would invoke simseq mode.
""")
        
def usage():
    sys.stdout.write("""simcov.py - Simulate short read coverage.

Usage: simcov.py [options]

Options:

  -l L | Set read length to L (defalt: {})
  -n N | Set number of reads to N (default: {})
  -g G | Set target genome size to G (default: {})
  -p   | Enable paired-end mode (default: single-end)
  -t T | Set insert size to T in paired-end mode (default: {})
  -f F | Set site frequency to F (default: no site analysis)
  -s S | Set simulation vector size to S (default: {})

Values for -l, -n, -g and -t can be followed by G (for billion) or M (for million).
The value for -f can be expressed as a float or a fraction (e.g. 1/8)

""".format(Simcov.readlen, Simcov.nreads, Simcov.genomeSize, Simcov.insertSize, Simcov.vectorSize))

def usage2():
    sys.stdout.write("""simreads.py - Generate random short reads from a sequence.

Usage: simreads.py [options] filename.fa

Write randomly-generated reads from the reference sequence in `filename.fa' to
file `outfile'.fastq.gz. In paired-end mode, writes to `outfile'_R1.fastq.gz
and `outfile'_R2.fastq.gz.

Options:

  -sn N | Base name for reads (default: {}).
  -nr R | Number of reads (default: {}).
  -rl L | Read length (default: {}).
  -i  I | Average insert size (default: {})
  -is S | Standard deviation of insert size (default: {}) 
  -e  E | Set sequencing error rate to E (default: {})
  -qs S | Average quality at first base (default: {})
  -qe E | Average quality at last base (default: {})
  -qv V | Quality standard deviation at last base (default: {})
  -o  O | Use O as base name for output files (default: {})
  -p    | Enable paired-end mode.
  -s P  | Simulate presence of P SNPs.

The value for -nr can be followed by G (for billion) or M (for million).

""".format(SimReads.seqname, SimReads.nreads, SimReads.readlen, SimReads.insertSize, SimReads.insertStdev, SimReads.errRate, SimReads.qstart, SimReads.qend, SimReads.qvend, SimReads.outfile))

def usage3():
    sys.stdout.write("""simseq.py - Generate random sequence.

Usage: simseq.py [options] filename.fa

Writes a random sequence in FASTA format to file `filename.fa'. Currently, all bases have equal probability.

Options:

  -sn S | Name of sequence (default: {}).
  -l L  | Set sequence length to L (default: {})

The Value for -l can be followed by G (for billion) or M (for million).

""".format(SimReads.seqname, SimReads.seqlen))

class Simcov(Script):
    readlen = 150
    paired = False
    insertSize = 400
    nreads = 10000000
    
    genomeSize = 3100000000
    vectorSize = 1000000
    vector = None
    scale = 1.0

    snpfreq = None
    snpstats = None
    
    def parseArgs(self, args):
        self.standardOpts(args)
        prev = ""
        for a in args:
            if prev == "-l":
                self.readlen = self.toInt(a, units=True)
                prev = ""
            elif prev == "-t":
                self.insertSize = self.toInt(a)
                prev = ""
            elif prev == "-n":
                self.nreads = self.toInt(a, units=True)
                prev = ""
            elif prev == "-g":
                self.genomeSize = self.toInt(a, units=True)
                prev = ""
            elif prev == "-s":
                self.vectorSize = self.toInt(a, units=True)
                prev = ""
            elif prev == "-f":
                self.snpfreq = self.toFloat(a)
                prev = ""
            elif a in ["-l", "-t", "-n", "-g", "-s", "-f"]:
                prev = a
            elif a == "-p":
                self.paired = True
            elif a == "-C":
                pass
            
    def run(self):
        self.vector = np.zeros(self.vectorSize, dtype=int)
        self.scale = 1.0 * self.vectorSize / self.genomeSize
        effReads = int(np.around(self.nreads * self.scale))
        self.writeConfiguration()
        if self.snpfreq:
            self.snpstats = SNPstats(self.vectorSize, self.snpfreq)
        if self.paired:
            self.simulatePaired(effReads)
        else:
            self.simulateUnpaired(effReads)
        CS = CovStats()
        CS.pass1(self.vector, self.vectorSize)
        CS.pass2(self.vector, self.vectorSize)
        self.report(CS)

    def simulateUnpaired(self, nreads):
        low = 1 - self.readlen
        high = self.vectorSize
        for i in range(nreads):
            pos = np.random.randint(low, high)
            start = max(pos, 0)
            end = min(pos + self.readlen, self.vectorSize)
            for p in range(start, end):
                self.vector[p] += 1
            if self.snpstats:
                self.snpstats.addread(start, end)

    def simulatePaired(self, nreads):
        low = 1 - self.insertSize
        high = self.vectorSize
        for i in range(nreads):
            pos = np.random.randint(low, high) # position of insert

            # left read
            start = max(pos - self.readlen, 0)
            end   = max(pos, 0)
            for p in range(start, end):
                self.vector[p] += 1
            if self.snpstats:
                self.snpstats.addread(start, end)

            # right read
            start = min(pos + self.insertSize, self.vectorSize)
            end   = min(start + self.readlen, self.vectorSize)
            for p in range(start, end):
                self.vector[p] += 1
            if self.snpstats:
                self.snpstats.addread(start, end)
            
    def writeConfiguration(self):
        sys.stdout.write("=== Configuration ===\n")
        sys.stdout.write("Genome size: {:,}bp\n".format(self.genomeSize))
        sys.stdout.write("Number of reads: {:,}\n".format(self.nreads))
        sys.stdout.write("Read length: {:,}bp\n".format(self.readlen))
        sys.stdout.write("Mode: {}\n".format("paired" if self.paired else "unpaired"))
        if self.paired:
            nbp = self.nreads * 2 * self.readlen
            sys.stdout.write("Insert size: {:,}bp\n".format(self.insertSize))
        else:
            nbp = self.nreads * self.readlen
        ecov = 1.0 * nbp / self.genomeSize
        sys.stdout.write("Expected average coverage: {:.1f}X\n".format(ecov))
        
    def report(self, CS):
        rcov = 1.0 * np.sum(self.vector) / self.vectorSize
        sys.stdout.write("\n=== Results ===\n")
        sys.stdout.write("Effective average coverage: {:.1f}X\n".format(rcov))
        sys.stdout.write("Max coverage: {:,}\n".format(CS.maxCov))
        pctcov = 1.0 * CS.numCovBases / self.vectorSize
        #sys.stdout.write("Genome covered at 1X: {:,}bp ({:.1f}%)\n".format(int(self.genomeSize * pctcov), pctcov * 100.0))
        CS.report(self.scale, self.genomeSize)

        if self.snpfreq:
            self.snpstats.report()
            #SS = SNPstats(self.vectorSize, self.snpfreq)
            #SS.pass1(self.vector)
            #print SS.counts

### Read simulator

class SeqSource(object):
    readlen = 100
    fpstart = 0
    fpend = 0
    linewidth = 0
    insertSize = 400
    insertStdev = 10
    filename = None
    stream = None

    def __init__(self, filename, readlen=100):
        self.filename = filename
        self.readlen = readlen
        self.getBounds()
    
    def __enter__(self):
        self.stream = open(self.filename, "r")
        return self.stream

    def __exit__(self, type, value, tb):
        self.stream.close()
    
    def getBounds(self):
        fl = os.path.getsize(self.filename)
        self.fpend = fl - self.readlen
        with open(self.filename, "r") as f:
            f.readline()
            self.fpstart = f.tell()
            r = f.readline()
            self.linewidth = len(r)

    def fpToPos(self, fp):
        nrows = (fp - self.fpstart) / self.linewidth # exploit integer division
        return fp - nrows - self.fpstart + 1

    def genPosition(self):
        return random.randint(self.fpstart, self.fpend)
    
    def genInsertPosition(self):
        """Returns a tuple containing the start positions of two mates in a read pair.
The positions are produced by selecting a start position at random, and adding to it
a random insert size, having mean insertSize and standard deviation insertStdev."""
        insize = np.random.normal(self.insertSize, self.insertStdev)
        while True:
            start = random.randint(self.fpstart, self.fpend)
            if start + insize < self.fpend:
                return (start, start + insize - self.readlen)

    def getOneRead(self, pos, probs):
        bases = []
        f = self.stream
        f.seek(pos)
        n = 0
        while True:
            b = f.read(1)
            if b == "\n":
                continue
            if random.random() < probs[n]:
                while True:
                    nb = random.choice('ACGT')
                    if nb != b:
                        b = nb
                        break
            bases.append(b)
            n += 1
            if n == self.readlen:
                break
        return bases
            
    def writeOneRead(self, pos, probs, out):
        f = self.stream
        f.seek(pos)
        n = 0
        while True:
            b = f.read(1)
            if b == "\n":
                continue
            if random.random() < probs[n]:
                while True:
                    nb = random.choice('ACGT')
                    if nb != b:
                        b = nb
                        break
            out.write(b)
            n += 1
            if n == self.readlen:
                break
            
class SimReads(Script, SeqSource):
    readlen = 150
    seqlen = 1000000
    nreads = 1000000
    paired = False
    seqname = "seq"
    errRate = 0.001
    qstart = 40
    qend = 30
    qvstart = 1
    qvend = 10
    nsnps = 0
    outfile = "reads"
    outfile2 = None
    snpfile = None

    qavgs = []
    qstdevs = []
    snps = {}
    
    def parseArgs(self, args):
        self.standardOpts(args)
        prev = ""
        for a in args:
            if prev == "-l":
                self.seqlen = self.toInt(a, units=True)
                prev = ""
            elif prev == "-sn":
                self.seqname = a
                prev = ""
            elif prev == "-nr":
                self.nreads = self.toInt(a, units=True)
                prev = ""
            elif prev == "-rl":
                self.readlen = self.toInt(a)
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-s":
                self.nsnps = self.toInt(a)
                prev = ""
            elif prev == "-qs":
                self.qstart = self.toFloat(a)
                prev = ""
            elif prev == "-qe":
                self.qend = self.toFloat(a)
                prev = ""
            elif prev == "-qv":
                self.qvend = self.toFloat(a)
                prev = ""
            elif prev == "-so":
                self.snpfile = a
                prev = ""
            elif prev == "-e":
                self.errRate = self.toFloat(a)
                prev = ""
            elif prev == "-i":
                self.insertSize = self.toInt(a)
                prev = ""
            elif prev == "-is":
                self.insertStdev = self.toInt(a)
                prev = ""
            elif a in ["-l", "-sn", "-nr", "-rl", "-o", "-so", "-s", "-qs", "-qe", "-qv", "-e", "-i", "-is"]:
                prev = a
            elif a == "-p":
                self.paired = True
            elif a in ["-R", "-S"]:
                pass
            elif self.filename is None:
                self.filename = a

        if self.paired:
            if self.outfile:
                self.outfile2 = self.outfile + "_R2.fastq.gz"
                self.outfile = self.outfile + "_R1.fastq.gz"
            else:
                self.errmsg(self.NOOUTFILE)
        else:
            self.outfile = self.outfile + ".fastq.gz"

    def initQuality(self):
        r = self.qstart - self.qend # range of quality scores
        x = np.power(r, 1.0 / self.readlen)
        q = 1
        for i in range(self.readlen):
            self.qavgs.append(41-q)
            q *= x
        r = self.qvend - self.qvstart
        x = np.power(r, 1.0 / self.readlen)
        q = self.qvstart
        for i in range(self.readlen):
            self.qstdevs.append(q)
            q *= x
        
    def initSNPs(self):
        if self.nsnps == 0:
            return
        sys.stderr.write("Simulating {} SNPs\n".format(self.nsnps))
        n = 0
        with open(self.filename, "r") as f:
            while True:
                sp = random.randint(self.fpstart, self.fpend)
                f.seek(sp)
                b = f.read(1)
                if b in "\r\n":
                    continue
                snp = SNPsite()
                snp.truepos = self.fpToPos(sp)
                snp.setRefAllele(b)
                self.snps[sp] = snp
                n += 1
                if n == self.nsnps:
                    break
        snpsfile = "snps.csv"
        sys.stderr.write("Writing SNPs to file {}\n".format(snpsfile))
        self.writeSNPs(snpsfile)

    def writeSNPs(self, filename):
        with open(filename, "w") as out:
            out.write("#Position\tRef\tAlt\tFreq\n")
            positions = sorted(self.snps.keys())
            for p in positions:
                snp = self.snps[p]
                out.write("{}\t{}\t{}\t{}\n".format(snp.truepos, snp.ref, snp.alt, snp.all1freq))

    def getAllele(self, b, pos):
        if pos in self.snps:
            snp = self.snps[pos]
            return snp.genAllele()
        else:
            return b
                
    def getOneRead(self, f, q, s):
        """Return the sequence of a read starting at position `s' using quality scores `q' from stream `f'."""
        probs = np.power(10, q / -10)
        bases = []
        f.seek(s)
        n = 0
        while True:
            b = f.read(1)
            if b == "\n":
                continue
            if random.random() < probs[n]:
                b = random.choice('ACGT')
            else:
                b = self.getAllele(b, f.tell() - 1)
            bases.append(b)
            n += 1
            if n == self.readlen:
                break
        return bases

    def getSingleRead(self, f, q):
        return self.getOneRead(f, q, random.randint(self.fpstart, self.fpend))

    def getPairedRead(self, f, q1, q2):
        (start1, start2) = self.genInsertPosition()
        r1 = self.getOneRead(f, q1, start1)
        r2 = self.getOneRead(f, q2, start2)
        r2 = reverseComplement(r2)
        return (r1, r2)
    
    def genQuality(self):
        """Returns a list of quality values sampling from the qavgs and qstdevs distribution."""
        return np.clip(np.random.normal(self.qavgs, self.qstdevs), 0, 40)
    
    def simSingleEnd(self):
        readname = self.seqname
        sys.stderr.write("Writing {} single-end reads to {}\n".format(self.nreads, self.outfile))
        with Output(self.outfile) as out:
            with open(self.filename, "r") as f:
                for i in range(1, self.nreads + 1):
                    q = self.genQuality()
                    r = self.getSingleRead(f, q)
                    self.writeRead(out, r, q, readname, i)

    def simPairedEnd(self):
        sys.stderr.write("Writing {} paired-end reads to {} and {}\n".format(self.nreads, self.outfile, self.outfile2))
        with genOpen(self.outfile, "w") as out1:
            with genOpen(self.outfile2, "w") as out2:
                with open(self.filename, "r") as f:
                    for i in range(1, self.nreads + 1):
                        q1 = self.genQuality()
                        q2 = self.genQuality()
                        np.flip(q2, axis=0)
                        (r1, r2) = self.getPairedRead(f, q1, q2)
                        self.writePairedRead(out1, out2, r1, q1, r2, q2, i)
                    
    def writeRead(self, out, r, q, name, i):
        out.write("@{}_{}\n".format(name, i))
        for b in r:
            out.write(b)
        out.write("\n+\n")
        for b in q:
            out.write(chr(int(b) + 33))
        out.write("\n")

    def writePairedRead(self, out1, out2, r1, q1, r2, q2, i):
        readname1 = self.seqname + "_1"
        readname2 = self.seqname + "_2"
        self.writeRead(out1, r1, q1, readname1, i)
        self.writeRead(out2, r2, q2, readname2, i)
        
    def run(self):
        self.getBounds()
        self.initQuality()
        self.initSNPs()
        if self.snpfile:
            self.writeSNPs(self.snpfile)
        if self.paired:
            self.simPairedEnd()
        else:
            self.simSingleEnd()

    def makeRandomSeq(self):
        with open(self.filename, "w") as out:
            out.write(">" + self.seqname + "\n")
            i = 0
            while True:
                out.write(random.choice('ACGT'))
                i += 1
                if i == self.seqlen:
                    break
                if i % 60 == 0:
                    out.write("\n")
            out.write("\n")
        sys.stderr.write("Random sequence of {}bp written to file {}\n".format(self.seqlen, self.filename))

    def run2(self):
        self.makeRandomSeq()

# Generate random reads with quality         
class SimFastq(Script):
    reference = None
    infile1 = None
    infile2 = None
    mode = 'single'
    outfile1 = "reads"
    outfile2 = None
    qualprob = []

    def init(self):
        self.qualprob = [0.0]*75
        idx = 33
        for i in range(0, 42):
            self.qualprob[idx] = math.pow(10, -(i/10.0))
            idx += 1

    def qualsToProbs(self, quals):
        return [ self.qualprob[ord(q)] for q in quals ]
            
    def parseArgs(self, args):
        for a in args:
            if not self.reference:
                self.reference = a
            elif not self.infile1:
                self.infile1 = a
            elif not self.infile2:
                self.infile2 = a
                self.mode = 'paired'

    def run(self):
        src = SeqSource(self.reference, 500) # To make sure we don't go out of bounds (since we don't know readlen yet)
        with src:
            with genOpen(self.infile1, "r") as f:
                while True:
                    readname = f.readline()
                    if not readname:
                        break
                    rl = len(f.readline()) - 1
                    src.readlen = rl
                    f.readline()
                    quals = f.readline().strip()
                    probs = self.qualsToProbs(quals)
                    sys.stdout.write(readname)
                    src.writeOneRead(src.genPosition(), probs, sys.stdout)
                    sys.stdout.write("\n+\n")
                    sys.stdout.write(quals)
                    sys.stdout.write("\n")
    
if __name__ == "__main__":
    prog = os.path.splitext(os.path.split(sys.argv[0])[1])[0]
    args = sys.argv[1:]
    if "-C" in args:
        prog = "simcov"
    elif "-R" in args:
        prog = "simreads"
    elif "-S" in args:
        prog = "simseq"
    elif "-F" in args:
        prog = "simfastq"
    if not args:
        usage0()
        sys.exit(1)
    if prog == "simcov":
        S = Simcov("simcov.py", version="1.0", usage=usage)
        S.parseArgs(args)
        S.run()
    elif prog == "simreads":
        S = SimReads("simreads.py", version="1.0", usage=usage2,
                     errors=[('NOOUTFILE', 'Missing output file name', 'The output file name must be specified in paired-end mode.')])
        S.parseArgs(args)
        if S.filename:
            S.run()
        else:
            S.errmsg(S.NOOUTFILE)
    elif prog == "simseq":
        S = SimReads("simseq.py", version="1.0", usage=usage3,
                     errors=[('NOFAFILE', 'Missing output file name', 'The name of the FASTA output file must be specified.')])
        S.parseArgs(args)
        if S.filename:
            S.run2()
        else:
            S.errmsg(S.NOFAFILE)
    elif prog == "simfastq":
        S = SimFastq("simfastq.py", version="1.0", usage=usage3)
        S.parseArgs(args)
        S.run()
