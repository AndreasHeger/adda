import sys, os, re, time

import alignlib
import ProfileLibrary

from AddaModule import AddaModule
import AddaIO
import IndexedFasta
import SegmentedFile

class FastaRecord:
    def __init__(self, title, sequence ):
        self.title = title
        self.pid = re.sub( "\s.*", "", title )
        self.sequence = sequence

def _iterate( infile ):

    h = infile.readline()[:-1]

    if not h: raise StopIteration

    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h: raise StopIteration
        continue

    h = h[1:]
    seq = []

    for line in infile:
        if line[0] == "#": continue
        line = line[:-1] # remove newline
        if not line: continue
        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            h = line[1:]
            seq = []
            continue

        seq.append(line)
    yield FastaRecord(h,''.join(seq))

class FastaIterator:

    def __init__(self, f, *args, **kwargs):
        self.mIterator = _iterate(f)
    def __iter__(self):
        return self
    def next(self):
        return self.mIterator.next()

class AddaSequences( AddaModule ):
    """output sequence information."""
    
    mName = "Sequences"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameNids = self.mConfig.get( "files", "output_nids", "adda.nids" )  
        self.mFilenameInputFasta = self.mConfig.get( "files", "input_fasta" )
        self.mFilenameOutputFasta = self.mConfig.get( "files", "output_fasta", "adda" )
        self.mMaxSequenceLength = self.mConfig.get( "segments", "max_sequence_length", 10000 )

        self.mFilenames = (self.mFilenameNids, )

    def startUp(self):
        if self.isComplete(): return
    
    def applyMethod(self ):

        self.mInput = 0
        self.mOutput = 0
        self.mRemoved = 0

        # use existing fasta file
        iterator = FastaIterator( open( self.mFilenameInputFasta, "r" ) )
        fasta = IndexedFasta.IndexedFasta( self.mFilenameOutputFasta, "w" )

        outfile = self.openOutputStream(self.mFilenameNids)
        outfile.write( "nid\tpid\thid\tlength\tsequence\n" )

        nid = 1
        hids = set()
        
        for seq in iterator:
            
            self.mInput += 1
            if len( seq.sequence ) > self.mMaxSequenceLength:
                self.mRemoved += 1
                continue

            hid = self.getHID( seq.sequence )
            if hid in hids:
                self.mDuplicates += 1
                continue
            
            hids.add(hid)
            outfile.write( "%s\t%s\t%s\t%i\t%s\n" % (nid, seq.pid, hid, len(seq.sequence), seq.sequence) )
            fasta.addSequence( nid, seq.sequence )
            nid += 1
            self.mOutput += 1

        fasta.close()
        outfile.close()

    def finish(self):
        
        self.info( "sequences: %i input, %i output, %i removed" %\
                   (self.mInput, self.mOutput, self.mRemoved ) )
        
        AddaModule.finish( self )
        
