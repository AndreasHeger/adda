# extract data from nrdb40.dump.gz and nrdb.dump.gz
# for nrdb40.fasta

import gzip, os, re, subprocess

TAG_NRDB40="INSERT INTO `nrdb40` VALUES ("
TAG_NRDB="INSERT INTO `nrdb` VALUES ("
TAG_SCOP="INSERT INTO `nrdb40_scop_domains_nr` VALUES ("
TAG_INTERPRO="INSERT INTO `nrdb40_interpro_domains_nr` VALUES ("

nids = {}
for line in gzip.open( "nrdb40.dump.gz"):
    if line.startswith(TAG_NRDB40):
        # remove trailing ");\n"
        data = line[len(TAG_NRDB40):-3].split("),(")
        for value in data:
            # remove brackets
            nid = int(value.strip())
            nids[nid] = 1

print "nrdb40: read %i nids" % len(nids)

if not os.path.exists( "nrdb40.fasta.gz"):
    outfile = gzip.open( "nrdb40.fasta.gz", "w" )
    noutput, nskipped = 0, 0
    for line in gzip.open( "nrdb.dump.gz" ):
        if line.startswith(TAG_NRDB):
            # remove trailing ");\n"
            data = line[len(TAG_NRDB):-3].split("),(")
            for value in data:
                fields = value.split(",")
                nid = int(fields[1])
                # remove quotes
                sequence = fields[3][1:-1]
                if nid in nids:
                    outfile.write( ">%i\n%s\n" % (nid, sequence) )
                    noutput += 1
                else:
                    nskipped += 1

    outfile.close()

    print "fasta: output %i sequences (%i skipped, %i total)" % (noutput, nskipped, nskipped + noutput)

# domains
def doDomains( outfilename, infilename, tag, map_old2new = None ):

    if os.path.exists(outfilename):
        return
    
    outfile = gzip.open( outfilename, "w" )
    outfile.write("nid\tstart\tend\tfamily\n" )
    noutput, nskipped = 0, 0
    for line in gzip.open( infilename ):
        if line.startswith(tag):
            # remove trailing ");\n"
            data = line[len(tag):-3].split("),(")
            for value in data:
                fields = value.split(",")
                nid = int(fields[0])
                start = int(fields[1]) - 1
                end = int(fields[2])
                # remove quotes
                family = fields[3][1:-1]
                if map_old2new:
                    try:
                        family = map_old2new[family]
                    except KeyError:
                        nskipped += 1
                        continue
                if nid in nids:
                    outfile.write( "%i\t%i\t%i\t%s\n" % (nid, start, end, family) )
                    noutput += 1
                else:
                    nskipped += 1

    print "%s: output %i domains (%i skipped, %i total)" % (outfilename, noutput, nskipped, nskipped + noutput)

    outfile.close()


doDomains( "reference.domains.gz", "nrdb40_scop_domains_nr.dump.gz", TAG_SCOP )

map_interpro2pfam = {}
with open( "map_interpro2pfam", "r") as infile:
    for line in infile:
        try:
            id1, description, id_pfam = line[:-1].split("\t")
        except ValueError,msg:
            print "parsing error in: ", line,
            continue
            
        # there is a one-to-many relationship, which is ignored
        map_interpro2pfam[id1] = id_pfam

doDomains( "benchmark.domains.gz", "nrdb40_interpro_domains_nr.dump.gz", TAG_INTERPRO, map_interpro2pfam )

## get PFAM domains from http://www.ebi.ac.uk/interpro/ISearch?mode=db&query=H
## and convert into a three-column table mapping old to new in the first two
## columns and an optional description afterwards.
