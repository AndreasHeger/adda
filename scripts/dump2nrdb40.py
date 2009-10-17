# extract data from nrdb40.dump.gz and nrdb.dump.gz
# for nrdb40.fasta

import gzip

TAG_NRDB40="INSERT INTO `nrdb40` VALUES ("
TAG_NRDB="INSERT INTO `nrdb` VALUES ("
TAG_SCOP="INSERT INTO `nrdb40_scop_domains_nr` VALUES ("


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
outfile = gzip.open( "reference.domains.gz", "w" )
outfile.write("nid\tstart\tend\tfamily\n" )
noutput, nskipped = 0, 0
for line in gzip.open( "nrdb40_scop_domains_nr.dump.gz" ):
    if line.startswith(TAG_SCOP):
        # remove trailing ");\n"
        data = line[len(TAG_SCOP):-3].split("),(")
        for value in data:
            fields = value.split(",")
            nid = int(fields[0])
            start = int(fields[1]) - 1
            end = int(fields[2])
            # remove quotes
            family = fields[3][1:-1]
            if nid in nids:
                outfile.write( "%i\t%i\t%i\t%s\n" % (nid, start, end, family) )
                noutput += 1
            else:
                nskipped += 1
    
print "domains: output %i domains (%i skipped, %i total)" % (noutput, nskipped, nskipped + noutput)

outfile.close()


