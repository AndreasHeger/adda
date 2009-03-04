USAGE = """python evaluate_domain_boundaries.py [OPTIONS]

Evaluate trees of domain partitions (or alternatively, parts or other benchmarking domains)

-s, --switch            switch between coverage of reference and size ratio if coverage is 1
-m,--skip_tms=          discard domains which contain transmembrane regions
-p, --parts=            table name of parts
-t, --trees=            table name of trees
-b, --bench=            table of domain table to be benchmarked (for example: nrdb40_domo_domains_nr)
-r, --reference=        table of reference table (for example: nrdb40_scop_domains_nr)
-o, --resolution=       resolution for scaling of domains
-q, --quality           take only sequences which are curated
--check_if_comparable   perform comparable check according to Islam95 (default level 85%)
--subset=               use only a subset of nids
"""

import string, getopt, sys, os, re, time, math, glob

import Experiment,Histogram

from Pairsdb import *
# from Table_nrdb90_masks import Table_nrdb90_masks
from Table_nrdb import Table_nrdb
from TableDomains import TableDomains

param_database = "pairsdb"
param_table_name_reference = None
param_table_name_trees = None
param_table_name_parts = None
param_table_name_bench = None
param_resolution = None
param_loglevel = 1
param_min_overlap = 1
param_switch = 0
param_combine_repeats = 1
param_skip_repeats = 0
param_skip_tms = 0
param_discard_full_length = 0
param_check_selection = 0
param_selection_threshold = 0.9
param_quality = None
param_no_full_length = None
param_only_full_length = None

## a full length domain should cover at least 90% of a sequence
param_min_length_ratio = 0.9

param_check_comparable = None
param_check_comparable_level = 0.85

param_bin_size = 1
param_subset = None

def BuildDomainArray( dbhandle, nid, length, table_name):
    
    statement = "SELECT domain_id, rep_from, rep_to FROM %s WHERE rep_nid = %i"\
                % (table_name, nid)
    
    domains = dbhandle.Execute(statement).fetchall()
    a = [0] * length

    max_domain = 0
    for domain_id,rep_from, rep_to in domains:
        a[rep_from-1:rep_to] = [domain_id] * (rep_to - rep_from + 1)
        max_domain = max(domain_id, max_domain)

    return a, max_domain

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "D:t:r::V:skb:meo:q:",
                                      ["Database=", "trees=", "reference=",
                                       "Verbose=", "switch", "parts=", "skip_repeats", "subset=",
                                       "bench=", "skip_tms", "check_selection", "resolution=", "quality",
                                       "no_full_length", "only_full_length", "check_if_comparable",
                                       "bin_size="])

    except getopt.error, msg:
        print USAGE
        print msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-t", "--trees" ):
            param_table_name_trees = a
        elif o in ("-D", "--Database"):
            param_database = a
        elif o in ("-r", "--parts"):
            param_table_name_parts = a
        elif o in ("-b", "--bench"):
            param_table_name_bench = a
        elif o in ("-r", "--reference"):
           param_table_name_reference = a
        elif o in ("-s", "--switch"):
            param_switch = 1
        elif o in ("-k", "--skip_repeats"):
            param_skip_repeats = 1
        elif o in ("-m", "--skip_tms"):
            param_skip_tms = 1
        elif o in ("e", "--check_selection"):
            param_check_selection = 1
        elif o in ("o", "--resolution"):
            param_resolution = string.atof(a)
        elif o in ("q", "--quality"):
            param_quality = 1
        elif o in ("-V", "--Verbose"):
            param_loglevel = string.atoi(a)
        elif o == "--no_full_length":
            param_no_full_length = 1
        elif o == "--only_full_length":
            param_only_full_length = 1
        elif o == "--check_if_comparable":
            param_check_comparable = 1
        elif o == "--bin_size":
            param_bin_size = string.atoi(a)
            
            
    dbhandle = Pairsdb()
    dbhandle.Connect()
    dbhandle.UseDatabase( param_database )

    print Experiment.GetHeader()
    print Experiment.GetParams()

    tbl_reference = TableDomains(dbhandle, "generic")
    tbl_reference.SetName(param_table_name_reference)
    
    # tbl_masks = Table_nrdb90_masks(dbhandle)
    tbl_nrdb = Table_nrdb(dbhandle)

    if param_table_name_trees:

        nids_statement = "SELECT DISTINCT s.nid FROM %s AS s, %s AS t %%s WHERE t.rep_nid = s.nid %%s" %\
                         (param_table_name_trees, param_table_name_reference)

        if param_quality:
            nids_statement = nids_statement % (", nrdb_quality AS q", "AND q.nid = s.rep_nid AND q.is_curated = 'T'")
        else:
            nids_statement = nids_statement % ("","")
            
        statement = """
        SELECT t.node, t.parent, t.level, t.xfrom, t.xto,
        ((LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i)) / (GREATEST( t.xto, %%i) - LEAST( t.xfrom, %%i))) AS ovl,
        ((LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i)) / (t.xto - t.xfrom)) AS cov_dom,
        ((LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i)) / (%%i - %%i)) AS cov_ref,
        ((t.xto - t.xfrom) / (%%i - %%i)) AS rat_ref
        FROM %s AS t
        WHERE t.nid = %%i
        AND (LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i) > %%i)
        ORDER BY ovl DESC
        LIMIT 1
        """ % (param_table_name_trees)
        
    elif param_table_name_parts or param_table_name_bench:

        if param_table_name_parts:
            table_name = param_table_name_parts
        else:
            table_name = param_table_name_bench

        if param_subset:
            nids_statement = "SELECT DISTINCT s.nid FROM %s AS s, %s AS t WHERE t.rep_nid = s.nid" % (param_subset, table_name)
        else:
            nids_statement = "SELECT DISTINCT s.rep_nid FROM %s AS s, %s AS t %%s WHERE t.rep_nid = s.rep_nid %%s" %\
                             (table_name, param_table_name_reference)

            if param_quality:
                nids_statement = nids_statement % (", nrdb_quality AS q", "AND q.nid = s.rep_nid AND q.is_curated = 'T'")
            else:
                nids_statement = nids_statement % ("","")

        statement = """
        SELECT 1, 0, 0, t.rep_from, t.rep_to,
        ((LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i)) / (GREATEST( t.rep_to, %%i) - LEAST( t.rep_from, %%i))) AS ovl,
        ((LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i)) / (t.rep_to - t.rep_from)) AS cov_dom,
        ((LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i)) / (%%i - %%i)) AS cov_ref,
        ((t.rep_to - t.rep_from) / (%%i - %%i)) AS rat_ref
        FROM %s AS t
        WHERE t.rep_nid = %%i
        AND (LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i) > %%i)
        ORDER BY ovl DESC
        LIMIT 1
        """ % (table_name)
    else:
        print "what shall I compare?"
        sys.exit(1)

    if param_check_selection:
        selection_statement = """
        SELECT t.domain_from, t.domain_to,
        ((LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i)) / (GREATEST( t.domain_to, %%i) - LEAST( t.domain_from, %%i))) AS ovl,
        ((LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i)) / (t.domain_to - t.domain_from)) AS cov_dom,
        ((LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i)) / (%%i - %%i)) AS cov_ref,
        ((t.domain_to - t.domain_from) / (%%i - %%i)) AS rat_ref
        FROM %s AS t
        WHERE t.domain_nid = %%i
        AND (LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i) > %%i)
        ORDER BY ovl DESC
        LIMIT 1
        """ % (param_table_name_parts)
        param_table_name_parts = None
        
        parts_same_as_trees, parts_larger_than_trees, parts_smaller_than_trees, parts_much_smaller_than_trees =  0,0,0,0

    
    nids = map(lambda x:x[0],dbhandle.Execute(nids_statement).fetchall())

    overlaps = []
    cov_doms = []
    cov_refs = []
    touched  = {}

    if param_check_selection:
        print "# NID\tDNODE\tDPARENT\tDLEVEL\tDFROM\tDTO\tRID\tRFROM\tRTO\tOVL\tDCOV\tRCOV\tRRCOV\tMRCOV"
    else:
        print "# NID\tDNODE\tDPARENT\tDLEVEL\tDFROM\tDTO\tRID\tRFROM\tRTO\tOVL\tDCOV\tRCOV\tRRCOV\tMRCOV"

    
    if param_loglevel >= 1:
        print "# --> processing %i nids" % len(nids)
        sys.stdout.flush()

    nskipped_no_assignments = 0
    nskipped_no_overlap = 0
    nskipped_wrong_domaintype = 0
    nfound = 0
    
    it = 0
    for nid in nids:

        it += 1
        # if it is > 1000: break

        if param_loglevel >= 2:
            print "# --> processing %i" % nid
            sys.stdout.flush()

        domains = tbl_reference.GetDomainBoundariesForNid( nid )

        length = tbl_nrdb.GetLength( nid )
        
        if not domains:
            nskipped_no_assignments +=1
            continue

        if param_no_full_length and len(domains) == 1:
            ## check if domain is actually full length, otherwise keep
            id, domain_from, domain_to = domains[0]
            if float(domain_to-domain_from+1) / float(length) >= param_min_length_ratio:
                nskipped_wrong_domaintype += 1
                continue
            
        if param_only_full_length:
            if len(domains) == 1:
                id, domain_from, domain_to = domains[0]
                if float(domain_to-domain_from+1) / float(length) <= param_min_length_ratio:
                    nskipped_wrong_domaintype += 1
                    continue
            else:
                nskipped_wrong_domaintype += 1                
                continue

        nfound += 1
        
        last_id = None
        x = 0

        while x < len(domains):
            
            id, domain_from, domain_to = domains[x]
                
            is_repeat = -1
            
            while x < len(domains) and domains[x][0] == id:
                domain_to = domains[x][2]
                x += 1
                is_repeat += 1

            if param_skip_repeats and is_repeat:
                continue

            # if param_skip_tms and tbl_masks.HasMask( nid, 2, domain_from, domain_to):
            #    continue

            if param_resolution:
                xdomain_from = int(float(domain_from-1)/param_resolution)
                xdomain_to   = int(float(domain_to-1)/param_resolution) + 1
            else:
                xdomain_from = domain_from
                xdomain_to   = domain_to

            if param_loglevel >= 2:
                print "--> processing domain %s (%i-%i) (%i-%i)" % ( id, domain_from, domain_to, xdomain_from, xdomain_to)
                sys.stdout.flush()

            s = statement % (xdomain_to, xdomain_from, xdomain_to, xdomain_from,
                             xdomain_to, xdomain_from,
                             xdomain_to, xdomain_from, xdomain_to, xdomain_from,
                             xdomain_to, xdomain_from,
                             nid,
                             xdomain_to, xdomain_from, param_min_overlap)

            if param_loglevel >= 4:
                print s
                
            result = dbhandle.Execute(s).fetchone()

            if not result:
                continue

            node, parent, level, xfrom, xto, overlap, cov_dom, cov_ref, rat_ref = result

            key = "%i-%s-%i-%i" % (nid, id, xfrom, xto)
            if touched.has_key(key):
                continue
            else:
                touched[key] = 1

            # discard full length domains
            if param_discard_full_length:
                if param_table_name_trees:            
                    if node == 0: continue
                else:
                    if length == xto - xfrom + 1: continue
            
            if param_switch and cov_ref == 1.0:
                xcov_ref = rat_ref
            else:
                xcov_ref = cov_ref
                
            # check, if selection did take a domain lower or further up
            if param_check_selection:
                yfrom = (xfrom * 10) + 1
                yto   = min(xto * 10 + 1, length)
                s = selection_statement % (yto, yfrom, yto, yfrom,
                                           yto, yfrom,
                                           yto, yfrom, yto, yfrom,
                                           yto, yfrom,
                                           nid,
                                           yto, yfrom, param_min_overlap)
                
                result = dbhandle.Execute(s).fetchone()
                if result:
                    parts_from, parts_to, ovl_parts, cov_parts, cov_tree, rat_parts = result


                    if rat_parts > 1.0:
                        parts_larger_than_trees += 1
                        token = ">"
                    elif rat_parts == 1.0:
                        parts_same_as_trees += 1
                        token = "="
                    else:
                        parts_smaller_than_trees += 1
                        token = "<"
                        if rat_parts < param_selection_threshold:
                            parts_much_smaller_than_trees += 1

                    sys.stdout.write(string.join(map(str, (nid,
                                                id, domain_from, domain_to,
                                                level,
                                                yfrom, yto,
                                                parts_from, parts_to,
                                                overlap, cov_dom, cov_ref, rat_ref, xcov_ref,
                                                ovl_parts, cov_parts, cov_tree, rat_parts,
                                                token)), "\t") + "\n")
                        
            else:
            # this is actually the default
                sys.stdout.write(string.join(map(str, (nid, node, parent, level, xfrom, xto,
                                            id,
                                            xdomain_from, xdomain_to,
                                            overlap, cov_dom, cov_ref, rat_ref, xcov_ref)), "\t") + "\n")
            
                overlaps.append( int(overlap * 100) )
                cov_doms.append( int(cov_dom * 100))
                cov_refs.append( int(xcov_ref * 100))            

        
    if param_loglevel >= 1:
        print "## skipped nids because of no overlap with reference: %i" % nskipped_no_overlap
        print "## skipped nids because of no assignments: %i" % nskipped_no_assignments
        print "## skipped nids because of wrong domain type: %i" % nskipped_wrong_domaintype
        print "## nids in comparison: %i" % nfound
        
    if param_check_selection:
        print "## parts larger than trees=", parts_larger_than_trees
        print "## parts like trees=", parts_same_as_trees
        print "## parts smaller than trees=", parts_smaller_than_trees                
        print "## parts much smaller than trees (<%f)=" % param_selection_threshold, parts_much_smaller_than_trees
        
        
    else:
        print "## histogram over overlaps"
        Histogram.Print(Histogram.Calculate( overlaps, min_value=0, increment=10, no_empty_bins = 1))

        average = float(reduce(lambda x,y: x+y, overlaps))/float(len(overlaps))
        print "min=%i\tmax=%i\tave=%f" % (min(overlaps), max(overlaps), average)

        print "## histogram over domain coverage"
        Histogram.Print(Histogram.Calculate( cov_doms, min_value=0, increment=10, no_empty_bins = 1))

        average = float(reduce(lambda x,y: x+y, cov_doms))/float(len(cov_doms))
        print "min=%i\tmax=%i\tave=%f" % (min(cov_doms), max(cov_doms), average)

        print "## histogram over reference coverage"
        Histogram.Print(Histogram.AddRelativeAndCumulativeDistributions(Histogram.Calculate( cov_refs, min_value=0, increment=param_bin_size, no_empty_bins = 1)))

        average = float(reduce(lambda x,y: x+y, cov_refs))/float(len(cov_refs))
        print "min=%i\tmax=%i\tave=%f" % (min(cov_refs), max(cov_refs), average)

    




















