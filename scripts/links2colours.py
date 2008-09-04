"""add colours to a list of nids together with labels.

-m, --multi_labels      if label is a comma-separated list, use grey for for colours.
-l, --legend            add legend to color list
"""


import sys,os,string,getopt

colors=( "red", "blue", "green", "yellow", "cyan", "magenta",
         "lightblue", "lightred", "lightgreen", "orange")

increment_colors= ("lightgrey", "darkgrey")

param_multi_labels = None

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "ml",
                                      ["multi_labels", "legend"])

    except getopt.error, msg:
        print USAGE
        sys.exit(2)

    for o,a in optlist:
        if o in ("m", "--multi_labels"):
            param_multi_labels = 1
        if o in ("l", "--legend"):
            param_legend = 1
    
    current_color = 0

    labels = {}
    labels["18"] = "black"
    labels["na"] = "white"
    while 1:

        line = sys.stdin.readline()
        if not line: break

        if line[0] == "#": continue
        
        id, label = string.split(line[:-1], "\t")[:2]

        color = None
        
        if param_multi_labels:
            n = string.count(label, ",")
            if n:
                n = min(n+1, len(increment_colors))
                color = increment_colors[n-1]
                print "%s\t%s" % (id, color)
                continue
            
        if not color:
            if labels.has_key(label):
                color = labels[label]
            else:
                labels[label] = colors[current_color]
                current_color += 1
                if current_color >= len(colors):
                    current_color = 0

        print "%s\t%s" % (id, labels[label])


    if param_legend:
        print "# legend:",
        for x in labels.keys():
            if x == "18": continue
            print "%s=%s;" % (x, labels[x]),
        print

