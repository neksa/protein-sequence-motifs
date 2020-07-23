#
# PyMOL plugin
# Smoothing Calpha coordinates in selection by a sliding window over 8 Calphas
#
#
#

from pymol import cmd, stored

def smoothca( selection, window = 4 ):
    print "Hello from Smoothing plugin\n"
    stored.alphaCarbons = []

    ca_selection = selection + " and n. CA"

    cmd.iterate_state(1, selector.process(ca_selection), "stored.alphaCarbons.append([x,y,z])")

    ca_list = stored.alphaCarbons[:]
    print "Count alpha carbons:", len(ca_list)
    for i in range(window, len(ca_list)-window):
        (x, y, z) = ca_list[i]
        subset = ca_list[i-window:i+window]
        sx, sy, sz = 0, 0 , 0
        #print "subset", len(subset)
        for (bx, by, bz) in subset:
            sx += bx
            sy += by
            sz += bz
        ax, ay, az = sx/len(subset), sy/len(subset), sz/len(subset)
        #print "old", x, y, z
        #print "new", ax, ay, az
        stored.alphaCarbons[i] = (ax, ay, az)

    #objName = cmd.identify(ca_selection, 1)[0][0]
    cmd.alter_state(1, ca_selection, "(x,y,z)=stored.alphaCarbons.pop(0)")

cmd.extend( "smoothca", smoothca );
