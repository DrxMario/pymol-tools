from pymol import cmd, stored

def zero_residues_sub(sel1, start=0, end=0, offset=0,chains=0):
        """
DESCRIPTION

    Renumbers the residues so that the given residue range starts at zero, or offset

USAGE

    zero_residues_sub selection, start, end [, offset [, chains ]]

EXAMPLES

    zero_residues_sub protName, 0, 10            # first residue is 0
    zero_residues_sub protName, 0, 10, 5         # first residue is 5
    zero_residues_sub protName, 0, 10, chains=1  # each chain starts at 0
    zero_residues_sub *
        """
        offset = int(offset)

        # variable to store the offset
        stored.first = None
        # get the names of the proteins in the selection

        names = ['(model %s and (%s))' % (p, sel1)
                        for p in cmd.get_object_list('(' + sel1 + ')')]

        if int (chains):
                names = ['(%s and chain %s)' % (p, chain)
                                for p in names
                                for chain in cmd.get_chains(p)]

        # for each name shown
        for p in names:
                # get this offset
                ok = cmd.iterate("first %s and polymer and n. CA" % p,"stored.first=resv")
                # don't waste time if we don't have to
                #if not ok or stored.first == offset:
                if not ok:
                        continue;
                # reassign the residue numbers
                p = p+ " and resi " + start + "-" + end
                cmd.alter("%s" % p, "resi=str(int(resi)-%s)" % str(int(start)-offset))
                # update pymol

        cmd.rebuild()

# let pymol know about the function
cmd.extend("zero_residues_sub", zero_residues_sub)
