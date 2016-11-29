#######################################################
# tmFRET Utility
# Mario Rosasco, 2016
# University of Washington, Dept. Physiol & Biophysics
#######################################################

from pymol import cmd, stored
import math

def distance(p1, p2):
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

def markAccessible(model, chain, resi):
        # this is an arbitrary threshold for solvent accessibility.
        # Average accessible surface areas of residues are ~70-200 A^2
        # (http://www.proteinsandproteomics.org/content/free/tables_1/table08.pdf)
        sasaThreshold = 12.0
        selStr = model + " and chain " + chain + " and resi " + resi 
        cmd.set('dot_solvent', 1)
        sasa = cmd.get_area(selStr)
        if (sasa > sasaThreshold):
                cmd.show("spheres", selStr + " and name ca")
                cmd.color("red", selStr + " and name ca")
                print selStr + " SASA: " + str(cmd.get_area(selStr))
        return
        


def tmFRET(obj, chain, r1, r2, dist, force=False, strict=False):
        """
DESCRIPTION

    Finds all residues within dist of the midpoint between n1 and n2, and calculates their solvent accessible radii.

USAGE

    tmFRET obj, chain, r1, r2, dist[, force=False[, strict=False]]

EXAMPLES

    tmFRET 3j5p, c, 378, 380, 20            # Finds calphas w/in 20 ang. of the MBS between 378 and 380 on chain c of object 3j5p
        """
        # Identify target metal binding site residues
        res1str = "chain " + chain + " and resi " + r1
        res2str = "chain " + chain + " and resi " + r2

        cmd.set('dot_solvent', 1)
        sasa1 = cmd.get_area(res1str)
        sasa2 = cmd.get_area(res2str)

        if ((sasa1 < 10.0) or (sasa2 < 10.0)):
                if (not force):
                        print "WARNING: One or both of the residues you selected are likely buried."
                        print res1str + " SASA: " + str(sasa1)
                        print res2str + " SASA: " + str(sasa2)
                        print "If you wish to run the tmFRET helper on this pair anyway, please use:"
                        print "tmFRET object, chain, residue1, residue2, distance, force=True"
                        return
                else:
                        print "WARNING: being forced to run the tmFRET helper on a poorly accessible site:"
                        print res1str + " SASA: " + str(sasa1)
                        print res2str + " SASA: " + str(sasa2)
        else:
                print res1str + " SASA: " + str(sasa1)
                print res2str + " SASA: " + str(sasa2)
        print "Beginning mutagenesis. Please be patient, this may take some time..."

        # mutate residues to his
        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")

        cmd.get_wizard().set_mode("HIS")
        

        # Set up lists of nitrogen positions
        # also track CE1 posns to prevent odd orientations
        ND1_1list = []
        NE2_1list = []
        CE1_1list = []
        ND1_2list = []
        NE2_2list = []
        CE1_2list = []

        # find all posns for nitrogens for first histidine
        cmd.get_wizard().do_select(res1str)
        nRot = cmd.count_states() # Find number of possible rotomers
        for i in range(1,nRot+1):
                cmd.get_wizard().do_select(res1str)
                cmd.frame(i)
                cmd.get_wizard().apply()
                ND1pos = cmd.get_coords(res1str + " and name ND1")[0]
                NE2pos = cmd.get_coords(res1str + " and name NE2")[0]
                CE1pos = cmd.get_coords(res1str + " and name CE1")[0]
                pos1String = "[%3.2f,%3.2f,%3.2f]" % (ND1pos[0], ND1pos[1], ND1pos[2])
                #print res1str + " rotamer " + str(i) + " ND1 at " + pos1String
                pos2String = "[%3.2f,%3.2f,%3.2f]" % (NE2pos[0], NE2pos[1], NE2pos[2])
                #print res1str + " rotamer " + str(i) + " NE2 at " + pos2String
                ND1_1list.append(ND1pos)
                NE2_1list.append(NE2pos)
                CE1_1list.append(CE1pos)

        # find all posns for nitrogens for second histidine
        cmd.get_wizard().do_select(res2str)
        nRot = cmd.count_states() # Find number of possible rotomers
        for i in range(1,nRot+1):
                cmd.get_wizard().do_select(res2str)
                cmd.frame(i)
                cmd.get_wizard().apply()
                ND1pos = cmd.get_coords(res2str + " and name ND1")[0]
                NE2pos = cmd.get_coords(res2str + " and name NE2")[0]
                CE1pos = cmd.get_coords(res2str + " and name CE1")[0]
                pos1String = "[%3.2f,%3.2f,%3.2f]" % (ND1pos[0], ND1pos[1], ND1pos[2])
                #print res2str + " rotamer " + str(i) + " ND1 at " + pos1String
                pos2String = "[%3.2f,%3.2f,%3.2f]" % (NE2pos[0], NE2pos[1], NE2pos[2])
                #print res2str + " rotamer " + str(i) + " NE2 at " + pos2String
                ND1_2list.append(ND1pos)
                NE2_2list.append(NE2pos)
                CE1_2list.append(CE1pos)
        
        # iterate over all pairs of nitrogens, determining distance for each
        # Compare to distance for an "ideal" metal binding site, and save best pair
        idealDistance = 3.0
        bestDistance = 3.6 # use default bestDist to gate what range of values we'll accept
        delta = math.fabs(idealDistance - bestDistance)
        besti = -1
        bestj = -1
        r1N = ""
        r2N = ""
        print "r1rot,r2rot,ND1-ND1,ND1-NE2,NE2-ND1,NE2-NE2"
        for i in range(len(ND1_1list)): # iterate over 1st H rotamers
                for j in range(len(ND1_2list)): # iterate over 2nd H rotamers
                        # Debug diagnostic, or can be used to generate a printout to make a histogram in Excel
                        dND11_ND12 = distance(ND1_1list[i],ND1_2list[j])
                        dND11_NE22 = distance(ND1_1list[i],NE2_2list[j])
                        dNE21_ND12 = distance(NE2_1list[i],ND1_2list[j])
                        dNE21_NE22 = distance(NE2_1list[i],NE2_2list[j])
                        print "%i,%i,%3.2f,%3.2f,%3.2f,%3.2f" % (i, j, dND11_ND12, dND11_NE22, dNE21_ND12, dNE21_NE22)
                        tmpDist = distance(ND1_1list[i],ND1_2list[j])
                        tmpDelta = math.fabs(tmpDist-idealDistance)
                        if (tmpDelta < delta):
                                #check other nitrogen pair, to make sure they're facing away from MBS (ie: not clashing)
                                tmpDist2 = distance(NE2_1list[i],NE2_2list[j])
                                if (tmpDist2 > tmpDist):
                                        # check apical carbons, to make sure they're not pointing into MBS
                                        tmpDist3 = distance(CE1_1list[i],CE1_2list[j])
                                        if strict and (tmpDist3 > tmpDist):
                                                delta = tmpDelta
                                                bestDistance = tmpDist
                                                besti = i
                                                bestj = j
                                                r1N = "ND1"
                                                r2N = "ND1"
                        tmpDist = distance(ND1_1list[i],NE2_2list[j])
                        tmpDelta = math.fabs(tmpDist-idealDistance)
                        if (tmpDelta < delta):
                                tmpDist2 = distance(NE2_1list[i],ND1_2list[j])
                                if (tmpDist2 > tmpDist):
                                        tmpDist3 = distance(CE1_1list[i],CE1_2list[j])
                                        if strict and (tmpDist3 > tmpDist):
                                                delta = tmpDelta
                                                bestDistance = tmpDist
                                                besti = i
                                                bestj = j
                                                r1N = "ND1"
                                                r2N = "NE2"
                        tmpDist = distance(NE2_1list[i],ND1_2list[j])
                        tmpDelta = math.fabs(tmpDist-idealDistance)
                        if (tmpDelta < delta):
                                tmpDist2 = distance(ND1_1list[i],NE2_2list[j])
                                if (tmpDist2 > tmpDist):
                                        tmpDist3 = distance(CE1_1list[i],CE1_2list[j])
                                        if strict and (tmpDist3 > tmpDist):
                                                delta = tmpDelta
                                                bestDistance = tmpDist
                                                besti = i
                                                bestj = j
                                                r1N = "NE2"
                                                r2N = "ND1"
                        tmpDist = distance(NE2_1list[i],NE2_2list[j])
                        tmpDelta = math.fabs(tmpDist-idealDistance)
                        if (tmpDelta < delta):
                                tmpDist2 = distance(ND1_1list[i],ND1_2list[j])
                                if (tmpDist2 > tmpDist):
                                        tmpDist3 = distance(CE1_1list[i],CE1_2list[j])
                                        if strict and (tmpDist3 > tmpDist):
                                                delta = tmpDelta
                                                bestDistance = tmpDist
                                                besti = i
                                                bestj = j
                                                r1N = "NE2"
                                                r2N = "NE2"
                                
        if (besti == -1 or bestj == -1): #couldn't find a rotamer pair in the accepted distance range
                print "No appropriate rotamers could be found for the selected metal binding site residue pair. Please try again with different residues."
                # recover original residues here somehow, should go back and build this in
                cmd.set_wizard("done") #cleanup
                return
        
        # set up the MBS
        cmd.get_wizard().do_select(res1str)
        cmd.frame(besti+1)
        cmd.get_wizard().apply()
        cmd.get_wizard().do_select(res2str)
        cmd.frame(bestj+1)
        cmd.get_wizard().apply()
        cmd.set_wizard("done") #cleanup

        # make the metal as a pseudoatom
        posn1 = cmd.get_coords(res1str + " and name " + r1N)
        posn2 = cmd.get_coords(res2str + " and name " + r2N)

        # midpoint between 2 coordinating nitrogens is the MBS
        MBS = (posn1 + posn2) / 2
        MBS = MBS[0] # necessary b/c get_coords returns a coord list: [[x,y,z]]

        # make a pseudoatom at the MBS
        posString = "[%3.2f,%3.2f,%3.2f]" % (MBS[0], MBS[1], MBS[2])
        print "Making metal pseudoatom at" + posString
        
        cmd.pseudoatom("XitionMetal", pos=posString) 
        cmd.show("spheres", "XitionMetal")

        # find the residues within dist of the MBS
        selString="XitionMetal around " + dist + " and name cb"
        cmd.select("potentialSites", selString)

        # iterate over the residues in the selection, coloring Calpha by the solvent accessible surface area
        myspace={'markAccessible': markAccessible}
        cmd.iterate('(potentialSites)', 'markAccessible(model, chain, resi)', space=myspace)
        

# let pymol know about the function
cmd.extend("tmFRET", tmFRET)
