#######################################################
# Color Residues by Correlation Values
# Mario Rosasco, 2016
# University of Washington, Dept. Physiol & Biophysics
#######################################################
from pymol import cmd, stored
import Tkinter
import tkFileDialog
import tkMessageBox
import pickle

#######################################################
# Registers plugin, allowing access to main Tk loop
#######################################################
def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'Color By Correlation', label='Color residues according to a correlation array', command = lambda s=self : CBC(s))

#######################################################
# Helper Function:
# Called by iterate to set color of resi in chain
# in model according to corr value
#######################################################
def ColorByCorr(model, chain, resi, corr, maxCorr, minCorr):
    # Set up the RGB scaling
    dCorr = maxCorr - minCorr
    # min = white = [1.0, 1.0, 1.0]
    # max = purple = [0.5, 0.1, 0.6]
    # dCol = [-0.5, -0.9, -0.4]
    RSlope = -0.5/dCorr
    GSlope = -0.9/dCorr
    BSlope = -0.4/dCorr
    
    # Compute the appropriate color value
    getcontext().prec = 3
    selStr = model + " and chain " + chain + " and resi " + resi
    Rval = 1 if (RSlope*corr + 1) > 1 else (RSlope*corr + 1)
    Gval = 1 if (GSlope*corr + 1) > 1 else (GSlope*corr + 1)
    Bval = 1 if (BSlope*corr + 1) > 1 else (BSlope*corr + 1)
    
    # Color the residue accordingly
    cstring = "corrcol" + str(resi)
    cval = [Rval, Gval, Bval]
    cmd.set_color(cstring, cval)
    try:
        cmd.color(cstring, selStr)
    except:
        print "Could not generate color for RGB val: " + str(cval)
        return
    print selStr + " Correlation: " + str(corr)
    return
    
#######################################################
# Main Program Driver
#######################################################
def CBC(app):
    # first check and make sure the user has made a selection
    if (not 'sele' in cmd.get_names('selections')):
        tkMessageBox.showwarning(
            title='Color by Correlation', 
            message='No selection detected. This plugin requires a selection named (sele) to apply the per-residue coloring to.'
        )
        return

    # ask the user to select the file with the correlations
    # NB - this file should be an iterable of numerical values
    # with length greater than or equal to the length of the protein
    # sequence, prepared/saved using pickle.
    fname = tkFileDialog.askopenfilename(
        title='Color by Correlation',
        # this doesn't work, even though askopenfilename should accept a '-message' option.
        # message='Please select the file where the correlations were pickled', 
        parent=app.root
    )
    print fname
    C = []
    if fname:
        f = open(fname, "rb")
        C = pickle.load(f)
        f.close()
    else:
        return
    if C == []:
        print "The correlation array could not be loaded.\nPlease confirm that the file was saved using pickle."
        return
        
    # iterate over the residues in the selection, coloring by the correlation values
    myspace={'ColorByCorr': ColorByCorr, 'C': C, 'maxCorr': max(C), 'minCorr': min(C)}
    cmd.iterate('(sele)', 'ColorByCorr(model, chain, resi, C[resv], maxCorr, minCorr)', space=myspace)