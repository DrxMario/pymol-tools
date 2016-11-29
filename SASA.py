#######################################################
# SASA Utility
# Mario Rosasco, 2016
# University of Washington, Dept. Physiol & Biophysics
#######################################################
import Tkinter
from Tkinter import *
import tkSimpleDialog
import tkMessageBox
from pymol import cmd, stored
import math


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'SASA', label='Compute Solvent Accessible Area', command = lambda s=self : SASA(s))

def SASA(app):
	# first check and make sure the user has made a selection
	if (not 'sele' in cmd.get_names('selections')):
		tkMessageBox.showwarning(
			'SASA Utility', 
			'No selection detected. This plugin requires a selection named (sele) to measure the surface area for.'
		)
		return

	# ask what size rolling ball the user would like to use for the solvent
	sol_rad = tkSimpleDialog.askstring(
		'SASA Utility', 
		'Please enter the solvent radius (in Angstrom):', 
		parent=app.root, 
		initialvalue='1.4'
	)
	if sol_rad == None:
		return
	else:
		sol_rad = float(sol_rad)
	
	
	dot_dens = tkSimpleDialog.askstring(
		'SASA Utility', 
		'Please enter the sampling density (1-4).\nHigher density is more accurate,\nbut slower for large selections:', 
		parent=app.root, initialvalue='3'
	)
	
	if dot_dens == None:
		return
	else:
		dot_dens = float(dot_dens)
	
	if dot_dens > 4:
		dot_dens = 4
		print 'Warning: selected sampling density greater than max allowed val (4). Defaulting to 4.'

	cmd.set('dot_solvent', 1)
	cmd.set('dot_density', dot_dens)
	cmd.set('solvent_radius', sol_rad)
	sasa_res = cmd.get_area('(sele)')
	
	# print out the results
	sasa_str = '%.2f' % sasa_res
	print 'Solvent accessible surface area for the given selection is ' + sasa_str + ' square Angstrom'
