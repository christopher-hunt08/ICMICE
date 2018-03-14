#!/usr/bin/env python

import MAUS

import math
import sys
import os


import ROOT
import array

import mickey
from mickey import tof_cuts
from mickey import scifi_cuts
from mickey import diffuser_cut
from mickey import banana_cut
from mickey import analysis_modules



if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  try :

    engine = mickey.Engine('emittance_reconstruction')

    engine.add_cut( tof_cuts.Cut_tof01_spacepoints() )
    engine.add_cut( tof_cuts.Cut_tof01_time() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_chisq_ndf() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_pt() )
    engine.add_cut( banana_cut.Cut_banana_plot_mass() )
#    engine.add_cut( diffuser_cut.Cut_diffuser_aperture() )

    engine.add_analysis( analysis_modules.EmittanceAnalysis() )

    namespace = engine.process_arguments()

    engine.go()

  except :
    raise

  print 
  print "Complete."
  print

