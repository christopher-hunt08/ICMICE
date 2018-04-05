#!/usr/bin/env python

import MAUS

import math
import sys
import os


import ROOT
import array

import mickey
from mickey import analysis_modules
from mickey import cuts_modules



if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  try :

    engine = mickey.Engine('emittance_reconstruction')

    engine.add_cut( cuts_modules.TOF01Spacepoints() )
    engine.add_cut( cuts_modules.TOF01Time() )
    engine.add_cut( cuts_modules.SciFiUpstreamChisqNDF() )
    engine.add_cut( cuts_modules.SciFiUpstreamPt() )
    engine.add_cut( cuts_modules.BananaPlot() )
    engine.add_cut( cuts_modules.SciFiRefitStatus() )
    engine.add_cut( cuts_modules.DiffuserAperture() )

    engine.add_analysis( analysis_modules.EmittanceAnalysis() )

    namespace = engine.process_arguments()

    engine.go()

  except :
    raise

  print 
  print "Complete."
  print

