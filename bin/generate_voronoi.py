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

    engine = mickey.Engine('generate_voronoi')

    engine.add_cut( cuts_modules.SciFiUpstreamChisqNDF() )
    engine.add_cut( cuts_modules.SciFiUpstreamMomentum() )
    engine.add_cut( cuts_modules.TOF01Time() )

    engine.beam_selection()

#    engine.add_analysis( analysis_modules.OpticsAnalysis() )
    engine.add_analysis( analysis_modules.VoronoiTessellation() )

    namespace = engine.process_arguments()

    engine.go()

  except :
    raise

  print 
  print "Complete."
  print


