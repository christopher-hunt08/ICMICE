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

    engine = mickey.Engine('emittance_systematic_error')

    engine.add_cut( scifi_cuts.Cut_scifi_upstream_chisq_ndf() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_pt() )

    engine.add_analysis( analysis_modules.EmittanceSystematicErrors() )

    namespace = engine.process_arguments()

    engine.go()

  except :
    raise

  print 
  print "Complete."
  print

