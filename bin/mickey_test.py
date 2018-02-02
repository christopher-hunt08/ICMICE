#!/usr/bin/env python

# This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
# MAUS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MAUS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MAUS.  If not, see <http://www.gnu.org/licenses/>.
#

import MAUS

import math
import sys
import os

import mickey

import ROOT

from mickey import tof_cuts
from mickey import scifi_cuts
from mickey import diffuser_cut
from mickey import analysis_modules

if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  print "\nInitialising..."

  try :

    engine = mickey.Engine('testing_mickey')

    engine.add_cut( tof_cuts.Cut_tof01_spacepoints() )
    engine.add_cut( tof_cuts.Cut_tof01_time() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_chisq_ndf() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_pt() )
    engine.add_cut( diffuser_cut.Cut_diffuser_aperture() )

    engine.add_analysis( analysis_modules.EmittanceAnalysis() )
#    engine.add_analysis( analysis_modules.MCEmittanceAnalysis() )

    namespace = engine.process_arguments()

    print
    print namespace
    print

    engine.go()

  except :
    raise

  print 
  print "Complete."
  print

