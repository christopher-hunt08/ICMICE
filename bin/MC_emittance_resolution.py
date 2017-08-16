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

import MPA
from MPA import scifi_analysis

import ROOT


if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

##### 1. Intialise plots ######################################################
  print "\nInitialising..."

  try :

    analysis_engine = MPA.analysis_engine('emittance_analysis', description=
      'Performs a simple emittance reconstruction for all planes specified '+\
      'in the arguments.')

    analysis_engine.add_processor(MPA.scifi_analysis.scifi_emittance_reconstruction(use_mc=True))

    namespace = analysis_engine.process_arguments()

  except :
    raise
  else :
##### 2. Load Events ##########################################################
    print "\nLoading Spills...\n"
    file_reader = analysis_engine.get_file_reader()

    try :
      while file_reader.next_event() :
        try :
          sys.stdout.write( 
              ' Spill ' + str(file_reader.get_current_spill_number()) + \
              ' of ' + str(file_reader.get_current_number_spills()) + \
              ' in File ' + str(file_reader.get_current_filenumber()) + \
              ' of ' + str(file_reader.get_number_files()) + '             \r')

          sys.stdout.flush()

##### 3. Perform Analysis  ####################################################

          analysis_engine.analyse_event()

        except ValueError as ex:
          print "An Error Occured. Skipping Spill: " + \
                str(file_reader.get_current_spill_number()) + \
                " In File: " + str(file_reader.get_current_filenumber())
          print "ERROR =", ex
          print
          continue
    except KeyboardInterrupt :
      print
      print "Keyboard Interrupt"
      print

##### 4. Perform Alignment Calculation #######################################
    print "All Spills Loaded                                                  "
    print "\nStarting Analysis"

    try :
      analysis_engine.conclude()
    except ValueError as ex :
      print "Analysis Failed:", ex
      print
      print "Stopping Execution"

##### 5. Save plots and data ##################################################
    print "\nSaving Plots and Data"
    analysis_engine.save_plots()
    analysis_engine.save_data()

  print 
  print "Complete."
  print

