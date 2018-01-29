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


if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  print "\nInitialising..."

  try :

    engine = mickey.Engine('testing_mickey')

    namespace = engine.process_arguments()

    print
    print namespace
    print

  except :
    raise
  else :
    print "\nLoading Spills...\n"

    try :
      while engine.next_event() :
        try :


          engine.analyse_event()


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


    try :
      engine.conclude()
    except ValueError as ex :
      print "Analysis Failed:", ex
      print
      print "Stopping Execution"

  print 
  print "Complete."
  print

