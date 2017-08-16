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
"""
  This analysis script will perform the first reconstrucion of straight trcks
  from the the MICE SciFi Trackers.

  It is critical to have a fast and simple way to check the tracker-tracker 
  alignment. 


  Script Aglorithm :

   - Load MAUS Configuration files. We need the correct geometry! - Or Do we?
   - Iterate through the SciFi Events
   - Locate single, straight track events, where we have one track in each of
     the upstream and downstream trackers.
   - Apply and sanity checks/cuts i.e. P-Value, Angle cuts, etc
   - Use upstream tracker to predict downstream tracker reconstruction.
   - Plot differences between the two reconstructions. Deviations will point to
     misalignments between the two trackers.
   
     [ AND FINALLY... ]
   - Measure the Earths Magnetic field local to the MICE Hall!
"""


# Import MAUS framework. 
import MAUS

# Generic Python imports
import math
import sys
import os
import argparse

# Third Party library import statements
import event_loader
import MPA
from MPA import scifi_analysis
from MPA import tof_analysis
from MPA import alignment_analysis

import ROOT
import json



# Useful Constants and configuration
STRAIGHT_ALGORITHM_ID = 0
STARTING_TRACKER = 0
RECON_STATION = 1
RECON_PLANE = 0
MIN_NUMBER_TRACKPOINTS = 0
TRACKER_LENGTH = 1101.0624

# Alignment Info
DELTA_X = 0.0
DELTA_Y = 0.0
DELTA_PHI_X = 0.0
DELTA_PHI_Y = 0.0
#DELTA_THETA = -2.0 * math.pi / 3.0
TRACKER_SEPARATION = -1.0

#IGNORE_PLANES = [ -7, -8, -9 ]
IGNORE_PLANES = []




if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

##### 1. Intialise plots ######################################################
  print "\nInitialising..."

  try :

    analysis_engine = MPA.analysis_engine("straight_track_alignment", description=
      'Performs a straight track reconstruction using the MICE SciFi '+\
      'Trackers in order to measure the gross misalignments between the '+\
      ' up- and downstream trackers')

    analysis_engine.add_processor(MPA.alignment_analysis.scifi_straight_track_alignment())
#    analysis_engine.add_processor(MPA.emr_scifi_analysis.emr_scifi_correlations())
    analysis_engine.add_processor(MPA.tof_scifi_analysis.tof_scifi_analyser())
    analysis_engine.add_processor(MPA.scifi_analysis.scifi_internal_alignment())

    namespace = analysis_engine.process_arguments()

  except :
    raise
  else :
##### 2. Load Events ##########################################################
    print "\nLoading Spills...\n"
    file_reader = analysis_engine.get_file_reader()
    file_reader.set_print_progress('spill')

    try :
      while analysis_engine.next_event() :
        try :
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




#def analyse_plots(plot_dict, data_dict) :
#  """
#    Use the plot dictionary to perform some analysis of the data
#  """
#
## Fixed Tracker Parameterisation
#  plot_dict['x_dist_profile'] = plot_dict['separation']['position_gradient_x'].ProfileX(\
#                                                              'x_dist_profile')
#  plot_dict['y_dist_profile'] = plot_dict['separation']['position_gradient_y'].ProfileX(\
#                                                              'y_dist_profile')
#
#  plot_dict['dist_profile'] = scifi.gaussian_profile_x( plot_dict['separation']['position_gradient'], 0.0, 150.0 )
#
#  x_fit_result = plot_dict['x_dist_profile'].Fit('pol1', 'QS', "", -0.01, 0.01)
#  y_fit_result = plot_dict['y_dist_profile'].Fit('pol1', 'QS', "", -0.01, 0.01)
#  fit_result = plot_dict['dist_profile'].Fit('pol1', 'QS', "", -0.0, 0.050)
#
#  x_fit_length = x_fit_result.Parameter(1)
#  x_fit_length_err = x_fit_result.ParError(1)
#  y_fit_length = y_fit_result.Parameter(1)
#  y_fit_length_err = y_fit_result.ParError(1)
#
#  fit_length = fit_result.Parameter(1)
#  fit_length_err = fit_result.ParError(1)
#
#  length = plot_dict['separation']['tracker_separation'].GetMean()
#  num_tracks = plot_dict['separation']['tracker_separation'].GetEntries()
#
#  x_delta_raw = plot_dict['residual']['xy'].GetMean(1)
#  y_delta_raw = plot_dict['residual']['xy'].GetMean(2)
#  x_delta_raw_err = plot_dict['residual']['xy'].GetRMS(1) / \
#                             math.sqrt( num_tracks )
#  y_delta_raw_err = plot_dict['residual']['xy'].GetRMS(2) / \
#                             math.sqrt( num_tracks )
#
#  x_delta_phi = plot_dict['residual']['txty'].GetMean(1)
#  y_delta_phi = plot_dict['residual']['txty'].GetMean(2)
#  x_delta_phi_err = plot_dict['residual']['txty'].GetRMS(1) / \
#                           math.sqrt( num_tracks )
#  y_delta_phi_err = plot_dict['residual']['txty'].GetRMS(2) / \
#                           math.sqrt( num_tracks )
#
#  tracker_separation = (x_fit_length + y_fit_length) / 2.0
#  tracker_separation_err = math.sqrt(x_fit_length_err**2 \
#                                                   + y_fit_length_err**2) / 2.0
#
#
#  alignment_dict = data_dict['alignment']
#  alignment_dict['distance'] = length
#
#  alignment_dict['x_delta'] += x_delta_raw
#  alignment_dict['y_delta'] += y_delta_raw
#  alignment_dict['x_delta_err'] = x_delta_raw_err
#  alignment_dict['y_delta_err'] = y_delta_raw_err
#
#  alignment_dict['x_delta_phi'] += x_delta_phi
#  alignment_dict['y_delta_phi'] += y_delta_phi
#  alignment_dict['x_delta_phi_err'] = x_delta_phi_err
#  alignment_dict['y_delta_phi_err'] = y_delta_phi_err
#
#  alignment_dict['tracker_separation'] = tracker_separation
#  alignment_dict['tracker_separation_err'] = tracker_separation_err
#
#  alignment_dict['x_fit_length'] = x_fit_length
#  alignment_dict['x_fit_length_err'] = x_fit_length_err
#  alignment_dict['y_fit_length'] = y_fit_length
#  alignment_dict['y_fit_length_err'] = y_fit_length_err
#
#
#  print
#  print "Downstream Tracker Misalignments:"
#  print
#  print "Assuming Center of tracker reference plane is stationary."
#  print
#  print "Distance Between Trackers Set To:", length, "mm"
#  print
#  print "Analysed {0:0.0f} Tracks".format(num_tracks)
#  print
#  print "Fitted Length, X  =", "{0:0.0f} +/- {1:0.0f}".format(\
#                                 x_fit_length + length, x_fit_length_err), "mm"
#  print "Fitted Length, Y  =", "{0:0.0f} +/- {1:0.0f}".format(\
#                                 y_fit_length + length, y_fit_length_err), "mm"
#  print
#  print "Fitted Length     =", "{0:0.0f} +/- {1:0.0f}".format(fit_length, \
#                                                          fit_length_err), "mm"
#  print
#  print "Ref Plane X Shift =", "{0:0.3e} +/- {1:0.3e}".format(x_delta_raw, \
#                                                         x_delta_raw_err), "mm"
#  print "Ref Plane Y Shift =", "{0:0.3e} +/- {1:0.3e}".format(y_delta_raw, \
#                                                         y_delta_raw_err), "mm"
#  print
#  print "Rotation in X-Z   =", "{0:0.3e} +/- {1:0.3e}".format(x_delta_phi, \
#                                                        x_delta_phi_err), "rad"
#  print "Rotation in Y-Z   =", "{0:0.3e} +/- {1:0.3e}".format(y_delta_phi, \
#                                                        y_delta_phi_err), "rad"
#
#  print
#  print "Recalculated Algnments:"
#  print
#  print "Ref Plane X Shift =", "{0:0.3e} +/- {1:0.3e}".format( \
#                alignment_dict['x_delta'], alignment_dict['x_delta_err']), "mm"
#  print "Ref Plane Y Shift =", "{0:0.3e} +/- {1:0.3e}".format( \
#                alignment_dict['y_delta'], alignment_dict['y_delta_err']), "mm"
#  print
#  print "Rotation in X-Z   =", "{0:0.3e} +/- {1:0.3e}".format( \
#        alignment_dict['x_delta_phi'], alignment_dict['x_delta_phi_err']), "mm"
#  print "Rotation in Y-Z   =", "{0:0.3e} +/- {1:0.3e}".format( \
#        alignment_dict['y_delta_phi'], alignment_dict['y_delta_phi_err']), "mm"
#  print 
#
#  return data_dict
#
#
#def load_data(filename) :
#  """
#    Save all the data to a json file
#  """
#  data_dict = None
#  with open(filename, 'r') as infile :
#    data_dict = json.load(infile)
#
## Reset Counters
#  for counter in data_dict['counters'] :
#    data_dict['counters'][counter] = 0
#
#  return data_dict
#
# 
#def save_data(data_dict, directory, filename) :
#  """
#    Save all the data to a json file
#  """
#  filename = os.path.join(directory, filename+".json")
#  with open(filename, 'w') as outfile :
#    json.dump(data_dict, outfile)

