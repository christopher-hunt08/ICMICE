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


# pylint: disable = W0311, E1101, W0613, C0111, W0621, C0103, W0702, W0611
# pylint: disable = W0603, R0911, R0912, R0914, R0915

# Import MAUS framework. 
import MAUS

# Generic Python imports
import math
import sys
import os
import argparse

# Third Party library import statements
import event_loader
import scifi_tools as scifi
import scifi_analysis
import ROOT
import json



# Useful Constants and configuration
STRAIGHT_ALGORITHM_ID = 0
STARTING_TRACKER = 0
RECON_STATION = 1
RECON_PLANE = 0
MIN_NUMBER_TRACKPOINTS = 0
TRACKER_LENGTH = 1101.0624
TOF_CUT_LOW = 0.0
TOF_CUT_HIGH = 100.0
MX_CUT = 1.05
MY_CUT = 1.05
GRAD_CUT = 1.05
RADIUS_CUT = 200.0
PROJECTED_RADIUS_CUT = 200.0
LOW_GRAD_CUT = 0.0
P_VALUE_CUT = 0.00

tof_ana = None
tof_scifi_ana = None

# Alignment Info
DELTA_X = 0.0
DELTA_Y = 0.0
DELTA_PHI_X = 0.0
DELTA_PHI_Y = 0.0
#DELTA_THETA = -2.0 * math.pi / 3.0
TRACKER_SEPARATION = -1.0

UP_DELTA_THETA = 2.0 * math.pi / 3.0
DOWN_DELTA_THETA = -2.0 * math.pi / 3.0
#UP_DELTA_THETA = 0.0
#DOWN_DELTA_THETA = 0.0

#IGNORE_PLANES = [ -7, -8, -9 ]
IGNORE_PLANES = []

def rotate_by_theta( vector, theta ) :
  new_vector = [ 0.0, 0.0 ]
  cos = math.cos( theta )
  sin = math.sin( theta )

  new_vector[0] = vector[0] * cos - vector[1] * sin
  new_vector[1] = vector[0] * sin + vector[1] * cos

  return new_vector


def up_tracker_correction( vector ) :
  new_vec = vector
  return new_vec
  new_vec = rotate_by_theta( vector, UP_DELTA_THETA )
  new_vec[0] = -new_vec[0]
  new_vec[1] = -new_vec[1]
  

def down_tracker_correction( vector ) :
  new_vec = vector
  return new_vec
  new_vec = rotate_by_theta( vector, DOWN_DELTA_THETA )
  new_vec[0] = -new_vec[0]
  new_vec[1] = -new_vec[1]


def init_plots() :
  """
    Initialises a dictionary of plots for later use
  """
  plot_dict = {}

  correlation_plots = {}


  correlation_plots['x_delta_x'] = ROOT.TH2F('x_delta_x', \
              "X Against Change in X", 200, -200.0, 200.0, 200, -200.0, 200.0 )
  correlation_plots['y_delta_y'] = ROOT.TH2F('y_delta_y', \
              "Y Against Change in Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
  correlation_plots['x_delta_y'] = ROOT.TH2F('x_delta_y', \
              "X Against Change in Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
  correlation_plots['y_delta_x'] = ROOT.TH2F('y_delta_x', \
              "Y Against Change in X", 200, -200.0, 200.0, 200, -200.0, 200.0 )

  correlation_plots['delta_x_residual_x'] = ROOT.TH2F('delta_x_residual_x', \
     "Delta X Against Residual X", 200, -200.0, 200.0, 200, -200.0, 200.0 )
  correlation_plots['delta_y_residual_y'] = ROOT.TH2F('delta_y_residual_y', \
     "Delta Y Against Residual Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
  correlation_plots['delta_x_residual_y'] = ROOT.TH2F('delta_x_residual_y', \
     "Delta X Against Residual Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
  correlation_plots['delta_y_residual_x'] = ROOT.TH2F('delta_y_residual_x', \
     "Delta Y Against Residual X", 200, -200.0, 200.0, 200, -200.0, 200.0 )


  separation_plots = {}

  separation_plots['tracker_separation'] = ROOT.TH1F('tracker_separation', \
                        "Separation of the two trackers", 10000, 0.0, 10000.0 )

  separation_plots['position_gradient_x'] = ROOT.TH2F('position_gradient_x', \
           "Change in X Against Gradient", 500, -1.0, 1.0, 500, -500.0, 500.0 )
  separation_plots['position_gradient_y'] = ROOT.TH2F('position_gradient_y', \
           "Change in Y Against Gradient", 500, -1.0, 1.0, 500, -500.0, 500.0 )
  separation_plots['position_gradient'] = ROOT.TH2F('position_gradient', \
     "Change in Position Against Gradient", 100, 0.0, 0.5, 100, -300.0, 300.0 )


  upstream_plots = {}

  upstream_plots['xy'] = ROOT.TH2F('upstream_xy', \
           'Upstream Beam Position', 1000, -200.0, 200.0, 1000, -200.0, 200.0 )
  upstream_plots['mxmy'] = ROOT.TH2F('upstream_mxmy', \
                   'Upstream Beam Gradient', 1000, -1.0, 1.0, 1000, -1.0, 1.0 )
  upstream_plots['theta'] = ROOT.TH1F('upstream_theta', \
                                     'Upstream Beam Rotation', 1000, -4.0, 4.0) 
  upstream_plots['xmy'] = ROOT.TH2F('upstream_xmy', \
                   'Upstream Beam X-My', 1000, -200.0, 200.0, 1000, -1.0, 1.0 )
  upstream_plots['ymx'] = ROOT.TH2F('upstream_ymx', \
                   'Upstream Beam Y-Mx', 1000, -200.0, 200.0, 1000, -1.0, 1.0 )

  downstream_plots = {}

  downstream_plots['xy'] = ROOT.TH2F('downstream_xy', \
         'Downstream Beam Position', 1000, -200.0, 200.0, 1000, -200.0, 200.0 )
  downstream_plots['mxmy'] = ROOT.TH2F('downstream_mxmy', \
                 'Downstream Beam Gradient', 1000, -1.0, 1.0, 1000, -1.0, 1.0 )
  downstream_plots['theta'] = ROOT.TH1F('downstream_theta', \
                                   'Downstream Beam Rotation', 1000, -4.0, 4.0) 
  downstream_plots['xmy'] = ROOT.TH2F('downstream_xmy', \
                 'Downstream Beam X-My', 1000, -200.0, 200.0, 1000, -1.0, 1.0 )
  downstream_plots['ymx'] = ROOT.TH2F('downstream_ymx', \
                 'Downstream Beam Y-Mx', 1000, -200.0, 200.0, 1000, -1.0, 1.0 )

  projected_plots = {}

  projected_plots['xy'] = ROOT.TH2F('projected_xy', \
           'Projected Beam Position', 1000, -200.0, 200.0, 1000, -200.0, 200.0) 
  projected_plots['mxmy'] = ROOT.TH2F('projected_mxmy', \
                   'Projected Beam Gradient', 1000, -1.0, 1.0, 1000, -1.0, 1.0) 
  projected_plots['txty'] = ROOT.TH2F('projected_txty', \
                      'Projected Beam Angle', 1000, -1.0, 1.0, 1000, -1.0, 1.0) 
  projected_plots['theta'] = ROOT.TH1F('projected_theta', \
                                    'Projected Beam Rotation', 1000, -4.0, 4.0) 

  residual_plots = {}

  residual_plots['xy'] = ROOT.TH2F('residual_xy', \
               'Residuals in Position', 100, -400.0, 400.0, 100, -400.0, 400.0)
  residual_plots['mxmy'] = ROOT.TH2F('residual_mxmy', \
                       'Residuals in Gradient', 100, -0.2, 0.2, 100, -0.2, 0.2)
  residual_plots['txty'] = ROOT.TH2F('residual_txty', \
                          'Residuals in Angle', 100, -0.2, 0.2, 100, -0.2, 0.2)
  residual_plots['theta'] = ROOT.TH1F('residual_theta', \
                                     'Residual Beam Rotation', 1000, -8.0, 8.0) 
  residual_plots['xmy'] = ROOT.TH2F('residual_xmy', \
                   'Residual Beam X-My', 1000, -200.0, 200.0, 1000, -0.2, 0.2 )
  residual_plots['ymx'] = ROOT.TH2F('residual_ymx', \
                   'Residual Beam Y-Mx', 1000, -200.0, 200.0, 1000, -0.2, 0.2 )


  global tof_ana
  global tof_scifi_ana
  tof_ana = scifi_analysis.tof_analyser()
  tof_scifi_ana = scifi_analysis.tof_scifi_analyser()

#  plot_dict['tof_correlations'] = tof_correlation_plots
  plot_dict['correlations'] = correlation_plots
  plot_dict['separation'] = separation_plots
  plot_dict['upstream'] = upstream_plots
  plot_dict['downstream'] = downstream_plots
  plot_dict['projected'] = projected_plots
  plot_dict['residual'] = residual_plots

  return plot_dict



def init_data() :
  data_dict = {}

  data_dict['counters'] = {}


  alignment_dict = {}

  alignment_dict['x_delta'] = 0.0
  alignment_dict['x_delta_err'] = 0.0
  alignment_dict['y_delta_err'] = 0.0
  alignment_dict['y_delta'] = 0.0
  alignment_dict['x_delta_phi'] = 0.0
  alignment_dict['x_delta_phi_err'] = 0.0
  alignment_dict['y_delta_phi'] = 0.0
  alignment_dict['y_delta_phi_err'] = 0.0
  alignment_dict['tracker_separation'] = -1.0
  alignment_dict['tracker_separation_err'] = 0.0

  data_dict['alignment'] = alignment_dict

  return data_dict


def set_alignment(data_dict) :
  global DELTA_X
  global DELTA_Y
  global DELTA_PHI_X
  global DELTA_PHI_Y
  global TRACKER_SEPARATION
  DELTA_X = data_dict['alignment']['x_delta']
  DELTA_Y = data_dict['alignment']['y_delta']
  DELTA_PHI_X = data_dict['alignment']['x_delta_phi']
  DELTA_PHI_Y = data_dict['alignment']['y_delta_phi']

  if TRACKER_SEPARATION < 0.0 :
    TRACKER_SEPARATION = data_dict['alignment']['tracker_separation']


def find_straight_tracks(data_dict, scifi_event) :
  """
    Extracts straight tracks from scifi events
  """
  scifi_tracks = scifi_event.scifitracks()
  upstream_tracks = []
  downstream_tracks = []
  for track in scifi_tracks :
    if track.GetAlgorithmUsed() != STRAIGHT_ALGORITHM_ID :
      continue

    if track.tracker() == 0 :
      upstream_tracks.append(track)
    elif track.tracker() == 1 :
      downstream_tracks.append(track)

# Only looking for single track events at present implementation
  track_list = []
  if len( upstream_tracks ) != 1 :
    return track_list
  if len( downstream_tracks) != 1 :
    return track_list

  track_list.append((upstream_tracks[0], downstream_tracks[0]))

  return track_list


def find_tof_spacepoints(data_dict, tof_event) :
  event_spacepoints = tof_event.GetTOFEventSpacePoint()

  tof0_sp_size = event_spacepoints.GetTOF0SpacePointArraySize()
  tof1_sp_size = event_spacepoints.GetTOF1SpacePointArraySize()
  tof2_sp_size = event_spacepoints.GetTOF2SpacePointArraySize()

  tof0_sp = []
  tof1_sp = []
  tof2_sp = []

  for tof0_i in range(tof0_sp_size) : 
    tof0_sp.append( event_spacepoints.GetTOF0SpacePointArrayElement(tof0_i) )
  for tof0_i in range(tof1_sp_size) : 
    tof1_sp.append(event_spacepoints.GetTOF1SpacePointArrayElement(0))
  for tof2_i in range(tof2_sp_size) : 
    tof2_sp.append(event_spacepoints.GetTOF2SpacePointArrayElement(0))

  spacepoints =  (tof0_sp, tof1_sp, tof2_sp)

  return spacepoints


def cut_scifi_event(data_dict, event) :
  """
    Examine scifi event to see if it should be vetoed
  """
  digits = event.digits()
  saturation_counter = 0

  for digit in digits :
    if digit.get_adc() == 255 :
      saturation_counter += 1

  if saturation_counter > 1000 :
    return True

  return False


def cut_tracks(data_dict, up_trk, down_trk) :
  """
    Examine a pair of tracks to see if they should be vetoed
  """
  up_counter = 0
  down_counter = 0

  for tp in up_trk.scifitrackpoints() :
    if tp.has_data() :
      up_counter += 1
    elif scifi.calculate_plane_id( tp.tracker(), tp.station(), tp.plane() ) in IGNORE_PLANES :
      up_counter += 1
  for tp in down_trk.scifitrackpoints() :
    if tp.has_data() :
      down_counter += 1
    elif scifi.calculate_plane_id( tp.tracker(), tp.station(), tp.plane() ) in IGNORE_PLANES :
      down_counter += 1

  if up_counter < MIN_NUMBER_TRACKPOINTS :
    return True
  if down_counter < MIN_NUMBER_TRACKPOINTS :
    return True

  if up_trk.P_value() < P_VALUE_CUT :
    return True
  if down_trk.P_value() < P_VALUE_CUT :
    return True

  up_ref = None
  down_ref = None
  for tp in up_trk.scifitrackpoints() :
    if tp.station() == RECON_STATION and tp.plane() == RECON_PLANE :
      up_ref = tp
  for tp in down_trk.scifitrackpoints() :
    if tp.station() == RECON_STATION and tp.plane() == RECON_PLANE :
      down_ref = tp

  if up_ref is None :
    return True
  if down_ref is None :
    return True


  length = TRACKER_SEPARATION

  up_pos = [ up_ref.pos().x(), up_ref.pos().y() ]
  up_gra = [ up_ref.mom().x() / up_ref.mom().z(), \
                                          up_ref.mom().y() / up_ref.mom().z() ]

  up_pos = up_tracker_correction( up_pos )
  up_gra = up_tracker_correction( up_gra )

  pro_pos = [ up_pos[0] + length*up_gra[0], up_pos[1] + length*up_gra[1] ]

  rad = math.sqrt( up_pos[0]**2 + up_pos[1]**2 )
  grad = math.sqrt( up_gra[0]**2 + up_gra[1]**2 )
  pro_rad = math.sqrt( pro_pos[0]**2 + pro_pos[1]**2 )

  if grad > GRAD_CUT :
    return True
  if grad < LOW_GRAD_CUT :
    return True
  if rad > RADIUS_CUT :
    return True
  if pro_rad > PROJECTED_RADIUS_CUT :
    return True

  return False


def fill_plots(plot_dict, data_dict, up_trk, down_trk, tof_spacepoints=None) :
  """
    Fill the plots in the plot dictionary with useful data
  """
  up_z = None
  down_z = None
  length = 0.0

  up_pos = None
  up_gra = None
  up_ang = None
  up_theta = None
  down_pos = None
  down_gra = None
  down_ang = None
  down_theta = None
  delta_theta = None

  projected_pos = [0.0, 0.0]
  projected_gra = [0.0, 0.0]
  projected_ang = [0.0, 0.0]

  for tp in up_trk.scifitrackpoints() :
    if tp.station() == RECON_STATION and tp.plane() == RECON_PLANE :
      up_z = tp.pos().z()
      up_pos = [ tp.pos().x(), tp.pos().y() ]
      up_ang = [ math.atan2(tp.mom().x(), tp.mom().z()), \
                                    math.atan2(tp.mom().y(), tp.mom().z()) ]
      break
  if up_z is None :
    raise ValueError("Could not find a trackpoint in Tracker 0, " + \
               "Station " + str(RECON_STATION) + ", Plane " + str(RECON_PLANE))

  up_pos = up_tracker_correction( up_pos )
  up_ang = up_tracker_correction( up_ang )

  up_gra = [ math.tan(up_ang[0]), math.tan(up_ang[1]) ]
  up_theta = math.atan2(up_gra[1], up_gra[0])


  for tp in down_trk.scifitrackpoints() :
    if tp.station() == RECON_STATION and tp.plane() == RECON_PLANE :
      down_z = tp.pos().z()
      down_pos = [ tp.pos().x(), tp.pos().y() ]
      down_ang = [ math.atan2(tp.mom().x(), tp.mom().z()), \
                                       math.atan2(tp.mom().y(), tp.mom().z()) ]
      break
  if down_z is None :
    raise ValueError("Could not find a trackpoint in Tracker 1, " + \
               "Station " + str(RECON_STATION) + ", Plane " + str(RECON_PLANE))

  down_pos = down_tracker_correction( down_pos )
  down_ang = down_tracker_correction( down_ang )

  down_pos = [ down_pos[0] + DELTA_X,  down_pos[1] + DELTA_Y ]
  down_ang = [ down_ang[0] + DELTA_PHI_X, down_ang[1] + DELTA_PHI_Y ]

  down_gra = [ math.tan(down_ang[0]), math.tan(down_ang[1]) ]
  down_theta = math.atan2(down_gra[1], down_gra[0])



  if TRACKER_SEPARATION < 0.0 : # Not already determined
    length = down_z - up_z
  else :
    length = TRACKER_SEPARATION

  if down_theta - up_theta > math.pi :
    delta_theta = down_theta - up_theta - 2.0*math.pi
  elif down_theta - up_theta < -math.pi :
    delta_theta = down_theta - up_theta + 2.0*math.pi
  else :
    delta_theta = down_theta - up_theta

  projected_pos[0] = up_pos[0] + length*up_gra[0]
  projected_pos[1] = up_pos[1] + length*up_gra[1]
  projected_gra[0] = up_gra[0]
  projected_gra[1] = up_gra[1]
  projected_ang[0] = up_ang[0]
  projected_ang[1] = up_ang[1]

  plot_dict['separation']['tracker_separation'].Fill(length)

# Calculate the difference in Z position from expected
  plot_dict['separation']['position_gradient_x'].Fill(up_gra[0], down_pos[0]-projected_pos[0])
  plot_dict['separation']['position_gradient_y'].Fill(up_gra[1], down_pos[1]-projected_pos[1])
# Calculate the Z separation
#  plot_dict['separation']['position_gradient_x'].Fill(up_gra[0], down_pos[0]-up_pos[0])
#  plot_dict['separation']['position_gradient_y'].Fill(up_gra[1], down_pos[1]-up_pos[1])

  grad = math.sqrt( up_gra[0]**2 + up_gra[1]**2 )
  delta_pos = [ down_pos[0]-up_pos[0], down_pos[1]-up_pos[1] ]

  delta = (delta_pos[0]*up_gra[0] + delta_pos[1]*up_gra[1]) / grad
  plot_dict['separation']['position_gradient'].Fill( grad, delta )

  plot_dict['correlations']['x_delta_x'].Fill( up_pos[0], down_pos[0]-up_pos[0])
  plot_dict['correlations']['y_delta_y'].Fill( up_pos[1], down_pos[1]-up_pos[1])
  plot_dict['correlations']['x_delta_y'].Fill( up_pos[0], down_pos[1]-up_pos[1])
  plot_dict['correlations']['y_delta_x'].Fill( up_pos[1], down_pos[0]-up_pos[0])

  plot_dict['correlations']['delta_x_residual_x'].Fill( down_pos[0]-up_pos[0], projected_pos[0]-down_pos[0] )
  plot_dict['correlations']['delta_y_residual_y'].Fill( down_pos[1]-up_pos[1], projected_pos[1]-down_pos[1] )
  plot_dict['correlations']['delta_x_residual_y'].Fill( down_pos[0]-up_pos[0], projected_pos[1]-down_pos[1] )
  plot_dict['correlations']['delta_y_residual_x'].Fill( down_pos[1]-up_pos[1], projected_pos[0]-down_pos[0] )

  plot_dict['upstream']['xy'].Fill(up_pos[0], up_pos[1])
  plot_dict['upstream']['mxmy'].Fill(up_gra[0], up_gra[1])
  plot_dict['upstream']['theta'].Fill(up_theta)
  plot_dict['upstream']['xmy'].Fill(up_pos[0], up_gra[1])
  plot_dict['upstream']['ymx'].Fill(up_pos[1], up_gra[0])

  plot_dict['downstream']['xy'].Fill(down_pos[0], down_pos[1])
  plot_dict['downstream']['mxmy'].Fill(down_gra[0], down_gra[1])
  plot_dict['downstream']['theta'].Fill(down_theta)
  plot_dict['downstream']['xmy'].Fill(up_pos[0], up_gra[1])
  plot_dict['downstream']['ymx'].Fill(up_pos[1], up_gra[0])

  plot_dict['projected']['xy'].Fill(projected_pos[0], projected_pos[1])
  plot_dict['projected']['mxmy'].Fill(projected_gra[0], projected_gra[1])
  plot_dict['projected']['txty'].Fill(projected_ang[0], projected_ang[1])
  plot_dict['projected']['theta'].Fill(up_theta)

  plot_dict['residual']['xy'].Fill(projected_pos[0] - down_pos[0], \
                                                projected_pos[1] - down_pos[1])
  plot_dict['residual']['mxmy'].Fill(projected_gra[0] - down_gra[0], \
                                                projected_gra[1] - down_gra[1])
  plot_dict['residual']['txty'].Fill(projected_ang[0] - down_ang[0], \
                                                projected_ang[1] - down_ang[1])
  plot_dict['residual']['theta'].Fill(delta_theta)
  plot_dict['residual']['xmy'].Fill(up_pos[0], up_gra[1])
  plot_dict['residual']['ymx'].Fill(up_pos[1], up_gra[0])


def analyse_plots(plot_dict, data_dict) :
  """
    Use the plot dictionary to perform some analysis of the data
  """

# Fixed Tracker Parameterisation
  plot_dict['x_dist_profile'] = plot_dict['separation']['position_gradient_x'].ProfileX(\
                                                              'x_dist_profile')
  plot_dict['y_dist_profile'] = plot_dict['separation']['position_gradient_y'].ProfileX(\
                                                              'y_dist_profile')

  plot_dict['dist_profile'] = scifi.gaussian_profile_x( plot_dict['separation']['position_gradient'], 0.0, 150.0 )

  x_fit_result = plot_dict['x_dist_profile'].Fit('pol1', 'QS', "", -0.01, 0.01)
  y_fit_result = plot_dict['y_dist_profile'].Fit('pol1', 'QS', "", -0.01, 0.01)
  fit_result = plot_dict['dist_profile'].Fit('pol1', 'QS', "", -0.0, 0.050)

  x_fit_length = x_fit_result.Parameter(1)
  x_fit_length_err = x_fit_result.ParError(1)
  y_fit_length = y_fit_result.Parameter(1)
  y_fit_length_err = y_fit_result.ParError(1)

  fit_length = fit_result.Parameter(1)
  fit_length_err = fit_result.ParError(1)

  length = plot_dict['separation']['tracker_separation'].GetMean()
  num_tracks = plot_dict['separation']['tracker_separation'].GetEntries()

  x_delta_raw = plot_dict['residual']['xy'].GetMean(1)
  y_delta_raw = plot_dict['residual']['xy'].GetMean(2)
  x_delta_raw_err = plot_dict['residual']['xy'].GetRMS(1) / \
                             math.sqrt( num_tracks )
  y_delta_raw_err = plot_dict['residual']['xy'].GetRMS(2) / \
                             math.sqrt( num_tracks )

  x_delta_phi = plot_dict['residual']['txty'].GetMean(1)
  y_delta_phi = plot_dict['residual']['txty'].GetMean(2)
  x_delta_phi_err = plot_dict['residual']['txty'].GetRMS(1) / \
                           math.sqrt( num_tracks )
  y_delta_phi_err = plot_dict['residual']['txty'].GetRMS(2) / \
                           math.sqrt( num_tracks )

  tracker_separation = (x_fit_length + y_fit_length) / 2.0
  tracker_separation_err = math.sqrt(x_fit_length_err**2 \
                                                   + y_fit_length_err**2) / 2.0


  alignment_dict = data_dict['alignment']
  alignment_dict['distance'] = length

  alignment_dict['x_delta'] += x_delta_raw
  alignment_dict['y_delta'] += y_delta_raw
  alignment_dict['x_delta_err'] = x_delta_raw_err
  alignment_dict['y_delta_err'] = y_delta_raw_err

  alignment_dict['x_delta_phi'] += x_delta_phi
  alignment_dict['y_delta_phi'] += y_delta_phi
  alignment_dict['x_delta_phi_err'] = x_delta_phi_err
  alignment_dict['y_delta_phi_err'] = y_delta_phi_err

  alignment_dict['tracker_separation'] = tracker_separation
  alignment_dict['tracker_separation_err'] = tracker_separation_err

  alignment_dict['x_fit_length'] = x_fit_length
  alignment_dict['x_fit_length_err'] = x_fit_length_err
  alignment_dict['y_fit_length'] = y_fit_length
  alignment_dict['y_fit_length_err'] = y_fit_length_err


  print
  print "Downstream Tracker Misalignments:"
  print
  print "Assuming Center of tracker reference plane is stationary."
  print
  print "Distance Between Trackers Set To:", length, "mm"
  print
  print "Analysed {0:0.0f} Tracks".format(num_tracks)
  print
  print "Fitted Length, X  =", "{0:0.0f} +/- {1:0.0f}".format(\
                                 x_fit_length + length, x_fit_length_err), "mm"
  print "Fitted Length, Y  =", "{0:0.0f} +/- {1:0.0f}".format(\
                                 y_fit_length + length, y_fit_length_err), "mm"
  print
  print "Fitted Length     =", "{0:0.0f} +/- {1:0.0f}".format(fit_length, \
                                                          fit_length_err), "mm"
  print
  print "Ref Plane X Shift =", "{0:0.3e} +/- {1:0.3e}".format(x_delta_raw, \
                                                         x_delta_raw_err), "mm"
  print "Ref Plane Y Shift =", "{0:0.3e} +/- {1:0.3e}".format(y_delta_raw, \
                                                         y_delta_raw_err), "mm"
  print
  print "Rotation in X-Z   =", "{0:0.3e} +/- {1:0.3e}".format(x_delta_phi, \
                                                        x_delta_phi_err), "rad"
  print "Rotation in Y-Z   =", "{0:0.3e} +/- {1:0.3e}".format(y_delta_phi, \
                                                        y_delta_phi_err), "rad"

  print
  print "Recalculated Algnments:"
  print
  print "Ref Plane X Shift =", "{0:0.3e} +/- {1:0.3e}".format( \
                alignment_dict['x_delta'], alignment_dict['x_delta_err']), "mm"
  print "Ref Plane Y Shift =", "{0:0.3e} +/- {1:0.3e}".format( \
                alignment_dict['y_delta'], alignment_dict['y_delta_err']), "mm"
  print
  print "Rotation in X-Z   =", "{0:0.3e} +/- {1:0.3e}".format( \
        alignment_dict['x_delta_phi'], alignment_dict['x_delta_phi_err']), "mm"
  print "Rotation in Y-Z   =", "{0:0.3e} +/- {1:0.3e}".format( \
        alignment_dict['y_delta_phi'], alignment_dict['y_delta_phi_err']), "mm"
  print 

  return data_dict


def load_data(filename) :
  """
    Save all the data to a json file
  """
  data_dict = None
  with open(filename, 'r') as infile :
    data_dict = json.load(infile)

# Reset Counters
  for counter in data_dict['counters'] :
    data_dict['counters'][counter] = 0

  return data_dict

 
def save_data(data_dict, directory, filename) :
  """
    Save all the data to a json file
  """
  filename = os.path.join(directory, filename+".json")
  with open(filename, 'w') as outfile :
    json.dump(data_dict, outfile)



if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  parser = argparse.ArgumentParser( description='Performs a straight track '+\
      'reconstruction using the MICE SciFi Trackers in order to measure the '+\
      'gross misalignments between the up- and downstream trackers' )

  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS '+\
      'output root files containing reconstructed straight tracks')

  parser.add_argument( '-N', '--max_num_events', type=int, \
                                   help='Maximum number of events to analyse.')

  parser.add_argument( '-O', '--output_filename', \
            default='straight_track_alignment', help='Set the output filename')

  parser.add_argument( '-D', '--output_directory', \
       default='./', help='Set the output directory')

  parser.add_argument( '-P', '--print_plots', action='store_true', \
                        help="Flag to save the plots as individual pdf files" )

  parser.add_argument( '-T', '--cut_tof', action='store_true', \
                                    help='Flag to cut on Muon Time of Flight' )

  parser.add_argument( '--tof_window', type=float, nargs=2, \
                                                         default=[0.0,100.0], \
                 help='Sepecify the upper and lower bounds of the TOF window' )
  parser.add_argument( '--cut_p_value', type=float, default=0.0, \
                              help='Set the cut on the tracker P-Value [0-1]' )
  parser.add_argument( '--cut_number_trackpoints', type=int, default=0, \
                    help='Set the cut on the number of trackpoints per track' )
  parser.add_argument( '--cut_gradient', type=float, default=1.0, \
                                     help='Set the cut on the track gradient' )
  parser.add_argument( '--cut_min_gradient', type=float, default=0.0, \
                             help='Set the cut on the minimum track gradient' )
  parser.add_argument( '--cut_radius', type=float, default=200.0, \
                              help='Set the cut on the track upstream radius' )
  parser.add_argument( '--cut_projected_radius', type=float, default=200.0, \
                  help='Set the cut on the track downstream projected radius' )

  parser.add_argument( '--alignment_data_file', type=str, default=None, \
                              help='File containing data from a previous run' )
  parser.add_argument( '--tracker_separation', type=float, default=-1.0, \
                      help='Manually set the separation between the trackers' )

#  parser.add_argument( '-C', '--configuration_file', help='Configuration '+\
#      'file for the reconstruction. I need the geometry information' )

  try :
    namespace = parser.parse_args()

    P_VALUE_CUT = namespace.cut_p_value
    MIN_NUMBER_TRACKPOINTS = namespace.cut_number_trackpoints
    GRAD_CUT = math.fabs( namespace.cut_gradient )
    LOW_GRAD_CUT = math.fabs( namespace.cut_min_gradient )
    RADIUS_CUT = namespace.cut_radius
    PROJECTED_RADIUS_CUT = namespace.cut_projected_radius
    TRACKER_SEPARATION = namespace.tracker_separation

    if namespace.cut_tof :
      TOF_CUT_LOW = namespace.tof_window[0]
      TOF_CUT_HIGH = namespace.tof_window[1]

  except :
    raise
  else :
##### 1. Load MAUS globals and geometry.
    # geom = load_tracker_geometry(namespace.configuration_file)

##### 2. Intialise plots ######################################################
    print "\nInitialising..."
    plot_dict = init_plots()
    data_dict = init_data()
    if namespace.alignment_data_file :
      print "\nLoading Alignment File:", namespace.alignment_data_file
      data_dict = load_data(namespace.alignment_data_file)
      set_alignment(data_dict)
    

##### 3. Load SciFi Events ####################################################
    print "\nLoading Spills...\n"
    file_reader = event_loader.maus_reader(namespace.maus_root_files)

    try :
      while file_reader.next_event() and \
               file_reader.get_total_num_events() != namespace.max_num_events :
        try :
          sys.stdout.write( 
              ' Spill ' + str(file_reader.get_current_spill_number()) + \
              ' of ' + str(file_reader.get_current_number_spills()) + \
              ' in File ' + str(file_reader.get_current_filenumber()) + \
              ' of ' + str(file_reader.get_number_files()) + '             \r')

          sys.stdout.flush()

          scifi_event = file_reader.get_event( 'scifi' )
          tof_event = file_reader.get_event( 'tof' )

##### 4. Extract potential tracks #############################################
          straights = find_straight_tracks(data_dict, scifi_event)
          tof_spacepoints = find_tof_spacepoints(data_dict, tof_event)

          if cut_scifi_event(data_dict, scifi_event) :
            continue

          tof_ana.analyse(file_reader)
          tof_scifi_ana.analyse(file_reader)
          if namespace.cut_tof and ( tof_ana.is_cut() or tof_scifi_ana.is_cut() ) :
            continue
          
          for up_str, down_str in straights :

##### 5. Apply Cuts ###########################################################
            if cut_tracks(data_dict, up_str, down_str) :
              continue

##### 6. Fill plots ###########################################################
            else :
              fill_plots(plot_dict, data_dict, up_str, down_str, tof_spacepoints)

        except ValueError as ex:
          print "An Error Occured. Skipping Spill: " + \
                str(file_reader.get_current_spill_number()) + \
                " In File: " + str(file_reader.get_current_filenumber())
          print "ERROR =", ex
          print
          continue

##### 7. Analysis Plots #######################################################
    except KeyboardInterrupt :
      print
      print "Keyboard Interrupt"
      print
    print "All Spills Loaded                                                  "
    print "\nStarting Analysis"
    out_file_name = os.path.join( namespace.output_directory, \
                                            namespace.output_filename+".root" )

    try :
      plot_dict = tof_ana.get_plot_dict(plot_dict)
      plot_dict = tof_scifi_ana.get_plot_dict(plot_dict)
      analyse_plots(plot_dict, data_dict)
    except ValueError as ex :
      print "Analysis Failed:", ex
      print
      print "Stopping Execution"
      scifi.save_plots(plot_dict, out_file_name)

##### 8. Save plots and data ##################################################
    print "\nSaving Plots and Data"
    scifi.save_plots(plot_dict, out_file_name)
    if namespace.print_plots :
      scifi.print_plots(plot_dict, namespace.output_directory)

    save_data(data_dict, namespace.output_directory, namespace.output_filename)

  print 
  print "Complete."
  print

