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
   - Plot overall distribution of track pairs to reconstruct the transported 
     beam centroid.
   - Use the reconstructed beam centroid parameters to calculate the tracker
     (and beam) alignments

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
import types
import array

# Third Party library import statements
import event_loader
import scifi_tools as scifi
import ROOT
import json



# Useful Constants and configuration
STRAIGHT_ALGORITHM_ID = 0
STARTING_TRACKER = 0
RECON_STATION = 1
RECON_PLANE = 0
MIN_NUMBER_TRACKPOINTS = 5
TRACKER_LENGTH = 1101.0624
TOF_CUT_LOW = 0.0
TOF_CUT_HIGH = 100.0
GRADIENT_CUT = 1.05
RADIUS_CUT = 200.0
PROJECTED_RADIUS_CUT = 200.0
P_VALUE_CUT = 0.00

IGNORE_STATIONS = [ -3 ]

UPSTREAM_REF_POSITION = -1350.0
DOWNSTREAM_REF_POSITION = 1350.0

# Alignment Info
DELTA_X = 0.0
DELTA_Y = 0.0
DELTA_PHI_X = 0.0
DELTA_PHI_Y = 0.0
TRACKER_SEPARATION = 3500.0

UP_DELTA_THETA = 2.0 * math.pi / 3.0
DOWN_DELTA_THETA = -2.0 * math.pi / 3.0
#UP_DELTA_THETA = 0.0
#DOWN_DELTA_THETA = 0.0

# Plot Options
PLOT_OPTIONS = {}

def rotate_by_theta( vector, theta ) :
  new_vector = [ 0.0, 0.0 ]
  cos = math.cos( theta )
  sin = math.sin( theta )

  new_vector[0] = vector[0] * cos - vector[1] * sin
  new_vector[1] = vector[0] * sin + vector[1] * cos

  return new_vector


def up_tracker_correction( vector ) :
  new_vec = rotate_by_theta( vector, UP_DELTA_THETA )
  new_vec[0] = -new_vec[0]
  new_vec[1] = -new_vec[1]
  return new_vec
  

def down_tracker_correction( vector ) :
  new_vec = rotate_by_theta( vector, DOWN_DELTA_THETA )
  new_vec[0] = -new_vec[0]
  new_vec[1] = -new_vec[1]
  return new_vec



def init_plots() :
  """
    Initialises a dictionary of plots for later use
  """
  plot_dict = {}

  station_dict = {}

  for st_id in [ -5, -4, -3, -2, -1, 1, 2, 3, 4, 5 ] :
    prefix = 'station_' + str( st_id ) + '_'
    station_dict[prefix+'spacepoints_xy'] = \
              ROOT.TH2D( prefix+'spacepoints_xy', "Spacepoint X-Y Positions", \
                                      1000, -200.0, 200.0, 1000, 200.0, 200.0 )

  plot_dict['station_plots'] = station_dict


  plot_dict['beam_positions_x'] = ROOT.TH2D( 'beam_positions_x', \
                              "Distribution of X Positions for each station", \
                                         11, -5.5, 5.5, 1000, -200.0, 200.0 )
  plot_dict['beam_positions_y'] = ROOT.TH2D( 'beam_positions_y', \
                              "Distribution of Y Positions for each station", \
                                         11, -5.5, 5.5, 1000, -200.0, 200.0 )
  plot_dict['beam_profile_x'] = None
  plot_dict['beam_profile_y'] = None
  plot_dict['beam_profile_x_up_fit'] = None
  plot_dict['beam_profile_y_up_fit'] = None
  plot_dict['beam_profile_x_down_fit'] = None
  plot_dict['beam_profile_y_down_fit'] = None

  plot_dict['tof_0_1'] = ROOT.TH1F( 'tof_0_1', 'Time TOF0 - TOF1', \
                                                             1000, 0.0, 100.0 )
  plot_dict['tof_1_2'] = ROOT.TH1F( 'tof_1_2', 'Time TOF1 - TOF2', \
                                                             1000, 0.0, 100.0 )
  plot_dict['tof_0_1_cut'] = ROOT.TH1F( 'tof_0_1_cut', 'Time TOF0 - TOF1', \
                                                             1000, 0.0, 100.0 )
  plot_dict['tof_1_2_cut'] = ROOT.TH1F( 'tof_1_2_cut', 'Time TOF1 - TOF2', \
                                                             1000, 0.0, 100.0 )

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

  data_dict['alignment'] = alignment_dict

  beam_dict = {}

  beam_dict['x_delta'] = 0.0
  beam_dict['x_delta_err'] = 0.0
  beam_dict['y_delta'] = 0.0
  beam_dict['y_delta_err'] = 0.0
  beam_dict['x_delta_phi'] = 0.0
  beam_dict['x_delta_phi_err'] = 0.0
  beam_dict['y_delta_phi'] = 0.0
  beam_dict['y_delta_phi_err'] = 0.0

  data_dict['beam_alignment'] = beam_dict

  positions_dict = {}

  positions_dict[-5] = None
  positions_dict[-4] = None
  positions_dict[-3] = None
  positions_dict[-2] = None
  positions_dict[-1] = None
  positions_dict[5] = None
  positions_dict[4] = None
  positions_dict[3] = None
  positions_dict[2] = None
  positions_dict[1] = None

  data_dict['station_positions'] = positions_dict


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


def find_station_positions(data_dict, scifi_event) :
  spacepoints = scifi_event.spacepoints()

  for spacepoint in spacepoints :
    tracker = spacepoint.get_tracker()
    station = spacepoint.get_station()

    st_id = station * ( -1.0 if tracker == 0 else 1.0 )
    position = spacepoint.get_position().z()
    if tracker == 0 :
      position = -position + UPSTREAM_REF_POSITION
    else :
      position += DOWNSTREAM_REF_POSITION

    data_dict['station_positions'][st_id] = position

  for station in data_dict['station_positions'] :
    if data_dict['station_positions'][station] is None :
      return False

  return True


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

def cut_tof_event(data_dict, plot_dict, event) :
  """
    Examine TOF event to see if it should be vetoed
  """
  event_spacepoints = event.GetTOFEventSpacePoint()

  tof0_sp_size = event_spacepoints.GetTOF0SpacePointArraySize()
  tof1_sp_size = event_spacepoints.GetTOF1SpacePointArraySize()
  tof2_sp_size = event_spacepoints.GetTOF2SpacePointArraySize()

  if tof0_sp_size < 1 or tof1_sp_size < 1 or tof2_sp_size < 1 :
    return True

  tof0_sp = event_spacepoints.GetTOF0SpacePointArrayElement(0)
  tof1_sp = event_spacepoints.GetTOF1SpacePointArrayElement(0)
  tof2_sp = event_spacepoints.GetTOF2SpacePointArrayElement(0)

  if tof1_sp_size != 1 or tof2_sp_size != 1 :
    return True

  diff_0_1 = tof1_sp.GetTime() - tof0_sp.GetTime()
  diff_1_2 = tof2_sp.GetTime() - tof1_sp.GetTime()

  plot_dict['tof_0_1'].Fill( diff_0_1 )
  plot_dict['tof_1_2'].Fill( diff_1_2 )

  if diff_1_2 < TOF_CUT_LOW or diff_1_2 > TOF_CUT_HIGH :
    return True

  plot_dict['tof_0_1_cut'].Fill( tof1_sp.GetTime() - tof0_sp.GetTime() )
  plot_dict['tof_1_2_cut'].Fill( tof2_sp.GetTime() - tof1_sp.GetTime() )

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
  for tp in down_trk.scifitrackpoints() :
    if tp.has_data() :
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

  if grad > GRADIENT_CUT :
    return True
  if rad > RADIUS_CUT :
    return True
  if pro_rad > PROJECTED_RADIUS_CUT :
    return True

  return False


def fill_plots(plot_dict, data_dict, up_trk, down_trk) :
  """
    Fill the plots in the plot dictionary with useful data
  """
  up_patrec = None
  down_patrec = None
  up_ref_pos = 0.0
  down_ref_pos = 0.0

  for trackpoint in up_trk.scifitrackpoints() :
    if trackpoint.station() == 1 and trackpoint.plane() == 0 :
      up_ref_pos = trackpoint.pos().z()

  for trackpoint in down_trk.scifitrackpoints() :
    if trackpoint.station() == 1 and trackpoint.plane() == 0 :
      down_ref_pos = trackpoint.pos().z()


  if up_trk.GetAlgorithmUsed() == 1 :
    up_patrec = up_trk.pr_track_pointer_straight()
  else :
    up_patrec = up_trk.pr_track_pointer_helical()

  if down_trk.GetAlgorithmUsed() == 1 :
    down_patrec = down_trk.pr_track_pointer_straight()
  else :
    down_patrec = down_trk.pr_track_pointer_helical()

  up_spacepoints = up_patrec.get_spacepoints_pointers()
  down_spacepoints = down_patrec.get_spacepoints_pointers()

  for sp_list in [ up_spacepoints, down_spacepoints ] :
    for spacepoint in sp_list :
      station = spacepoint.get_station()
      pos = spacepoint.get_position()
      if spacepoint.get_tracker() == 0 :
        st_id = station * -1
        pos.setX( pos.x() * -1.0 )
        pos.setZ( up_ref_pos - pos.z() )

        pos = up_tracker_correction( [ pos.x(), pos.y() ] )
      else :
        st_id = station
        pos.setZ( up_ref_pos + pos.z() )

        pos = down_tracker_correction( [ pos.x(), pos.y() ] )

      if st_id in IGNORE_STATIONS :
        continue

      plot_name = 'station_' + str( st_id ) + '_spacepoints_xy'
      plot_dict['station_plots'][plot_name].Fill( pos[0], pos[1] )

      st_id = float(st_id)

      plot_dict['beam_positions_x'].Fill( st_id, pos[0] )
      plot_dict['beam_positions_y'].Fill( st_id, pos[1] )

      if st_id < 0.499999999 and st_id > -0.4999999999 : 
        print st_id


def analyse_plots(plot_dict, data_dict) :
  """
    Use the plot dictionary to perform some analysis of the data
  """
  for component in [ '_x', '_y' ] :
    z_pos = array.array( 'd' )
    trans_pos = array.array( 'd' )
    errors = array.array( 'd' )
    zeros = array.array( 'd' )

    plot = plot_dict['beam_positions'+component]

    for i in range( plot.GetXaxis().GetNbins()+2 ) :
      projection = plot.ProjectionY( \
                     'profile'+component+'_pro_'+str(i), i, i )
      if projection.GetEntries() == 0 :
        continue

      pro_mean, pro_mean_err, pro_std, pro_std_err = \
                                               scifi.fit_gaussian( projection )

      errors.append( pro_mean_err )
      trans_pos.append( pro_mean )
      z_pos.append( data_dict['station_positions'][ i-6 ] )
      zeros.append(0.0)

    position_graph = ROOT.TGraphErrors( len(zeros), z_pos, trans_pos, \
                                                                zeros, errors )
    position_graph.SetName('beam_profile'+component)
    plot_dict['beam_profile'+component] = position_graph

  profile_x = plot_dict['beam_profile_x']
  profile_y = plot_dict['beam_profile_y']

  up_x_func = ROOT.TF1( "up_fit_x", "pol1", -5000.0, 0.0 )
  up_y_func = ROOT.TF1( "up_fit_y", "pol1", -5000.0, 0.0 )
  down_x_func = ROOT.TF1( "down_fit_x", "pol1", 0.0, 5000.0 )
  down_y_func = ROOT.TF1( "down_fit_y", "pol1", 0.0, 5000.0 )

  up_fit_x = profile_x.Fit( 'up_fit_x', "QSR" )
  up_fit_y = profile_y.Fit( 'up_fit_y', "QSR" )
  down_fit_x = profile_x.Fit( 'down_fit_x', "QSR" )
  down_fit_y = profile_y.Fit( 'down_fit_y', "QSR" )

  plot_dict['beam_profile_x_up_fit'] = up_x_func
  plot_dict['beam_profile_y_up_fit'] = up_y_func
  plot_dict['beam_profile_x_down_fit'] = down_x_func
  plot_dict['beam_profile_y_down_fit'] = down_y_func


  up_beam_gra_x = up_x_func.GetParameter(1)
  up_beam_gra_x_err = up_x_func.GetParError(1)
  up_beam_gra_y = up_y_func.GetParameter(1)
  up_beam_gra_y_err = up_y_func.GetParError(1)

  up_beam_pos_x = data_dict['station_positions'][-1]*up_beam_gra_x + up_x_func.GetParameter(0)
  up_beam_pos_x_err = up_x_func.GetParError(0)
  up_beam_pos_y = data_dict['station_positions'][-1]*up_beam_gra_y + up_y_func.GetParameter(0)
  up_beam_pos_y_err = up_y_func.GetParError(0)

  up_beam_rot_x = math.atan( up_beam_gra_x )
  up_beam_rot_x_err = up_beam_gra_x_err # Approx linear
  up_beam_rot_y = math.atan( up_beam_gra_y )
  up_beam_rot_y_err = up_beam_gra_y_err # Approx linear



  down_beam_gra_x = down_x_func.GetParameter(1)
  down_beam_gra_x_err = down_x_func.GetParError(1)
  down_beam_gra_y = down_y_func.GetParameter(1)
  down_beam_gra_y_err = down_y_func.GetParError(1)

  down_beam_pos_x = data_dict['station_positions'][1]*down_beam_gra_x + down_x_func.GetParameter(0)
  down_beam_pos_x_err = down_x_func.GetParError(0)
  down_beam_pos_y = data_dict['station_positions'][1]*down_beam_gra_y + down_y_func.GetParameter(0)
  down_beam_pos_y_err = down_y_func.GetParError(0)

  down_beam_rot_x = math.atan( down_beam_gra_x )
  down_beam_rot_x_err = down_beam_gra_x_err # Approx linear
  down_beam_rot_y = math.atan( down_beam_gra_y )
  down_beam_rot_y_err = down_beam_gra_y_err # Approx linear


#  down_pos_x = down_beam_pos_x - data_dict['station_positions'][1]*up_beam_gra_x + up_x_func.GetParameter(0)
#  down_pos_x_err = math.sqrt( up_x_func.GetParError(0)**2 + down_beam_pos_x_err**2 )
#  down_pos_y = down_beam_pos_y - data_dict['station_positions'][1]*up_beam_gra_y + up_y_func.GetParameter(0)
#  down_pos_y_err = math.sqrt( up_y_func.GetParError(0)**2 + down_beam_pos_y_err**2 )

  length = TRACKER_SEPARATION
  down_pos_x = down_beam_pos_x - ( up_beam_pos_x + length*up_beam_gra_x )
  down_pos_x_err = math.sqrt( up_beam_pos_x_err**2 + down_beam_pos_x_err**2 + (length*up_beam_gra_x_err)**2 )
  down_pos_y = down_beam_pos_y - ( up_beam_pos_y + length*up_beam_gra_y )
  down_pos_y_err = math.sqrt( up_beam_pos_y_err**2 + down_beam_pos_y_err**2 + (length*up_beam_gra_y_err)**2 )

  down_rot_x = down_beam_rot_x - up_beam_rot_x
  down_rot_x_err = math.sqrt( down_beam_rot_x_err**2 + up_beam_rot_x_err**2 )
  down_rot_y = down_beam_rot_y - up_beam_rot_y
  down_rot_y_err = math.sqrt( down_beam_rot_y_err**2 + up_beam_rot_y_err**2 )


  print
  print "Incoming Beam Misalignments:"
  print
  print "Displacement and rotation of beam with respect to upstream tracker:"
  print
  print "X Position       =  {0:0.3f} +/- {1:0.3f} mm".format( up_beam_pos_x, up_beam_pos_x_err )
  print "Y Position       =  {0:0.3f} +/- {1:0.3f} mm".format( up_beam_pos_y, up_beam_pos_y_err )
  print
  print "X Rotation       =  {0:0.3f} +/- {1:0.3f} mrad".format( up_beam_rot_x*1000.0, up_beam_rot_x_err*1000.0 )
  print "Y Rotation       =  {0:0.3f} +/- {1:0.3f} mrad".format( up_beam_rot_y*1000.0, up_beam_rot_y_err*1000.0 )
  print

  print
  print "Downstream Tracker Beam Misalignments:"
  print
  print "Displacement and rotation of beam with respect to downstream tracker:"
  print
  print "X Position       =  {0:0.3f} +/- {1:0.3f} mm".format( down_beam_pos_x, down_beam_pos_x_err )
  print "Y Position       =  {0:0.3f} +/- {1:0.3f} mm".format( down_beam_pos_y, down_beam_pos_y_err )
  print
  print "X Rotation       =  {0:0.3f} +/- {1:0.3f} mrad".format( down_beam_rot_x*1000.0, down_beam_rot_x_err*1000.0 )
  print "Y Rotation       =  {0:0.3f} +/- {1:0.3f} mrad".format( down_beam_rot_y*1000.0, down_beam_rot_y_err*1000.0 )
  print

  print
  print "Downstream Tracker Alignment:"
  print
  print "Displacement and rotation of between the two trackers:"
  print
  print "X Position       =  {0:0.3f} +/- {1:0.3f} mm".format( down_pos_x, down_pos_x_err )
  print "Y Position       =  {0:0.3f} +/- {1:0.3f} mm".format( down_pos_y, down_pos_y_err )
  print
  print "X Rotation       =  {0:0.3f} +/- {1:0.3f} mrad".format( down_rot_x*1000.0, down_rot_x_err*1000.0 )
  print "Y Rotation       =  {0:0.3f} +/- {1:0.3f} mrad".format( down_rot_y*1000.0, down_rot_y_err*1000.0 )
  print



#  print
#  print "Downstream Tracker Misalignments:"
#  print
#  print "Assuming Center of tracker reference plane is stationary."
#  print
#  print "Distance Between Trackers Set To:", length, "mm"
#  print
#  print "Analysed {0:0.0f} Tracks".format(num_tracks)
#  print
#  print "Fitted Length, X  =", "{0:0.0f} +/- {1:0.0f}".format(x_fit_length, \
#                                                        x_fit_length_err), "mm"
#  print "Fitted Length, Y  =", "{0:0.0f} +/- {1:0.0f}".format(y_fit_length, \
#                                                        y_fit_length_err), "mm"
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


def save_plots(plot_dict, filename) :
  """
    Save all the plots to file. Assumes a directory like strucutre of plots to
    recursively save them.
  """
  outfile = ROOT.TFile(filename, "RECREATE")

  for key in sorted(plot_dict) :
    if type( plot_dict[key] ) is types.DictType :
      directory = outfile.mkdir( key )
      directory.cd()
      save_plot( plot_dict[key], directory )
      outfile.cd()
    elif plot_dict[key] is not None :
      plot_dict[key].Write()
      
  outfile.Close()


def save_plot(plot_dict, outfile) :
  """
    The recursive saving function for the plot dictionary
  """
  for key in sorted(plot_dict) :
    if type( plot_dict[key] ) is types.DictType :
      directory = outfile.mkdir( key )
      directory.cd()
      save_plot( plot_dict[key], directory )
      outfile.cd()
    else :
      plot_dict[key].Write()


def print_plots(plot_dict, location) :
  """
    Print the plots to PDF files rather than a root file. Assumes a recursive
    structure of directories and plots to save then to file
  """
  if not os.path.exists( location ) :
    os.makedirs( location )

  for key in sorted(plot_dict) :
    if type( plot_dict[key] ) is types.DictType :
      new_location = os.path.join( location, key )
      print_plots( plot_dict[key], new_location )
      
    else :
      canvas = ROOT.TCanvas( key+'_canvas' )
      plot_dict[key].Draw()
      if key in PLOT_OPTIONS :
        apply_options( canvas, plot_dict[key], PLOT_OPTIONS[key] )
      canvas.SaveAs( os.path.join(location, key ) + ".pdf", "pdf" )


 
def save_data(data_dict, directory, filename) :
  """
    Save all the data to a json file
  """
  filename = os.path.join(directory, filename+".json")
  with open(filename, 'w') as outfile :
    json.dump(data_dict, outfile)



if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )

  parser = argparse.ArgumentParser( description='Performs a straight track '+\
      'reconstruction using the MICE SciFi Trackers in order to measure the '+\
      'gross misalignments between the up- and downstream trackers' )

  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS '+\
      'output root files containing reconstructed straight tracks')

  parser.add_argument( '-N', '--max_num_events', type=int, \
                                   help='Maximum number of events to analyse.')

  parser.add_argument( '-O', '--output_filename', \
       default='beam_centroid_track_alignment', help='Set the output filename')

  parser.add_argument( '-D', '--output_directory', \
       default='./', help='Set the output directory')

  parser.add_argument( '-P', '--print_plots', action='store_true', \
                        help="Flag to save the plots as individual pdf files" )

  parser.add_argument( '-T', '--cut_tof', action='store_true', \
                                    help='Flag to cut on Muon Time of Flight' )
  parser.add_argument( '--tof_window', type=float, nargs=2, \
                                                         default=[0.0, 100.0],\
                 help='Sepecify the upper and lower bounds of the TOF window' )

  parser.add_argument( '--cut_p_value', type=float, default=0.0, \
                              help='Set the cut on the tracker P-Value [0-1]' )
  parser.add_argument( '--cut_number_trackpoints', type=int, default=0, \
                    help='Set the cut on the number of trackpoints per track' )
  parser.add_argument( '--cut_gradient', type=float, default=1.0, \
                                   help='Cut tracks above a certain gradient' )
  parser.add_argument( '--cut_radius', type=float, default=200.0, \
                                     help='Cut tracks above a certain radius' )
  parser.add_argument( '--cut_projected_radius', type=float, default=200.0, \
                help='Cut tracks above a certain radius projected downstream' )


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
    TRACKER_SEPARATION = namespace.tracker_separation
    GRADIENT_CUT = namespace.cut_gradient
    RADIUS_CUT = namespace.cut_radius
    PROJECTED_RADIUS_CUT = namespace.cut_projected_radius

    if namespace.cut_tof :
      TOF_CUT_LOW = namespace.tof_window[0]
      TOF_CUT_HIGH = namespace.tof_window[1]

  except :
    raise
  else :

##### 1. Intialise plots ######################################################
    print "\nInitialising..."
    plot_dict = init_plots()
    data_dict = init_data()
    if namespace.alignment_data_file :
      print "\nLoading Alignment File:", namespace.alignment_data_file
      data_dict = load_data(namespace.alignment_data_file)
      set_alignment(data_dict)


##### 2. Load MAUS globals and geometry #######################################
    print "\nLocating all Tracker Stations..."

    file_reader = event_loader.maus_reader(namespace.maus_root_files)

    while file_reader.next_event() :
      scifi_event = file_reader.get_event( 'scifi' )
      if find_station_positions(data_dict, scifi_event) :
        break
    else :
        raise ValueError("Could not location all stations")
    
    file_reader.reset()

##### 3. Load SciFi Events ####################################################
    print "\nLoading Spills...\n"

    try :
      while file_reader.next_event() and \
               file_reader.get_total_num_events() != namespace.max_num_events :
        try :
          sys.stdout.write( 
              '  Spill ' + str(file_reader.get_current_spill_number()) + \
              ' of ' + str(file_reader.get_current_number_spills()) + \
              ' in File ' + str(file_reader.get_current_filenumber()) + \
              ' of ' + str(file_reader.get_number_files()) + '             \r')
          sys.stdout.flush()

          scifi_event = file_reader.get_event( 'scifi' )
          tof_event = file_reader.get_event( 'tof' )

##### 4. Extract potential tracks #############################################
          straights = find_straight_tracks(data_dict, scifi_event)

          if cut_scifi_event(data_dict, scifi_event) :
            continue

          if namespace.cut_tof and \
                               cut_tof_event(data_dict, plot_dict, tof_event) :
            continue
          
          for up_str, down_str in straights :

##### 5. Apply Cuts ###########################################################
            if cut_tracks(data_dict, up_str, down_str) :
              continue

##### 6. Fill plots ###########################################################
            else :
              fill_plots(plot_dict, data_dict, up_str, down_str)

        except ValueError :
          print "An Error Occured. Skipping Spill: " + \
                str(file_reader.get_current_spill_number()) + \
                " In File: " + str(file_reader.get_current_filenumber()) + "\n"
          continue

##### 7. Analysis Plots #######################################################
    except KeyboardInterrupt :
      print
      print "Keyboard Interrupt"
      print
    print "All Spills Loaded                                                  "
    print "\nStarting Analysis"
    analyse_plots(plot_dict, data_dict)

##### 8. Save plots and data ##################################################
    print "\nSaving Plots and Data"

    outdir = namespace.output_directory
    if not os.path.exists( outdir ) :
      os.makedirs( outdir )
    outfile = os.path.join( outdir, namespace.output_filename+".root" )

    save_plots(plot_dict, outfile)

    if namespace.print_plots :
      print_plots(plot_dict, outdir)

    save_data(data_dict, namespace.output_directory, namespace.output_filename)


  print 
  print "Complete."
  print

