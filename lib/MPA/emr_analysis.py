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

# pylint: disable = W0311, E1101, W0613, C0111, R0911, W0621, C0103, R0902

import ROOT
import os
import math

import framework
import analysis.tools as tools


"""
  EMR Analysis Classes are stored here.
"""

class emr_entry_positions(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "emr_entry_positions")

    self.__emr_tracks = []
    self.__emr_trackpoints = []
    self.__emr_spacepoints = [] 
    self.__count_hits = 0
    self.__count_primaries = 0

    self.__entry_positions = ROOT.TH2F( 'emr_entry_positions', \
               "Track Entry To EMR", 1000, -500.0, 500.0, 1000, -500.0, 500.0 )


  def get_trackpoints(self) :
    return self.__emr_trackpoints


  def get_spacepoints(self) :
    return self.__emr_spacepoints


  def get_tracks(self) :
    return self.__emr_tracks


  def get_dependencies(self, inserter) :
    pass


  def _reset(self) :
    self.__emr_tracks = []
    self.__emr_spacepoints = [] 


  def _process(self, file_reader) :
    emr_event = file_reader.get_event('emr')

    emr_tracks = emr_event.GetEMREventTrackArray()
    self.__count_primaries += len(emr_tracks)

    for track_event in emr_tracks :
      track = track_event.GetEMRTrack()
      origin = track.GetOrigin()

      for trackpoint in track.GetEMRTrackPointArray() :
        self.__emr_trackpoints.append(trackpoint)

      self.__emr_tracks.append(track)
      self.__entry_positions.Fill(origin.x(), origin.y())
    
#    plane_hits = emr_event.GetEMRPlaneHitArray()
#    for plane_hit_i in range(len(plane_hits)) :
#      plane_hit = plane_hits[plane_hit_i]
#      if plane_hit.GetPlane() != 0 :
#        continue
#      bars = plane_hit.GetEMRBarArrayPrimary()
#      self.__count_primaries += len(bars)
#
#      for bars_i in range(len(bars)) :
#        bar_hits = bars[bars_i].GetEMRBarHitArray()
#
#        for hits_i in range(len(bar_hits)) :
#          hit = bar_hits[hits_i]
#
#          self.__emr_hits.append(hit)
#          self.__entry_positions.Fill(hit.GetX(), hit.GetY())

#          print plane_hit.GetPlane(), bars[bars_i].GetBar(), " | ", hit.GetX(), hit.GetErrorX(), hit.GetY(), hit.GetErrorY(), hit.GetZ(), hit.GetErrorZ()


  def _store_data(self, data_dict) :
    emr_entry = {}
    emr_entry['primaries'] = self.__count_primaries

    data_dict['emr_entry'] = emr_entry
    return data_dict


  def _store_plots(self, plot_dict) :
    emr_entry_plots = {}
    emr_entry_plots['entry_positions'] = self.__entry_positions

    plot_dict['emr_entry'] = emr_entry_plots
    return plot_dict


