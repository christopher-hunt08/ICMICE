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
import scifi_extractors
import emr_analysis


"""
  Classes that compare the EMR to the Scifi Recon are stored here.
"""

class emr_scifi_correlations(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "emr_scifi_correlations")

    self.__emr_z_pos = None

    self.__correlation_xx = ROOT.TH2F("scifi_x_emr_x", \
                                     "Residual In EMR From SciFI Projection", \
                                      1000, -500.0, 500.0, 1000, -500.0, 500.0)
    self.__correlation_yy = ROOT.TH2F("scifi_y_emr_y", \
                                     "Residual In EMR From SciFI Projection", \
                                      1000, -500.0, 500.0, 1000, -500.0, 500.0)
    self.__correlation_xy = ROOT.TH2F("scifi_x_emr_y", \
                                     "Residual In EMR From SciFI Projection", \
                                      1000, -500.0, 500.0, 1000, -500.0, 500.0)
    self.__correlation_yx = ROOT.TH2F("scifi_y_emr_x", \
                                     "Residual In EMR From SciFI Projection", \
                                      1000, -500.0, 500.0, 1000, -500.0, 500.0)
    self.__count_mismatch = 0


  def get_args(self, parser) :
    parser.add_argument('--emr_z_pos', type=float, default=None,
                                          help="Manually set the EMR Position")


  def process_args(self, namespace) :
    self.__emr_z_pos = namespace.emr_z_pos


  def get_dependencies(self, inserter) :
    self.__emr_entry = inserter(emr_analysis.emr_entry_positions())
    self.__scifi_tracks = inserter(scifi_extractors.scifi_straight_track_extractor())


  def _reset(self) :
    self.__count_mismatch = 0


  def _process(self, file_reader) :
    emr_tracks = self.__emr_entry.get_tracks()
    downstream_tracks = self.__scifi_tracks.get_downstream_tracks()

    if len(downstream_tracks) != 1 or len(emr_tracks) != 1 : 
      self.__count_mismatch += 1
      return True
    else :
      emr_track = emr_tracks[0]
      emr_origin = emr_track.GetOrigin()

      for tp in downstream_tracks[0].scifitrackpoints() :
        if tp.station() == 5 and tp.plane() == 2 :
          trackpoint = tp
          break

      if self.__emr_z_pos is not None :
        separation = self.__emr_z_pos - trackpoint.pos().z()
      else :
        separation = emr_origin.GetZ() - trackpoint.pos().z()
    
      grad_x = trackpoint.mom().x() / trackpoint.mom().z()
      grad_y = trackpoint.mom().y() / trackpoint.mom().z()

      projected_x = trackpoint.pos().x() + grad_x*separation
      projected_y = trackpoint.pos().y() + grad_y*separation

      self.__correlation_xx.Fill(projected_x, emr_origin.GetX())
      self.__correlation_yy.Fill(projected_y, emr_origin.GetY())
      self.__correlation_xy.Fill(projected_x, emr_origin.GetY())
      self.__correlation_yx.Fill(projected_y, emr_origin.GetX())
      return False


  def _store_data(self, data_dict) :
    emr_scifi_correlation = {}
    emr_scifi_correlation['mismatched'] = self.__count_mismatch

    data_dict['emr_scifi_correlation'] = emr_scifi_correlation
    return data_dict


  def _store_plots(self, plot_dict) :
    emr_scifi_correlation_plots = {}
    emr_scifi_correlation_plots['scifi_x_emr_x'] = self.__correlation_xx
    emr_scifi_correlation_plots['scifi_y_emr_y'] = self.__correlation_yy
    emr_scifi_correlation_plots['scifi_x_emr_y'] = self.__correlation_xy
    emr_scifi_correlation_plots['scifi_y_emr_x'] = self.__correlation_yx

    plot_dict['emr_scifi_correlation'] = emr_scifi_correlation_plots
    return plot_dict


