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
import tof_analysis


class tof_scifi_analyser(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "tof_scifi_correlations")

    self.__tof_residual_upstream = ROOT.TH2F('upstream_tof_residual', \
        "Projected Residuals at TOF1", 200, -200.0, 200.0, 200, -200.0, 200.0 )

    self.__tof_residual_downstream = ROOT.TH2F('downstream_tof_residual', \
        "Projected Residuals at TOF2", 200, -200.0, 200.0, 200, -200.0, 200.0 )

    self.__tof_corr_upstream_x = ROOT.TH2F('upstream_tof_x', \
           "Upstream X Against TOF1 X", 400, -400.0, 400.0, 7, -210.0, 210.0 )
    self.__tof_corr_upstream_y = ROOT.TH2F('upstream_tof_y', \
           "Upstream Y Against TOF1 Y", 400, -400.0, 400.0, 7, -210.0, 210.0 )

    self.__tof_corr_upstream_x_y = ROOT.TH2F('upstream_x_tof_y', \
           "Upstream X Against TOF1 Y", 400, -400.0, 400.0, 7, -210.0, 210.0 )
    self.__tof_corr_upstream_y_x = ROOT.TH2F('upstream_y_tof_x', \
           "Upstream Y Against TOF1 X", 400, -400.0, 400.0, 7, -210.0, 210.0 )


    self.__tof_corr_downstream_x = ROOT.TH2F('downstream_tof_x', \
         "Downstream X Against TOF2 X", 400, -400.0, 400.0, 7, -210.0, 210.0 )
    self.__tof_corr_downstream_y = ROOT.TH2F('downstream_tof_y', \
         "Downstream Y Against TOF2 Y", 400, -400.0, 400.0, 7, -210.0, 210.0 )

    self.__tof_corr_downstream_x_y = ROOT.TH2F('downstream_x_tof_y', \
         "Downstream X Against TOF2 Y", 400, -400.0, 400.0, 7, -210.0, 210.0 )
    self.__tof_corr_downstream_y_x = ROOT.TH2F('downstream_y_tof_x', \
         "Downstream Y Against TOF2 X", 400, -400.0, 400.0, 7, -210.0, 210.0 )


  def get_dependencies(self, inserter) :
    self.__spacepoints = inserter(tof_analysis.tof_extract_spacepoints())
    self.__tof = inserter(tof_analysis.tof_analyser())
    self.__straight_tracks = inserter(scifi_extractors.scifi_straight_track_extractor())


  def _reset(self) :
    pass


  def _process(self, file_reader) :

    if self.__tof.is_cut() :
      return True

    upstream_tracks = self.__straight_tracks.get_upstream_tracks()
    downstream_tracks = self.__straight_tracks.get_downstream_tracks()

    tof1_sp = self.__spacepoints.get_tof1_sp()
    tof2_sp = self.__spacepoints.get_tof2_sp()

# Only looking for single track events at present implementation
    if len( upstream_tracks ) == 1 :
      up_trk = upstream_tracks[0]

      for tp in up_trk.scifitrackpoints() :
        if tp.station() == 5 and tp.plane() == 2 :
          gra = [ tp.mom().x() / tp.mom().z(), tp.mom().y() / tp.mom().z() ]
          tof_dist = tof1_sp[0].GetGlobalPosZ() - tp.pos().z()

          pos = [ tp.pos().x() + gra[0]*tof_dist, tp.pos().y() + gra[1]*tof_dist ]

          self.__tof_residual_upstream.Fill( pos[0] - tof1_sp[0].GetGlobalPosX(), 
                                                pos[1] - tof1_sp[0].GetGlobalPosY() )

          self.__tof_corr_upstream_x.Fill( pos[0], tof1_sp[0].GetGlobalPosX() )
          self.__tof_corr_upstream_y.Fill( pos[1], tof1_sp[0].GetGlobalPosY() )
          self.__tof_corr_upstream_x_y.Fill( pos[0], tof1_sp[0].GetGlobalPosY() )
          self.__tof_corr_upstream_y_x.Fill( pos[1], tof1_sp[0].GetGlobalPosX() )


    if len( downstream_tracks) == 1 :
      down_trk = downstream_tracks[0]

      for tp in down_trk.scifitrackpoints() :
        if tp.station() == 5 and tp.plane() == 2 :
          gra = [ tp.mom().x() / tp.mom().z(), tp.mom().y() / tp.mom().z() ]
          tof_dist = tof2_sp[0].GetGlobalPosZ() - tp.pos().z()

          pos = [ tp.pos().x() + gra[0]*tof_dist, tp.pos().y() + gra[1]*tof_dist ]

          self.__tof_residual_downstream.Fill( pos[0] - tof2_sp[0].GetGlobalPosX(), 
                                                pos[1] - tof2_sp[0].GetGlobalPosY() )

          self.__tof_corr_downstream_x.Fill( pos[0], tof2_sp[0].GetGlobalPosX() )
          self.__tof_corr_downstream_y.Fill( pos[1], tof2_sp[0].GetGlobalPosY() )
          self.__tof_corr_downstream_x_y.Fill( pos[0], tof2_sp[0].GetGlobalPosY() )
          self.__tof_corr_downstream_y_x.Fill( pos[1], tof2_sp[0].GetGlobalPosX() )



  def _store_plots(self, plot_dict) :
    tof_corr_upstream = {}

    tof_corr_upstream['residuals'] = self.__tof_residual_upstream

    tof_corr_upstream['x'] = self.__tof_corr_upstream_x
    tof_corr_upstream['y'] = self.__tof_corr_upstream_y
    tof_corr_upstream['x_y'] = self.__tof_corr_upstream_x_y
    tof_corr_upstream['y_x'] = self.__tof_corr_upstream_y_x


    tof_corr_downstream = {}

    tof_corr_downstream['residuals'] = self.__tof_residual_downstream

    tof_corr_downstream['x'] = self.__tof_corr_downstream_x
    tof_corr_downstream['y'] = self.__tof_corr_downstream_y
    tof_corr_downstream['x_y'] = self.__tof_corr_downstream_x_y
    tof_corr_downstream['y_x'] = self.__tof_corr_downstream_y_x

    plot_dict['tof_correlations'] = {}

    plot_dict['tof_correlations']['upstream'] = tof_corr_upstream
    plot_dict['tof_correlations']['downstream'] = tof_corr_downstream

    return plot_dict


  def _store_data(self, data_dict) :
    return data_dict

