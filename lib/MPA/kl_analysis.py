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
  KL Analysis Classes are stored here.
"""

class kl_spacepoint_extractor(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "kl_spacepoint_extractor")

    self.__quality_cut = False
    self.__counter_hits = 0
    self.__counter_cut = 0
    self._reset()

    self.__hit_positions = ROOT.TH2F('kl_hit_positions', \
               "Positions of KL Hits", 200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__hit_positions_cut = ROOT.TH2F('kl_hit_positions_cut', \
         "Positions of KL Hits (CUT)", 200, -100.0, 100.0, 200, -100.0, 100.0 )


  def set_quality_cut(self, cut) :
    self.__quality_cut = cut
  

  def get_dependencies(self, inserter) :
    pass


  def get_args(self, parser) :
    parser.add_argument('--cut_kl_quality', type=bool, default=False, \
                             help='Set to cut poor quality reconstructed hits')


  def process_args(self, namespace) :
    self.__quality_cut = namespace.cut_kl_quality


  def _reset(self) :
    self.__pos = []
    self.__charge = 0.0
    self.__charge_product = 0.0
    self.__cell = -1
    self.__quality_flag = False

    self.__hits = []


  def _process(self, file_reader) :
    kl_event = file_reader.get_event('kl')

    cell_hit_container = kl_event.GetKLEventCellHit()
    number_hits = cell_hit_container.GetKLCellHitArraySize()

    self.__counter_hits += number_hits

    for hit_i in range(number_hits) :
      cell_hit = cell_hit_container.GetKLCellHitArrayElement(hit_i)

      self.__hit_positions.Fill(cell_hit.GetGlobalPosX(), \
                                                    cell_hit.GetGlobalPosY())

      if self.__quality_cut and not cell_hit.GetFlag() :
        self.__counter_cut += 1
      else :
        self.__hits.append(cell_hit)
        self.__hit_positions_cut.Fill(cell_hit.GetGlobalPosX(), \
                                                    cell_hit.GetGlobalPosY())


  def _store_plots(self, plot_dict) :
    kl_dict = {}

    kl_dict['hit_positions'] = self.__hit_positions
    kl_dict['hit_positions_cut'] = self.__hit_positions_cut

    plot_dict['kl'] = kl_dict
    return plot_dict


  def _store_data(self, data_dict) :
    kl_dict = {}

    kl_dict['total_hits'] = self.__counter_hits
    kl_dict['total_events'] = self.get_number_events()
    kl_dict['hits_event'] = self.__counter_hits / self.get_number_events()
    kl_dict['cut_hits'] = self.__counter_cut

    data_dict['kl_extractor'] = kl_dict
    return data_dict

