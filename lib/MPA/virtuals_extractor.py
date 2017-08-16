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
  Virtual Hits Extractor Classes are stored here.
"""

class virtual_hit_extractor(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "virtual_hit_extractor")
    
    self.__virtual_hits = []


  def get_virtual_hits(self) :
    return self.__virtual_hits


  def _reset(self) :
    self.__virtual_hits = []


  def get_dependencies(self, inserter) :
    pass


  def _process(self, file_reader) :
    mc_event = file_reader.get_event('mc')
    virtual_hits_count = mc_event.GetVirtualHitsSize()

    for virt_i in range(virtual_hits_count) :
      virt = mc_event.GetAVirtualHit(virt_i)
      self.__virtual_hits.append(virt)

  def _store_plots(self, plot_dict) :
    return plot_dict


  def _store_data(self, data_dict) :
    return data_dict


