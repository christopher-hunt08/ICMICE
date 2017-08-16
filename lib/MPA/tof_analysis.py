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


class tof_extract_spacepoints(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "tof_extract_spacepoints")

    self.__tof0_sp = []
    self.__tof1_sp = []
    self.__tof2_sp = []


  def get_tof0_sp(self) :
    return self.__tof0_sp


  def get_tof1_sp(self) :
    return self.__tof1_sp


  def get_tof2_sp(self) :
    return self.__tof2_sp


  def get_dependencies(self, inserter) :
    pass


  def _reset(self) :
    self.__tof0_sp = []
    self.__tof1_sp = []
    self.__tof2_sp = []


  def _process(self, file_reader) :
    tof_event = file_reader.get_event('tof')
    event_spacepoints = tof_event.GetTOFEventSpacePoint()

    tof0_sp_size = event_spacepoints.GetTOF0SpacePointArraySize()
    tof1_sp_size = event_spacepoints.GetTOF1SpacePointArraySize()
    tof2_sp_size = event_spacepoints.GetTOF2SpacePointArraySize()

    for tof0_i in range(tof0_sp_size) : 
      self.__tof0_sp.append( event_spacepoints.GetTOF0SpacePointArrayElement(tof0_i) )
    for tof1_i in range(tof1_sp_size) : 
      self.__tof1_sp.append(event_spacepoints.GetTOF1SpacePointArrayElement(tof1_i))
    for tof2_i in range(tof2_sp_size) : 
      self.__tof2_sp.append(event_spacepoints.GetTOF2SpacePointArrayElement(tof2_i))

    return False


################################################################################


class tof_analyser(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "tof")

    self.__plot_tof_0_1 = ROOT.TH1F( 'tof_0_1', 'Time TOF0 - TOF1', \
                                                       1000, 0.0, 100.0 )
    self.__plot_tof_1_2 = ROOT.TH1F( 'tof_1_2', 'Time TOF1 - TOF2', \
                                                       1000, 0.0, 100.0 )
    self.__plot_tof_0_2 = ROOT.TH1F( 'tof_0_2', 'Time TOF0 - TOF2', \
                                                       1000, 0.0, 200.0 )
    self.__plot_tof_0_1_cut = ROOT.TH1F( 'tof_0_1_cut', 'Time TOF0 - TOF1', \
                                                       1000, 0.0, 100.0 )
    self.__plot_tof_1_2_cut = ROOT.TH1F( 'tof_1_2_cut', 'Time TOF1 - TOF2', \
                                                       1000, 0.0, 100.0 )
    self.__plot_tof_0_2_cut = ROOT.TH1F( 'tof_0_2_cut', 'Time TOF0 - TOF2', \
                                                       1000, 0.0, 200.0 )
    self.__plot_tof_0_spacepoints = ROOT.TH2F( 'tof_0_spacepoints', \
                   "TOF0 Spacepoints", 20, -200.0, 200.0, 20, -200.0, 200.0 )
    self.__plot_tof_1_spacepoints = ROOT.TH2F( 'tof_1_spacepoints', \
                   "TOF1 Spacepoints", 20, -200.0, 200.0, 20, -200.0, 200.0 )
    self.__plot_tof_2_spacepoints = ROOT.TH2F( 'tof_2_spacepoints', \
                   "TOF2 Spacepoints", 20, -200.0, 200.0, 20, -200.0, 200.0 )

    self.__tof_cut_low = 0.0
    self.__tof_cut_high = 1000.0
    self.__window = False
    self.__tof_window_a = -1
    self.__tof_window_b = -1
    self.__require_tof0 = False
    self.__require_tof1 = False
    self.__require_tof2 = False


  def set_cut_tof_window(self, low_time, high_time) :
    self.__tof_cut_low = low_time
    self.__tof_cut_high = high_time


  def get_dependencies(self, inserter) :
    self.__tof_sp = inserter(tof_extract_spacepoints())


  def get_args(self, parser) :
    parser.add_argument( '--tof_window', nargs=4, default=None, \
                 help='Sepecify the <TOFa> <TOFb> <UPPER> <LOWER> for the TOF window' )
    parser.add_argument( '--require_tofs', nargs='+', type=int, default=[], help='Specify which TOF detectors should requrie a spacepoint' )


  def process_args(self, namespace) :
    if 0 in namespace.require_tofs :
      self.require_tof0
    if 1 in namespace.require_tofs :
      self.require_tof1
    if 2 in namespace.require_tofs :
      self.require_tof2

    if namespace.tof_window is not None :
      self.__window = True

      if '0' in namespace.tof_window[0:2] :
        self.__require_tof0 = True
      if '1' in namespace.tof_window[0:2] :
        self.__require_tof1 = True
      if '2' in namespace.tof_window[0:2] :
        self.__require_tof2 = True

      self.__tof_window_a = int(namespace.tof_window[0])
      self.__tof_window_b = int(namespace.tof_window[1])

      self.__tof_cut_low = float(namespace.tof_window[2])
      self.__tof_cut_high = float(namespace.tof_window[3])


  def _reset(self) :
    pass


  def _process(self, file_reader) :

    tof_sp = []
    tof_sp.append(self.__tof_sp.get_tof0_sp())
    tof_sp.append(self.__tof_sp.get_tof1_sp())
    tof_sp.append(self.__tof_sp.get_tof2_sp())

    tof0 = (len(tof_sp[0]) > 0)
    tof1 = (len(tof_sp[1]) > 0)
    tof2 = (len(tof_sp[2]) > 0)

#    tof0_sp = self.__tof_sp.get_tof0_sp()
#    tof1_sp = self.__tof_sp.get_tof1_sp()
#    tof2_sp = self.__tof_sp.get_tof2_sp()
#
#    tof0 = (len(tof0_sp) > 0)
#    tof1 = (len(tof1_sp) > 0)
#    tof2 = (len(tof2_sp) > 0)
  
    if not tof0 and self.__require_tof0 :
      self._cut()
      return True
    if not tof1 and self.__require_tof1 :
      self._cut()
      return True
    if not tof2 and self.__require_tof2 :
      self._cut()
      return True

    if tof0 :
      self.__plot_tof_0_spacepoints.Fill( tof_sp[0][0].GetGlobalPosX(), \
                                                   tof_sp[0][0].GetGlobalPosY() )
    if tof1 :
      self.__plot_tof_1_spacepoints.Fill( tof_sp[1][0].GetGlobalPosX(), \
                                                   tof_sp[1][0].GetGlobalPosY() )
    if tof2 :
      self.__plot_tof_2_spacepoints.Fill( tof_sp[2][0].GetGlobalPosX(), \
                                                   tof_sp[2][0].GetGlobalPosY() )

    if tof0 and tof1 :
      self.__plot_tof_0_1.Fill( tof_sp[1][0].GetTime() - tof_sp[0][0].GetTime() )
    if tof1 and tof2 :
      self.__plot_tof_1_2.Fill( tof_sp[2][0].GetTime() - tof_sp[1][0].GetTime() )
    if tof0 and tof2 :
      self.__plot_tof_0_2.Fill( tof_sp[2][0].GetTime() - tof_sp[0][0].GetTime() )

    if self.__window :
      diff = tof_sp[self.__tof_window_b][0].GetTime() - tof_sp[self.__tof_window_a][0].GetTime()

      if (diff < self.__tof_cut_low) or (diff > self.__tof_cut_high) :
        self._cut()
        return True

    if tof0 and tof1 :
      self.__plot_tof_0_1_cut.Fill( tof_sp[1][0].GetTime() - tof_sp[0][0].GetTime() )
    if tof1 and tof2 :
      self.__plot_tof_1_2_cut.Fill( tof_sp[2][0].GetTime() - tof_sp[1][0].GetTime() )
    if tof0 and tof2 :
      self.__plot_tof_0_2_cut.Fill( tof_sp[2][0].GetTime() - tof_sp[0][0].GetTime() )

    return False


  def _store_plots(self, plot_dict) :
    tof_plots = {}

    tof_plots['tof_0_1'] = self.__plot_tof_0_1
    tof_plots['tof_1_2'] = self.__plot_tof_1_2
    tof_plots['tof_0_2'] = self.__plot_tof_0_2
    tof_plots['tof_0_1_cut'] = self.__plot_tof_0_1_cut
    tof_plots['tof_1_2_cut'] = self.__plot_tof_1_2_cut
    tof_plots['tof_0_2_cut'] = self.__plot_tof_0_2_cut
    tof_plots['tof_0_spacepoints'] = self.__plot_tof_0_spacepoints
    tof_plots['tof_1_spacepoints'] = self.__plot_tof_1_spacepoints
    tof_plots['tof_2_spacepoints'] = self.__plot_tof_2_spacepoints

    plot_dict['tof'] = tof_plots

    return plot_dict


  def _store_data(self, data_dict) :
    return data_dict


