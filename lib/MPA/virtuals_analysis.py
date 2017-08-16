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
import array
from sets import Set

import framework
import virtuals_extractor

import analysis
from analysis import covariances
from analysis import tools
from analysis import hit_types
from analysis import inspectors

"""
  Virtual Hit Analysis Classes are stored here.
"""

################################################################################


class plane_data() :

  def __init__(self, plane_id) :
    self.plane_id = plane_id
    self.plane_position = 0.0
    self.field_sum = [0.0, 0.0, 0.0]
    self.inspector = None

    self.covariance = covariances.CovarianceMatrix()

  def get_mean_field(self) :
    num = self.covariance.length()
    return [self.field_sum[0]/num, self.field_sum[1]/num, self.field_sum[2]/num]


  def add_hit(self, hit) :
    self.plane_position = hit.GetPosition().z()
    field = hit.GetBField()
    self.field_sum[0] += field.x()*1000.0
    self.field_sum[1] += field.y()*1000.0
    self.field_sum[2] += field.z()*1000.0
    self.covariance.add_hit(hit_types.AnalysisHit(virtual_track_point=hit))

    if self.inspector is not None :
      self.inspector.add_hit(hit_types.AnalysisHit(virtual_track_point=hit))


  def add_primary(self, hit) :
    self.plane_position = hit.GetPosition().z()
    pos = hit.GetPosition()
    mom = hit.GetMomentum()
    time = hit.GetTime()
    station = -1
    pid = hit.GetParticleId()

    ana_hit = hit_types.AnalysisHit(x=pos.x(), y=pos.y(), z=pos.z(), px=mom.x(), py=mom.y(), pz=mom.z(), pid=pid, station=station)

    self.covariance.add_hit(ana_hit)

    if self.inspector is not None :
      self.inspector.add_hit(ana_hit)


################################################################################


class virtual_cut() :

  def __init__(self, plane_id) :
    self.__plane_id = plane_id

    self.__radius = 1.0E+30
    self.__pt_window = (0.0, 1.0E+30)
    self.__pz_window = (0.0, 1.0E+30)
    self.__momentum_window = (0.0, 1.0E+30)


  def __repr__(self) :
    return "{0} : {7}, ({1} : {2}), ({3} : {4}), ({5} : {6})".format(self.__plane_id, self.__pt_window[0], self.__pt_window[1], self.__pz_window[0], self.__pz_window[1], self.__momentum_window[0], self.__momentum_window[1], self.__radius)

  def set_radius(self, radius) :
    self.__radius = radius

  def set_pt_window(self, pt_low, pt_high) :
    self.__pt_window = (pt_low, pt_high)

  def set_momentum_window(self, p_low, p_high) :
    self.__momentum_window = (p_low, p_high)

  def set_pz_window(self, pz_low, pz_high) :
    self.__pz_window = (pz_low, pz_high)


  def get_plane_id(self) :
    return self.__plane_id


  def is_cut(self, vhit) :
    rad = math.sqrt(vhit.GetPosition().x()**2 + vhit.GetPosition().y()**2)
    pt = math.sqrt(vhit.GetMomentum().x()**2 + vhit.GetMomentum().y()**2)
    p = math.sqrt(vhit.GetMomentum().x()**2 + vhit.GetMomentum().y()**2 + vhit.GetMomentum().z()**2)

    if rad > self.__radius :
      return True
    elif pt < self.__pt_window[0] :
      return True
    elif pt > self.__pt_window[1] :
      return True
    elif vhit.GetMomentum().z() < self.__pz_window[0] :
      return True
    elif vhit.GetMomentum().z() > self.__pz_window[1] :
      return True
    elif p < self.__momentum_window[0] :
      return True
    elif p > self.__momentum_window[1] :
      return True
  
    return False



################################################################################


class virtual_beam_properties(framework.processor_base) :

  def __init__(self, inspector_position_resolution=0.0001) :
    framework.processor_base.__init__(self, "virtual_beam_properties")

    self.__current_track = None
    self.__global_cuts = True

    self.__cut_transmision = [0.0, 0.0]
    self.__cut_pid = [0]
    self.__save_data_file = None

    self.__inspectors = {}
    self.__plane_cuts = {}

    self.__primaries = plane_data("primaries")
    self.__primaries.inspector = inspectors.PhaseSpace2DInspector("primaries", 0)
    self.__analyse_primaries = False

    self.__plane_list = []
    self.__plane_list_length = 0
    self.__virt_ensemble_size = 0


  def get_number_planes(self) :
    return len(self.__plane_list)


  def get_beam_position(self, plane_id) :
    plot = self.__plane_list[plane_id].position_plot
    return plot.GetMean(0), plot.GetMean(1)


  def get_beam_momentum(self, plane_id) :
    plot = self.__plane_list[plane_id].momentum_plot
    return plot.GetMean(0), plot.GetMean(1)


  def is_cut(self) :
    if self.__current_track is None :
      return True
    else :
      return False


  def get_virtual_hits(self, plane_id) :
    if plane_id in self.__current_track :
      return self.__current_track[plane_id]
    else :
      return None


  def interpolate_beam_position(self, z_position) :
    pass


  def interpolate_beam_momentum(self, z_position) :
    pass


  def interpolate_beam_emittance(self, z_position) :
    pass


  def interpolate_beam_momentum(self, z_position) :
    pass


  def interpolate_beam_alpha(self, z_position) :
    pass


  def interpolate_beam_beta(self, z_position) :
    pass


  def _add_plane(self, plane_id) :

    while self.__plane_list_length <= plane_id :
      self.__plane_list.append(plane_data(self.__plane_list_length))
      if self.__plane_list_length in self.__inspections :
        self.__plane_list[self.__plane_list_length].inspector = inspectors.PhaseSpace2DInspector(self.__plane_list_length+1, self.__virt_ensemble_size)

      self.__plane_list_length += 1


  def get_dependencies(self, inserter) :
    self.__virt_extractor = inserter( \
                                    virtuals_extractor.virtual_hit_extractor())


  def get_args(self, parser) :
    parser.add_argument('--virt_cut_pid', nargs='+', type=int, \
         default=[13,-13], help='Cut on the PIDs of the particles in the beam')
    parser.add_argument('--virt_cut_transmission', \
                                           type=float, nargs=2, default=None, \
                        help='Expect particles to pass through virtual planes')
    parser.add_argument('--virt_inspections', type=int, nargs='+', \
 default=[], help="Specify virtual plane ids to analyse the MC Beam in detail")
    parser.add_argument('--virt_ensemble_size', type=int, default=0, \
              help ='Set the size of the ensemble of particles to calculate '+\
                               'emittance from. Zero uses the entire dataset.')
    parser.add_argument('--save_data_file', default=None, \
          help='Save the results to a plain text data file')

    parser.add_argument('--virt_cut_radius', nargs='+', metavar='cut_info', \
            action='append',\
            default=None, help='Cut on the radius of beam at multiple planes.')
    parser.add_argument('--virt_pt_window', nargs='+', metavar='cut_info', \
        action='append',\
        default=None, help='Cut on a window of Pt of the particles in the ' + \
                                              'beam at multple virtual planes')
    parser.add_argument('--virt_pz_window', nargs='+', action='append',\
                                                          metavar='cut_info', \
   default=None, help='Cut on a window of Pz of the particles in the beam ' + \
                                                  'at multiple virtual planes')
    parser.add_argument('--virt_momentum_window', nargs='+', action='append',\
                                            metavar='cut_info', default=None, \
     help='Cut on a window of total momentum of the particles in the beam ' + \
                                                  'at multiple virtual planes')
    parser.add_argument('--not_global_cuts', action='store_true',\
        help='Select whether a plane-based cut globally removes the track from the analysis.' )

    parser.add_argument('--analyse_primaries', action='store_true', help='Analyse the primary beam distribution')


  def process_args(self, namespace) :
    self.__cut_pid = namespace.virt_cut_pid
    self.__cut_transmission = namespace.virt_cut_transmission
    self.__inspections = [ ins - 1 for ins in namespace.virt_inspections ]
    self.__virt_ensemble_size = namespace.virt_ensemble_size
    self.__save_data_file = namespace.save_data_file
    self.__analyse_primaries = namespace.analyse_primaries

    if namespace.not_global_cuts :
      self.__global_cuts = False

    if namespace.virt_pt_window is not None :
      for option in namespace.virt_pt_window :
        cut_low = float(option[0])
        cut_high = float(option[1])
        for plane_id_str in option[2:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = virtual_cut(plane_id)
          self.__plane_cuts[plane_id].set_pt_window(cut_low, cut_high)

    if namespace.virt_pz_window is not None :
      for option in namespace.virt_pz_window :
        cut_low = float(option[0])
        cut_high = float(option[1])
        for plane_id_str in option[2:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = virtual_cut(plane_id)
          self.__plane_cuts[plane_id].set_pz_window(cut_low, cut_high)

    if namespace.virt_momentum_window is not None :
      for option in namespace.virt_momentum_window :
        cut_low = float(option[0])
        cut_high = float(option[1])
        for plane_id_str in option[2:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = virtual_cut(plane_id)
          self.__plane_cuts[plane_id].set_momentum_window(cut_low, cut_high)

    if namespace.virt_cut_radius is not None :
      for option in namespace.virt_cut_radius :
        cut = float(option[0])
        for plane_id_str in option[1:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = virtual_cut(plane_id)
          self.__plane_cuts[plane_id].set_radius(cut)


  def _analyse_primaries(self, primary) :
    self.__primaries.add_primary(primary)
    self.__primaries.plane_position = primary.GetPosition().z()



  def _reset(self) :
    pass


  def _process(self, file_reader) :
    if self.__analyse_primaries :
      self._analyse_primaries(file_reader.get_event('mc').GetPrimary())

    virtual_hits = self.__virt_extractor.get_virtual_hits()
    self.__current_track = None

    transmission_high = Set([])
    transmission_low = Set([])
    bad_tracks = Set([])

    for virt in virtual_hits :
      plane_id = virt.GetStationId()

      if plane_id >= self.__plane_list_length :
        pos = virt.GetPosition() 
        self._add_plane(plane_id)
        self.__plane_list[plane_id].plane_position = pos.z()

    
      if virt.GetParticleId() not in self.__cut_pid :
        bad_tracks.add(virt.GetTrackId())
        continue
      if self.__global_cuts :
        if plane_id in self.__plane_cuts and self.__plane_cuts[plane_id].is_cut(virt) :
          bad_tracks.add(virt.GetTrackId())
          continue
      if self.__cut_transmission is not None :
        if plane_id == self.__cut_transmission[0] :
          transmission_low.add(virt.GetTrackId())
        if plane_id == self.__cut_transmission[1] :
          transmission_high.add(virt.GetTrackId())
        continue
      else :
        transmission_high.add(virt.GetTrackId())
        transmission_low.add(virt.GetTrackId())


    allowed_tracks = (transmission_low & transmission_high) - bad_tracks
#    print allowed_tracks, transmission_low, transmission_high, bad_tracks
    if len(allowed_tracks) > 0 :
      self.__current_track = {}
      for virt in virtual_hits :
        if virt.GetTrackId() in allowed_tracks :
          station_id = virt.GetStationId() 

          if station_id in self.__plane_cuts and self.__plane_cuts[station_id].is_cut(virt) :
            allowed_tracks.remove(virt.GetTrackId())
            continue
          self.__plane_list[station_id].add_hit(virt)

          if station_id in self.__current_track and virt.GetTrackId() == 1 :
            self.__current_track[station_id][virt.GetTrackId()] = virt
          else :
            self.__current_track[station_id] = { virt.GetTrackId() : virt }


  def conclude(self) :
    for datum in self.__plane_list :
      if datum.covariance.length() <= 1 :
        continue
      if datum.inspector is not None :
        datum.inspector.fill_plots()

    if self.__save_data_file is not None :
      with open(self.__save_data_file+'.dat', 'w') as outfile :
        outfile.write("# ID  Num  X  Y  Z  Px  Py  Pz  P  Bz  Emittance Alpha Beta\n")
        for datum in self.__plane_list :
          if datum.covariance.length() <= 1 :
            continue
          else :
            means = datum.covariance.get_means(['x', 'y', 'z', 'px', 'py', 'pz'])
            P = datum.covariance.get_momentum()
            outfile.write(str(datum.plane_id) + ' ' +
                          str(datum.covariance.length()) + ' ' + 
                          str(means['x']) + ' ' +
                          str(means['y']) + ' ' +
                          str(means['z']) + ' ' +
                          str(means['px']) + ' ' +
                          str(means['py']) + ' ' +
                          str(means['pz']) + ' ' +
                          str(P) + ' ' +
                          str(datum.get_mean_field()[2]) + ' ' +
                          str(datum.covariance.get_emittance(['x', 'px', 'y', 'py'])) + ' ' +
                          str(datum.covariance.get_alpha(['x', 'y'])) + ' ' + 
                          str(datum.covariance.get_beta(['x', 'y'])) + '\n')



  def _store_plots(self, plot_dict) :
    inspection_plots = {}
    virt_plots = {}
    z_pos = array.array('d')
    beam_rms = array.array('d')
    alpha = array.array('d')
    beta = array.array('d')
    momentum = array.array('d')
    emittance = array.array('d')
    emittance_x = array.array('d')
    emittance_y = array.array('d')
    am = array.array('d')
#    conam = array.array('d')
    field = array.array('d')
    number = array.array('d')

    if self.__analyse_primaries :
      virt_plots['primaries'] = self.__primaries.inspector.get_plot_dictionary()

      z_pos.append(self.__primaries.plane_position)
      beam_rms.append(self.__primaries.inspector.covariance.get_rms(['x', 'y']))
      alpha.append(self.__primaries.inspector.covariance.get_alpha(['x', 'y']))
      beta.append(self.__primaries.inspector.covariance.get_beta(['x', 'y']))
      momentum.append(self.__primaries.inspector.covariance.get_momentum())
      emittance.append(self.__primaries.inspector.covariance.get_emittance(['x', 'px', 'y', 'py']))
      emittance_x.append(self.__primaries.inspector.covariance.get_emittance(['x', 'px']))
      emittance_y.append(self.__primaries.inspector.covariance.get_emittance(['y', 'py']))
      am.append(self.__primaries.inspector.covariance.get_angular_momentum())
#      conam.append(self.__primaries.inspector.covariance.get_canonical_angular_momentum(datum.get_mean_field()[2]))
      number.append(self.__primaries.inspector.covariance.length())


    for datum in self.__plane_list :
      if datum.covariance.length() <= 1 :
        continue
      z_pos.append(datum.plane_position)
      beam_rms.append(datum.covariance.get_rms(['x', 'y']))
      alpha.append(datum.covariance.get_alpha(['x', 'y']))
      beta.append(datum.covariance.get_beta(['x', 'y']))
      momentum.append(datum.covariance.get_momentum())
      emittance.append(datum.covariance.get_emittance(['x', 'px', 'y', 'py']))
      emittance_x.append(datum.covariance.get_emittance(['x', 'px']))
      emittance_y.append(datum.covariance.get_emittance(['y', 'py']))
      am.append(datum.covariance.get_angular_momentum())
#      conam.append(datum.covariance.get_canonical_angular_momentum(datum.get_mean_field()[2]))
      field.append(datum.get_mean_field()[2])
      number.append(datum.covariance.length())

      if datum.inspector is not None :
        inspection_plots[str(datum.inspector.plane)] = datum.inspector.get_plot_dictionary()

    virt_plots['inspections'] = inspection_plots
    if len(z_pos) == 0 :
      virt_plots['alpha_function'] = None
      virt_plots['beta_function'] = None
      virt_plots['emittance'] = None
      virt_plots['emittance_x'] = None
      virt_plots['emittance_y'] = None
      virt_plots['angular_momentum'] = None
      virt_plots['canonical_angular_momentum'] = None
      virt_plots['momentum'] = None
      virt_plots['beam_rms'] = None
      virt_plots['field'] = None
      virt_plots['number'] = None
      virt_plots['beam_rms'] = None
    else :
      virt_plots['alpha_function'] = ROOT.TGraph(len(z_pos), z_pos, alpha)
      virt_plots['beta_function'] = ROOT.TGraph(len(z_pos), z_pos, beta)
      virt_plots['emittance'] = ROOT.TGraph(len(z_pos), z_pos, emittance)
      virt_plots['emittance_x'] = ROOT.TGraph(len(z_pos), z_pos, emittance_x)
      virt_plots['emittance_y'] = ROOT.TGraph(len(z_pos), z_pos, emittance_y)
#      virt_plots['angular_momentum'] = ROOT.TGraph(len(z_pos), z_pos, am)
#      virt_plots['canonical_angular_momentum'] = ROOT.TGraph(len(z_pos), z_pos, conam)
      virt_plots['momentum'] = ROOT.TGraph(len(z_pos), z_pos, momentum)
      virt_plots['beam_rms'] = ROOT.TGraph(len(z_pos), z_pos, beam_rms)
      virt_plots['field'] = ROOT.TGraph(len(z_pos), z_pos, field)
      virt_plots['number'] = ROOT.TGraph(len(z_pos), z_pos, number)
      virt_plots['beam_rms'] = ROOT.TGraph(len(z_pos), z_pos, beam_rms)

    plot_dict['virtual_plots'] = virt_plots

    return plot_dict


  def _store_data(self, data_dict) :

    inspector_counter = 0
    inspector_data = []

    if self.__analyse_primaries :
      data_dict['primaries'] = self.__primaries.inspector.get_data_dictionary()

    for plane_i in range(len(self.__plane_list)) :
      datum = self.__plane_list[plane_i]
      num = datum.covariance.length()
      if num <= 1 :
        continue

      if datum.inspector is not None :
        inspector_data.append(datum.inspector.get_data_dictionary())

    data_dict['inspections'] = inspector_data
    data_dict['number_planes_analysed'] = len(self.__plane_list)


    return data_dict



