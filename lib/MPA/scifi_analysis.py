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
import json
import os
import math
import array
import numpy
import random

import framework
import scifi_extractors
import tof_analysis
import virtuals_analysis
import analysis
from analysis import tools
from analysis import covariances
from analysis import hit_types
from analysis import beam_sampling

"""
  SciFi Analysis Classes are stored here.
"""

################################################################################


class scifi_internal_alignment(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_internal_alignment")

    self.__upstream_1_5_x = ROOT.TH2F('upstream_1_5_x', \
                                      "Upstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__upstream_1_5_y = ROOT.TH2F('upstream_1_5_y', \
                                      "Upstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__downstream_1_5_x = ROOT.TH2F('downstream_1_5_x', \
                                    "Downstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__downstream_1_5_y = ROOT.TH2F('downstream_1_5_y', \
                                    "Downstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )

    self.__upstream_5_1_x = ROOT.TH2F('upstream_5_1_x', \
                                      "Upstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__upstream_5_1_y = ROOT.TH2F('upstream_5_1_y', \
                                      "Upstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__downstream_5_1_x = ROOT.TH2F('downstream_5_1_x', \
                                    "Downstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )
    self.__downstream_5_1_y = ROOT.TH2F('downstream_5_1_y', \
                                    "Downstream Track End to End Comparison", \
                                       200, -100.0, 100.0, 200, -100.0, 100.0 )

    self.__upstream_patrec_kalman_x = ROOT.TH2F('upstream_patrec_kalman_x', \
                        "Upstream pattern recognition-kalman residuals in x", \
                                           15, 0.5, 15.5, 201, -100.5, 100.5 )
    self.__upstream_patrec_kalman_y = ROOT.TH2F('upstream_patrec_kalman_y', \
                        "Upstream pattern recognition-kalman residuals in y", \
                                           15, 0.5, 15.5, 201, -100.5, 100.5 )

    self.__downstream_patrec_kalman_x = ROOT.TH2F('downstream_patrec_kalman_x', \
                        "Downstream pattern recognition-kalman residuals in x", \
                                           15, 0.5, 15.5, 201, -100.5, 100.5 )
    self.__downstream_patrec_kalman_y = ROOT.TH2F('downstream_patrec_kalman_y', \
                        "Downstream pattern recognition-kalman residuals in y", \
                                           15, 0.5, 15.5, 201, -100.5, 100.5 )

    self.__up_station_align = []
    self.__down_station_align = []
    for station_i in range(1, 6) :
      self.__up_station_align.append(ROOT.TH2F(
                                       'station_alignment_0_'+str(station_i), \
                            "Spacepoint distance from fit 0."+str(station_i), \
                                          200, -10.0, 10.0, 200, -10.0, 10.0 ))
      self.__down_station_align.append(ROOT.TH2F(
                                       'station_alignment_1_'+str(station_i), \
                            "Spacepoint distance from fit 1."+str(station_i), \
                                          200, -10.0, 10.0, 200, -10.0, 10.0 ))

    self.__track_extractor = None


  def get_dependencies(self, inserter) :
    self.__track_extractor = inserter(scifi_extractors.scifi_straight_track_candidates())


  def _reset(self) :
    pass


  def _process(self, file_reader) :
    if not self.__track_extractor.is_cut() :
      up_tp_1 = self.__track_extractor.get_upstream_trackpoint(1)
      up_tp_15 = self.__track_extractor.get_upstream_trackpoint(15)

      down_tp_1 = self.__track_extractor.get_downstream_trackpoint(1)
      down_tp_15 = self.__track_extractor.get_downstream_trackpoint(15)

      up_delta_z = up_tp_15.pos().z() - up_tp_1.pos().z()
      down_delta_z = down_tp_15.pos().z() - down_tp_1.pos().z()

      up_1_grad = [ up_tp_1.mom().x() / up_tp_1.mom().z(), \
                                        up_tp_1.mom().y() / up_tp_1.mom().z() ]

      down_1_grad = [ down_tp_1.mom().x() / down_tp_1.mom().z(), \
                                    down_tp_1.mom().y() / down_tp_1.mom().z() ]


      self.__upstream_1_5_x.Fill((up_tp_1.pos().x() + up_1_grad[0]*up_delta_z), \
                                                              up_tp_15.pos().x())
      self.__upstream_1_5_y.Fill((up_tp_1.pos().y() + up_1_grad[1]*up_delta_z), \
                                                              up_tp_15.pos().y())

      self.__downstream_1_5_x.Fill((down_tp_1.pos().x() + down_1_grad[0]*down_delta_z), \
                                                              down_tp_15.pos().x())
      self.__downstream_1_5_y.Fill((down_tp_1.pos().y() + down_1_grad[1]*down_delta_z), \
                                                              down_tp_15.pos().y())

      up_15_grad = [ up_tp_15.mom().x() / up_tp_15.mom().z(), \
                                      up_tp_15.mom().y() / up_tp_15.mom().z() ]

      down_15_grad = [ down_tp_15.mom().x() / down_tp_15.mom().z(), \
                                  down_tp_15.mom().y() / down_tp_15.mom().z() ]

      self.__upstream_5_1_x.Fill((up_tp_15.pos().x() - up_15_grad[0]*up_delta_z), \
                                                              up_tp_1.pos().x())
      self.__upstream_5_1_y.Fill((up_tp_15.pos().y() - up_15_grad[1]*up_delta_z), \
                                                              up_tp_1.pos().y())

      self.__downstream_5_1_x.Fill((down_tp_15.pos().x() - down_15_grad[0]*down_delta_z), \
                                                              down_tp_1.pos().x())
      self.__downstream_5_1_y.Fill((down_tp_15.pos().y() - down_15_grad[1]*down_delta_z), \
                                                              down_tp_1.pos().y())



      up_ref_pos = self.__track_extractor.get_upstream_track().\
                           pr_track_pointer_straight().get_reference_position()
      up_ref_gra = self.__track_extractor.get_upstream_track().\
                           pr_track_pointer_straight().get_reference_momentum()
      down_ref_pos = self.__track_extractor.get_downstream_track().\
                           pr_track_pointer_straight().get_reference_position()
      down_ref_gra = self.__track_extractor.get_downstream_track().\
                           pr_track_pointer_straight().get_reference_momentum()

      up_ref_gra.SetX(up_ref_gra.x() / up_ref_gra.z())
      up_ref_gra.SetY(up_ref_gra.y() / up_ref_gra.z())

      down_ref_gra.SetX(down_ref_gra.x() / down_ref_gra.z())
      down_ref_gra.SetY(down_ref_gra.y() / down_ref_gra.z())

      up_spacepoints = self.__track_extractor.get_upstream_track().\
                         pr_track_pointer_straight().get_spacepoints_pointers()
      down_spacepoints = self.__track_extractor.get_downstream_track().\
                         pr_track_pointer_straight().get_spacepoints_pointers()

      for sp in up_spacepoints :
        station = sp.get_station()
      
        sp_pos = sp.get_position()
        diff_z = - sp_pos.z()

        ex_x = -(up_ref_pos.x() + diff_z * up_ref_gra.x())
        ex_y =   up_ref_pos.y() + diff_z * up_ref_gra.y()

        diff_x = ex_x - sp_pos.x()
        diff_y = ex_y - sp_pos.y()

        self.__up_station_align[station-1].Fill(diff_x, diff_y)

      for sp in down_spacepoints :
        station = sp.get_station()
      
        sp_pos = sp.get_position()
        diff_z = sp_pos.z()

        ex_x = down_ref_pos.x() + diff_z * down_ref_gra.x()
        ex_y = down_ref_pos.y() + diff_z * down_ref_gra.y()

        diff_x = ex_x - sp_pos.x()
        diff_y = ex_y - sp_pos.y()

        self.__down_station_align[station-1].Fill(diff_x, diff_y)

      for plane in range( 1, 16 ) :
        tp = self.__track_extractor.get_upstream_trackpoint(plane)
        sp = None
        for test_sp in up_spacepoints :
          if test_sp.get_station() == int( (plane-1) / 3 ) + 1 :
            sp = test_sp
            break
        else :
          continue

        self.__upstream_patrec_kalman_x.Fill( plane, tp.pos().x() - sp.get_global_position().x() )
        self.__upstream_patrec_kalman_y.Fill( plane, tp.pos().y() - sp.get_global_position().y() )


      for plane in range( 1, 16 ) :
        tp = self.__track_extractor.get_downstream_trackpoint(plane)
        sp = None
        for test_sp in down_spacepoints :
          if test_sp.get_station() == int( (plane-1) / 3 ) + 1 :
            sp = test_sp
            break
        else :
          continue

        self.__downstream_patrec_kalman_x.Fill( plane, tp.pos().x() - sp.get_global_position().x() )
        self.__downstream_patrec_kalman_y.Fill( plane, tp.pos().y() - sp.get_global_position().y() )

      
    return False


  def _store_plots(self, plot_dict) :
    internal_alignment = {}
    internal_alignment['upstream_1_5_x'] = self.__upstream_1_5_x
    internal_alignment['upstream_1_5_y'] = self.__upstream_1_5_y
    internal_alignment['downstream_1_5_x'] = self.__downstream_1_5_x
    internal_alignment['downstream_1_5_y'] = self.__downstream_1_5_y

    internal_alignment['upstream_5_1_x'] = self.__upstream_5_1_x
    internal_alignment['upstream_5_1_y'] = self.__upstream_5_1_y
    internal_alignment['downstream_5_1_x'] = self.__downstream_5_1_x
    internal_alignment['downstream_5_1_y'] = self.__downstream_5_1_y

    internal_alignment['upstream_patrec_kalman_x'] = self.__upstream_patrec_kalman_x
    internal_alignment['upstream_patrec_kalman_y'] = self.__upstream_patrec_kalman_y
    internal_alignment['downstream_patrec_kalman_x'] = self.__downstream_patrec_kalman_x
    internal_alignment['downstream_patrec_kalman_y'] = self.__downstream_patrec_kalman_y

    for station_i in range(1, 6) :
      internal_alignment['upstream_station_alignment_'+str(station_i)] = \
                                           self.__up_station_align[station_i-1]
      internal_alignment['downstream_station_alignment_'+str(station_i)] = \
                                         self.__down_station_align[station_i-1]

    plot_dict['scifi_internal_alignment'] = internal_alignment
    return plot_dict


  def _store_data(self, data_dict) :
    return data_dict



################################################################################


class scifi_beam_selection(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_beam_selection")

    self.__mass = None
    self.__samplers = []
    self.__sampler_names = []
    self.__vector_function = None
    self.__do_weights = False

    self.__selection_tracker = None
    self.__selection_station = None
    self.__selection_plane = None

    self.__misses = analysis.inspectors.PhaseSpace2DInspector("selection_misses", -1)
    self.__number_parent = 0
    self.__number_pass = 0
    self.__number_fail = 0



  def get_dependencies(self, inserter) :

    self.__tracks_extract = [ inserter(scifi_extractors.scifi_helical_track_processor(0)), inserter(scifi_extractors.scifi_helical_track_processor(1, False)) ]


  def get_args(self, parser) :
    parser.add_argument( '--select_with_weights', action="store_true", help='Forces the selection algorithm to apply statistical weights - not cuts' )
    parser.add_argument( '--select_plane', type=int, nargs=3, default=None, help='Specify the <tracker> <station> <plane> to act as the selection plane' )
    parser.add_argument( '--select_momentum', type=float, nargs=2, default=None, help='Specify the mean momentum and standard deviation to select' )
    parser.add_argument( '--select_momentum_window', type=float, nargs=2, default=None, help='Specify the mean momentum and window size' )
    parser.add_argument( '--select_beam', nargs=3, default=None, help='Specify the <emittance> <beta> <p> to define the daughter distribution.')
    parser.add_argument( '--select_phasespace', nargs=4, default=None, help='Specify the <emittance> <beta> <alpha> <p> to define the daughter distribution in uncorrelated separate X and Y phase-space distributions.')
    parser.add_argument( '--select_covariance_matrix', default=None, help='Specify a file containing a covariance matrix to define the daughter distribution.')
    parser.add_argument( '--select_beta', nargs=1, default=None, help='Specify the <beta> to define the daughter distribution.')
    parser.add_argument( '--select_amplitude', nargs='?', const=0, default=None, type=float, help='Specify the emittance to define the daughter distribution.')
    parser.add_argument( '--select_uniform_amplitude', nargs=1, default=None, type=float, help='Specify the max emittance to define the daughter distribution.')


  def process_args(self, namespace) :
    self.__mass = namespace.mass_assumption
    self.__do_weights = namespace.select_with_weights

    if framework.get_last_analysis_root() is None :
      if (namespace.select_momentum is not None) or (namespace.select_momentum_window is not None) or (namespace.select_beam is not None) or (namespace.select_phasespace is not None) \
        or (namespace.select_covariance_matrix is not None) or (namespace.select_beta is not None) or (namespace.select_amplitude is not None) or (namespace.select_uniform_amplitude is not None) :
        raise ValueError("Last analysis MUST be specified for beam selection to be possible. Parent distributions are required.")

    if namespace.select_plane is not None :
      self.__selection_tracker = namespace.select_plane[0]
      self.__selection_station = namespace.select_plane[1]
      self.__selection_plane = namespace.select_plane[2]
    else :
      return

    tr = "tracker_"+str(namespace.select_plane[0])
    st = "station_"+str(namespace.select_plane[1])
    pl = "plane_"+str(namespace.select_plane[2])

    if namespace.select_momentum is not None :
      parent_dist = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/p")
      mean = namespace.select_momentum[0]
      std = namespace.select_momentum[1]

      self.__samplers.append(beam_sampling.GaussianMomentumSampler(parent_dist, mean, std))
      self.__sampler_names.append("momentum_sampler")
      

    if namespace.select_momentum_window is not None :
      parent_dist = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/p")
      mean = namespace.select_momentum_window[0]
      window = namespace.select_momentum_window[1]

      self.__samplers.append(beam_sampling.MomentumWindowSampler(parent_dist, mean, window))
      self.__sampler_names.append("momentum_window_sampler")
      

    if namespace.select_beam is not None :
      parent_covariance = numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tr][st][pl]['covariance_matrix'])

      emittance = float(namespace.select_beam[0])
      beta = float(namespace.select_beam[1])
      momentum = float(namespace.select_beam[2])


      cov_xx = beta*emittance*self.__mass / momentum
      cov_pxpx = (1.0 / beta)*emittance*momentum*self.__mass

      selection_covariance = numpy.diag([cov_xx, cov_pxpx, cov_xx, cov_pxpx])
      means = [0.0, 0.0, 0.0, 0.0]
      self.__samplers.append(beam_sampling.Gaussian4DSampler(means, parent_covariance, means, selection_covariance))
      self.__sampler_names.append("analytic_4D_guassian")


    if namespace.select_beta is not None :
      temp_covariance = numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tr][st][pl]['covariance_matrix'])
      parent_covariance = [ [ temp_covariance[0][0], temp_covariance[0][2] ], \
                            [ temp_covariance[2][0], temp_covariance[2][2] ] ]

      cov_xx = float(namespace.select_beam[1])

      selection_covariance = numpy.diag([cov_xx, cov_xx])
      means = [0.0, 0.0]
      self.__samplers.append(beam_sampling.GaussianSampler(means, parent_covariance, means, selection_covariance))
      self.__sampler_names.append("analytic_xy_gaussian")


    if namespace.select_phasespace is not None :
      parent_dist_x_px = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/x_px")
      parent_dist_x_py = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/x_py")
      parent_dist_y_px = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/y_px")
      parent_dist_y_py = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/y_py")
      parent_dist_x_y = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/x_y")
      parent_dist_px_py = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/px_py")

      emittance = float(namespace.select_phasespace[0])
      beta = float(namespace.select_phasespace[1])
      alpha = float(namespace.select_phasespace[2])
      momentum = float(namespace.select_phasespace[3])


      cov_xx = beta*emittance*self.__mass / momentum
      cov_xpx = -1.0*alpha*emittance*self.__mass
      cov_pxpx = ((1.0 + alpha**2) / beta)*emittance*momentum*self.__mass

      selection_covariance = numpy.array([[cov_xx,  cov_xpx,  0.0,     0.0     ], \
                                          [cov_xpx, cov_pxpx, 0.0,     0.0     ], \
                                          [0.0,     0.0,      cov_xx,  cov_xpx ], \
                                          [0.0,     0.0,      cov_xpx, cov_pxpx]])
      means = numpy.array([0.0, 0.0, 0.0, 0.0])
      self.__samplers.append(beam_sampling.XY4DPhaseSpaceSampler(parent_dist_x_px, parent_dist_x_py, parent_dist_y_px, parent_dist_y_py, parent_dist_x_y, parent_dist_px_py, means, selection_covariance))
      self.__sampler_names.append("binned_4D_projections")


    if namespace.select_covariance_matrix is not None :
      parent_dist_x_px = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/x_px")
      parent_dist_x_py = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/x_py")
      parent_dist_y_px = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/y_px")
      parent_dist_y_py = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/y_py")
      parent_dist_x_y = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/x_y")
      parent_dist_px_py = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/px_py")


      selection_covariance = numpy.loadtxt(namespace.select_covariance_matrix)
      means = numpy.array([0.0, 0.0, 0.0, 0.0])
      self.__samplers.append(beam_sampling.XY4DPhaseSpaceSampler(parent_dist_x_px, parent_dist_x_py, parent_dist_y_px, parent_dist_y_py, parent_dist_x_y, parent_dist_px_py, means, selection_covariance))
      self.__sampler_names.append("binned_4D_projections_from_matrix")


    if namespace.select_amplitude is not None :
      parent_dist = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/amplitude")

      if namespace.select_amplitude == 0 :
        emittance = parent_dist.GetMean() * 4
      else :
        emittance = float(namespace.select_amplitude)

      max_val = emittance*20.0
      covariance = numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tr][st][pl]['covariance_matrix'])

      self.__samplers.append(beam_sampling.Amplitude4DSampler(parent_dist, covariance, emittance, max_x=max_val))
      self.__sampler_names.append("chi-squared_4D_amplitude")


    if namespace.select_uniform_amplitude is not None :
      parent_dist = framework.get_last_analysis_root().Get("emittance_reconstruction/"+tr+"/"+st+"/"+pl+"/amplitude")

      max_amplitude = float(namespace.select_uniform_amplitude[0])

      covariance = numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tr][st][pl]['covariance_matrix'])

      self.__samplers.append(beam_sampling.UniformAmplitude4DSampler(parent_dist, covariance, max_amplitude))
      self.__sampler_names.append("uniform_4D_amplitude")



  def _process(self, file_reader) :
    # If no sampling rely on individual tracker cuts
    if self.__selection_tracker is None :
      framework.set_event_statistical_weight(1.0)
      return False

    # Cut if the specific selection tracker is cut
    if self.__tracks_extract[self.__selection_tracker].is_cut() :
      return True

    tp = self.__tracks_extract[self.__selection_tracker].get_trackpoint_byplane(self.__selection_station, self.__selection_plane)
    if tp is None : return True

    if len(self.__samplers) != 0 :
      hit = hit_types.AnalysisHit(scifi_track_point=tp)

      self.__number_parent += 1

      if self.__do_weights :
        weight = framework.get_event_statistical_weight()
        for sampler in self.__samplers :
          sampler.weight(hit)
          if new_weight < 1.0e-9 :
            self.__number_fail += 1
            self.__misses.add_hit(hit)
            return True
          else :
            weight *= new_weight
        else :
          framework.set_event_statistical_weight(weight)
          self.__number_pass += 1
          return False
      else :
        for sampler in self.__samplers :
          if not sampler.accept(hit) :
            framework.set_event_statistical_weight(0.0)
            self.__number_fail += 1
            self.__misses.add_hit(hit)
            return True
        else :
          framework.set_event_statistical_weight(1.0)
          self.__number_pass += 1
          return False

    else :
      return False


  def get_track_extractor(self, tracker) :
    return self.__tracks_extract[tracker]


  def conclude(self) :
    if self.__selection_tracker is not None :
      print
      print "Beam Selection Tested {0} Events.".format(self.__number_parent)
      print "  Accepted {0}".format(self.__number_pass)
      print "  Tossed   {0}".format(self.__number_fail)
      print


  def _store_data(self, data_dict) :
    data_dict['beam_sampling'] = {'selection_misses' : self.__misses.get_data_dictionary()}
    return data_dict


  def _store_plots(self, plot_dict) :
    sampler_dict = {}

    for sampler, name in zip(self.__samplers, self.__sampler_names) :
      sampler_dict[name] = sampler.get_plots()

    plot_dict['beam_sampling'] = sampler_dict
    plot_dict['beam_sampling']['selection_misses'] = self.__misses.get_plot_dictionary()
    return plot_dict



################################################################################


class scifi_emittance_reconstruction(framework.processor_base) :

  def __init__(self, use_mc=False) :
    framework.processor_base.__init__(self, "scifi_emittance_recon")

    self.__output_filename = None
    self.__output_directory = None
    self.__use_mc = use_mc
    self.__require_upstream = False
    self.__require_downstream = False

    self.__recon_trackers = [0,1]
    self.__recon_stations = [1]
    self.__recon_planes = [0]
    self.__analyse_cooling = False
    self.__reference_station = None
    self.__reference_plane = None
    self.__ensemble_size = 0
    self.__save_data_file = None
    self.__save_covariances = False
    self.__mc_track_extractor = None
    self.__mc_scifi_dictionary = None
    self.__graphs = {'data':{}, 'mc':{}, 'residuals':{}}
    self.__cooling_inspector = None

    self.__correction_upstream = None
    self.__correction_downstream = None
    self.__p_correction_upstream = [0.0, 0.0]
    self.__p_correction_downstream = [0.0, 0.0]

    self.__inspectors = []
    self.__mc_inspectors = []

    self.__emittance_graph = None
    self.__beta_graph = None
    self.__alpha_graph = None
    self.__momentum_graph = None

    for tr in [ 0, 1 ] :
      self.__inspectors.append([])
      self.__mc_inspectors.append([])

      for st in range(5) :
        self.__inspectors[tr].append([])
        self.__mc_inspectors[tr].append([])

        for pl in range(15) :
          self.__inspectors[tr][st].append(None)
          self.__mc_inspectors[tr][st].append(None)


  def get_dependencies(self, inserter) :
    self.__beam_selector = inserter(scifi_beam_selection())
#    self.__tracks_extract = [ inserter(scifi_extractors.scifi_helical_track_processor(0)), inserter(scifi_extractors.scifi_helical_track_processor(1, False)) ]
    self.__tof_analysis = inserter(tof_analysis.tof_analyser())
    if self.__use_mc :
      self.__mc_track_extractor = inserter(virtuals_analysis.virtual_beam_properties())


  def get_args(self, parser) :
    parser.add_argument( '--ensemble_size', type=int, default=0, help='Set the size of the ensemble of particles to calculate emittance from. Zero uses the entire dataset.' )

    parser.add_argument( '--recon_planes', type=int, nargs='+', default=[0], help='Choose the reconstruction plane in the trackers. This is the location that the emittance will be reconstructed for cooling calculations.', metavar='PLANE' )
    parser.add_argument( '--recon_stations', type=int, nargs='+', default=[1], help='Choose the reconstruction station in the trackers. This is the location that the emittance will be reconstructed for cooling calculations.', metavar='STATION' )
    parser.add_argument( '--recon_trackers', type=int, nargs='+', default=[0,1], help='Choose the reconstruction tracker. This is the location that the emittance will be reconstructed for cooling calculations.', metavar='TRACKER' )
    parser.add_argument( '--unbias_results', action='store_true', help='Only measure one plane for each particle track. Means the number of particles is reduced by the number of measurement planes, but protectes the results from unintentional bias' )
    if self.__use_mc :
      parser.add_argument( '--scifi_virt_dict', type=str, default='scifi_virtual_dict.json', help='Specify the json file name to load the ScifiPlane:VirtualPlane Dictionary.')

    parser.add_argument( '--correction_upstream', help='Enter the location of the upstream correction matrix to use in emittance calculations.', metavar='FILENAME', default=None )
    parser.add_argument( '--correction_downstream', help='Enter the location of the downstream correction matrix to use in emittance calculations.', metavar='FILENAME', default=None )
#    parser.add_argument( '--correction_directory', default=None, help='Specify a directory in which the up- and downstream correction matrices can be found. They are expected to be named "measure_emittance_correction_matrix_up.dat, etc.' )

    parser.add_argument( '--reference_plane', default=[-1, -1], nargs=2, type=int, help='Enter the Station and Plane of the trackers that are the nominated reference planes' )

    parser.add_argument( '--p_correction_upstream', nargs=2, help='Enter the linear correction (const + gradient) to the upstream total momentum.', default=[0.0, 0.0] )
    parser.add_argument( '--p_correction_downstream', nargs=2, help='Enter the linear correction (const + gradient) to the downstream total momentum.', default=[0.0, 0.0] )

    parser.add_argument( '--require_upstream', action='store_true', help='Requires there the be precisely 1 track in the upstream tracker' )
    parser.add_argument( '--require_downstream', action='store_true', help='Requires there the be precisely 1 track in the downstream tracker' )

    parser.add_argument('--save_data_file', default=None, \
          help='Select whether to save the analysis to a plain text data file')

    parser.add_argument('--save_covariances', action="store_true", \
          help='Select whether to save the covariance matrices in text files')


  def process_args(self, namespace) :
    self.__output_filename = namespace.output_filename
    self.__output_directory = namespace.output_directory

    self.__recon_trackers = namespace.recon_trackers
    self.__recon_stations = [st-1 for st in namespace.recon_stations]
    self.__recon_planes = namespace.recon_planes

    self.__reference_station = namespace.reference_plane[0]
    self.__reference_plane = namespace.reference_plane[1]

    self.__unbias_results = namespace.unbias_results
    self.__ensemble_size = namespace.ensemble_size
    self.__save_data_file = namespace.save_data_file
    self.__save_covariances = namespace.save_covariances
    if namespace.require_upstream :
      self.__require_upstream = True
    if namespace.require_downstream :
      self.__require_downstream = True

    if self.__use_mc :
      with open(namespace.scifi_virt_dict, 'r') as infile :
        self.__mc_scifi_dictionary = json.load(infile)

    if namespace.correction_upstream is not None :
      self.__correction_upstream = numpy.loadtxt(namespace.correction_upstream)
      if self.__correction_upstream.shape != (4, 4) :
        raise ValueError('Upstream Covariance Correction must be a 4x4 matrix')
    if namespace.correction_downstream is not None :
      self.__correction_downstream = numpy.loadtxt(namespace.correction_downstream)
      if self.__correction_downstream.shape != (4, 4) :
        raise ValueError('Downstream Covariance Correction must be a 4x4 matrix')

    self.__p_correction_upstream = namespace.p_correction_upstream
    self.__p_correction_downstream = namespace.p_correction_downstream

    for tr in self.__recon_trackers :
      for st in self.__recon_stations :
        for pl in self.__recon_planes :
          plane_id = str(tr)+"."+str(st)+"."+str(pl)
          self.__inspectors[tr][st][pl] = analysis.inspectors.PhaseSpace2DInspector(plane_id, self.__ensemble_size)
          if self.__use_mc :
            self.__mc_inspectors[tr][st][pl] = analysis.inspectors.PhaseSpace2DInspector(plane_id, self.__ensemble_size)

          if tr == 0 and self.__correction_upstream is not None :
            self.__inspectors[tr][st][pl].set_covariance_correction(self.__correction_upstream)
          elif tr == 1 and self.__correction_downstream is not None :
            self.__inspectors[tr][st][pl].set_covariance_correction(self.__correction_downstream)

          if tr == 0 :
            self.__inspectors[tr][st][pl].set_momentum_correction(float(self.__p_correction_upstream[0]), float(self.__p_correction_upstream[1]))
        
          elif tr == 1 :
            self.__inspectors[tr][st][pl].set_momentum_correction(float(self.__p_correction_downstream[0]), float(self.__p_correction_downstream[1]))

    if framework.get_last_analysis_json() is not None :
      if (self.__reference_station in self.__recon_stations) and (self.__reference_plane in self.__recon_planes) :
        self.__analyse_cooling = True
        self.__cooling_inspector = analysis.inspectors.CoolingInspector()

        tracker = "tracker_0"
        station = "station_"+str(self.__reference_station)
        plane = "plane_"+str(self.__reference_plane)
        self.__cooling_inspector.set_upstream_covariance(numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tracker][station][plane]['covariance_matrix']))
        tracker = "tracker_1"
        self.__cooling_inspector.set_downstream_covariance(numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tracker][station][plane]['covariance_matrix']))

      for tr in self.__recon_trackers :
        tracker = "tracker_"+str(tr)
        for st in self.__recon_stations :
          station = "station_"+str(st+1)
          for pl in self.__recon_planes :
            plane = "plane_"+str(pl)
            
            covariance = numpy.array(framework.get_last_analysis_json()['emittance_reconstruction'][tracker][station][plane]['covariance_matrix'])
            self.__inspectors[tr][st][pl].set_parent_covariance(covariance)
              

  def _reset(self) :
    pass


  def _process(self, file_reader) :

    if self.__beam_selector.is_cut() or self.__tof_analysis.is_cut() :
      self._cut()
      return True

    if self.__unbias_results :
      random_plane = random.choice( self.__recon_planes )
      random_station = random.choice( self.__recon_stations )

    if self.__require_upstream :
      if self.__beam_selector.get_track_extractor(0).get_track() is None :
        return True
    if self.__require_downstream :
      if self.__beam_selector.get_track_extractor(1).get_track() is None :
        return True

    weight = framework.get_event_statistical_weight()

    if self.__use_mc :
      for tr in self.__recon_trackers :
        for st in self.__recon_stations :
          for pl in self.__recon_planes :
            tp = self.__beam_selector.get_track_extractor(tr).get_trackpoint_byplane(st+1, pl)
            if tp is None : continue

            hit = hit_types.AnalysisHit(scifi_track_point=tp)
            if numpy.isnan(hit.get_x()) or numpy.isinf(hit.get_x()) :
              print hit
              continue

            virt_plane = self.__mc_scifi_dictionary[str(tools.calculate_plane_id(tr, st+1, pl))][0]
            virt = self.__mc_track_extractor.get_virtual_hits(virt_plane)
            if virt is None :
              print "No Virt!"
              continue
            vhit = hit_types.AnalysisHit(virtual_track_point=virt[1])

            hit.set_weight(weight)
            self.__inspectors[tr][st][pl].add_hit(hit)

            vhit.set_weight(weight)
            self.__mc_inspectors[tr][st][pl].add_hit(vhit)

    else :
      for tr in self.__recon_trackers :
        for st in self.__recon_stations :
          for pl in self.__recon_planes :
            tp = self.__beam_selector.get_track_extractor(tr).get_trackpoint_byplane(st+1, pl)
            if tp is None : continue

            hit = hit_types.AnalysisHit(scifi_track_point=tp)
            hit.set_weight(weight)
            if numpy.isnan(hit.get_x()) or numpy.isinf(hit.get_x()) :
              print hit
            self.__inspectors[tr][st][pl].add_hit(hit)


    if self.__analyse_cooling :
      up_tp = self.__beam_selector.get_track_extractor(0).get_trackpoint_byplane(self.__reference_station, self.__reference_plane)
      down_tp = self.__beam_selector.get_track_extractor(1).get_trackpoint_byplane(self.__reference_station, self.__reference_plane)

      if up_tp is None :
        up_hit = None 
      else :
        up_hit = hit_types.AnalysisHit(scifi_track_point=up_tp)
        up_hit.set_weight(weight)
      if down_tp is None :
        down_hit = None 
      else :
        down_hit = hit_types.AnalysisHit(scifi_track_point=down_tp)
        down_hit.set_weight(weight)

      self.__cooling_inspector.add_hits(up_hit, down_hit)

    return False


  def conclude(self) :
    inspectors = [ (self.__inspectors, "data") ]
    if self.__use_mc :
      inspectors.append( (self.__mc_inspectors, "mc") )

    if self.__save_covariances :
      for tr in self.__recon_trackers :
        for st in self.__recon_stations :
          for pl in self.__recon_planes :
            name = "_"+str(tr)+"."+str(st)+"."+str(pl)+".txt"
            out_file = os.path.join( self.__output_directory, self.__output_filename+name )
            self.__inspectors[tr][st][pl].covariance.save_covariance(out_file, [ 'x', 'px', 'y', 'py' ])


    for inspector_list, name in inspectors :
      zeros = array.array('d')
      plane_ids = array.array('d')
      positions = array.array('d')
      emittance_x = array.array('d')
      emittance_x_err = array.array('d')
      emittance_y = array.array('d')
      emittance_y_err = array.array('d')
      emittance = array.array('d')
      emittance_err = array.array('d')
      alpha_x = array.array('d')
      alpha_x_err = array.array('d')
      alpha_y = array.array('d')
      alpha_y_err = array.array('d')
      alpha = array.array('d')
      alpha_err = array.array('d')
      beta_x = array.array('d')
      beta_x_err = array.array('d')
      beta_y = array.array('d')
      beta_y_err = array.array('d')
      beta = array.array('d')
      beta_err = array.array('d')
      momentum = array.array('d')
      momentum_err = array.array('d')

      for tr in self.__recon_trackers :
        for st in self.__recon_stations :
          for pl in self.__recon_planes :
            inspector_list[tr][st][pl].fill_plots() # Include the last ensemble

            num = inspector_list[tr][st][pl].get_number_ensembles()
#            if num > 0 :
#              err = math.sqrt(num) 
#            else :
#              err = 1.0

            ins_data = inspector_list[tr][st][pl].get_data_dictionary()

            zeros.append(0.0)
            plane_ids.append(tools.calculate_plane_id(tr, st+1, pl))
            positions.append(ins_data['position'])

            emittance_x.append(ins_data['emittance_x'])
            emittance_x_err.append(ins_data['emittance_x_error'])

            emittance_y.append(ins_data['emittance_y'])
            emittance_y_err.append(ins_data['emittance_y_error'])

            emittance.append(ins_data['emittance'])
            emittance_err.append(ins_data['emittance_error'])

            alpha_x.append(ins_data['alpha_x'])
            alpha_x_err.append(ins_data['alpha_x_error'])

            alpha_y.append(ins_data['alpha_y'])
            alpha_y_err.append(ins_data['alpha_y_error'])

            alpha.append(ins_data['alpha'])
            alpha_err.append(ins_data['alpha_error'])

            beta_x.append(ins_data['beta_x'])
            beta_x_err.append(ins_data['beta_x_error'])

            beta_y.append(ins_data['beta_y'])
            beta_y_err.append(ins_data['beta_y_error'])

            beta.append(ins_data['beta'])
            beta_err.append(ins_data['beta_error'])

            momentum.append(ins_data['momentum'])
            momentum_err.append(ins_data['momentum_error'])

      if len(plane_ids) > 0 :

        pos_array, emittance_x, emittance_x_err = tools.sort_arrays([positions, emittance_x, emittance_x_err])
        self.__graphs[name]['emittance_x'] = ROOT.TGraphErrors(len(zeros), pos_array, emittance_x, zeros, emittance_x_err)
        self.__graphs[name]['emittance_x'].SetTitle("Reconstructed Normalised 2D X Emittance")
        self.__graphs[name]['emittance_x'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['emittance_x'].GetYaxis().SetTitle("#epsilon_{x}^{N}  [mm]")

        pos_array, emittance_y, emittance_y_err = tools.sort_arrays([positions, emittance_y, emittance_y_err])
        self.__graphs[name]['emittance_y'] = ROOT.TGraphErrors(len(zeros), pos_array, emittance_y, zeros, emittance_y_err)
        self.__graphs[name]['emittance_y'].SetTitle("Reconstructed Normalised 2D Y Emittance")
        self.__graphs[name]['emittance_y'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['emittance_y'].GetYaxis().SetTitle("#epsilon_{y}^{N}  [mm]")

        pos_array, emittance, emittance_err = tools.sort_arrays([positions, emittance, emittance_err])
        self.__graphs[name]['emittance'] = ROOT.TGraphErrors(len(zeros), pos_array, emittance, zeros, emittance_err)
        self.__graphs[name]['emittance'].SetTitle("Reconstructed Normalised 4D Emittance")
        self.__graphs[name]['emittance'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['emittance'].GetYaxis().SetTitle("#epsilon_{4D}^{N}  [mm]")

        pos_array, alpha_x, alpha_x_err = tools.sort_arrays([positions, alpha_x, alpha_x_err])
        self.__graphs[name]['alpha_x'] = ROOT.TGraphErrors(len(zeros), pos_array, alpha_x, zeros, alpha_x_err)
        self.__graphs[name]['alpha_x'].SetTitle("Reconstructed Normalised 2D X Alpha")
        self.__graphs[name]['alpha_x'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['alpha_x'].GetYaxis().SetTitle("#alpha_{x}^{N}  [mm]")

        pos_array, alpha_y, alpha_y_err = tools.sort_arrays([positions, alpha_y, alpha_y_err])
        self.__graphs[name]['alpha_y'] = ROOT.TGraphErrors(len(zeros), pos_array, alpha_y, zeros, alpha_y_err)
        self.__graphs[name]['alpha_y'].SetTitle("Reconstructed Normalised 2D Y Alpha")
        self.__graphs[name]['alpha_y'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['alpha_y'].GetYaxis().SetTitle("#alpha_{y}^{N}  [mm]")

        pos_array, alpha, alpha_err = tools.sort_arrays([positions, alpha, alpha_err])
        self.__graphs[name]['alpha'] = ROOT.TGraphErrors(len(zeros), pos_array, alpha, zeros, alpha_err)
        self.__graphs[name]['alpha'].SetTitle("Reconstructed Normalised 4D Alpha")
        self.__graphs[name]['alpha'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['alpha'].GetYaxis().SetTitle("#alpha_{4D}^{N}  [mm]")

        pos_array, beta_x, beta_x_err = tools.sort_arrays([positions, beta_x, beta_x_err])
        self.__graphs[name]['beta_x'] = ROOT.TGraphErrors(len(zeros), pos_array, beta_x, zeros, beta_x_err)
        self.__graphs[name]['beta_x'].SetTitle("Reconstructed Normalised 2D X Beta")
        self.__graphs[name]['beta_x'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['beta_x'].GetYaxis().SetTitle("#beta_{x}^{N}  [mm]")

        pos_array, beta_y, beta_y_err = tools.sort_arrays([positions, beta_y, beta_y_err])
        self.__graphs[name]['beta_y'] = ROOT.TGraphErrors(len(zeros), pos_array, beta_y, zeros, beta_y_err)
        self.__graphs[name]['beta_y'].SetTitle("Reconstructed Normalised 2D Y Beta")
        self.__graphs[name]['beta_y'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['beta_y'].GetYaxis().SetTitle("#beta_{y}^{N}  [mm]")

        pos_array, beta, beta_err = tools.sort_arrays([positions, beta, beta_err])
        self.__graphs[name]['beta'] = ROOT.TGraphErrors(len(zeros), pos_array, beta, zeros, beta_err)
        self.__graphs[name]['beta'].SetTitle("Reconstructed Normalised 4D Beta")
        self.__graphs[name]['beta'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['beta'].GetYaxis().SetTitle("#beta_{4D}^{N}  [mm]")

        pos_array, momentum, momentum_err = tools.sort_arrays([positions, momentum, momentum_err])
        self.__graphs[name]['momentum'] = ROOT.TGraphErrors(len(zeros), pos_array, momentum, zeros, momentum_err)
        self.__graphs[name]['momentum'].SetTitle("Reconstructed Total Momentum")
        self.__graphs[name]['momentum'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs[name]['momentum'].GetYaxis().SetTitle("p  [MeV]")
        

    if self.__use_mc :
      zeros = array.array('d')
      plane_ids = array.array('d')
      positions = array.array('d')
      emittance_x = array.array('d')
      emittance_x_err = array.array('d')
      emittance_y = array.array('d')
      emittance_y_err = array.array('d')
      emittance = array.array('d')
      emittance_err = array.array('d')
      alpha_x = array.array('d')
      alpha_x_err = array.array('d')
      alpha_y = array.array('d')
      alpha_y_err = array.array('d')
      alpha = array.array('d')
      alpha_err = array.array('d')
      beta_x = array.array('d')
      beta_x_err = array.array('d')
      beta_y = array.array('d')
      beta_y_err = array.array('d')
      beta = array.array('d')
      beta_err = array.array('d')
      momentum = array.array('d')
      momentum_err = array.array('d')

      for tr in self.__recon_trackers :
        for st in self.__recon_stations :
          for pl in self.__recon_planes :

            data = self.__inspectors[tr][st][pl].get_data_dictionary()
            mc = self.__mc_inspectors[tr][st][pl].get_data_dictionary()

            zeros.append(0.0)
            plane_ids.append(tools.calculate_plane_id(tr, st+1, pl))
            positions.append(data['position'])

            emittance_x.append(data['emittance_x'] - mc['emittance_x'])
            emittance_x_err.append(data['emittance_x_error'])

            emittance_y.append(data['emittance_y'] - mc['emittance_y'])
            emittance_y_err.append(data['emittance_y_error'])

            emittance.append(data['emittance'] - mc['emittance'])
            emittance_err.append(data['emittance_error'])

            alpha_x.append(data['alpha_x'] - mc['alpha_x'])
            alpha_x_err.append(data['alpha_x_error'])

            alpha_y.append(data['alpha_y'] - mc['alpha_y'])
            alpha_y_err.append(data['alpha_y_error'])

            alpha.append(data['alpha'] - mc['alpha'])
            alpha_err.append(data['alpha_error'])

            beta_x.append(data['beta_x'] - mc['beta_x'])
            beta_x_err.append(data['beta_x_error'])

            beta_y.append(data['beta_y'] - mc['beta_y'])
            beta_y_err.append(data['beta_y_error'])

            beta.append(data['beta'] - mc['beta'])
            beta_err.append(data['beta_error'])

            momentum.append(data['momentum'] - mc['momentum'])
            momentum_err.append(data['momentum_error'])

      if len(plane_ids) > 0 :
        pos_array, emittance_x, emittance_x_err = tools.sort_arrays([positions, emittance_x, emittance_x_err])
        self.__graphs["residuals"]['emittance_x'] = ROOT.TGraphErrors(len(zeros), pos_array, emittance_x, zeros, emittance_x_err)
        self.__graphs["residuals"]['emittance_x'].SetTitle("Reconstructed Normalised 2D X Emittance")
        self.__graphs["residuals"]['emittance_x'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['emittance_x'].GetYaxis().SetTitle("#epsilon_{x}^{N}  [mm]")

        pos_array, emittance_y, emittance_y_err = tools.sort_arrays([positions, emittance_y, emittance_y_err])
        self.__graphs["residuals"]['emittance_y'] = ROOT.TGraphErrors(len(zeros), pos_array, emittance_y, zeros, emittance_y_err)
        self.__graphs["residuals"]['emittance_y'].SetTitle("Reconstructed Normalised 2D Y Emittance")
        self.__graphs["residuals"]['emittance_y'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['emittance_y'].GetYaxis().SetTitle("#epsilon_{y}^{N}  [mm]")

        pos_array, emittance, emittance_err = tools.sort_arrays([positions, emittance, emittance_err])
        self.__graphs["residuals"]['emittance'] = ROOT.TGraphErrors(len(zeros), pos_array, emittance, zeros, emittance_err)
        self.__graphs["residuals"]['emittance'].SetTitle("Reconstructed Normalised 4D Emittance")
        self.__graphs["residuals"]['emittance'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['emittance'].GetYaxis().SetTitle("#epsilon_{4D}^{N}  [mm]")

        pos_array, alpha_x, alpha_x_err = tools.sort_arrays([positions, alpha_x, alpha_x_err])
        self.__graphs["residuals"]['alpha_x'] = ROOT.TGraphErrors(len(zeros), pos_array, alpha_x, zeros, alpha_x_err)
        self.__graphs["residuals"]['alpha_x'].SetTitle("Reconstructed Normalised 2D X Alpha")
        self.__graphs["residuals"]['alpha_x'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['alpha_x'].GetYaxis().SetTitle("#alpha_{x}^{N}  [mm]")

        pos_array, alpha_y, alpha_y_err = tools.sort_arrays([positions, alpha_y, alpha_y_err])
        self.__graphs["residuals"]['alpha_y'] = ROOT.TGraphErrors(len(zeros), pos_array, alpha_y, zeros, alpha_y_err)
        self.__graphs["residuals"]['alpha_y'].SetTitle("Reconstructed Normalised 2D Y Alpha")
        self.__graphs["residuals"]['alpha_y'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['alpha_y'].GetYaxis().SetTitle("#alpha_{y}^{N}  [mm]")

        pos_array, alpha, alpha_err = tools.sort_arrays([positions, alpha, alpha_err])
        self.__graphs["residuals"]['alpha'] = ROOT.TGraphErrors(len(zeros), pos_array, alpha, zeros, alpha_err)
        self.__graphs["residuals"]['alpha'].SetTitle("Reconstructed Normalised 4D Alpha")
        self.__graphs["residuals"]['alpha'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['alpha'].GetYaxis().SetTitle("#alpha_{4D}^{N}  [mm]")

        pos_array, beta_x, beta_x_err = tools.sort_arrays([positions, beta_x, beta_x_err])
        self.__graphs["residuals"]['beta_x'] = ROOT.TGraphErrors(len(zeros), pos_array, beta_x, zeros, beta_x_err)
        self.__graphs["residuals"]['beta_x'].SetTitle("Reconstructed Normalised 2D X Beta")
        self.__graphs["residuals"]['beta_x'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['beta_x'].GetYaxis().SetTitle("#beta_{x}^{N}  [mm]")

        pos_array, beta_y, beta_y_err = tools.sort_arrays([positions, beta_y, beta_y_err])
        self.__graphs["residuals"]['beta_y'] = ROOT.TGraphErrors(len(zeros), pos_array, beta_y, zeros, beta_y_err)
        self.__graphs["residuals"]['beta_y'].SetTitle("Reconstructed Normalised 2D Y Beta")
        self.__graphs["residuals"]['beta_y'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['beta_y'].GetYaxis().SetTitle("#beta_{y}^{N}  [mm]")

        pos_array, beta, beta_err = tools.sort_arrays([positions, beta, beta_err])
        self.__graphs["residuals"]['beta'] = ROOT.TGraphErrors(len(zeros), pos_array, beta, zeros, beta_err)
        self.__graphs["residuals"]['beta'].SetTitle("Reconstructed Normalised 4D Beta")
        self.__graphs["residuals"]['beta'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['beta'].GetYaxis().SetTitle("#beta_{4D}^{N}  [mm]")

        pos_array, momentum, momentum_err = tools.sort_arrays([positions, momentum, momentum_err])
        self.__graphs["residuals"]['momentum'] = ROOT.TGraphErrors(len(zeros), pos_array, momentum, zeros, momentum_err)
        self.__graphs["residuals"]['momentum'].SetTitle("Reconstructed Total Momentum")
        self.__graphs["residuals"]['momentum'].GetXaxis().SetTitle("z  [mm]")
        self.__graphs["residuals"]['momentum'].GetYaxis().SetTitle("p  [MeV]")



  def _store_plots(self, plot_dict) :
    tracker_plots = {}
    for tr in self.__recon_trackers :
      station_plots = {}
      for st in self.__recon_stations :
        plane_plots = {}
        for pl in self.__recon_planes :
          plane_plots['plane_'+str(pl)] = self.__inspectors[tr][st][pl].get_plot_dictionary()
        station_plots['station_'+str(st+1)] = plane_plots
      tracker_plots['tracker_'+str(tr)] = station_plots

    if self.__graphs['data']['emittance'] is not None :
      tracker_plots['emittance_x'] = self.__graphs['data']['emittance_x']
      tracker_plots['emittance_y'] = self.__graphs['data']['emittance_y']
      tracker_plots['emittance'] = self.__graphs['data']['emittance']
      tracker_plots['alpha_x'] = self.__graphs['data']['alpha_x']
      tracker_plots['alpha_y'] = self.__graphs['data']['alpha_y']
      tracker_plots['alpha'] = self.__graphs['data']['alpha']
      tracker_plots['beta_x'] = self.__graphs['data']['beta_x']
      tracker_plots['beta_y'] = self.__graphs['data']['beta_y']
      tracker_plots['beta'] = self.__graphs['data']['beta']
      tracker_plots['momentum'] = self.__graphs['data']['momentum']


    if self.__use_mc :
      mc_tracker_plots = {}
      mc_residual_plots = {}
      for tr in self.__recon_trackers :
        station_plots = {}
        for st in self.__recon_stations :
          plane_plots = {}
          for pl in self.__recon_planes :
            plane_plots['plane_'+str(pl)] = self.__mc_inspectors[tr][st][pl].get_plot_dictionary()
          station_plots['station_'+str(st+1)] = plane_plots
        mc_tracker_plots['tracker_'+str(tr)] = station_plots

      if self.__graphs['mc']['emittance'] is not None :
        mc_tracker_plots['emittance_x'] = self.__graphs['mc']['emittance_x']
        mc_tracker_plots['emittance_y'] = self.__graphs['mc']['emittance_y']
        mc_tracker_plots['emittance'] = self.__graphs['mc']['emittance']
        mc_tracker_plots['alpha_x'] = self.__graphs['mc']['alpha_x']
        mc_tracker_plots['alpha_y'] = self.__graphs['mc']['alpha_y']
        mc_tracker_plots['alpha'] = self.__graphs['mc']['alpha']
        mc_tracker_plots['beta_x'] = self.__graphs['mc']['beta_x']
        mc_tracker_plots['beta_y'] = self.__graphs['mc']['beta_y']
        mc_tracker_plots['beta'] = self.__graphs['mc']['beta']
        mc_tracker_plots['momentum'] = self.__graphs['mc']['momentum']

        mc_residual_plots['emittance_x'] = self.__graphs['residuals']['emittance_x']
        mc_residual_plots['emittance_y'] = self.__graphs['residuals']['emittance_y']
        mc_residual_plots['emittance'] = self.__graphs['residuals']['emittance']
        mc_residual_plots['alpha_x'] = self.__graphs['residuals']['alpha_x']
        mc_residual_plots['alpha_y'] = self.__graphs['residuals']['alpha_y']
        mc_residual_plots['alpha'] = self.__graphs['residuals']['alpha']
        mc_residual_plots['beta_x'] = self.__graphs['residuals']['beta_x']
        mc_residual_plots['beta_y'] = self.__graphs['residuals']['beta_y']
        mc_residual_plots['beta'] = self.__graphs['residuals']['beta']
        mc_residual_plots['momentum'] = self.__graphs['residuals']['momentum']



    plot_dict['emittance_reconstruction'] = tracker_plots
    if self.__use_mc :
      plot_dict['mc_emittance_reconstruction'] = mc_tracker_plots
      plot_dict['mc_recon_comparisons'] = mc_residual_plots

    if self.__analyse_cooling :
      plot_dict['cooling_analysis'] = self.__cooling_inspector.get_plot_dictionary()

    return plot_dict


  def _store_data(self, data_dict) :
    tracker_data = {}
    for tr in self.__recon_trackers :
      station_data = {}
      for st in self.__recon_stations :
        plane_data = {}
        for pl in self.__recon_planes :
          data = self.__inspectors[tr][st][pl].get_data_dictionary()

          plane_data['plane_'+str(pl)] = data
        station_data['station_'+str(st+1)] = plane_data
      tracker_data['tracker_'+str(tr)] = station_data

    data_dict['emittance_reconstruction'] = tracker_data


    if self.__use_mc :
      tracker_data = {}
      for tr in self.__recon_trackers :
        station_data = {}
        for st in self.__recon_stations :
          plane_data = {}
          for pl in self.__recon_planes :
            data = self.__mc_inspectors[tr][st][pl].get_data_dictionary()

            plane_data['plane_'+str(pl)] = data
          station_data['station_'+str(st+1)] = plane_data
        tracker_data['tracker_'+str(tr)] = station_data

      data_dict['mc_emittance_reconstruction'] = tracker_data


################################################################################

