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
Alignment Analysis Classes in Here
"""

import ROOT
import math
import json

import framework
import scifi_extractors
import scifi_analysis
import tof_analysis
import analysis.tools as tools


class scifi_straight_track_alignment(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_straight_track_alignment")

    self.__number_tracks = 0
    self.__non_linear = False

    self.__track_candidates = None
    self.__tof_analysis = None

    self.__default_tracker_separation = -1.0
    self.__default_tracker_alignment = [0.0, 0.0, 0.0, 0.0]

    self.__projection_radius = 400.0
    self.__projection_cut = 0

    self.__x_delta_x = ROOT.TH2F('x_delta_x', \
              "X Against Change in X", 200, -200.0, 200.0, 200, -200.0, 200.0 )
    self.__y_delta_y = ROOT.TH2F('y_delta_y', \
              "Y Against Change in Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
    self.__x_delta_y = ROOT.TH2F('x_delta_y', \
              "X Against Change in Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
    self.__y_delta_x = ROOT.TH2F('y_delta_x', \
              "Y Against Change in X", 200, -200.0, 200.0, 200, -200.0, 200.0 )

    self.__projected_x_x = ROOT.TH2F('projected_x_x', \
         "Delta X Against Residual X", 200, -200.0, 200.0, 200, -200.0, 200.0 )
    self.__projected_y_y = ROOT.TH2F('projected_y_y', \
         "Delta Y Against Residual Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
    self.__projected_x_y = ROOT.TH2F('projected_x_y', \
         "Delta X Against Residual Y", 200, -200.0, 200.0, 200, -200.0, 200.0 )
    self.__projected_y_x = ROOT.TH2F('projected_y_x', \
         "Delta Y Against Residual X", 200, -200.0, 200.0, 200, -200.0, 200.0 )

    self.__phi_projected_x_x = ROOT.TH2F('phi_projected_x_x', \
                 "Delta X Against Residual X", 200, -0.2, 0.2, 200, -0.2, 0.2 )
    self.__phi_projected_y_y = ROOT.TH2F('phi_projected_y_y', \
                 "Delta Y Against Residual Y", 200, -0.2, 0.2, 200, -0.2, 0.2 )
    self.__phi_projected_x_y = ROOT.TH2F('phi_projected_x_y', \
                 "Delta X Against Residual Y", 200, -0.2, 0.2, 200, -0.2, 0.2 )
    self.__phi_projected_y_x = ROOT.TH2F('phi_projected_y_x', \
                 "Delta Y Against Residual X", 200, -0.2, 0.2, 200, -0.2, 0.2 )

    self.__tracker_separation = ROOT.TH1F('tracker_separation', \
                        "Separation of the two trackers", 10000, 0.0, 10000.0 )
    self.__position_gradient_x = ROOT.TH2F('position_gradient_x', \
           "Change in X Against Gradient", 500, -1.0, 1.0, 500, -500.0, 500.0 )
    self.__position_gradient_y = ROOT.TH2F('position_gradient_y', \
           "Change in Y Against Gradient", 500, -1.0, 1.0, 500, -500.0, 500.0 )
    self.__position_gradient = ROOT.TH2F('position_gradient', \
     "Change in Position Against Gradient", 100, 0.0, 0.5, 100, -300.0, 300.0 )


    self.__residual_xy = ROOT.TH2F('residual_xy', \
#               'Residuals in Position', 100, -400.0, 400.0, 100, -400.0, 400.0)
               'Residuals in Position', 200, -200.0, 200.0, 200, -200.0, 200.0)
    self.__residual_mxmy = ROOT.TH2F('residual_mxmy', \
#                       'Residuals in Gradient', 100, -0.2, 0.2, 100, -0.2, 0.2)
                       'Residuals in Gradient', 200, -0.2, 0.2, 200, -0.2, 0.2)
    self.__residual_txty = ROOT.TH2F('residual_txty', \
#                          'Residuals in Angle', 100, -0.2, 0.2, 100, -0.2, 0.2)
                          'Residuals in Angle', 200, -0.2, 0.2, 200, -0.2, 0.2)
    self.__residual_theta = ROOT.TH1F('residual_theta', \
                                     'Residual Beam Rotation', 320, -8.0, 8.0) 
    self.__residual_xmy = ROOT.TH2F('residual_xmy', \
                   'Residual Beam X-My', 200, -200.0, 200.0, 200, -0.2, 0.2 )
    self.__residual_ymx = ROOT.TH2F('residual_ymx', \
                   'Residual Beam Y-Mx', 200, -200.0, 200.0, 200, -0.2, 0.2 )
   
    self.__x_dist_profile = None
    self.__y_dist_profile = None
    self.__dist_profile = None

    self.__alignment_dict = {}
    self.__alignment_dict['x_delta'] = 0.0
    self.__alignment_dict['x_delta_err'] = 0.0
    self.__alignment_dict['y_delta'] = 0.0
    self.__alignment_dict['y_delta_err'] = -1.0
    self.__alignment_dict['z_delta'] = 0.0
    self.__alignment_dict['z_delta_err'] = -1.0
    self.__alignment_dict['x_delta_phi'] = 0.0
    self.__alignment_dict['x_delta_phi_err'] = 0.0
    self.__alignment_dict['y_delta_phi'] = 0.0
    self.__alignment_dict['y_delta_phi_err'] = 0.0
    self.__alignment_dict['delta_theta'] = 0.0
    self.__alignment_dict['delta_theta_err'] = -1.0
    self.__alignment_dict['tracker_separation'] = -1.0
    self.__alignment_dict['tracker_separation_err'] = 0.0

    self.__up_corr_x = 0.0
    self.__up_corr_y = 0.0
    self.__up_corr_z = 0.0
    self.__up_corr_tx = 0.0
    self.__up_corr_ty = 0.0
    self.__up_corr_theta = 0.0

    self.__down_corr_x = 0.0
    self.__down_corr_y = 0.0
    self.__down_corr_z = 0.0
    self.__down_corr_tx = 0.0
    self.__down_corr_ty = 0.0
    self.__down_corr_theta = 0.0

    self.__corr_sine_x = 0.0
    self.__corr_sine_y = 0.0
    self.__corr_cosine_x = 1.0
    self.__corr_cosine_y = 1.0


  def add_correction_up(self, pos_x, pos_y, pos_z, angle_x, angle_y, theta) :
    self.__up_corr_x += pos_x
    self.__up_corr_y += pos_y
    self.__up_corr_z += pos_z
    self.__up_corr_tx += angle_x
    self.__up_corr_ty += angle_y
    self.__up_corr_theta += theta


  def add_correction_down(self, pos_x, pos_y, pos_z, angle_x, angle_y, theta) :
    self.__down_corr_x += pos_x
    self.__down_corr_y += pos_y
    self.__down_corr_z += pos_z
    self.__down_corr_tx += angle_x
    self.__down_corr_ty += angle_y
    self.__down_corr_theta += theta

    self.__corr_sine_x = math.sin(self.__down_corr_tx)
    self.__corr_sine_y = math.sin(self.__down_corr_ty)
    self.__corr_cosine_x = math.cos(self.__down_corr_tx)
    self.__corr_cosine_y = math.cos(self.__down_corr_ty)


  def set_correction_up(self, pos_x, pos_y, pos_z, angle_x, angle_y, theta) :
    self.__up_corr_x = pos_x
    self.__up_corr_y = pos_y
    self.__up_corr_z = pos_z
    self.__up_corr_tx = angle_x
    self.__up_corr_ty = angle_y
    self.__up_corr_theta = theta


  def set_correction_down(self, pos_x, pos_y, pos_z, angle_x, angle_y, theta) :
    self.__down_corr_x = pos_x
    self.__down_corr_y = pos_y
    self.__down_corr_z = pos_z
    self.__down_corr_tx = angle_x
    self.__down_corr_ty = angle_y
    self.__down_corr_theta = theta

    self.__corr_sine_x = math.sin(self.__down_corr_tx)
    self.__corr_sine_y = math.sin(self.__down_corr_ty)
    self.__corr_cosine_x = math.cos(self.__down_corr_tx)
    self.__corr_cosine_y = math.cos(self.__down_corr_ty)


  def get_dependencies(self, inserter) :
    self.__track_candidates = inserter(scifi_extractors.scifi_straight_track_candidates())
    self.__tof_analysis = inserter(tof_analysis.tof_analyser())


  def get_args(self, parser) :
    parser.add_argument( '--tracker_separation', type=float, default=-1.0, \
                      help='Manually set the separation between the trackers' )
    parser.add_argument( '--alignment_data_file', type=str, default=None, \
                              help='File containing data from a previous run' )
    parser.add_argument( '--non_linear_alignment', action='store_true', \
               help='Turn on the second order corrections to the calculation' )
    parser.add_argument( '--projection_cut', type=float, default=400.0, \
                      help='Impose a cut on the radius of the projected beam' )
    parser.add_argument( '--tracker_placement', type=float, nargs=4, \
                                                default=[0.0, 0.0, 0.0, 0.0], \
     help='Specify the known offset in transverse separation and rotation ' + \
                'of the trackers, to be included in the alignemnt values. ' + \
                             'Format: delta_x delta_y delta_phi_x delta_phi_y')


  def process_args(self, namespace) :
    self.__default_tracker_separation = namespace.tracker_separation
    self.__default_tracker_alignment = namespace.tracker_placement
    self.__non_linear = namespace.non_linear_alignment

    self.__projection_radius = namespace.projection_cut
    
    alignment_data = namespace.alignment_data_file

    if alignment_data is not None :
      with open(alignment_data, 'r') as infile :
        data_dict = json.load(infile)

        self.__alignment_dict = data_dict['alignment']
        self.set_correction_down(self.__alignment_dict['x_delta'], \
            self.__alignment_dict['y_delta'], \
            self.__alignment_dict['z_delta'], \
            self.__alignment_dict['x_delta_phi'], \
            self.__alignment_dict['y_delta_phi'], \
            self.__alignment_dict['delta_theta'] )


  def _reset(self) :
    pass


  def _process(self, file_reader) :

    if self.__track_candidates.is_cut() :
      return True

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

    up_tp = self.__track_candidates.get_upstream_reference()
    up_z = up_tp.pos().z()
    up_pos = [ (up_tp.pos().x() + self.__up_corr_x), \
                                         (up_tp.pos().y() + self.__up_corr_y) ]
#    up_ang = [ \
#          (math.atan2(-up_tp.mom().x(), up_tp.mom().z()) + self.__up_corr_tx), \
#          (math.atan2(-up_tp.mom().y(), up_tp.mom().z()) + self.__up_corr_ty) ]
    up_ang = [ \
          (math.atan(up_tp.gradient().x()) + self.__up_corr_tx), \
          (math.atan(up_tp.gradient().y()) + self.__up_corr_ty) ]

    up_gra = [ math.tan(up_ang[0]), math.tan(up_ang[1]) ]
    up_theta = math.atan2(up_gra[1], up_gra[0])


    down_tp = self.__track_candidates.get_downstream_reference()
    down_z = down_tp.pos().z()
    raw_down_pos = [down_tp.pos().x(), down_tp.pos().y()]

    if self.__non_linear :
      correction_x = raw_down_pos[0]*self.__corr_sine_x*math.tan(up_ang[0])
      correction_y = raw_down_pos[1]*self.__corr_sine_y*math.tan(up_ang[1])

      down_pos = [ ((raw_down_pos[0]*self.__corr_cosine_x) + self.__down_corr_x), \
                   ((raw_down_pos[1]*self.__corr_cosine_y) + self.__down_corr_y) ]
#      down_pos = [ ((raw_down_pos[0]*self.__corr_cosine_x) + self.__down_corr_x) + correction_x, \
#                   ((raw_down_pos[1]*self.__corr_cosine_y) + self.__down_corr_y) + correction_y ]
    else :
      down_pos = [ (down_tp.pos().x() + self.__down_corr_x), \
                                     (down_tp.pos().y() + self.__down_corr_y) ]

    down_ang = [ \
    (math.atan2(down_tp.mom().x(), down_tp.mom().z()) + self.__down_corr_tx), \
    (math.atan2(down_tp.mom().y(), down_tp.mom().z()) + self.__down_corr_ty) ]

#    down_pos[0] = down_pos[0] + self.__default_tracker_alignment[0]
#    down_pos[1] = down_pos[1] + self.__default_tracker_alignment[1]
#    down_ang[0] = down_ang[0] + self.__default_tracker_alignment[2]
#    down_ang[1] = down_ang[1] + self.__default_tracker_alignment[3]

    down_gra = [ math.tan(down_ang[0]), math.tan(down_ang[1]) ]
    down_theta = math.atan2(down_gra[1], down_gra[0])


    if self.__default_tracker_separation > 0.0 :
      length = self.__default_tracker_separation
    else :
      length = down_z - up_z

    if down_theta - up_theta > math.pi :
      delta_theta = down_theta - up_theta - 2.0*math.pi
    elif down_theta - up_theta < -math.pi :
      delta_theta = down_theta - up_theta + 2.0*math.pi
    else :
      delta_theta = down_theta - up_theta

    projected_pos = [ (up_pos[0] + length*up_gra[0]), \
                                               (up_pos[1] + length*up_gra[1]) ]
# WARNING HACK!
#    projected_pos = [ (up_pos[0] - length*up_gra[0]), \
#                                               (up_pos[1] - length*up_gra[1]) ]
    projected_gra = [ up_gra[0], up_gra[1] ]
    projected_ang = [ up_ang[0], up_ang[1] ]

    proj_radius = math.sqrt( projected_pos[0]**2 + projected_pos[1]**2)
    if proj_radius > self.__projection_radius :
      self.__projection_cut += 1
      return True



    self.__tracker_separation.Fill(length)
    self.__position_gradient_x.Fill(up_gra[0], down_pos[0] - projected_pos[0])
    self.__position_gradient_y.Fill(up_gra[1], down_pos[1] - projected_pos[1])

    grad = math.sqrt( up_gra[0]**2 + up_gra[1]**2 )
    delta_pos = [ down_pos[0]-up_pos[0], down_pos[1]-up_pos[1] ]

    delta = (delta_pos[0]*up_gra[0] + delta_pos[1]*up_gra[1]) / grad
    self.__position_gradient.Fill(grad, delta)

    self.__x_delta_x.Fill(up_pos[0], down_pos[0]-up_pos[0])
    self.__y_delta_y.Fill(up_pos[1], down_pos[1]-up_pos[1])
    self.__x_delta_y.Fill(up_pos[0], down_pos[1]-up_pos[1])
    self.__y_delta_x.Fill(up_pos[1], down_pos[0]-up_pos[0])

    self.__projected_x_x.Fill( projected_pos[0], down_pos[0] )
    self.__projected_y_y.Fill( projected_pos[1], down_pos[1] )
    self.__projected_x_y.Fill( projected_pos[0], down_pos[1] )
    self.__projected_y_x.Fill( projected_pos[1], down_pos[0] )

    self.__phi_projected_x_x.Fill( projected_ang[0], down_ang[0] )
    self.__phi_projected_y_y.Fill( projected_ang[1], down_ang[1] )
    self.__phi_projected_x_y.Fill( projected_ang[0], down_ang[1] )
    self.__phi_projected_y_x.Fill( projected_ang[1], down_ang[0] )

    self.__residual_xy.Fill(projected_pos[0] - down_pos[0], \
                                                projected_pos[1] - down_pos[1])
    self.__residual_mxmy.Fill(projected_gra[0] - down_gra[0], \
                                                projected_gra[1] - down_gra[1])
    self.__residual_txty.Fill(projected_ang[0] - down_ang[0], \
                                                projected_ang[1] - down_ang[1])
    self.__residual_theta.Fill(delta_theta)
    self.__residual_xmy.Fill(up_pos[0], up_gra[1])
    self.__residual_ymx.Fill(up_pos[1], up_gra[0])

    self.__number_tracks += 1

    return False


  def conclude(self) :
    self.__x_dist_profile = self.__position_gradient_x.ProfileX('x_dist_profile')
    self.__y_dist_profile = self.__position_gradient_y.ProfileX('y_dist_profile')

    self.__dist_profile = tools.gaussian_profile_x(self.__position_gradient, 0.0, 150.0 )

#    x_fit_result = self.__x_dist_profile.Fit('pol1', 'OS', "", -0.01, 0.01)
#    y_fit_result = self.__y_dist_profile.Fit('pol1', 'OS', "", -0.01, 0.01)
#    fit_result = self.__dist_profile.Fit('pol1', 'OS', "", -0.0, 0.050)
#
#    x_fit_length = x_fit_result.Parameter(1)
#    x_fit_length_err = x_fit_result.ParError(1)
#    y_fit_length = y_fit_result.Parameter(1)
#    y_fit_length_err = y_fit_result.ParError(1)
#
#    fit_length = fit_result.Parameter(1)
#    fit_length_err = fit_result.ParError(1)

    self.__x_dist_profile.Fit('pol1', 'QOS', "", -0.01, 0.01)
    self.__y_dist_profile.Fit('pol1', 'QOS', "", -0.01, 0.01)
    self.__dist_profile.Fit('pol1', 'QOS', "", -0.0, 0.050)

    x_fit_result = self.__x_dist_profile.GetFunction('pol1')
    y_fit_result = self.__y_dist_profile.GetFunction('pol1')
    fit_result = self.__dist_profile.GetFunction('pol1')

    if not x_fit_result : 
      print "ERROR Could not perform x-fit"
      x_fit_length = 0.0
      x_fit_length_err = 0.0
    else :
      x_fit_length = x_fit_result.GetParameter(1)
      x_fit_length_err = x_fit_result.GetParError(1)

    if not y_fit_result : 
      print "ERROR Could not perform y-fit"
      y_fit_length = 0.0
      y_fit_length_err = 0.0
    else :
      y_fit_length = y_fit_result.GetParameter(1)
      y_fit_length_err = y_fit_result.GetParError(1)

    if not fit_result : 
      print "ERROR Could not perform fit"
      fit_length = 0.0
      fit_length_err = 0.0
    else :
      fit_length = fit_result.GetParameter(1)
      fit_length_err = fit_result.GetParError(1)


    length = self.__tracker_separation.GetMean()
    num_tracks = self.__tracker_separation.GetEntries()

    self.__x_residual_hist = self.__residual_xy.ProjectionX("x_residual_profile")
    self.__y_residual_hist = self.__residual_xy.ProjectionY("y_residual_profile")

    x_delta_raw, x_delta_raw_err, std, std_err = tools.fit_gaussian(self.__x_residual_hist)
    y_delta_raw, y_delta_raw_err, std, std_err = tools.fit_gaussian(self.__y_residual_hist)

    self.__tx_residual_hist = self.__residual_txty.ProjectionX("tx_residual_profile")
    self.__ty_residual_hist = self.__residual_txty.ProjectionY("ty_residual_profile")

    x_delta_phi, x_delta_phi_err, std, std_err = tools.fit_gaussian(self.__tx_residual_hist)
    y_delta_phi, y_delta_phi_err, std, std_err = tools.fit_gaussian(self.__ty_residual_hist)

    tracker_separation = (x_fit_length + y_fit_length) / 2.0
    tracker_separation_err = math.sqrt(x_fit_length_err**2 \
                                                   + y_fit_length_err**2) / 2.0

    self.__alignment_dict['distance'] = length

    self.__alignment_dict['x_delta'] += x_delta_raw
    self.__alignment_dict['y_delta'] += y_delta_raw
    self.__alignment_dict['x_delta_err'] = x_delta_raw_err
    self.__alignment_dict['y_delta_err'] = y_delta_raw_err

    self.__alignment_dict['x_delta_phi'] += x_delta_phi
    self.__alignment_dict['y_delta_phi'] += y_delta_phi
    self.__alignment_dict['x_delta_phi_err'] = x_delta_phi_err
    self.__alignment_dict['y_delta_phi_err'] = y_delta_phi_err

    self.__alignment_dict['tracker_separation'] = tracker_separation
    self.__alignment_dict['tracker_separation_err'] = tracker_separation_err

    self.__alignment_dict['x_fit_length'] = x_fit_length
    self.__alignment_dict['x_fit_length_err'] = x_fit_length_err
    self.__alignment_dict['y_fit_length'] = y_fit_length
    self.__alignment_dict['y_fit_length_err'] = y_fit_length_err


  def _store_plots(self, plot_dict) :
    correlation_plots = {}
    separation_plots = {}
    residual_plots = {}

    correlation_plots['x_delta_x'] = self.__x_delta_x
    correlation_plots['y_delta_y'] = self.__y_delta_y
    correlation_plots['x_delta_y'] = self.__x_delta_y
    correlation_plots['y_delta_x'] = self.__y_delta_x

    correlation_plots['projected_x_x'] = self.__projected_x_x
    correlation_plots['projected_y_y'] = self.__projected_y_y
    correlation_plots['projected_x_y'] = self.__projected_x_y
    correlation_plots['projected_y_x'] = self.__projected_y_x

    correlation_plots['phi_projected_x_x'] = self.__phi_projected_x_x
    correlation_plots['phi_projected_y_y'] = self.__phi_projected_y_y
    correlation_plots['phi_projected_x_y'] = self.__phi_projected_x_y
    correlation_plots['phi_projected_y_x'] = self.__phi_projected_y_x

    separation_plots['tracker_separation'] = self.__tracker_separation
    separation_plots['position_gradient_x'] = self.__position_gradient_x
    separation_plots['position_gradient_y'] = self.__position_gradient_y
    separation_plots['position_gradient'] = self.__position_gradient

    residual_plots['xy'] = self.__residual_xy
    residual_plots['mxmy'] = self.__residual_mxmy
    residual_plots['txty'] = self.__residual_txty
    residual_plots['theta'] = self.__residual_theta
    residual_plots['xmy'] = self.__residual_xmy
    residual_plots['ymx'] = self.__residual_ymx
    residual_plots['x'] = self.__x_residual_hist
    residual_plots['y'] = self.__y_residual_hist
    residual_plots['tx'] = self.__tx_residual_hist
    residual_plots['ty'] = self.__ty_residual_hist

    plot_dict['correlations'] = correlation_plots
    plot_dict['tracker_separation'] = separation_plots
    plot_dict['residuals'] = residual_plots
    return plot_dict


  def _store_data(self, data_dict) :
    data_dict['alignment'] = self.__alignment_dict
    data_dict['alignment']['projection_cut'] = self.__projection_cut
    data_dict['alignment']['number_tracks'] = self.__number_tracks

    return data_dict


