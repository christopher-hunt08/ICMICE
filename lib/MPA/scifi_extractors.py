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
  SciFi Extractor Classes are stored here.
"""


class scifi_straight_track_extractor(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_straight_track_extractor")

    self.__upstream_tracks = []
    self.__downstream_tracks = []

    self.__counter_upstream_tracks = 0
    self.__counter_downstream_tracks = 0

    self.__track_distribution = ROOT.TH2F('ste_track_distribution', \
                 "Number of tracks in eac tracker", 6, -0.5, 5.5, 6, -0.5, 5.5)

  def get_upstream_tracks(self) :
    return self.__upstream_tracks


  def get_downstream_tracks(self) :
    return self.__downstream_tracks


  def get_dependencies(self, inserter) :
    pass


  def _reset(self) :
    self.__upstream_tracks = []
    self.__downstream_tracks = []


  def _process(self, file_reader) :
    scifi_event = file_reader.get_event('scifi')
    scifi_tracks = scifi_event.scifitracks()

    event_up_tracks = 0
    event_down_tracks = 0

    for track in scifi_tracks :
      if track.GetAlgorithmUsed() != tools.STRAIGHT_ALGORITHM_ID :
        continue

      if track.tracker() == 0 :
        self.__upstream_tracks.append(track)
        self.__counter_upstream_tracks += 1
        event_up_tracks += 1
      elif track.tracker() == 1 :
        self.__downstream_tracks.append(track)
        self.__counter_downstream_tracks += 1
        event_down_tracks += 1

    self.__track_distribution.Fill(event_up_tracks, event_down_tracks)
    return False


  def _store_plots(self, plot_dict) :
    temp_dict = {}
    temp_dict['track_distribution'] = self.__track_distribution

    plot_dict['scifi_straight_track_extractor'] = temp_dict
    return plot_dict


  def _store_data(self, data_dict) :
    temp_dict = {}
    temp_dict['number_events'] = self.get_number_events()
    temp_dict['upstream_tracks'] = self.__counter_upstream_tracks
    temp_dict['downstream_tracks'] = self.__counter_downstream_tracks

    data_dict['scifi_straight_track_extractor'] = temp_dict
    return data_dict


################################################################################


class scifi_helical_track_extractor(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_helical_track_extractor")

    self.__upstream_tracks = []
    self.__downstream_tracks = []

    self.__counter_upstream_tracks = 0
    self.__counter_downstream_tracks = 0


  def get_upstream_tracks(self) :
    return self.__upstream_tracks


  def get_downstream_tracks(self) :
    return self.__downstream_tracks


  def get_dependencies(self, inserter) :
    pass


  def _reset(self) :
    self.__upstream_tracks = []
    self.__downstream_tracks = []


  def _process(self, file_reader) :
    scifi_event = file_reader.get_event('scifi')
    scifi_tracks = scifi_event.scifitracks()

    for track in scifi_tracks :
      if track.GetAlgorithmUsed() != tools.HELICAL_ALGORITHM_ID :
        continue

      if track.tracker() == 0 :
        self.__upstream_tracks.append(track)
        self.__counter_upstream_tracks += 1
      elif track.tracker() == 1 :
        self.__downstream_tracks.append(track)
        self.__counter_downstream_tracks += 1

    return False


  def _store_data(self, data_dict) :
    temp_dict = {}
    temp_dict['number_events'] = self.get_number_events()
    temp_dict['upstream_tracks'] = self.__counter_upstream_tracks
    temp_dict['downstream_tracks'] = self.__counter_downstream_tracks

    data_dict['scifi_helical_track_extractor'] = temp_dict
    return data_dict


################################################################################


class scifi_straight_track_candidates(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_straight_track_candidates")

    self.__upstream_track = None
    self.__downstream_track = None

    self.__upstream_trackpoints = []
    self.__downstream_trackpoints = []

    self.__ref_station = 1
    self.__ref_plane = 0

    self.__selection_tracker = 1
    self.__selection_station = 1
    self.__selection_plane = 0

    self.__count_mismatch = 0
    self.__count_chisqndf = 0
    self.__count_no_tps = 0
    self.__count_grad = 0
    self.__count_rad = 0
    self.__count_success = 0

    self.__pval = 0.0
    self.__chisqndf = 0.0
    self.__no_tps = 0
    self.__grad_low = 0.0
    self.__grad_high = 100.0
    self.__rad_low = 0.0
    self.__rad_high = 1000.0

    self.__up_n_tracks = ROOT.TH1F('ste_up_number_tracks', \
                                            "Number of Tracks", 21, -0.5, 20.5)
    self.__up_pval = ROOT.TH1F('ste_up_p_value', \
                                               "Track P Values", 500, 0.0, 1.0)
    self.__up_chisqndf = ROOT.TH1F('ste_up_chi2ndf', \
                                               "Track P Values", 500, 0.0, 1.0)
    self.__up_chisqndf_cut = ROOT.TH1F('ste_up_chi2ndf_cut', \
                                     "Track Chi2 per NDF (Cut)", 500, 0.0, 1.0)
    self.__up_ntps = ROOT.TH1F('ste_up_no_trackpoints', \
                                       "Number of Trackpoints", 17, -0.5, 16.5)
    self.__up_ntps_cut = ROOT.TH1F('ste_up_no_trackpoints_cut', \
                                 "Number of Trackpoints (Cut)", 17, -0.5, 16.5)
    self.__up_grad = ROOT.TH2F('ste_up_grad', "Track Gradient", \
                                                100, -0.5, 0.5, 100, -0.5, 0.5)
    self.__up_grad_cut = ROOT.TH2F('ste_up_grad_cut', "Track Gradient (Cut)", \
                                                100, -0.5, 0.5, 100, -0.5, 0.5)
    self.__up_pos = ROOT.TH2F('ste_up_pos', "Track Radius", \
                                        500, -200.0, 200.0, 500, -200.0, 200.0)
    self.__up_pos_cut = ROOT.TH2F('ste_up_pos_cut', "Track Radius (Cut)", \
                                        500, -200.0, 200.0, 500, -200.0, 200.0)

    self.__down_n_tracks = ROOT.TH1F('ste_down_number_tracks', \
                                            "Number of Tracks", 21, -0.5, 20.5)
    self.__down_pval = ROOT.TH1F('ste_down_p_value', "Track P Values", \
                                                                 500, 0.0, 1.0)
    self.__down_chisqndf = ROOT.TH1F('ste_down_chi2ndf', \
                                           "Track Chi2 per NDF", 500, 0.0, 1.0)
    self.__down_chisqndf_cut = ROOT.TH1F('ste_down_chi2ndf_cut', \
                                     "Track Chi2 per NDF (Cut)", 500, 0.0, 1.0)
    self.__down_ntps = ROOT.TH1F('ste_down_no_trackpoints', \
                                       "Number of Trackpoints", 17, -0.5, 16.5)
    self.__down_ntps_cut = ROOT.TH1F('ste_down_no_trackpoints_cut', \
                                 "Number of Trackpoints (Cut)", 17, -0.5, 16.5)
    self.__down_grad = ROOT.TH2F('ste_down_grad', "Track Gradient", \
                                                100, -0.5, 0.5, 100, -0.5, 0.5)
    self.__down_grad_cut = ROOT.TH2F('ste_down_grad_cut', \
                        "Track Gradient (Cut)", 100, -0.5, 0.5, 100, -0.5, 0.5)
    self.__down_pos = ROOT.TH2F('ste_down_pos', "Track Radius", \
                                        500, -200.0, 200.0, 500, -200.0, 200.0)
    self.__down_pos_cut = ROOT.TH2F('ste_down_pos_cut', "Track Radius (Cut)", \
                                        500, -200.0, 200.0, 500, -200.0, 200.0)


  def set_cut_pvalue(self, pval) :
    self.__pval = pval


  def set_cut_number_trackpoints(self, num) :
    self.__no_tps = num


  def set_cut_gradient_window(self, low, high) :
    self.__grad_low = low
    self.__grad_high = high


  def set_cut_radius_window(self, low, high) :
    self.__rad_low = low
    self.__rad_high = high


  def set_reference_station(self, station) :
    self.__ref_station = station


  def set_reference_plane(self, plane) :
    self.__ref_plane = plane


  def get_trackpoint(self, tracker, station, plane) :
    if tracker == 0 :
      return self.get_upstream_trackpoint(int((1 + plane + (station - 1) * 3)))
    else :
      return self.get_downstream_trackpoint(int((1 + plane + \
                                                           (station - 1) * 3)))


  def get_upstream_track(self) :
    return self.__upstream_track


  def get_upstream_trackpoint(self, id_num) :
    return self.__upstream_trackpoints[id_num]


  def get_upstream_reference(self) :
    return self.__upstream_trackpoints[int((1 + self.__ref_plane + \
                                                (self.__ref_station - 1) * 3))]


  def get_downstream_track(self) :
    return self.__downstream_track


  def get_downstream_trackpoint(self, id_num) :
    return self.__downstream_trackpoints[id_num]


  def get_downstream_reference(self) :
    return self.__downstream_trackpoints[int((1 + self.__ref_plane + \
                                                (self.__ref_station - 1) * 3))]


  def get_dependencies(self, inserter) :
    self.__straight_tracks = inserter(scifi_straight_track_extractor())


  def get_args(self, parser) :
    parser.add_argument( '--cut_chisq_ndf', type=float, default=100.0, \
          help='Set the cut on the tracker chi-squared per degree of freedom' )
#   parser.add_argument( '--cut_p_value', type=float, default=0.0, \
#                             help='Set the cut on the tracker P-Value [0-1]' )
    parser.add_argument( '--cut_number_trackpoints', type=int, default=0, \
                    help='Set the cut on the number of trackpoints per track' )
    parser.add_argument( '--cut_gradient', type=float, default=1.0, \
                                     help='Set the cut on the track gradient' )
    parser.add_argument( '--cut_min_gradient', type=float, default=0.0, \
                             help='Set the cut on the minimum track gradient' )
    parser.add_argument( '--cut_radius', type=float, default=200.0, \
                              help='Set the cut on the track upstream radius' )
    parser.add_argument( '--cut_min_radius', type=float, default=0.0, \
                      help='Set the cut on the track upstream minimum radius' )


  def process_args(self, namespace) :
#    self.__pval = namespace.cut_p_value
    self.__chisqndf = namespace.cut_chisq_ndf
    self.__no_tps = namespace.cut_number_trackpoints
    self.__grad_low = namespace.cut_min_gradient
    self.__grad_high = namespace.cut_gradient
    self.__rad_low = namespace.cut_min_radius
    self.__rad_high = namespace.cut_radius


  def _reset(self) :
    self.__upstream_track = None
    self.__downstream_track = None

    self.__upstream_trackpoints = [ None for i in range(16) ]
    self.__downstream_trackpoints = [ None for i in range(16) ]


  def _process(self, file_reader) :
    upstream_tracks = self.__straight_tracks.get_upstream_tracks()
    downstream_tracks = self.__straight_tracks.get_downstream_tracks()

    self.__up_n_tracks.Fill(len(upstream_tracks))
    self.__down_n_tracks.Fill(len(downstream_tracks))

    best_upstream_track = None
    up_p_val = -1.0e100
    for trk in upstream_tracks :
      p_val = trk.P_value()
      if p_val > up_p_val :
        best_upstream_track = trk
        up_p_val = p_val
      
    best_downstream_track = None
    down_p_val = -1.0e100
    for trk in downstream_tracks :
      p_val = trk.P_value()
      if p_val > down_p_val :
        best_downstream_track = trk
        down_p_val = p_val

    if best_upstream_track is None or best_downstream_track is None :
      self._cut()
      self.__count_mismatch += 1
      return True # Can't continue otherwise

    up_trk = best_upstream_track
    down_trk = best_downstream_track

#    if len(good_upstream_tracks) != 1 or len(good_downstream_tracks) != 1 :
#      self._cut()
#      self.__count_mismatch += 1
#      return True # Can't continue otherwise

#    up_trk = upstream_tracks[0]
#    down_trk = downstream_tracks[0]


    up_p_val = up_trk.P_value()
    self.__up_pval.Fill(up_p_val)
#    if up_p_val < self.__pval :
#      self._cut()
#      self.__count_pval += 1
#    else :
#      self.__up_pval_cut.Fill(up_p_val)

    up_chisqn = up_trk.chi2() / up_trk.ndf()
    self.__up_chisqndf.Fill( up_chisqn )
    if up_chisqn > self.__chisqndf :
      self._cut()
      self.__count_chisqndf += 1
    else :
      self.__up_chisqndf_cut.Fill(up_chisqn)
    

    down_p_val = down_trk.P_value()
    self.__down_pval.Fill(down_p_val)
#    if down_p_val < self.__pval :
#      self._cut()
#      self.__count_pval += 1
#    else :
#      self.__down_pval_cut.Fill(down_p_val)

    down_chisqn = down_trk.chi2() / down_trk.ndf()
    self.__down_chisqndf.Fill( down_chisqn )
    if down_chisqn > self.__chisqndf :
      self._cut()
      self.__count_chisqndf += 1
    else :
      self.__down_chisqndf_cut.Fill(down_chisqn)


    up_tp_count = 0
    down_tp_count = 0

    for tp in up_trk.scifitrackpoints() :
      if tp.has_data() :
        up_tp_count += 1
      if tp.station() == self.__ref_station and tp.plane() == self.__ref_plane :
        up_ref = tp

    self.__up_ntps.Fill(up_tp_count)
    if up_tp_count < self.__no_tps :
      self._cut()
      self.__count_no_tps += 1
    else :
      self.__up_ntps_cut.Fill(up_tp_count)


    for tp in down_trk.scifitrackpoints() :
      if tp.has_data() :
        down_tp_count += 1
      if tp.station() == self.__ref_station and tp.plane() == self.__ref_plane :
        down_ref = tp

    self.__down_ntps.Fill(down_tp_count)
    if down_tp_count < self.__no_tps :
      self._cut()
      self.__count_no_tps += 1
    else :
      self.__down_ntps_cut.Fill(down_tp_count)


    up_pos = [ up_ref.pos().x(), up_ref.pos().y() ]
    up_gra = [ up_ref.mom().x() / up_ref.mom().z(), \
                                          up_ref.mom().y() / up_ref.mom().z() ]

    down_pos = [ down_ref.pos().x(), down_ref.pos().y() ]
    down_gra = [ down_ref.mom().x() / down_ref.mom().z(), \
                                      down_ref.mom().y() / down_ref.mom().z() ]

#    up_pos = up_tracker_correction( up_pos )
#    up_gra = up_tracker_correction( up_gra )

    up_rad = math.sqrt( up_pos[0]**2 + up_pos[1]**2 )
    up_grad = math.sqrt( up_gra[0]**2 + up_gra[1]**2 )

    self.__up_grad.Fill(up_gra[0], up_gra[1])
    if up_grad > self.__grad_high or up_grad < self.__grad_low :
      self._cut()
      self.__count_grad += 1
    else :
      self.__up_grad_cut.Fill(up_gra[0], up_gra[1])

    self.__up_pos.Fill(up_pos[0], up_pos[1])
    if up_rad > self.__rad_high or up_rad < self.__rad_low :
      self._cut()
      self.__count_rad += 1
    else :
      self.__up_pos_cut.Fill(up_pos[0], up_pos[1])


    down_rad = math.sqrt( down_pos[0]**2 + down_pos[1]**2 )
    down_grad = math.sqrt( down_gra[0]**2 + down_gra[1]**2 )

    self.__down_grad.Fill(down_gra[0], down_gra[1])
    if down_grad > self.__grad_high or down_grad < self.__grad_low :
      self._cut()
      self.__count_grad += 1
    else :
      self.__down_grad_cut.Fill(down_gra[0], down_gra[1])

    self.__down_pos.Fill(down_pos[0], down_pos[1])
    if down_rad > self.__rad_high or down_rad < self.__rad_low :
      self._cut()
      self.__count_rad += 1
    else :
      self.__down_pos_cut.Fill(down_pos[0], down_pos[1])


    for tp in up_trk.scifitrackpoints() :
      id_n = abs(tools.calculate_plane_id(0, tp.station(), tp.plane()))
      self.__upstream_trackpoints[id_n] = tp
    self.__upstream_track = up_trk

    for tp in down_trk.scifitrackpoints() :
      id_n = abs(tools.calculate_plane_id(1, tp.station(), tp.plane()))
      self.__downstream_trackpoints[id_n] = tp
    self.__downstream_track = down_trk


    if self.is_cut() :
      return True
    else :
      self.__count_success += 1
      return False


  def _store_plots(self, plot_dict) :

    upstream_plots = {}

    upstream_plots['number_tracks'] = self.__up_n_tracks
    upstream_plots['p_value'] = self.__up_pval
    upstream_plots['chisq_ndf'] = self.__up_chisqndf
    upstream_plots['chisq_ndf_cut'] = self.__up_chisqndf_cut
    upstream_plots['number_trackpoints'] = self.__up_ntps
    upstream_plots['number_trackpoints_cut'] = self.__up_ntps_cut
    upstream_plots['gradient'] = self.__up_grad
    upstream_plots['gradient_cut'] = self.__up_grad_cut
    upstream_plots['position'] = self.__up_pos
    upstream_plots['position_cut'] = self.__up_pos_cut

    downstream_plots = {}

    downstream_plots['number_tracks'] = self.__down_n_tracks
    downstream_plots['p_value'] = self.__down_pval
    downstream_plots['chisq_ndf'] = self.__down_chisqndf
    downstream_plots['chisq_ndf_cut'] = self.__down_chisqndf_cut
    downstream_plots['number_trackpoints'] = self.__down_ntps
    downstream_plots['number_trackpoints_cut'] = self.__down_ntps_cut
    downstream_plots['gradient'] = self.__down_grad
    downstream_plots['gradient_cut'] = self.__down_grad_cut
    downstream_plots['position'] = self.__down_pos
    downstream_plots['position_cut'] = self.__down_pos_cut

    plot_dict['tracks'] = {}

    plot_dict['tracks']['upstream'] = upstream_plots
    plot_dict['tracks']['downstream'] = downstream_plots

    return plot_dict


  def _store_data(self, data_dict) :

    track_cuts = {}

    track_cuts['mismatch'] = self.__count_mismatch
    track_cuts['chisqndf'] = self.__count_chisqndf
    track_cuts['number_trackpoints'] = self.__count_no_tps
    track_cuts['gradient'] = self.__count_grad 
    track_cuts['radius'] = self.__count_rad
    track_cuts['passed_pairs'] = self.__count_success
    
    data_dict['track_cuts'] = track_cuts


################################################################################



class scifi_helical_track_candidates(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_helical_track_candidates")

    self.__upstream_track = None
    self.__downstream_track = None

    self.__upstream_trackpoints = []
    self.__downstream_trackpoints = []

    self.__ref_station = 1
    self.__ref_plane = 0

    self.__count_mismatch = 0
    self.__count_chisqndf = 0
    self.__count_no_tps = 0
    self.__count_pt = 0
    self.__count_rad = 0
    self.__count_success = 0

    self.__chisqndf = 0.0
    self.__no_tps = 0
    self.__pt_low = 0.0
    self.__pt_high = 500.0
    self.__rad_low = 0.0
    self.__rad_high = 1000.0

    self.__up_n_tracks = ROOT.TH1F('hte_up_number_tracks', \
                                            "Number of Tracks", 21, -0.5, 20.5)
    self.__up_pval = ROOT.TH1F('hte_up_p_value', \
                                               "Track P Values", 500, 0.0, 1.0)
    self.__up_chisqndf = ROOT.TH1F('hte_up_chisqndf', \
                                             "Track P Chi2 NDF", 500, 0.0, 1.0)
    self.__up_chisqndf_cut = ROOT.TH1F('hte_up_chisqndf_cut', \
                                       "Track P Chi2 NDF (Cut)", 500, 0.0, 1.0)
    self.__up_ntps = ROOT.TH1F('hte_up_no_trackpoints', \
                                       "Number of Trackpoints", 17, -0.5, 16.5)
    self.__up_ntps_cut = ROOT.TH1F('hte_up_no_trackpoints_cut', \
                                 "Number of Trackpoints (Cut)", 17, -0.5, 16.5)
    self.__up_pt = ROOT.TH1F('hte_up_pt', "Track Pt", 200, 0.0, 200.0)
    self.__up_pt_cut = ROOT.TH1F('hte_up_pt_cut', "Track Pt (Cut)", \
                                                               200, 0.0, 200.0)

    self.__down_n_tracks = ROOT.TH1F('hte_down_number_tracks', \
                                            "Number of Tracks", 21, -0.5, 20.5)
    self.__down_pval = ROOT.TH1F('hte_down_p_value', "Track P Values", \
                                                                 500, 0.0, 1.0)
    self.__down_chisqndf = ROOT.TH1F('hte_down_chisqndf', \
                                             "Track P Chi2 NDF", 500, 0.0, 1.0)
    self.__down_chisqndf_cut = ROOT.TH1F('hte_down_chisqndf_cut', \
                                       "Track P Chi2 NDF (Cut)", 500, 0.0, 1.0)
    self.__down_ntps = ROOT.TH1F('hte_down_no_trackpoints', \
                                       "Number of Trackpoints", 17, -0.5, 16.5)
    self.__down_ntps_cut = ROOT.TH1F('hte_down_no_trackpoints_cut', \
                                 "Number of Trackpoints (Cut)", 17, -0.5, 16.5)
    self.__down_pt = ROOT.TH1F('hte_down_pt', "Track Pt", 200, 0.0, 200.0)
    self.__down_pt_cut = ROOT.TH1F('hte_down_pt_cut', "Track Pt (Cut)", \
                                                               200, 0.0, 200.0)


  def set_cut_number_trackpoints(self, num) :
    self.__no_tps = num


  def set_cut_pt_window(self, low, high) :
    self.__pt_low = low
    self.__pt_high = high


  def set_reference_station(self, station) :
    self.__ref_station = station


  def set_reference_plane(self, plane) :
    self.__ref_plane = plane


  def get_trackpoint(self, tracker, station, plane) :
    if tracker == 0 :
      return self.get_upstream_trackpoint(int((1 + plane + (station - 1) * 3)))
    else :
      return self.get_downstream_trackpoint(int((1 + plane + \
                                                           (station - 1) * 3)))


  def get_upstream_track(self) :
    return self.__upstream_track


  def get_upstream_trackpoint(self, id_num) :
    return self.__upstream_trackpoints[id_num]


  def get_upstream_reference(self) :
    return self.__upstream_trackpoints[int((1 + self.__ref_plane + \
                                                (self.__ref_station - 1) * 3))]


  def get_downstream_track(self) :
    return self.__downstream_track


  def get_downstream_trackpoint(self, id_num) :
    return self.__downstream_trackpoints[id_num]


  def get_downstream_reference(self) :
    return self.__downstream_trackpoints[int((1 + self.__ref_plane + \
                                                (self.__ref_station - 1) * 3))]


  def get_dependencies(self, inserter) :
    self.__helical_tracks = inserter(scifi_helical_track_extractor())


  def get_args(self, parser) :
    parser.add_argument( '--helix_cut_chisq_ndf', type=float, default=100.0, \
          help='Set the cut on the tracker chi-squared per degree of freedom' )
    parser.add_argument( '--helix_cut_number_trackpoints', type=int, default=0, \
                    help='Set the cut on the number of trackpoints per track' )
    parser.add_argument( '--helix_cut_pt', type=float, default=500.0, \
                                     help='Set the cut on the track gradient' )
    parser.add_argument( '--helix_cut_min_pt', type=float, default=0.0, \
                             help='Set the cut on the minimum track gradient' )


  def process_args(self, namespace) :
    self.__chisqndf = namespace.helix_cut_chisq_ndf
    self.__no_tps = namespace.helix_cut_number_trackpoints
    self.__pt_low = namespace.helix_cut_min_pt
    self.__pt_high = namespace.helix_cut_pt


  def _reset(self) :
    self.__upstream_track = None
    self.__downstream_track = None

    self.__upstream_trackpoints = [ None for i in range(16) ]
    self.__downstream_trackpoints = [ None for i in range(16) ]


  def _process(self, file_reader) :
    upstream_tracks = self.__helical_tracks.get_upstream_tracks()
    downstream_tracks = self.__helical_tracks.get_downstream_tracks()

    self.__up_n_tracks.Fill(len(upstream_tracks))
    self.__down_n_tracks.Fill(len(downstream_tracks))

    if len(upstream_tracks) != 1 or len(downstream_tracks) != 1 :
      self._cut()
      self.__count_mismatch += 1
      return True # Can't continue otherwise


    up_trk = upstream_tracks[0]
    down_trk = downstream_tracks[0]


    up_p_val = up_trk.P_value()
    self.__up_pval.Fill(up_p_val)
    down_p_val = down_trk.P_value()
    self.__down_pval.Fill(down_p_val)

    up_chisqndf = up_trk.chi2() / up_trk.ndf()
    self.__up_chisqndf.Fill( up_chisqndf )
    if up_chisqndf > self.__chisqndf :
      self._cut()
      self.__count_chisqndf += 1
    else :
      self.__up_chisqndf_cut.Fill( up_chisqndf )

    down_chisqndf = down_trk.chi2() / down_trk.ndf()
    self.__down_chisqndf.Fill( down_chisqndf )
    if down_chisqndf > self.__chisqndf :
      self._cut()
      self.__count_chisqndf += 1
    else :
      self.__down_chisqndf_cut.Fill( down_chisqndf )

    up_tp_count = 0
    down_tp_count = 0

    for tp in up_trk.scifitrackpoints() :
      if tp.has_data() :
        up_tp_count += 1
      if tp.station() == self.__ref_station and tp.plane() == self.__ref_plane :
        up_ref = tp

    self.__up_ntps.Fill(up_tp_count)
    if up_tp_count < self.__no_tps :
      self._cut()
      self.__count_no_tps += 1
    else :
      self.__up_ntps_cut.Fill(up_tp_count)


    for tp in down_trk.scifitrackpoints() :
      if tp.has_data() :
        down_tp_count += 1
      if tp.station() == self.__ref_station and tp.plane() == self.__ref_plane :
        down_ref = tp

    self.__down_ntps.Fill(down_tp_count)
    if down_tp_count < self.__no_tps :
      self._cut()
      self.__count_no_tps += 1
    else :
      self.__down_ntps_cut.Fill(down_tp_count)


    up_mom = [ up_ref.mom().x(), up_ref.mom().y() ]
    down_mom = [ down_ref.mom().x(), down_ref.mom().y() ]

    up_pt = math.sqrt( up_mom[0]**2 + up_mom[1]**2 )
    down_pt = math.sqrt( down_mom[0]**2 + down_mom[1]**2 )

    self.__up_pt.Fill(up_pt)
    if up_pt > self.__pt_high or up_pt < self.__pt_low :
      self._cut()
      self.__count_pt += 1
    else :
      self.__up_pt_cut.Fill(up_pt)


    self.__down_pt.Fill(down_pt)
    if down_pt > self.__pt_high or down_pt < self.__pt_low :
      self._cut()
      self.__count_pt += 1
    else :
      self.__down_pt_cut.Fill(down_pt)


    for tp in up_trk.scifitrackpoints() :
      id_n = abs(tools.calculate_plane_id(0, tp.station(), tp.plane()))
      self.__upstream_trackpoints[id_n] = tp
    self.__upstream_track = up_trk

    for tp in down_trk.scifitrackpoints() :
      id_n = abs(tools.calculate_plane_id(1, tp.station(), tp.plane()))
      self.__downstream_trackpoints[id_n] = tp
    self.__downstream_track = down_trk


    if not self.is_cut() :
      return True
    else:
      self.__count_success += 1
      return False


  def _store_plots(self, plot_dict) :

    upstream_plots = {}

    upstream_plots['number_tracks'] = self.__up_n_tracks
    upstream_plots['p_value'] = self.__up_pval
    upstream_plots['chisq_ndf'] = self.__up_chisqndf
    upstream_plots['chisq_ndf_cut'] = self.__up_chisqndf_cut
    upstream_plots['number_trackpoints'] = self.__up_ntps
    upstream_plots['number_trackpoints_cut'] = self.__up_ntps_cut
    upstream_plots['pt'] = self.__up_pt
    upstream_plots['pt_cut'] = self.__up_pt_cut

    downstream_plots = {}

    downstream_plots['number_tracks'] = self.__down_n_tracks
    downstream_plots['p_value'] = self.__down_pval
    downstream_plots['chisq_ndf'] = self.__down_chisqndf
    downstream_plots['chisq_ndf_cut'] = self.__down_chisqndf_cut
    downstream_plots['number_trackpoints'] = self.__down_ntps
    downstream_plots['number_trackpoints_cut'] = self.__down_ntps_cut
    downstream_plots['pt'] = self.__down_pt
    downstream_plots['pt_cut'] = self.__down_pt_cut

    plot_dict['tracks'] = {}

    plot_dict['tracks']['upstream'] = upstream_plots
    plot_dict['tracks']['downstream'] = downstream_plots

    return plot_dict


  def _store_data(self, data_dict) :

    track_cuts = {}

    track_cuts['mismatch'] = self.__count_mismatch
    track_cuts['number_trackpoints'] = self.__count_no_tps
    track_cuts['pt'] = self.__count_pt 
    track_cuts['passed_pairs'] = self.__count_success
    
    data_dict['track_cuts'] = track_cuts


################################################################################


class scifi_helical_track_processor(framework.processor_base) :

  def __init__(self, tracker, apply_cuts=True) :
    framework.processor_base.__init__(self, "scifi_helical_track_processor_"+str(tracker))

    self.__tracker = tracker
    self.__apply_cuts = apply_cuts
    self.__track = None

    self.__trackpoints = []

    self.__ref_station = 1
    self.__ref_plane = 0

    self.__count_mismatch = 0
    self.__count_chisqndf = 0
    self.__count_no_tps = 0
    self.__count_pt = 0
    self.__count_rad = 0
    self.__count_momentum_window = 0
    self.__count_success = 0

    self.__chisqndf = 0.0
    self.__no_tps = 0
    self.__pt_low = 0.0
    self.__pt_high = 500.0
    self.__momentum_window = [0.0, 500.0]
    self.__rad_low = 0.0
    self.__rad_high = 1000.0

    self.__plot_n_tracks = ROOT.TH1F('hte_number_tracks', \
                                                 "Number of Tracks", 21, -0.5, 20.5)
    self.__plot_pval = ROOT.TH1F('hte_p_value', \
                                                    "Track P Values", 500, 0.0, 1.0)
    self.__plot_chisqndf = ROOT.TH1F('hte_chisqndf', \
                                               "Track ChiSq per DoF", 500, 0.0, 1.0)
    self.__plot_chisqndf_cut = ROOT.TH1F('hte_chisqndf_cut', \
                                         "Track ChiSq per DoF (Cut)", 500, 0.0, 1.0)
    self.__plot_ntps = ROOT.TH1F('hte_no_trackpoints', \
                                            "Number of Trackpoints", 17, -0.5, 16.5)
    self.__plot_ntps_cut = ROOT.TH1F('hte_no_trackpoints_cut', \
                                      "Number of Trackpoints (Cut)", 17, -0.5, 16.5)
    self.__plot_pt = ROOT.TH1F('hte_pt', "Track Pt", 200, 0.0, 200.0)
    self.__plot_pt_cut = ROOT.TH1F('hte_pt_cut', "Track Pt (Cut)", \
                                                               200, 0.0, 200.0)
    self.__plot_p = ROOT.TH1F('hte_p', "Track P", 500, 0.0, 500.0)
    self.__plot_p_cut = ROOT.TH1F('hte_p_cut', "Track P (Cut)", \
                                                               500, 0.0, 500.0)


  def set_cut_number_trackpoints(self, num) :
    self.__no_tps = num


  def set_cut_pt_window(self, low, high) :
    self.__pt_low = low
    self.__pt_high = high


  def set_reference_station(self, station) :
    self.__ref_station = station


  def set_reference_plane(self, plane) :
    self.__ref_plane = plane


  def get_trackpoint_byplane(self, station, plane) :
#    if self.get_trackpoint(int((1 + plane + (station - 1) * 3))) is None :
#      print
#      print self.__tracker, self.is_cut(), station, plane, len(self.__trackpoints), int((1 + plane + (station - 1) * 3))
#      print self.__trackpoints
    return self.get_trackpoint(int((1 + plane + (station - 1) * 3)))


  def get_track(self) :
    return self.__track


  def get_trackpoint(self, id_num) :
    return self.__trackpoints[id_num]


  def get_reference(self) :
    return self.__trackpoints[int((1 + self.__ref_plane + \
                                                (self.__ref_station - 1) * 3))]


  def get_dependencies(self, inserter) :
    self.__helical_tracks = inserter(scifi_helical_track_extractor())


  def get_args(self, parser) :
    parser.add_argument( '--select_plane', type=int, nargs=3, default=None, help='Specify the <tracker> <station> <plane> to act as the selection plane' )
    parser.add_argument( '--helix_cut_chisq_ndf', type=float, default=0.0, \
           help='Set the cut on the tracker Chi Square per degree of freedom' )
    parser.add_argument( '--helix_cut_number_trackpoints', type=int, default=0, \
                    help='Set the cut on the number of trackpoints per track' )
    parser.add_argument( '--helix_cut_pt', type=float, default=500.0, \
                                     help='Set the cut on the track gradient' )
    parser.add_argument( '--helix_cut_min_pt', type=float, default=0.0, \
                             help='Set the cut on the minimum track gradient' )
    parser.add_argument( '--helix_momentum_window', type=float, nargs=2, \
                      default=[0.0, 500.0], help='Set the total momentum cut')


  def process_args(self, namespace) :
    if namespace.select_plane is None :
      self.__selection_tracker = 0
      self.__selection_station = 1
      self.__selection_plane = 0
    else :
      self.__selection_tracker = namespace.select_plane[0]
      self.__selection_station = namespace.select_plane[1]
      self.__selection_plane = namespace.select_plane[2]
    self.__chisqndf = namespace.helix_cut_chisq_ndf
    self.__no_tps = namespace.helix_cut_number_trackpoints
    self.__pt_low = namespace.helix_cut_min_pt
    self.__pt_high = namespace.helix_cut_pt
    self.__momentum_window = namespace.helix_momentum_window


  def _reset(self) :
    self.__track = None

    self.__trackpoints = [ None for i in range(16) ]


  def _process(self, file_reader) :
    if self.__tracker == 0 :
      tracks = self.__helical_tracks.get_upstream_tracks()
    else :
      tracks = self.__helical_tracks.get_downstream_tracks()

    self.__plot_n_tracks.Fill(len(tracks))

    if len(tracks) < 1 :
      self._cut()
      return True

    trk = tracks[0]
    tp_count = 0

    for tp in trk.scifitrackpoints() :
      if tp.has_data() :
        tp_count += 1
      if tp.station() == self.__ref_station and tp.plane() == self.__ref_plane :
        ref = tp

    mom = [ ref.mom().x(), ref.mom().y(), ref.mom().z() ]
    pt = math.sqrt( mom[0]**2 + mom[1]**2 )
    p = math.sqrt( mom[0]**2 + mom[1]**2 + mom[2]**2 )
    p_val = trk.P_value()
    chisqndf = trk.chi2() / trk.ndf()

    self.__plot_pt.Fill(pt)
    self.__plot_chisqndf.Fill(chisqndf)
    self.__plot_p.Fill(p)
    self.__plot_pval.Fill(p_val)
    self.__plot_ntps.Fill(tp_count)

    if chisqndf > self.__chisqndf :
      self._cut()
      self.__count_chisqndf += 1
    else :
      self.__plot_chisqndf_cut.Fill(chisqndf)

    if tp_count < self.__no_tps :
      self._cut()
      self.__count_no_tps += 1
    else :
      self.__plot_ntps_cut.Fill(tp_count)

    if self.__apply_cuts :
      if pt > self.__pt_high or pt < self.__pt_low :
        self._cut()
        self.__count_pt += 1
      else :
        self.__plot_pt_cut.Fill(pt)
  
      if p < self.__momentum_window[0] or p > self.__momentum_window[1] :
        self._cut()
        self.__count_momentum_window += 1
      else :
        self.__plot_p_cut.Fill(p)

    for tp in trk.scifitrackpoints() :
      id_n = abs(tools.calculate_plane_id(0, tp.station(), tp.plane()))
      self.__trackpoints[id_n] = tp
    self.__track = trk


    if self.is_cut() :
      return True
    else :
      self.__count_success += 1
      return False


  def _store_plots(self, plot_dict) :

    plots = {}

    plots['number_tracks'] = self.__plot_n_tracks
    plots['p_value'] = self.__plot_pval
    plots['chisq_ndf'] = self.__plot_chisqndf
    plots['chisq_ndf_cut'] = self.__plot_chisqndf_cut
    plots['number_trackpoints'] = self.__plot_ntps
    plots['number_trackpoints_cut'] = self.__plot_ntps_cut
    plots['pt'] = self.__plot_pt
    plots['pt_cut'] = self.__plot_pt_cut
    plots['p'] = self.__plot_p
    plots['p_cut'] = self.__plot_p_cut

    if 'tracks' not in plot_dict : 
      plot_dict['tracks'] = {}

    plot_dict['tracks']['tracker_'+str(self.__tracker)] = plots

    return plot_dict


  def _store_data(self, data_dict) :

    track_cuts = {}

    track_cuts['mismatch'] = self.__count_mismatch
    track_cuts['chisqndf'] = self.__count_chisqndf
    track_cuts['number_trackpoints'] = self.__count_no_tps
    track_cuts['pt'] = self.__count_pt 
    track_cuts['passed_pairs'] = self.__count_success
    
    data_dict['tracker_'+str(self.__tracker)+'_cuts'] = track_cuts


