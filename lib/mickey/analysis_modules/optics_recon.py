
from _analysis_base import Analysis_Base
from analysis import inspectors
from analysis import tools

import ROOT
import math
import array
import copy


class OpticsAnalysis(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "optics_reconstruction")

    self.__inspectors = { 0: {}, 1: {} }
    self.__data = { 0: {}, 1: {} }
    self.__do_upstream = False
    self.__do_downstream = False
    self.__graphs = {}


  def analyse_event(self, analysis_event, weight=1.0) :
    if self.__do_upstream and (analysis_event.num_upstream_tracks()) :
      track = analysis_event.upstream_track()
      for plane_id in self.__inspectors[0] :
        hit = track[plane_id]
        hit.set_weight(weight)
        self.__inspectors[0][plane_id].add_hit( hit )

    if self.__do_downstream and (analysis_event.num_downstream_tracks()):
      track = analysis_event.downstream_track()
      for plane_id in self.__inspectors[1] :
        hit = track[plane_id]
        hit.set_weight(weight)
        self.__inspectors[1][plane_id].add_hit( hit )


  def _get_plots(self, plot_dict) :
    for graph in self.__graphs :
      plot_dict[graph] = self.__graphs[graph]


  def _get_data(self, data_dict) :
    data_dict['upstream'] = self.__data[0]
    data_dict['downstream'] = self.__data[1]


  def configure_arguments(self, parser) :
    parser.add_argument( "--trackers", nargs='+', type=int, default=[0, 1], help="Specify which trackers to analyse" )
    parser.add_argument( "--planes", nargs='+', type=int, default=[1, 4, 7, 10, 13], help="Specify which tracker planes to analyse (1-15)" )


  def parse_arguments(self, namespace) :
    if 0 in namespace.trackers :
      self.__do_upstream = True
    if 1 in namespace.trackers :
      self.__do_downstream = True

    for i in namespace.trackers :
      for j in namespace.planes :
        self.__inspectors[i][j] = inspectors.PhaseSpace2DInspector(i, 0)
        self.__data[i][j] = None


  def conclude(self) :
    for i in self.__inspectors :
      for j in self.__inspectors[i] :
        self.__data[i][j] = self.__inspectors[i][j].get_data_dictionary()

    emittance = array.array("d")
    emittance_x = array.array("d")
    emittance_y = array.array("d")
    beta = array.array("d")
    beta_x = array.array("d")
    beta_y = array.array("d")
    alpha = array.array("d")
    alpha_x = array.array("d")
    alpha_y = array.array("d")
    momentum = array.array("d")
    number_particles = array.array("d")

    error_emittance = array.array("d")
    error_emittance_x = array.array("d")
    error_emittance_y = array.array("d")
    error_beta = array.array("d")
    error_beta_x = array.array("d")
    error_beta_y = array.array("d")
    error_alpha = array.array("d")
    error_alpha_x = array.array("d")
    error_alpha_y = array.array("d")
    error_momentum = array.array("d")
    zeros = array.array("d")
    position = array.array("d")

    for i in self.__inspectors :
      for j in self.__inspectors[i] :
        emittance.append(self.__data[i][j]['emittance'])
        emittance_x.append(self.__data[i][j]['emittance_x'])
        emittance_y.append(self.__data[i][j]['emittance_y'])
        beta.append(self.__data[i][j]['beta'])
        beta_x.append(self.__data[i][j]['beta_x'])
        beta_y.append(self.__data[i][j]['beta_y'])
        alpha.append(self.__data[i][j]['alpha'])
        alpha_x.append(self.__data[i][j]['alpha_x'])
        alpha_y.append(self.__data[i][j]['alpha_y'])
        momentum.append(self.__data[i][j]['momentum'])
        number_particles.append(self.__data[i][j]['number_particles'])

        error_emittance.append(self.__data[i][j]['emittance_error'])
        error_emittance_x.append(self.__data[i][j]['emittance_x_error'])
        error_emittance_y.append(self.__data[i][j]['emittance_y_error'])
        error_beta.append(self.__data[i][j]['beta_error'])
        error_beta_x.append(self.__data[i][j]['beta_x_error'])
        error_beta_y.append(self.__data[i][j]['beta_y_error'])
        error_alpha.append(self.__data[i][j]['alpha_error'])
        error_alpha_x.append(self.__data[i][j]['alpha_x_error'])
        error_alpha_y.append(self.__data[i][j]['alpha_y_error'])
        error_momentum.append(self.__data[i][j]['momentum_error'])

        zeros.append(0.0)
        position.append(self.__data[i][j]['position'])


    position_copy = copy.copy(position)

    position, emittance, emittance_x, emittance_y, alpha, alpha_x, alpha_y, beta, beta_x, beta_y, momentum, number_particles = tools.sort_arrays(\
        [position, emittance, emittance_x, emittance_y, alpha, alpha_x, alpha_y, beta, beta_x, beta_y, momentum, number_particles], 0)

    position_copy, error_emittance, error_emittance_x, error_emittance_y, error_alpha, error_alpha_x, error_alpha_y, error_beta, error_beta_x, error_beta_y, error_momentum = tools.sort_arrays(\
        [position_copy, error_emittance, error_emittance_x, error_emittance_y, error_alpha, error_alpha_x, error_alpha_y, error_beta, error_beta_x, error_beta_y, error_momentum], 0)


    self.__graphs['emittance'] = ROOT.TGraphErrors(len(position), position, emittance, zeros, error_emittance)
    self.__graphs['emittance_x'] = ROOT.TGraphErrors(len(position), position, emittance_x, zeros, error_emittance_x)
    self.__graphs['emittance_y'] = ROOT.TGraphErrors(len(position), position, emittance_y, zeros, error_emittance_y)

    self.__graphs['beta'] = ROOT.TGraphErrors(len(position), position, beta, zeros, error_beta)
    self.__graphs['beta_x'] = ROOT.TGraphErrors(len(position), position, beta_x, zeros, error_beta_x)
    self.__graphs['beta_y'] = ROOT.TGraphErrors(len(position), position, beta_y, zeros, error_beta_y)

    self.__graphs['alpha'] = ROOT.TGraphErrors(len(position), position, alpha, zeros, error_alpha)
    self.__graphs['alpha_x'] = ROOT.TGraphErrors(len(position), position, alpha_x, zeros, error_alpha_x)
    self.__graphs['alpha_y'] = ROOT.TGraphErrors(len(position), position, alpha_y, zeros, error_alpha_y)
      
    self.__graphs['momentum'] = ROOT.TGraphErrors(len(position), position, momentum, zeros, error_momentum)
    self.__graphs['number_particles'] = ROOT.TGraphErrors(len(position), position, number_particles, zeros, zeros)

    

