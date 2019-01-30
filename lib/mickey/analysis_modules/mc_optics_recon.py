
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array


class MCTruthOpticsRecon(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "mc_truth_optics_reconstruction", require_mc=True)

    self.__inspectors = []
    self.__data = []
    self.__graphs = {}
    self.__PID = -13
    self.__mass 


  def analyse_event(self, analysis_event, weight=1.0) :
    virtual_plane_ids = analysis_event.virtual_plane_ids()

    for plane_id in virtual_plane_ids :
      while plane_id > len(self.__inspectors) :
        self.__inspectors.append(inspectors.PhaseSpace2DInspector(i, 0))

      hit = analysis_event.mc_virtual_hit(plane_id)
      if hit is None :
        continue
      if self.__PID is not None :
        if hit.get_pid() != self.__PID :
          continue
      hit.set_weight(weight)
      self.__inspectors[plane_id].add_hit(hit)


  def _get_plots(self, plot_dict) :
    for graph in self.__graphs :
      plot_dict[graph] = self.__graphs[graph]


  def _get_data(self, data_dict) :
    data_dict['virtual_plane_data'] = self.__data


  def configure_arguments(self, parser) :
    parser.add_argument( "--pid", type=int, default=self.__PID, help="Particle ID to analyse" )


  def parse_arguments(self, namespace) :
    self.__PID = namespace.pid


  def conclude(self) :
    for i in range(len(self.__inspectors)) :
      self.__data[i] = self.__inspectors[i].get_data_dictionary()

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
    kinetic = array.array("d")
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
    error_kinetic = array.array("d")
    error_momentum = array.array("d")
    zeros = array.array("d")
    position = array.array("d")

    for i in self.__inspectors :
      emittance.append(self.__data[i]['emittance'])
      emittance_x.append(self.__data[i]['emittance_x'])
      emittance_y.append(self.__data[i]['emittance_y'])
      beta.append(self.__data[i]['beta'])
      beta_x.append(self.__data[i]['beta_x'])
      beta_y.append(self.__data[i]['beta_y'])
      alpha.append(self.__data[i]['alpha'])
      alpha_x.append(self.__data[i]['alpha_x'])
      alpha_y.append(self.__data[i]['alpha_y'])
      momentum.append(self.__data[i]['momentum'])
      kinetic.append(self.__data[i]['energy']-self.__mass)
      number_particles.append(self.__data[i]['number_particles'])

      error_emittance.append(self.__data[i]['emittance_error'])
      error_emittance_x.append(self.__data[i]['emittance_x_error'])
      error_emittance_y.append(self.__data[i]['emittance_y_error'])
      error_beta.append(self.__data[i]['beta_error'])
      error_beta_x.append(self.__data[i]['beta_x_error'])
      error_beta_y.append(self.__data[i]['beta_y_error'])
      error_alpha.append(self.__data[i]['alpha_error'])
      error_alpha_x.append(self.__data[i]['alpha_x_error'])
      error_alpha_y.append(self.__data[i]['alpha_y_error'])
      error_momentum.append(self.__data[i]['momentum_error'])
      error_kinetic.append(self.__data[i]['energy_error'])

      zeros.append(0.0)
      position.append(self.__data[i]['position'])


    position_copy = copy.copy(position)

    position, emittance, emittance_x, emittance_y, alpha, alpha_x, alpha_y, beta, beta_x, beta_y, momentum, kinetic, number_particles = tools.sort_arrays(\
        [position, emittance, emittance_x, emittance_y, alpha, alpha_x, alpha_y, beta, beta_x, beta_y, momentum, kinetic, number_particles], 0)

    position_copy, error_emittance, error_emittance_x, error_emittance_y, error_alpha, error_alpha_x, error_alpha_y, error_beta, error_beta_x, error_beta_y, error_momentum, error_kinetic = tools.sort_arrays(\
        [position_copy, error_emittance, error_emittance_x, error_emittance_y, error_alpha, error_alpha_x, error_alpha_y, error_beta, error_beta_x, error_beta_y, error_momentum, error_kinetic], 0)


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
    self.__graphs['kinetic'] = ROOT.TGraphErrors(len(position), position, kinetic, zeros, error_kinetic)
    self.__graphs['number_particles'] = ROOT.TGraphErrors(len(position), position, number_particles, zeros, zeros)

    

