

from .. import LastAnalysis

from _selection_base import Selection_Base
from analysis import covariances

import ROOT
import math
import array
import numpy


class ParentAnalysis(object) :

  def __init__(self) :
    self.__oosrt = 1.0/math.sqrt(2.0) #One Over Square Root Two
    self.__field = 3.0e-3 # Kilotesla?
    self.__charge = 1.0

    self.__covariance = covariances.CovarianceMatrix()

    self.__position_plot = ROOT.TH2F('position_parent', 'Beam Position', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__momentum_plot = ROOT.TH2F('momentum_parent', 'Beam Momentum', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__x_phasespace_plot = ROOT.TH2F('x_phasespace_parent', 'X-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__y_phasespace_plot = ROOT.TH2F('y_phasespace_parent', 'Y-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__xy_phasespace_plot = ROOT.TH2F('xpy_phasespace_parent', 'X-Py-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__yx_phasespace_plot = ROOT.TH2F('ypx_phasespace_parent', 'Y-Px-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__rpt_phasespace_plot = ROOT.TH2F('rpt_phasespace_parent', 'r-Pt-Phasespace', 100, 0.0, 400.0, 100, 0.0, 400.0)
    self.__phi_phasespace_plot = ROOT.TH1F('phi_phasespace_parent', '#phi-Phasespace', 100, -4.0, 4.0 )
    self.__theta_phasespace_plot = ROOT.TH1F('theta_phasespace_parent', '#theta-Phasespace', 100, -4.0, 4.0 )
    self.__pz_plot = ROOT.TH1F('pz_parent', 'p_{z}', 400, 0.0, 400.0 )
    self.__p_plot = ROOT.TH1F('p_parent', 'p', 400, 0.0, 400.0 )
    self.__L_plot = ROOT.TH1F('angular_momentum_parent', 'L', 2000, -1000.0, 1000.0)
    self.__L_canon_plot = ROOT.TH1F('canonical_angular_momentum_parent', 'L_{canon}', 2000, -1000.0, 1000.0)

    self.__parent_covariance = None
    self.__parent_covariance_inv = None
    self.__parent_emittance = 0.0

    self.__amplitude_plot = ROOT.TH1F('single_particle_amplitudes_parent', 'Amplitude', 1000, 0.0, 200.0)
    self.__amplitude_momentum_plot = ROOT.TH2F('A_p_phasespace_parent', 'A-p-Phasespace', 200, 0.0, 100.0, 260, 130.0, 260.0 )


    if LastAnalysis.LastData is not None :
      self.__parent_covariance = numpy.array(LastAnalysis.LastData['beam_selection']['parent_analysis']['covariance_matrix'])
      self.__parent_covariance_inv = numpy.linalg.inv(self.__parent_covariance)
      self.__parent_emittance = covariances.emittance_from_matrix(self.__parent_covariance)


  def fill_plots(self, event, event_weight) :
    hit = event.selection_trackpoint()
    hit.set_weight(event_weight)
    self.__covariance.add_hit(hit)

    self.__position_plot.Fill(hit.get_x(), hit.get_y(), event_weight)
    self.__momentum_plot.Fill(hit.get_px(), hit.get_py(), event_weight)
    self.__x_phasespace_plot.Fill(hit.get_x(), hit.get_px(), event_weight)
    self.__y_phasespace_plot.Fill(hit.get_y(), hit.get_py(), event_weight)
    self.__xy_phasespace_plot.Fill(hit.get_x(), hit.get_py(), event_weight)
    self.__yx_phasespace_plot.Fill(hit.get_y(), hit.get_px(), event_weight)
    self.__rpt_phasespace_plot.Fill(hit.get_r(), hit.get_pt(), event_weight)
    self.__phi_phasespace_plot.Fill(hit.get_phi(), event_weight)
    self.__theta_phasespace_plot.Fill(hit.get_theta(), event_weight)
    self.__pz_plot.Fill(hit.get_pz(), event_weight)
    self.__p_plot.Fill(hit.get_p(), event_weight)
    self.__L_plot.Fill((hit.get_x()*hit.get_py() - hit.get_y()*hit.get_px())/hit.get_pz(), event_weight)
    self.__L_canon_plot.Fill((hit.get_x()*hit.get_py() - hit.get_y()*hit.get_px() + 0.5*self.__charge*self.__field*hit.get_x()**2*hit.get_y()**2)/hit.get_pz(), event_weight)

    if self.__parent_covariance is not None :
      vector = numpy.array(hit.get_as_vector()[2:6]) # Just the x, px, y, py components
      amplitude = self.__parent_emittance*vector.transpose().dot(self.__parent_covariance_inv.dot(vector))
      self.__amplitude_plot.Fill(amplitude, event_weight)
      self.__amplitude_momentum_plot.Fill(amplitude, hit.get_p(), event_weight)


  def get_plots(self) :
    plot_dict = {}
    plot_dict['x_y'] = self.__position_plot
    plot_dict['px_py'] = self.__momentum_plot
    plot_dict['x_px'] = self.__x_phasespace_plot
    plot_dict['y_py'] = self.__y_phasespace_plot
    plot_dict['x_py'] = self.__xy_phasespace_plot
    plot_dict['y_px'] = self.__yx_phasespace_plot
    plot_dict['r_pt'] = self.__rpt_phasespace_plot
    plot_dict['phi'] = self.__phi_phasespace_plot
    plot_dict['theta'] = self.__theta_phasespace_plot
    plot_dict['pz'] = self.__pz_plot
    plot_dict['p'] = self.__p_plot
    plot_dict['L'] = self.__L_plot
    plot_dict['L_canon'] = self.__L_canon_plot
    plot_dict['amplitude'] = self.__amplitude_plot
    plot_dict['amplitude_momentum'] = self.__amplitude_momentum_plot

    return plot_dict


  def get_data(self) :
    data_dict = {}
    data_dict['x_mean'] = self.__position_plot.GetMean(1)
    data_dict['y_mean'] = self.__position_plot.GetMean(2)
    data_dict['px_mean'] = self.__momentum_plot.GetMean(1)
    data_dict['py_mean'] = self.__momentum_plot.GetMean(2)
    data_dict['x_rms'] = self.__position_plot.GetRMS(1)
    data_dict['y_rms'] = self.__position_plot.GetRMS(2)
    data_dict['px_rms'] = self.__momentum_plot.GetRMS(1)
    data_dict['py_rms'] = self.__momentum_plot.GetRMS(2)

    number = self.__covariance.length()
    if number > 1 :
      cov = self.__covariance.get_covariance_matrix(['x', 'px', 'y', 'py'])
      data_dict['covariance_matrix'] = [ [ cov[i][j] for i in range(4) ] for j in range(4) ]

      data_dict['emittance'] = self.__covariance.get_emittance(['x', 'px', 'y', 'py'])
      data_dict['emittance_x'] = self.__covariance.get_emittance(['x', 'px'])
      data_dict['emittance_y'] = self.__covariance.get_emittance(['y', 'py'])
      data_dict['beta'] = self.__covariance.get_beta(['x', 'y'])
      data_dict['beta_x'] = self.__covariance.get_beta(['x'])
      data_dict['beta_y'] = self.__covariance.get_beta(['y'])
      data_dict['alpha'] = self.__covariance.get_alpha(['x', 'y'])
      data_dict['alpha_x'] = self.__covariance.get_alpha(['x'])
      data_dict['alpha_y'] = self.__covariance.get_alpha(['y'])
      data_dict['momentum'] = self.__p_plot.GetMean()
      data_dict['number_particles'] = number
    else :
      data_dict['covariance_matrix'] = [ [ 0.0 for i in range(4) ] for j in range(4) ]

      data_dict['emittance'] = 0.0
      data_dict['emittance_x'] = 0.0
      data_dict['emittance_y'] = 0.0
      data_dict['beta'] = 0.0
      data_dict['beta_x'] = 0.0
      data_dict['beta_y'] = 0.0
      data_dict['alpha'] = 0.0
      data_dict['alpha_x'] = 0.0
      data_dict['alpha_y'] = 0.0
      data_dict['momentum'] = 0.0
      data_dict['number_particles'] = 0

    return data_dict

