
from _selection_base import Selection_Base
#from _selection_base import parent_plots
#from _selection_base import parent_data
from analysis import covariances

from .. import LastAnalysis

import ROOT
import math
import array
import numpy


class ParentAnalysis(Selection_Base) :

  def __init__(self) :
    Selection_Base.__init__(self, "parent_analysis", requires_parent=False)

    self.__covariance = covariances.CovarianceMatrix()

    self.__position_plot = ROOT.TH2F('inspected_position_parent', 'Beam Position', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__momentum_plot = ROOT.TH2F('inspected_momentum_parent', 'Beam Momentum', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__x_phasespace_plot = ROOT.TH2F('inspected_x_phasespace_parent', 'X-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__y_phasespace_plot = ROOT.TH2F('inspected_y_phasespace_parent', 'Y-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__xy_phasespace_plot = ROOT.TH2F('inspected_xpy_phasespace_parent', 'X-Py-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__yx_phasespace_plot = ROOT.TH2F('inspected_ypx_phasespace_parent', 'Y-Px-Phasespace', 100, -400.0, 400.0, 100, -400.0, 400.0)
    self.__rpt_phasespace_plot = ROOT.TH2F('inspected_rpt_phasespace_parent', 'r-Pt-Phasespace', 100, 0.0, 400.0, 100, 0.0, 400.0)
    self.__phi_phasespace_plot = ROOT.TH1F('inspected_phi_phasespace_parent', '#phi-Phasespace', 100, -4.0, 4.0 )
    self.__theta_phasespace_plot = ROOT.TH1F('inspected_theta_phasespace_parent', '#theta-Phasespace', 100, -4.0, 4.0 )
    self.__pz_plot = ROOT.TH1F('inspected_pz_parent', 'p_{z}', 400, 0.0, 400.0 )
    self.__p_plot = ROOT.TH1F('inspected_p_parent', 'p', 400, 0.0, 400.0 )

    self.__parent_covariance = None
    self.__parent_covariance_inv = None
    self.__parent_emittance = 0.0

    self.__amplitude_plot = ROOT.TH1F('single_particle_amplitudes_parent', 'Amplitude', 500, 0.0, 100.0)
    self.__amplitude_momentum_plot = ROOT.TH2F('inspected_A_p_phasespace_parent', 'A-p-Phasespace', 200, 0.0, 100.0, 260, 130.0, 260.0 )


    if LastAnalysis.LastData is not None :
      self.__parent_covariance = numpy.array(LastAnalysis.LastData['beam_selection']['parent_analysis']['covariance_matrix'])
      self.__parent_covariance_inv = numpy.linalg.inv(self.__parent_covariance)
      self.__parent_emittance = covariances.emittance_from_matrix(self.__parent_covariance)


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    self.__covariance.add_hit(hit)

    self.__position_plot.Fill(hit.get_x(), hit.get_y())
    self.__momentum_plot.Fill(hit.get_px(), hit.get_py())
    self.__x_phasespace_plot.Fill(hit.get_x(), hit.get_px())
    self.__y_phasespace_plot.Fill(hit.get_y(), hit.get_py())
    self.__xy_phasespace_plot.Fill(hit.get_x(), hit.get_py())
    self.__yx_phasespace_plot.Fill(hit.get_y(), hit.get_px())
    self.__rpt_phasespace_plot.Fill(hit.get_r(), hit.get_pt())
    self.__phi_phasespace_plot.Fill(hit.get_phi())
    self.__theta_phasespace_plot.Fill(hit.get_theta())
    self.__pz_plot.Fill(hit.get_pz())
    self.__p_plot.Fill(hit.get_p())

    if self.__parent_covariance is not None :
      vector = numpy.array(hit.get_as_vector()[2:6]) # Just the x, px, y, py components
      amplitude = self.__parent_emittance*vector.transpose().dot(self.__parent_covariance_inv.dot(vector))
      self.__amplitude_plot.Fill(amplitude)
      self.__amplitude_momentum_plot.Fill(amplitude, hit.get_p())

    return 1.0


  def _get_plots(self, plot_dict) :
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
    plot_dict['amplitude'] = self.__amplitude_plot
    plot_dict['amplitude_momentum'] = self.__amplitude_momentum_plot


  def _get_data(self, data_dict) :
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

