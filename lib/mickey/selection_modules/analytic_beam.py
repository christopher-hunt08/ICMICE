
from _selection_base import Selection_Base

from .. import LastAnalysis

import analysis
from analysis import beam_sampling

import ROOT
import numpy


class SelectAnalyticBeam(Selection_Base) :

  def __init__(self, emittance, alpha, beta, momentum, mass=analysis.tools.MUON_MASS) :
    Selection_Base.__init__(self, "analytic_beam_selection")

    self.__parent_means = numpy.array([LastAnalysis.LastData['beam_selection']['parent_analysis']['x_mean'], \
                           LastAnalysis.LastData['beam_selection']['parent_analysis']['px_mean'], \
                           LastAnalysis.LastData['beam_selection']['parent_analysis']['y_mean'], \
                           LastAnalysis.LastData['beam_selection']['parent_analysis']['py_mean'] ])

    self.__parent_covariance = numpy.array(LastAnalysis.LastData['beam_selection']['parent_analysis']['covariance_matrix'])

    cov_xx = beta*emittance*mass/momentum
    cov_xp = -1.0*alpha*emittance*mass
    cov_pp = ((1.0 + alpha**2) / beta) * emittance * momentum * mass

    self.__means = numpy.array( [ 0.0, 0.0, 0.0, 0.0 ] )
    self.__covariance = numpy.array( [[ cov_xx, cov_xp, 0.0, 0.0 ], \
                                      [ cov_xp, cov_pp, 0.0, 0.0 ], \
                                      [ 0.0, 0.0, cov_xx, cov_xp ], \
                                      [ 0.0, 0.0, cov_xp, cov_pp ]] )

    self.__sampler = beam_sampling.Gaussian4DSampler(self.__parent_means, self.__parent_covariance, \
                                                     self.__means, self.__covariance)


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    weight = self.__sampler.weight(hit)

    return weight


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass

