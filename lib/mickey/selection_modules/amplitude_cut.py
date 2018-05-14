
from _selection_base import Selection_Base

from .. import LastAnalysis

from analysis import covariances

import ROOT
import math
import array
import numpy


class CutAmplitude(Selection_Base) :

  def __init__(self, amp_cut) :
    Selection_Base.__init__(self, "amplitude_cut")

    self.__parent = LastAnalysis.LastPlots.Get('beam_selection/parent_analysis/amplitude')
    covariance = numpy.array(LastAnalysis.LastData['beam_selection']['parent_analysis']['covariance_matrix'])

    self.__emittance = covariances.emittance_from_matrix(covariance)
    self.__cov_inv = numpy.linalg.inv(covariance)
    self.__amp_cut = amp_cut


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    vec = numpy.array(hit.get_as_vector()[2:6])

    amplitude = self.__emittance * vec.transpose().dot(self.__cov_inv.dot(vec))

    if amplitude > self.__amp_cut :
      return 0.0
    else :
      return 1.0


  def _get_plots(self, plot_dict) :
    pass

  def _get_data(self, data_dict) :
    pass

