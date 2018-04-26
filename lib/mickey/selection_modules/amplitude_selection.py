
from _selection_base import Selection_Base

from .. import LastAnalysis

from analysis import beam_sampling

import ROOT
import math
import array
import numpy


class SelectAmplitude(Selection_Base) :

  def __init__(self, emittance) :
    Selection_Base.__init__(self, "amplitude_selection")

    self.__parent = LastAnalysis.LastPlots.Get('beam_selection/parent_analysis/amplitude')
    covariance = numpy.array(LastAnalysis.LastData['beam_selection']['parent_analysis']['covariance_matrix'])
    self.__sampler = beam_sampling.Amplitude4DSampler(self.__parent, covariance, emittance, max_x=16*emittance)


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    weight = self.__sampler.weight(hit)

    return weight


  def _get_plots(self, plot_dict) :
    pass

  def _get_data(self, data_dict) :
    pass

