
from _selection_base import Selection_Base

from .. import LastAnalysis

from analysis import beam_sampling

import ROOT
import math
import array
import numpy


class SelectUncorrelated4D(Selection_Base) :

  def __init__(self, alpha, beta) :
    Selection_Base.__init__(self, "uncorrelated4D_selection")

    self.__parent_x = LastAnalysis.LastPlots.Get('beam_selection/parent_analysis/x_px')
    self.__parent_y = LastAnalysis.LastPlots.Get('beam_selection/parent_analysis/y_py')

    emittance = LastAnalysis.LastData['beam_selection']['parent_analysis']['emittance']
    p = LastAnalysis.LastData['beam_selection']['parent_analysis']['momentum']

    means = numpy.array( [0.0, 0.0] )
    covariance = emittance*numpy.array( [[beta, -alpha], [-alpha, p*p*(1.0+alpha*alpha)/beta]] )

    print "Configuring Uncorrelated 4D selection"
    print means
    print
    print covariance
    print
    print emittance

    self.__sampler = beam_sampling.XYPhaseSpaceSampler(self.__parent_x, self.__parent_y, means, covariance)


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    weight = self.__sampler.weight(hit)

    return weight


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass

