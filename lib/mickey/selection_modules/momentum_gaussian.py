
from _selection_base import Selection_Base

from .. import LastAnalysis

from analysis import beam_sampling

import ROOT
import math
import array
import numpy


class SelectMomentum(Selection_Base) :

  def __init__(self, mean_p, rms_p) :
    Selection_Base.__init__(self, "momentum_selection")

    self.__p_plot = ROOT.TH1F('gaussian_momentum_p', 'p', 400, 0.0, 400.0 )

    self.__parent = LastAnalysis.LastPlots.Get('beam_selection/parent_analysis/p')
    self.__sampler = beam_sampling.GaussianMomentumSampler(self.__parent, mean_p, rms_p)
    self.__normalisation = self.__sampler.get_selection_normalisation() / self.__sampler.get_weight_normalisation()


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    weight = self.__sampler.weight(hit)
    self.__p_plot.Fill(hit.get_p(), weight)

    return weight

  def get_normalisation(self) :
    return self.__normalisation


  def _get_plots(self, plot_dict) :
    plot_dict['p'] = self.__p_plot


  def _get_data(self, data_dict) :
    pass



class SelectLongMomentum(Selection_Base) :

  def __init__(self, mean_pz, rms_pz) :
    Selection_Base.__init__(self, "momentum_selection")

    self.__pz_plot = ROOT.TH1F('gaussian_momentum_pz', 'pz', 400, 0.0, 400.0 )

    self.__parent = LastAnalysis.LastPlots.Get('beam_selection/parent_analysis/pz')
    self.__sampler = beam_sampling.GaussianSampler(lambda hit : hit.get_pz(), self.__parent, mean_pz, rms_pz)
    self.__normalisation = self.__sampler.get_selection_normalisation() / self.__sampler.get_weight_normalisation()


  def weigh_event(self, event) :
    hit = event.selection_trackpoint()
    weight = self.__sampler.weight(hit)
    self.__pz_plot.Fill(hit.get_pz(), weight)

    return weight

  def get_normalisation(self) :
    return self.__normalisation


  def _get_plots(self, plot_dict) :
    plot_dict['pz'] = self.__pz_plot


  def _get_data(self, data_dict) :
    pass

