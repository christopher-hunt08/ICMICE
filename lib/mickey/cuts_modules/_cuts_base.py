

import ROOT
import argparse


class Cut_Base(object) :

  def __init__(self, cut_name, require_mc=False) :
    self.__name = cut_name
    self.__cut_counter = 0
    self.__event_counter = 0
    self.__require_mc = require_mc


  def require_mc(self) :
    return self.__require_mc


  def is_cut(self, analysis_event) :
    self.__event_counter += 1

    if self._is_cut( analysis_event ) :
      self.__cut_counter += 1
      return True
    else :
      return False


  def get_plots(self) :
    plot_dict = {}

    self._get_plots(plot_dict)

    return self.__name, plot_dict


  def get_data(self) :
    data_dict = {}
    data_dict["number_events"] = self.__event_counter
    data_dict["number_cut"] = self.__cut_counter

    self._get_data( data_dict )

    return self.__name, data_dict


  def _is_cut(self, analysis_event) :
    raise NotImplementedError()


  def fill_histograms(self, analysis_event) :
    raise NotImplementedError()


  def _get_plots(self, plot_dict) :
    raise NotImplementedError()


  def _get_data(self, data_dict) :
    raise NotImplementedError()


  def configure_arguments(self, parser) :
    raise NotImplementedError()


  def parse_arguments(self, namespace) :
    raise NotImplementedError()

