
import ROOT
import argparse

class Analysis_Base(object) :

  def __init__(self, analysis_name, require_mc=False) :
    self.__name = analysis_name
    self.__require_mc = require_mc


  def require_mc(self) :
    return self.__require_mc


  def get_plots(self) :
    plot_dict = {}

    self._get_plots(plot_dict)

    return self.__name, plot_dict


  def get_data(self) :
    data_dict = {}

    self._get_data( data_dict )

    return self.__name, data_dict


  def analyse_event(self, analysis_event) :
    raise NotImplementedError()


  def _get_plots(self, plot_dict) :
    raise NotImplementedError()


  def _get_data(self, data_dict) :
    raise NotImplementedError()


  def configure_arguments(self, parser) :
    raise NotImplementedError()


  def parse_arguments(self, namespace) :
    raise NotImplementedError()

