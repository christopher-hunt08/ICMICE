

import ROOT
import argparse


####################################################################################################
class Selection_Base(object) :

  def __init__(self, selection_name, requires_parent=True) :
    self.__name = selection_name
    self.__requires_parent = requires_parent


  def requires_parent(self) :
    return self.__requires_parent


  def get_plots(self) :
    plot_dict = {}
    self._get_plots(plot_dict)
    return self.__name, plot_dict


  def get_data(self) :
    data_dict = {}
    self._get_data( data_dict )
    return self.__name, data_dict


  def weigh_event(self, analysis_event) :
    raise NotImplementedError()


  def _get_plots(self, plot_dict) :
    raise NotImplementedError()


  def _get_data(self, data_dict) :
    raise NotImplementedError()

