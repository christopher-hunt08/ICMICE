
from _selection_base import Selection_Base

from .. import LastAnalysis

import ROOT
import json
import analysis
import scipy.spatial
from analysis.beam_sampling import multivariate_gaussian
import numpy
import array
import os


class KDESelection(Selection_Base) :

  def __init__(self) :
    Selection_Base.__init__(self, "kde_phasespace_selection")

    self.__data_file_number = 0
    self.__kde_data_file = LastAnalysis.LastData['analysis']['KDE_Calculation']['kde_data']

    self.__event_counter = 0
    self.__current_file_events = 0

    self.advance_kde_data_file()



  def get_normalisation(self) :
    return self.__normalisation


  def weigh_event(self, event) :
    if self.__event_counter >= self.__current_file_events :
      self.advance_kde_data_file()

    weight = self.__weights[self.__event_counter]
    self.__event_counter += 1

    if weight < 1.0e-5 : 
      return 0.0
    else :
      return weight


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass


  def advance_kde_data_file(self) :
    kde_data = None
    filename = self.__kde_data_file+".{0:05d}".format(self.__data_file_number)
    if not os.path.exists(filename) :
      raise StopIteration
    else :
      with open(filename, 'r') as infile :
        kde_data = json.load(infile)

      self.__weights = kde_data['weights']
  #    self.__densities = kde_data['densities']
      self.__normalisation = kde_data['normalisation']
      
      self.__current_file_events = len(self.__weights)
      self.__data_file_number += 1
      self.__event_counter = 0


