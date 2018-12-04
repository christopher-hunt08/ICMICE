
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


class VoronoiPhaseSpaceSelection(Selection_Base) :

  def __init__(self) :
    Selection_Base.__init__(self, "voronoi_phasespace_selection")

    self.__data_file_number = 0
    self.__voronoi_data_file = LastAnalysis.LastData['analysis']['Voronoi_Tessellation']['voronoi_data']

    self.__event_counter = 0
    self.__current_file_events = 0

    self.advance_voronoi_data_file()



  def get_normalisation(self) :
    return self.__normalisation


  def weigh_event(self, event) :
    if self.__event_counter >= self.__current_file_events :
      self.advance_voronoi_data_file()

    weight = self.__weights[self.__event_counter]
    self.__event_counter += 1

    if weight < 1.0e-3 : 
      return 0.0
    else :
      return weight


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass


  def advance_voronoi_data_file(self) :
    voronoi_data = None
    filename = self.__voronoi_data_file+".{0:05d}".format(self.__data_file_number)
    if not os.path.exists(filename) :
      raise StopIteration
    else :
      with open(filename, 'r') as infile :
        voronoi_data = json.load(infile)

      self.__weights = voronoi_data['weights']
  #    self.__densities = voronoi_data['densities']
  #    self.__vertices = voronoi_data['vertices']
  #    self.__regions = voronoi_data['regions']
  #    self.__point_regions = voronoi_data['point_regions']
      self.__normalisation = voronoi_data['normalisation']
      
      self.__current_file_events = len(self.__weights)
      self.__data_file_number += 1
      self.__event_counter = 0


