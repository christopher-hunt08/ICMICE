
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array


class MCSciFiEfficiency(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "mc_scifi_efficiency")

    self.__inspectors = []

    self.__required_number_stations = 0
    self.__track_algorithm = -1

    self.__upstream_found = 0
    self.__upstream_expected = 0
    self.__downstream_found = 0
    self.__downstream_expected = 0

    self.__upstream_found_inspector = inspectors.PhaseSpace2DInspector(-1000, 0)
    self.__upstream_expected_inspector = inspectors.PhaseSpace2DInspector(-1001, 0)
    self.__downstream_found_inspector = inspectors.PhaseSpace2DInspector(1000, 0)
    self.__downstream_expected_inspector = inspectors.PhaseSpace2DInspector(1001, 0)

    self.__plots = {'upstream' : {}, 'downstream' : {}}
    self.__data = {'upstream' : {}, 'downstream' : {}}


  def analyse_event(self, analysis_event, weight=1.0) :
    pass


  def _get_plots(self, plot_dict) :
    plot_dict['upstream'] = self.__plots['upstream']
    plot_dict['downstream'] = self.__plots['downstream']


  def _get_data(self, data_dict) :
    data_dict['upstream'] = self.__data['upstream']
    data_dict['downstream'] = self.__data['downstream']


  def configure_arguments(self, parser) :
    parser.add_argument( '--required_number_stations', default=4, type=int, help="Number of stations an MC track must cross in order to require a track" )
    parser.add_argument( '--track_algorithm', default=1, choices=[0,1], help="Tracker reconstruction algorithm to require" )


  def parse_arguments(self, namespace) :
    self.__required_number_stations = namespace.required_number_stations
    self.__track_algorithm = namespace.track_algorithm


  def conclude(self) :
    self.__plots['upstream'] = { 'expected': self.__upstream_expected_inspector.get_plot_dictionary(), 'found': self.__upstream_found_inspector.get_plot_dictionary() }
    self.__plots['downstream'] = { 'expected': self.__downstream_expected_inspector.get_plot_dictionary(), 'found': self.__downstream_found_inspector.get_plot_dictionary() }

    self.__data['upstream'] = { 'expected': self.__upstream_expected_inspector.get_data_dictionary(), 'found': self.__upstream_found_inspector.get_data_dictionary() }
    self.__data['downstream'] = { 'expected': self.__downstream_expected_inspector.get_data_dictionary(), 'found': self.__downstream_found_inspector.get_data_dictionary() }

