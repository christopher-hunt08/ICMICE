
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array
import os
import numpy


class CoolingExtractor(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "cooling_data_extractor")

    self.__extract_plane = None

    self.__output_filename = None
    self.__the_upstream_data = []
    self.__the_downstream_data = []
    self.__the_upstream_weights = []
    self.__the_downstream_weights = []


  def analyse_event(self, analysis_event, weight=1.0) :
    if analysis_event.num_upstream_tracks() == 1 :
      upstream_hit = self._extract_upstream(self.__extract_plane, analysis_event)
      self.__the_upstream_data.append( [ upstream_hit.get_x(), upstream_hit.get_px(), upstream_hit.get_y(), upstream_hit.get_py(), upstream_hit.get_z(), upstream_hit.get_pz() ] )
      self.__the_upstream_weights.append( weight )

    if analysis_event.num_downstream_tracks() == 1 :
      downstream_hit = self._extract_downstream(self.__extract_plane, analysis_event)
      self.__the_downstream_data.append( [ downstream_hit.get_x(), downstream_hit.get_px(), downstream_hit.get_y(), downstream_hit.get_py(), downstream_hit.get_z(), downstream_hit.get_pz() ] )
      self.__the_downstream_weights.append( weight )


  def _extract_upstream(self, plane, event) :
    return event.upstream_track()[plane]


  def _extract_downstream(self, plane, event) :
    return event.downstream_track()[plane]


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--scifi_extraction_plane", type=int, default=-1, help="Specify the tracker plane at which to extract the muon information" )
    parser.add_argument( "--scifi_extraction_filename", default="cooling_data", help="Output file name for the extracted beam" )


  def parse_arguments(self, namespace) :
    self.__extract_plane = abs(namespace.scifi_extraction_plane)
    self.__output_filename = os.path.join(namespace.output_directory, namespace.scifi_extraction_filename)+".npz"


  def conclude(self) :
    with open(self.__output_filename, 'w') as outfile :
      numpy.savez(outfile, upstream_data=self.__the_upstream_data, upstream_weights=self.__the_upstream_weights, downstream_data=self.__the_downstream_data, downstream_weights=self.__the_downstream_weights)


