
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array
import os


class SciFiExtractor(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "scifi_particle_extraction")

    self.__extract_plane = None
    self.__extract_function = None

    self.__outfile = None
    self.__counter = 0


  def analyse_event(self, analysis_event, weight=1.0) :
    hit = self.__extract_function(self.__extract_plane, analysis_event)
    hit.set_weight(weight)

    string = str(self.__counter) + " 1 " + "-2" + " 0 " + str(hit.get_time()) + " 0.0 " + \
             str(1.0e-3*hit.get_x()) + " " + str(1.0e-3*hit.get_y()) + " " + str(1.0e-3*hit.get_z()) + " " + \
             str(1.0e-3*hit.get_px()) + " " + str(1.0e-3*hit.get_py()) + " " + str(1.0e-3*hit.get_pz()) + " 0 0 0\n"
## 0  : Event Number
#value_1 = 1 # Particle Number
#value_2 = -2 # Particle ID
#value_3 = 0 # Status
## 4  : Time
#value_5 = 0.0
## 6  : x
## 7  : y
## 8  : z
## 9  : px
## 10 : py
## 11 : pz
#value_12 = 0
#value_13 = 0
#value_14 = 0
#
    self.__outfile.write(string)
    self.__counter = self.__counter + 1



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
    parser.add_argument( "--scifi_extraction_filename", default="scifi_extracted_beam", help="Output file name for the extracted beam" )


  def parse_arguments(self, namespace) :
    if namespace.scifi_extraction_plane > 0 :
      self.__extract_function = self._extract_downstream
    else :
      self.__extract_function = self._extract_upstream
      
    self.__extract_plane = abs(namespace.scifi_extraction_plane)

    self.__outfile = open(os.path.join(namespace.output_directory, namespace.scifi_extraction_filename+".txt"), 'w')
    self.__outfile.write( "Custom Input file\n" )
    self.__outfile.write( "0. 0. 0. 0. 0. 0. 0. 0.\n" )


  def conclude(self) :
    self.__outfile.close()

