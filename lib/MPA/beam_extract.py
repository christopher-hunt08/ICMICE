
import ROOT
import json
import os
import math
import array
import numpy
import random

import framework
import scifi_analysis
import analysis
from analysis import hit_types



class scifi_beam_extractor(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_beam_extractor")

    self.__extract_planes = []
    self.__outfiles = []
    self.__counter = 0

    self.__beam_selector = None

    self.__format_function = None


  def get_dependencies(self, inserter) :
    self.__beam_selector = inserter(scifi_analysis.scifi_beam_selection())

  def get_args(self, parser) :
    parser.add_argument('--extraction_plane', type=int, nargs=3, action='append', help='Specify the <tracker> <station> <plane> at which point the data is extracted.')
    parser.add_argument('--extraction_format', help='Specify the format for the saved beam file. (Only icool003 available at present!')

  def process_args(self, namespace) :
    self.__extract_planes = namespace.extraction_plane

    self.__format_function = self._format_ecal9

    for plane in self.__extract_planes :
      tr = plane[0]
      st = plane[1]
      pl = plane[2]

      iden = "_" + str(tr) + "." + str(st) + "." + str(pl)
      outfile = open(os.path.join(namespace.output_directory, namespace.output_filename+iden+".txt"), 'w')
      outfile.write( "Custom Input file\n" )
      outfile.write( "0. 0. 0. 0. 0. 0. 0. 0.\n" )
      self.__outfiles.append(outfile)


  def _format_ecal9(self, hit) :

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

    return string



  def _process(self, file_reader) :
    if self.__beam_selector.is_cut() :
      self._cut()
      return True


    for outfile, plane in zip( self.__outfiles, self.__extract_planes ) :
      tr = plane[0]
      st = plane[1]
      pl = plane[2]

      tp = self.__beam_selector.get_track_extractor(tr).get_trackpoint_byplane(st, pl)
      if tp is None : continue

      hit = hit_types.AnalysisHit(scifi_track_point=tp)
      outfile.write(self.__format_function(hit))

      self.__counter = self.__counter + 1

    return False


  def conclude(self) :
    for outfile in self.__outfiles :
      outfile.close()


