
import ROOT
import json
import os
import math
import array
import numpy
import random

import framework
import scifi_extractors
import tof_analysis
import virtuals_analysis
import analysis
from analysis import tools
from analysis import covariances
from analysis import hit_types
from analysis import beam_sampling


class scifi_produce_corrections(framework.processor_base) :

  def __init__(self) :
    framework.processor_base.__init__(self, "scifi_correction_production")

    self.__output_filename = None
    self.__output_directory = None

    self.__track_extractors = []
    self.__mc_track_extractor = None
    self.__mc_scifi_dictionary = None

    self.__corrections = {}



  def get_dependencies(self, inserter) :
    self.__track_extractors = [ inserter(scifi_extractors.scifi_helical_track_processor(0, False)), inserter(scifi_extractors.scifi_helical_track_processor(1, False)) ]
    self.__mc_track_extractor = inserter(virtuals_analysis.virtual_beam_properties())


  def get_args(self, parser) :

    parser.add_argument( '--scifi_virt_dict', type=str, default=None, help='Specify the json file name to load the ScifiPlane:VirtualPlane Dictionary.')

    parser.add_argument( '--save_correction_plane', type=int, nargs=3, action='append', help="The <tracker> <station> <plane> to perform the analysis" )


  def process_args(self, namespace) :
    self.__output_filename = namespace.output_filename
    self.__output_directory = namespace.output_directory

    if namespace.scifi_virt_dict is not None :
      with open(namespace.scifi_virt_dict, 'r') as infile :
        self.__mc_scifi_dictionary = json.load(infile)
    else :
      print 
      print "ERROR. No scifi-virtual dictionary found."
      raise ValueError("Please specifiy a SciFi-Virtual plane lookup dictionary to use the correction matrix functionality")


    for tracker, station, plane in namespace.save_correction_plane :
      if tracker not in self.__corrections :
        self.__corrections[tracker] = {}
      if station not in self.__corrections[tracker] :
        self.__corrections[tracker][station] = {}

      self.__corrections[tracker][station][plane] = covariances.CorrectionMatrix()


  def _reset(self) :
    pass


  def _process(self, file_reader) :

    if 0 in self.__corrections :
      for station in self.__corrections[0] :
        for plane in self.__corrections[0][station] :

          trackpoint = self.__track_extractors[0].get_trackpoint_byplane(station, plane)
          if trackpoint is None :
            print "No SciFi"
            continue

          virt_plane = self.__mc_scifi_dictionary[str(tools.calculate_plane_id(0, station, plane))][0]
          virt = self.__mc_track_extractor.get_virtual_hits(virt_plane)
#          print virt
          if virt is None :
            print "No Virt!"
            continue

          recon = analysis.hit_types.AnalysisHit(scifi_track_point=trackpoint)
          mc = analysis.hit_types.AnalysisHit(virtual_track_point=virt[1])

#          if numpy.isinf(recon.get_x()) or numpy.isnan(recon.get_x()) :
#            print "FOUND RECON NAN/INF"
#            continue
#          if numpy.isinf(mc.get_x()) or numpy.isnan(mc.get_x()) :
#            print "FOUND MC NAN/INF"
#            continue

          self.__corrections[0][station][plane].add_hit(recon, mc)


    if 1 in self.__corrections :
      for station in self.__corrections[1] :
        for plane in self.__corrections[1][station] :

          trackpoint = self.__track_extractors[1].get_trackpoint_byplane(station, plane)
          if trackpoint is None :
            print "No SciFi"
            continue

          virt_plane = self.__mc_scifi_dictionary[str(tools.calculate_plane_id(1, station, plane))][0]
          virt = self.__mc_track_extractor.get_virtual_hits(virt_plane)
          if virt is None :
            print "No Virt!"
            continue

          recon = analysis.hit_types.AnalysisHit(scifi_track_point=trackpoint)
          mc = analysis.hit_types.AnalysisHit(virtual_track_point=virt[1])

          self.__corrections[1][station][plane].add_hit(recon, mc)

    return False


  def conclude(self) :
    for tracker in self.__corrections :
      for station in self.__corrections[tracker] :
        for plane in self.__corrections[tracker][station] :

          print
          print "PLANE:", tracker, station, plane
          print "Number of Hits: ", self.__corrections[tracker][station][plane].number_hits()
          matrix = self.__corrections[tracker][station][plane].get_full_correction(['x', 'px', 'y', 'py'])

          filename = os.path.join(self.__output_directory, self.__output_filename+"_{0:d}.{1:d}.{2:d}.dat".format(tracker, station, plane))

          numpy.savetxt(filename, matrix)


  def _store_plots(self, plot_dict) :
    return plot_dict


  def _store_data(self, data_dict) :
    return data_dict


