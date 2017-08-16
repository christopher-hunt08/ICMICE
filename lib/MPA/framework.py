#!/usr/bin/env python

# This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
# MAUS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MAUS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MAUS.  If not, see <http://www.gnu.org/licenses/>.
#

# pylint: disable = W0311, E1101, W0613, C0111, R0911, W0621, C0103, R0902

import ROOT
import os
import math
import json
import glob
import argparse
import event_loader

import analysis.tools as tools

"""
MPA : MAUS Python Analysis Framework Modules
"""

_last_json = None
_last_root = None
_event_weight = 0.0

def get_last_analysis_json() :
  return _last_json


def get_last_analysis_root() :
  return _last_root


def set_event_statistical_weight(weight) :
  global _event_weight
  _event_weight = weight


def get_event_statistical_weight() :
  return _event_weight


class analysis_engine(object) :
  """
    Container class to store all the selected processor objects, cache 
    pre-processed data for the classes to share, and combine cuts and plots
    for easy user access.

    This essentially the analysis engine.
  """
  def __init__(self, job_name, description="") :
    self.__analysis_name = job_name
    self.__description = description
    self.__processors = []
    self.__parser = argparse.ArgumentParser(description=description, conflict_handler='resolve')
    self.__file_reader = None
    self.__select_events = False
    self.__save_good_events = False

    self.__max_num_events = 0
    self.__event_counter = 0
    self.__event_weight = 0.0

    self.__output_directory = None
    self.__output_filename = None
    self.__print_plots = False

    self.__last_json = None
    self.__last_root = None

    self.__parser.add_argument( 'maus_root_files', nargs='+', \
                            help='List of MAUS output root files containing '+\
                                               'reconstructed data.')

    self.__parser.add_argument( '-q', '--quiet', action='store_true', help="Reduce the verbosity of the output")

    self.__parser.add_argument( '-f', '--find_files', action="store_true", \
                                   help='Assume the maus_root_files variables are directories and go looking inside for the individual root files')

    self.__parser.add_argument( '-N', '--max_num_events', type=int, \
                                   help='Maximum number of events to analyse.')

    self.__parser.add_argument( '-O', '--output_filename', \
            default=job_name, help='Set the output filename')

    self.__parser.add_argument( '-D', '--output_directory', \
                                 default='./', help='Set the output directory')
    self.__parser.add_argument( '-P', '--print_plots', action='store_true', \
                        help="Flag to save the plots as individual pdf files" )
    self.__parser.add_argument( '--selection_file', default=None, nargs='+', help='JSON files with a list of spill and event numbers to include in the analysis' )
    self.__parser.add_argument( '-S', '--save_good_events', type=str, default=None, help='Save the good events from the specified analysis module to file' )
    self.__parser.add_argument( '--list_modules', action='store_true', help='Prints out a list of all the modules used in this analysis' )
    self.__parser.add_argument( '--mass_assumption', default=tools.MUON_MASS, type=float, help='Default mass to assume for all tracks' )
    self.__parser.add_argument( '--last_analysis', default=None, help='Base name of the previous analysis run - some modules require this knowledge.' )
#    parser.add_argument( '-D', '--split_by_directory', action='store_true', \
#        help='Tells the script to split the supplied list of root files into "\
#              different data series, depending on the directory they are in.' )


  def get_argparser(self) :
    return self.__parser


  def process_arguments(self) :
    self.__namespace = self.__parser.parse_args()


    if self.__namespace.last_analysis is not None :
      global _last_json
      global _last_root
      with open(self.__namespace.last_analysis+'.json', 'r') as infile :
        _last_json = json.load(infile)
      _last_root = ROOT.TFile(self.__namespace.last_analysis+'.root', 'READ')


    if self.__namespace.list_modules :
      self.print_modules()
      raise SystemExit(0)

    for proc in self.__processors :
      proc.process_args(self.__namespace)

    load_events = None
    if self.__namespace.selection_file is not None :
      self.__select_events = True
      load_events = {}
      for filename in self.__namespace.selection_file :
        with open(filename, 'r') as infile :
          events = json.load(infile)
        load_events.update(events)

    self.__good_events_module = self.__namespace.save_good_events
    if self.__good_events_module is not None :
      self.__save_good_events = True

    self.__output_filename = self.__namespace.output_filename
    self.__output_directory = self.__namespace.output_directory
    self.__mass_assumption = self.__namespace.mass_assumption

    self.__print_plots = self.__namespace.print_plots

    self.__quiet_out = self.__namespace.quiet


    if self.__namespace.find_files :
      self.__maus_root_files = []
      for directory in self.__namespace.maus_root_files :
        self.__maus_root_files.extend(glob.glob(os.path.join(directory, "*.root")))
    else :
      self.__maus_root_files = self.__namespace.maus_root_files

    if self.__quiet_out :
      self.__file_reader = event_loader.maus_reader(self.__maus_root_files, print_progress='file')
    else :
      self.__file_reader = event_loader.maus_reader(self.__maus_root_files, print_progress='spill')
    self.__file_reader.set_max_num_events(self.__namespace.max_num_events)
    self.__file_reader.select_events(load_events)
    return self.__namespace


  def print_modules(self) :
    print
    print "The modules used in this analysis are:"
    print
    for proc in self.__processors :
      print " o ", proc.get_name()


  def set_output_filename(self, filename) :
    self.__output_filename = filename


  def set_output_directory(self, directory) :
    self.__output_directory = directory


  def get_file_reader(self) :
    return self.__file_reader


  def next_event(self) :
    return self.__file_reader.next_selected_event()


  def _prepend_processor(self, processor) :
    for proc in self.__processors :
      if proc.get_name() == processor.get_name() :
        return proc
    else :
      self.__processors.insert(0, processor)
      processor.get_dependencies(self._prepend_processor)
      processor.get_args(self.__parser)
      return proc


  def add_processor(self, processor) :
    for proc in self.__processors :
      if proc.get_name() == processor.get_name() :
        return proc
    else :
      processor.get_dependencies(self.add_processor)
      self.__processors.append(processor)
      proc = self.__processors[-1]
      proc.set_analysis_name(self.__analysis_name)
      proc.get_args(self.__parser)
      return proc


  def analyse_event(self, file_reader = None) :
    reader = None
    self.__event_counter += 1
    if file_reader is not None :
      reader = file_reader
    else :
      reader = self.__file_reader

    set_event_statistical_weight(reader.get_current_statistical_weight())

    for proc in self.__processors :
      proc.reset()

    for proc in self.__processors :
      proc.process(reader)
      if proc.is_cut() :
        set_event_statistical_weight(0.0)

    if self.__save_good_events :
      result = True
      for proc in self.__processors :
        if proc.get_name() == self.__good_events_module :
          if not proc.is_cut() and (get_event_statistical_weight() > 1.0e-9) :
            reader.save_event(get_event_statistical_weight())
          break

 
  def conclude(self) :
    for proc in self.__processors :
      proc.conclude()

    if self.__save_good_events :
      filename = os.path.join(self.__output_directory, self.__output_filename)
      with open(filename+'-good_events.json', 'w') as outfile :
        json.dump(self.__file_reader.get_saved_events(), outfile)


#  def print_processors(self) :
#    for count, proc in enumerate(self.__processors) :
#      print count, ":", proc.get_name()


  def get_plot_dict(self, plot_dict = None) :
    if plot_dict is None :
      plot_dict = {}

    for proc in self.__processors :
      proc.get_plot_dict(plot_dict)

    return plot_dict

  def save_plots(self, outfilename = None, print_plots = False) :
    if outfilename is not None :
      filename = outfilename
    else :
      filename = os.path.join(self.__output_directory, \
                                                self.__output_filename+".root")

    plot_dict = self.get_plot_dict()
    tools.save_plots(plot_dict, filename)

    if print_plots or self.__print_plots :
      tools.print_plots(plot_dict, self.__output_directory)


  def get_data_dict(self, data_dict = None) :
    if data_dict is None :
      data_dict = {}

    data_dict['events_analysed'] = self.__event_counter
    data_dict['arguments'] = vars(self.__namespace)

    for proc in self.__processors :
      proc.get_data_dict(data_dict)

    return data_dict


  def save_data(self, outfilename = None) :
    if outfilename is not None :
      filename = outfilename
    else :
      filename = os.path.join(self.__output_directory, \
                                                self.__output_filename+".json")

    with open(filename, 'w') as outfile :
      json.dump(self.get_data_dict(), outfile)


  def find_processor(self, proc_name) :
    found = None
    for proc in self.__processors :
      if proc.get_name() == proc_name :
        found = proc
        break
    return found


################################################################################


class processor_base(object) :
  """
    Small classes (essentially functors) to process little pieces of the 
    events and perform individual analysis routines.
  """
  def __init__(self, name) :
    self.__name = name
    self.__analysis_name = ""
    self.__counter_cut = 0
    self.__counter_events = 0
    self.__is_cut = False
    self.__processed = False


  def set_analysis_name(self, ana_name) :
    self.__analysis_name = ana_name


  def get_analysis_name(self) :
    return self.__analysis_name


  def get_name(self) :
    return self.__name


  def get_number_events(self) :
    return self.__counter_events


  def is_cut(self) :
    return self.__is_cut


  def get_number_cut(self) :
    return self.__counter_cut


  def reset(self) :
    self.__is_cut = False
    self.__processed = False
    self._reset()


  def process(self, file_reader) :
    self.__counter_events += 1
    if not self.__processed : 
      result = self._process(file_reader)
      self.__is_cut = result
      self.__processed = True


  def is_processed(self) :
    return self.__processed


  def get_plot_dict(self, plot_dict=None) :
    if plot_dict is None :
      plot_dict = {}

    return self._store_plots(plot_dict)


  def get_data_dict(self, data_dict=None) :
    if data_dict is None :
      data_dict = {}

    return self._store_data(data_dict)


  def _cut(self) :
    if not self.__is_cut :
      self.__is_cut = True
      self.__counter_cut += 1


  def get_dependencies(self, inserter) :
    raise NotImplementedError(\
                       "This function should be overloaded in a derived class")


  def get_args(self, parser) :
    pass


  def process_args(self, parser) :
    pass


  def _reset(self) :
    pass


  def _process(self, file_reader) :
    raise NotImplementedError(\
                       "This function should be overloaded in a derived class")
    return False


  def conclude(self) :
    pass


  def _store_plots(self, plot_dict) :
    return plot_dict


  def _store_data(self, data_dict) :
    return data_dict



