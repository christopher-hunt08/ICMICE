

import ROOT
import os
import json
import glob
import argparse
import event_loader

import analysis.tools as tools

import _parsing

## Some global variables to keep track of things between modules
_last_json = None
_last_root = None
_QUIET = False


## Define the engine object
class Engine(object) :

####################################################################################################
  def __init__(self, job_name) :
    self.__analysis_name = job_name

    self.__parser = argparse.ArgumentParser(description=description, conflict_handler='resolve')
    self.__file_reader = None
    self.__select_events = False
    self.__save_good_events = False
    self.__do_cuts = True
    self.__do_analysis = True

    self.__max_num_events = 0
    self.__event_counter = 0
    self.__event_weight = 0.0

    self.__output_directory = None
    self.__output_filename = None
    self.__print_plots = False

    self.__last_json = None
    self.__last_root = None


####################################################################################################
  def get_argparser(self) :
    return self.__parser


####################################################################################################
  def analyse_event(self) :
    self.__event_counter += 1
    set_event_statistical_weight(self.__file_reader.get_current_statistical_weight())


    event = event.build_event(self.__file_reader.get_event)

    if self.__do_cuts :


    if self.__do_analysis :


   ## MAGIC GOES HERE 

      # Build event
      # Analyse pre-cuts
      # Perform cuts
      # Analyse Post Cuts
      # Perform Analyses...





####################################################################################################
  def set_output_filename(self, filename) :
    self.__output_filename = filename


####################################################################################################
  def set_output_directory(self, directory) :
    self.__output_directory = directory


####################################################################################################
  def get_file_reader(self) :
    return self.__file_reader


####################################################################################################
  def next_event(self) :
    return self.__file_reader.next_selected_event()


####################################################################################################
  def conclude(self) :
    for proc in self.__processors :
      proc.conclude()

    if self.__save_good_events :
      filename = os.path.join(self.__output_directory, self.__output_filename)
      with open(filename+'-good_events.json', 'w') as outfile :
        json.dump(self.__file_reader.get_saved_events(), outfile)


####################################################################################################
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


####################################################################################################
  def get_plot_dict(self, plot_dict = None) :
    if plot_dict is None :
      plot_dict = {}

    for proc in self.__processors :
      proc.get_plot_dict(plot_dict)

    return plot_dict


####################################################################################################
  def save_plots(self, outfilename = None, print_plots = False) :
    if outfilename is not None :
      filename = outfilename
    else :
      filename = os.path.join(self.__output_directory, self.__output_filename+".root")

    plot_dict = self.get_plot_dict()
    tools.save_plots(plot_dict, filename)

    if print_plots or self.__print_plots :
      tools.print_plots(plot_dict, self.__output_directory)


####################################################################################################
  def get_data_dict(self, data_dict = None) :
    if data_dict is None :
      data_dict = {}

    data_dict['events_analysed'] = self.__event_counter
    data_dict['arguments'] = vars(self.__namespace)

    for proc in self.__processors :
      proc.get_data_dict(data_dict)

    return data_dict


####################################################################################################
  def save_data(self, outfilename = None) :
    if outfilename is not None :
      filename = outfilename
    else :
      filename = os.path.join(self.__output_directory, self.__output_filename+".json")

    with open(filename, 'w') as outfile :
      json.dump(self.get_data_dict(), outfile)


####################################################################################################
  def process_arguments(self) :
    """
       Analyse the argparse arguments
    """
    self.__namespace = self.__parser.parse_args()
    global _last_json
    global _last_root
    global _QUIET


    if self.__namespace.last_analysis is not None :
      with open(self.__namespace.last_analysis+'.json', 'r') as infile :
        _last_json = json.load(infile)
      _last_root = ROOT.TFile(self.__namespace.last_analysis+'.root', 'READ')

    if namespace.no_cuts :
      self.__do_cuts = False

    if namespace.no_analysis :
      self.__do_analysis = False


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

    _QUIET = self.__namespace.quiet


    if self.__namespace.find_files :
      self.__maus_root_files = []
      for directory in self.__namespace.maus_root_files :
        self.__maus_root_files.extend(glob.glob(os.path.join(directory, "*.root")))
    else :
      self.__maus_root_files = self.__namespace.maus_root_files

    if _QUIET :
      self.__file_reader = event_loader.maus_reader(self.__maus_root_files, print_progress='file')
    else :
      self.__file_reader = event_loader.maus_reader(self.__maus_root_files, print_progress='spill')

    self.__file_reader.set_max_num_events(self.__namespace.max_num_events)
    self.__file_reader.select_events(load_events)
    return self.__namespace


