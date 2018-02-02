

import ROOT
import os
import json
import glob
import argparse
import event_loader

import analysis.tools as tools

import _parsing
import event

## Some global variables to keep track of things between modules
_last_json = None
_last_root = None
_QUIET = False


## Define the engine object
class Engine(object) :

####################################################################################################
  def __init__(self, job_name, description="") :
    self.__analysis_name = job_name

    self.__parser = argparse.ArgumentParser(description=description, conflict_handler='resolve')
    _parsing.get_engine_parser(self.__parser, job_name)
    self.__file_reader = None
    self.__select_events = False
    self.__save_good_events = False
    self.__do_cuts = True
    self.__do_analysis = True
    self.__use_mc = False
    self.__mc_lookup = None

    self.__cuts = []
    self.__analyses = []

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
  def add_cut(self, cut) :
    cut.configure_arguments(self.__parser)
    if cut.require_mc() :
      self.__use_mc = True
    self.__cuts.append(cut)


####################################################################################################
  def add_analysis(self, analysis) :
    analysis.configure_arguments(self.__parser)
    if analysis.require_mc() :
      self.__use_mc = True
    self.__analyses.append(analysis)


####################################################################################################
  def go(self, conclude=True, save_plots=True, save_data=True) :
    """
      Run through all the events in the file reader, before running the conclusion if required.
    """

    if self.__use_mc :
      if self.__mc_lookup is None :
        raise ValueError("MC Lookup File not specified")

    try :
      while self.next_event() :
        try :

          self.analyse_event()

        except ValueError as ex:
          print
          print "An Error Occured. Skipping Spill..."
          print "ERROR =", ex
          print
          continue
    except KeyboardInterrupt :
      print
      print "Keyboard Interrupt"
      print

    if conclude :
      try :
        self.conclude()
      except ValueError as ex :
        print "Analysis Failed:", ex
        print
        print "Stopping Execution"

    if save_plots :
      self.save_plots()

    if save_data :
      self.save_data()


####################################################################################################
  def analyse_event(self) :
    self.__event_counter += 1
#    set_event_statistical_weight(self.__file_reader.get_current_statistical_weight())

    if self.__use_mc :
      maus_event = event.build_event(self.__file_reader.get_event, self.__mc_lookup)
    else :
      maus_event = event.build_event(self.__file_reader.get_event)

    is_ok = True

    if self.__do_cuts :
      for cut in self.__cuts :
        cut.fill_histograms(maus_event)
        if cut.is_cut(maus_event) :
          is_ok = False

    if not is_ok :
      return

    if self.__save_good_events :
      self.__file_reader.save_event()


#    maus_event.print_me()

    if self.__do_analysis :
      for analysis in self.__analyses :
        analysis.analyse_event(maus_event)

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
    for proc in self.__analyses :
      proc.conclude()

    if self.__save_good_events :
      filename = os.path.join(self.__output_directory, self.__output_filename)
      with open(filename+'-good_events.json', 'w') as outfile :
        json.dump(self.__file_reader.get_saved_events(), outfile)


####################################################################################################
  def get_plot_dict(self, plot_dict = None) :
    if plot_dict is None :
      plot_dict = {}

    cut_dict = {}
    analysis_dict = {}

    for proc in self.__cuts :
      name, plots = proc.get_plots()
      cut_dict[name] = plots
    for ana in self.__analyses :
      name, plots = ana.get_plots()
      analysis_dict[name] = plots

    plot_dict["cuts"] = cut_dict
    plot_dict["analysis"] = analysis_dict

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

    cut_dict = {}
    analysis_dict = {}

    data_dict['events_analysed'] = self.__event_counter
    data_dict['arguments'] = vars(self.__namespace)

    for proc in self.__cuts :
      name, data = proc.get_data()
      cut_dict[name] = data
    for ana in self.__analyses :
      name, data = ana.get_data()
      analysis_dict[name] = data

    data_dict["cuts"] = cut_dict
    data_dict["analysis"] = analysis_dict

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


#    if self.__namespace.last_analysis is not None :
#      with open(self.__namespace.last_analysis+'.json', 'r') as infile :
#        _last_json = json.load(infile)
#      _last_root = ROOT.TFile(self.__namespace.last_analysis+'.root', 'READ')


    if self.__namespace.no_cuts :
      self.__do_cuts = False
    else :
      for cut in self.__cuts :
        cut.parse_arguments(self.__namespace)

    if self.__namespace.no_analysis :
      self.__do_analysis = False
    else :
      for analysis in self.__analyses :
        analysis.parse_arguments(self.__namespace)

    if self.__namespace.virtual_plane_lookup is not None :
      with open(self.__namespace.virtual_plane_lookup, "r") as infile :
        temp_dict = json.load(infile)
      self.__mc_lookup = {}
      for key, item in temp_dict.iteritems() :
        self.__mc_lookup[int(key)] = int(item)

    load_events = None
    if self.__namespace.selection_file is not None :
      self.__select_events = True
      load_events = {}
      for filename in self.__namespace.selection_file :
        with open(filename, 'r') as infile :
          events = json.load(infile)
        load_events.update(events)

    self.__save_good_events = self.__namespace.save_good_events
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


