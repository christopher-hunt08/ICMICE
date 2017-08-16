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
"""
  This scripts runs over a single or multiple files to create a summary of the
  data that is saved. Things like the number of events, number of digits in 
  each detector, etc. Just some useful scalars that could show the location of
  any issues.
"""

import MAUS

import ROOT
import sys

import event_loader


# EVENT_TYPES = {'scifi' : process_scifi, 'emr' : process_emr, \
#        'tof' : process_tof, 'kl' : process_kl, 'ckov' : process_ckov, \
#        'trigger' process_trigger, 'global' : process_global, 'mc' : process_mc}



def init_data() :
  summary = {}

  detectors = {}
  detectors['scifi'] = { 'Summary' : { 'digits': {"string" : '# Digits', 'value' : 0}, \
                                       'clusters' : { 'string' : '# Clusters', 'value' : 0 }, \
                                       'spacepoints' : { 'string' : '# Spacepoints', 'value' : 0 }, \
                                       'helix_tracks' : { 'string' : '# Helix Tracks', 'value' : 0 }, \
                                       'straight_tracks' : { 'string' : '# Straight Tracks', 'value' : 0 }, \
                                       'tracks' : { 'string' : '# Tracks', 'value' : 0 } } }

  data =  { 'summary' : summary, 'detectors' : detectors }

  return data


def process_scifi( data, event ) :

  data['Summary']['digits']['value'] += event.digits().size()
  data['Summary']['clusters']['value'] += event.clusters().size()
  data['Summary']['spacepoints']['value'] += event.spacepoints().size()
  data['Summary']['helix_tracks']['value'] += event.helicalprtracks().size()
  data['Summary']['straight_tracks']['value'] += event.straightprtracks().size()
  data['Summary']['tracks']['value'] += event.scifitracks().size()


def process_files( filenames ) :
  loader = event_loader.maus_reader( filenames )
  data = init_data()

  print data
  print

  while loader.next_event() :
# Some global counters

# Run detector-specific counters
    for event_type in EVENT_TYPES :
# Load Detector Event
      event = loader.get_event( event_type )
# Run Respective summary function
      data[event_type] = EVENT_TYPES[event_type]( data['detectors'][event_type], event )
  

  print_summary( data )
      


EVENT_TYPES = {'scifi' : process_scifi }

def print_summary( data ) :
  print 
  print "Data Summary"
  print
  print "  Overall Summary:"
  for counter in data['summary'] :
    print "     {0:>30s}  =  {1:<12}".format( \
                                           counter, data['summary'][counter] )

  print
  print
  print "  Detector Summaries"
  print
  for detector in EVENT_TYPES :
    print "  Detector : {0}".format( detector.upper() )
    for section in data['detectors'][detector] :
      section_data = data['detectors'][detector][section]
      print
      print "        {0}".format( section )
      for counter in section_data :
        print "     {0:>30s}  =  {1:<12}". format( \
               section_data[counter]['string'], section_data[counter]['value'] )
  print
  print
  print 'Complete'
  print


if __name__ == "__main__" :
  ROOT.gROOT.SetBatch( True )
  
  filenames = sys.argv[1:]

  try :
    process_files( filenames )

  except :
    raise




