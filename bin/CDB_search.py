#! /bin/env python

import sys
import os
import argparse
import pprint
from cdb import Beamline
from cdb import CoolingChannel
CDB_SERVER = 'http://cdb.mice.rl.ac.uk'

START_RUN = 7400
SEARCH = False
ALL_RUNS = []


def compare(this_one, that_one) :
  if this_one == that_one :
    return True
  else :
    return False

def find_channel(channel) :
  channel_interface = CoolingChannel(CDB_SERVER)
  found = []

  for run in ALL_RUNS :
    run_channel = channel_interface.get_coolingchannel_for_run(run)
    if compare(run_channel, channel) :
      found.append(run)

  return found


def find_beamline(beamline) :
  beamline_interface = Beamline(CDB_SERVER)
  found = []

  for run in ALL_RUNS :
    run_beamline = beamline_interface.get_beamline_for_run(run)
    if compare(run_beamline, beamline) :
      found.append(run)

  return found




if __name__ == "__main__" :
  parser = argparse.ArgumentParser()

  parser.add_argument('--channel', help="Channel Name")
  parser.add_argument('--beamline', help="Beamline Config")
  parser.add_argument('--diffuser', help="Diffuser Setting")
  parser.add_argument('--comment', help="Begin-of-run OR End-of-run comment")
  parser.add_argument('--date', help="Runs that happened on that date")
  parser.add_argument('--run_number', help="Detail a specific run number")
  parser.add_argument('-S', '--search', action="store_true", help="Assume arguments are search terms")
  parser.add_argument('--starting_run', type=int, default=7400, help="Starting run number to search from")
  parser.add_argument('-Q', '--quiet', action="store_true", help="Minimal outout for piping")

  namespace = parser.parse_args()

  SEARCH = namespace.search
  START_RUN = namespace.starting_run

  CHANNEL = CoolingChannel(CDB_SERVER)
  BEAMLINE = Beamline(CDB_SERVER)

  all_runs = sorted( Beamline(CDB_SERVER).get_run_numbers() )
  for run in all_runs :
    if int(run) > START_RUN :
      ALL_RUNS.append(run)

  all_found = []

  if not namespace.quiet :
    print
    sys.stdout.write("Searching: [                                                                                                    ]\rSearching: [")
    sys.stdout.flush()
    total = len(ALL_RUNS)
    percent = int(-1)

  for number, run in enumerate(ALL_RUNS) :
    int_run = int(run)
    found = True
    run_channel = CHANNEL.get_coolingchannel_for_run(run)
    run_beamline = BEAMLINE.get_beamline_for_run(run)

    if not namespace.quiet :
      new_percent = int((float(number) / float(total)) * 100.0)
      if new_percent != percent :
        sys.stdout.write("-")
        sys.stdout.flush()
        percent = new_percent

    if namespace.channel is not None :
      if compare(run_channel['tag'], namespace.channel) :
        found = True
      else :
        found = False

    if found and ( namespace.beamline is not None ) :
      if compare(run_beamline[int_run]['optics'], namespace.beamline) :
        found = True
      else :
        found = False


    if found :
      all_found.append( {"run_number": run, "channel": run_channel['tag'], "beamline": run_beamline[int_run]['optics'], 'diffuser':run_beamline[int_run]['diffuser_thickness']} )

  if namespace.quiet :
    for run in all_found :
      print run["run_number"],

  else :
    print
    print
    print "Found {0} Runs".format(len(all_found))
    print

    for run in all_found :
      print run["run_number"], run["channel"], run["beamline"], run["diffuser"]

