#!/usr/bin/env python

# Import MAUS Framework (Required!)
import MAUS

# Generic Python imports
import sys
import os
import json
import math
import argparse

# Third Party library import statements
import analysis
import ROOT
import event_loader

ALIGNMENT_TOLERANCE=1.0E-3


def find_virtual_planes(file_reader, max_num_events) :
  virtual_plane_dict = {}

  try :
    while file_reader.next_event() and file_reader.get_total_num_events() <= max_num_events :
      mc_event = file_reader.get_event( 'mc' )

      for vhit_num in xrange(mc_event.GetVirtualHitsSize()) :
        vhit = mc_event.GetAVirtualHit(vhit_num)
        plane_id = vhit.GetStationId()
        virtual_plane_dict[ plane_id ] = vhit.GetPosition().z()
  except KeyboardInterrupt :
    pass

  return virtual_plane_dict



def print_virtual_plane_dict(virtual_plane_dict) :
  print
  for station_id, position in sorted(virtual_plane_dict.iteritems()) :
    print position, " : ", station_id



def create_virtual_plane_dict(file_reader, max_num_events) :
  """
    Matches up scifitrackpoints to virtual planes to make a lookup dictionary 
  """
  virtual_plane_dict = {}
  virtual_plane_list = {}
  for num in range( -15, 0, 1 ) :
    virtual_plane_dict[ num ] = ( -1, (ALIGNMENT_TOLERANCE * 100.0) )
  for num in range( 1, 16, 1 ) :
    virtual_plane_dict[ num ] = ( -1, (ALIGNMENT_TOLERANCE * 100.0) )

  while file_reader.next_event() and file_reader.get_total_num_events() <= max_num_events :
    scifi_event = file_reader.get_event( 'scifi' )
    mc_event = file_reader.get_event( 'mc' )

    for vhit_num in xrange(mc_event.GetVirtualHitsSize()) :
      vhit = mc_event.GetAVirtualHit(vhit_num)
      plane_id = vhit.GetStationId()
      virtual_plane_list[ plane_id ] = vhit.GetPosition().z()


    tracks = scifi_event.scifitracks()
    for track in tracks :
      trackpoints = track.scifitrackpoints()
      for trkpt in trackpoints :
        z_pos = trkpt.pos().z()
        plane_id = analysis.tools.calculate_plane_id(\
                               trkpt.tracker(), trkpt.station(), trkpt.plane())

        for vhit_num in xrange(mc_event.GetVirtualHitsSize()) :
          vhit = mc_event.GetAVirtualHit(vhit_num)
          diff = math.fabs(vhit.GetPosition().z() - z_pos)

          if diff < virtual_plane_dict[ plane_id ][1] :
            virtual_plane_dict[ plane_id ] = ( vhit.GetStationId(), diff )

    done = True
    for plane in virtual_plane_dict :
      if virtual_plane_dict[plane][1] > ALIGNMENT_TOLERANCE :
#        print plane, virtual_plane_dict[plane]
        done = False
    if done :
      break

  file_reader.reset()
  return virtual_plane_list, virtual_plane_dict



def inverse_virtual_plane_dict(virtual_plane_dict) :
  """
    Create the inverse lookup.
  """
  inverse_dict = {}
  for num in range( -15, 0, 1 ) :
    inverse_dict[virtual_plane_dict[num][0]] = num
  for num in range( 1, 16, 1 ) :
    inverse_dict[virtual_plane_dict[num][0]] = num

  return inverse_dict


def trim_virtual_plane_dict(virtual_plane_dict) :
  """
    Create the inverse lookup.
  """
  new_dict = {}
  for station, plane in virtual_plane_dict.iteritems() :
    new_dict[station] = plane[0]

  return new_dict



if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  parser = argparse.ArgumentParser( description='Prints out a list of virtual planes and creates the scifi-virtual plane dictionary' )

  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS output root files containing reconstructed straight tracks')

  parser.add_argument( '-O', '--output_filename', default='virtual_plane_lookup', help='Set the output filename')

  parser.add_argument( '-D', '--output_directory', default='./', help='Set the output directory')

  parser.add_argument( '-N', '--max_num_events', type=int, default=100, help='Maximum number of events to analyse.')


  try :
    namespace = parser.parse_args()
  except BaseException as ex:
    raise
  else :
    print
    file_reader = event_loader.maus_reader(namespace.maus_root_files)
##### 3. Initialise Plane Dictionary ##########################################
    print "Locating all virtual planes."
    print "Please Wait."
    print
#    virtual_plane_dict = find_virtual_planes(file_reader, namespace.max_num_events)
    virtual_plane_dict, scifi_virtual_dict = create_virtual_plane_dict(file_reader, namespace.max_num_events)
#    scifi_virtual_dict = inverse_virtual_plane_dict(scifi_virtual_dict)
    scifi_virtual_dict = trim_virtual_plane_dict(scifi_virtual_dict)

    outfile = os.path.join(namespace.output_directory, namespace.output_filename) + '.json'
    with open(outfile, 'w') as output :
      json.dump(scifi_virtual_dict, output)

  print

