#!/bin/env python

import MAUS
import ROOT

import numpy
import argparse
import sys
import os
import math
from math import sqrt
import array
import event_loader
import glob

from analysis import tools



SPACEPOINT_POSITIONS = array.array('d')
SPACEPOINT_POSITIONS.append( 0.0 )
SPACEPOINT_POSITIONS.append( 200.0 )
SPACEPOINT_POSITIONS.append( 450.0 )
SPACEPOINT_POSITIONS.append( 750.0 )
SPACEPOINT_POSITIONS.append( 1100.0 )



def get_tracks(scifi_event) :
  spacepoints = scifi_event.spacepoints()

  upstream_track = {}
  downstream_track = {}

  bad_up = False
  bad_down = False

  for spacepoint in spacepoints :
    if ( spacepoint.get_tracker() == 0 ) and ( bad_up == False ) :
      station = spacepoint.get_station() - 1
      if station in upstream_track :
        bad_up = True
        continue
      upstream_track[station] = spacepoint

    if ( spacepoint.get_tracker() == 1 ) and ( bad_down == False ) :
      station = spacepoint.get_station() - 1
      if station in downstream_track :
        bad_down = True
        continue
      downstream_track[station] = spacepoint

  if bad_up or len( upstream_track ) != 5 :
    upstream_track = None
  if bad_down or len( downstream_track ) != 5 :
    downstream_track = None

  return upstream_track, downstream_track



Z_MATRIX = numpy.array( [ [ 1.0, 0.0 ], [ 1.0, 200.0 ], [ 1.0, 450.0 ], [ 1.0, 750.0 ], [ 1.0, 1100.0 ] ] )
z_T = Z_MATRIX.transpose()
FIT_MATRIX = numpy.linalg.inv( z_T.dot(Z_MATRIX) ).dot( z_T )
def fit_track( track ) :
  x_matrix = []
  y_matrix = []

  for station in range(5) :
    spacepoint = track[station]
    pos = spacepoint.get_position()

    x_matrix.append( pos.X() )
    y_matrix.append( pos.Y() )

  x_matrix = numpy.array( x_matrix )
  y_matrix = numpy.array( y_matrix )

  beta_x = FIT_MATRIX.dot( x_matrix )
  beta_y = FIT_MATRIX.dot( y_matrix )

  solve_x = Z_MATRIX.dot( beta_x ) - x_matrix
  solve_y = Z_MATRIX.dot( beta_y ) - y_matrix
  chisq = solve_x.dot(solve_x) + solve_y.dot(solve_y)

  return beta_x, beta_y, ( chisq / 4.0 ) / 4.0




class Cuts(object) :
  def __init__(self, name="") :
    self.__tof_0_1_low = 0.0
    self.__tof_0_1_high = 100.0

    self.__chi_squared_max = 100.0

    self.__min_num_trackpoints = 15

    self.__tof_0_1 = ROOT.TH1F( name+"tof_0_1_time", "", 500, 0.0, 100.0 )
    self.__tof_0_1_cut = ROOT.TH1F( name+"tof_0_1_time_cut", "", 500, 0.0, 100.0 )

    self.__chi_squared = ROOT.TH1F( name+"chi_squared", "", 200, 0.0, 100.0 )
    self.__chi_squared_cut = ROOT.TH1F( name+"chi_squared_cut", "", 200, 0.0, 100.0 )


  def get_plots(self) :
    plot_dict = {}
    plot_dict['tof_0_1_time'] = self.__tof_0_1
    plot_dict['tof_0_1_time_cut'] = self.__tof_0_1_cut
    plot_dict['chi_squared'] = self.__chi_squared
    plot_dict['chi_squared_cut'] = self.__chi_squared_cut
    return plot_dict


  def set_tof_cut(self, tof_low, tof_high) :
    self.__tof_0_1_low = float(tof_low)
    self.__tof_0_1_high = float(tof_high)


  def set_chi_squared_cut(self, chi) :
    self.__chi_squared_max = float(chi)


  def is_cut_tof(self, tof_event) :
    event_spacepoints = tof_event.GetTOFEventSpacePoint()

    tof0_sp_size = event_spacepoints.GetTOF0SpacePointArraySize()
    tof1_sp_size = event_spacepoints.GetTOF1SpacePointArraySize()

    if tof0_sp_size == 0 or tof1_sp_size == 0 :
      return False # Must be MC

    if tof0_sp_size != 1 or tof1_sp_size != 1 :
      return True

    tof0_sp = event_spacepoints.GetTOF0SpacePointArrayElement(0)
    tof1_sp = event_spacepoints.GetTOF1SpacePointArrayElement(0)

    diff_0_1 = tof1_sp.GetTime() - tof0_sp.GetTime()

    self.__tof_0_1.Fill( diff_0_1 )

    if ( diff_0_1 < self.__tof_0_1_low ) or ( diff_0_1 > self.__tof_0_1_high ) :
      return True

    self.__tof_0_1_cut.Fill( tof1_sp.GetTime() - tof0_sp.GetTime() )
    return False


  def is_cut_scifi(self, track) :
    _, _, chisq = track

    self.__chi_squared.Fill( chisq )

    if ( chisq > self.__chi_squared_max ) :
      return True

    self.__chi_squared_cut.Fill( chisq )

    return False





class Analysis(object) :
  def __init__(self, name="") :
    self.__name = name
    self.__gradient_reconstruction = ROOT.TH2F( name+"coaxial_gradients", "", 200, -0.05, 0.05, 200, -0.05, 0.05 )
    self.__station_distributions = [ ROOT.TH2F( name+"_station_{0}".format( i ), "", 500, -100., 100.0, 500, -100.0, 100.0 ) for i in range(5) ]

    self.__x_graph = None
    self.__y_graph = None
    self.__r_graph = None


  def analyse_track(self, track) :
    beta_x, beta_y, _ = track

    self.__gradient_reconstruction.Fill( beta_x[1], beta_y[1] )


  def analyse_trackpoints(self, track) :
    for i in range(5) :
      sp = track[i].get_position()

      self.__station_distributions[i].Fill(sp.X(), sp.Y())


  def conclude(self) :
    a_x = array.array('d')
    a_ex = array.array('d')
    a_y = array.array('d')
    a_ey = array.array('d')
    zeros = array.array('d')
    for i in range(5) :
      zeros.append( 0.0 )
      a_x.append( self.__station_distributions[i].GetMean(1) )
      a_ex.append( self.__station_distributions[i].GetMeanError(1) )
      a_y.append( self.__station_distributions[i].GetMean(2) )
      a_ey.append( self.__station_distributions[i].GetMeanError(2) )

    self.__x_graph = ROOT.TGraphErrors( len(a_x), SPACEPOINT_POSITIONS, a_x, zeros, a_ex )
    self.__y_graph = ROOT.TGraphErrors( len(a_y), SPACEPOINT_POSITIONS, a_y, zeros, a_ey )


    y_projection = self.__gradient_reconstruction.ProjectionY( self.__name+'_pro_Y', 1, self.__gradient_reconstruction.GetYaxis().GetNbins() )
    y_pro_mean, y_pro_mean_err, y_pro_std, y_pro_std_err = tools.fit_gaussian( y_projection )
    y_pro_mean, y_pro_mean_err, y_pro_std, y_pro_std_err = tools.fit_gaussian( y_projection, y_pro_mean-y_pro_std, y_pro_mean+y_pro_std )

    x_projection = self.__gradient_reconstruction.ProjectionX( self.__name+'_pro_X', 1, self.__gradient_reconstruction.GetXaxis().GetNbins() )
    x_pro_mean, x_pro_mean_err, x_pro_std, x_pro_std_err = tools.fit_gaussian( x_projection )
    x_pro_mean, x_pro_mean_err, x_pro_std, x_pro_std_err = tools.fit_gaussian( x_projection, x_pro_mean-x_pro_std, x_pro_mean+x_pro_std )

    results_string = "\nFit Results\nName : {0}\n\n Mean x = {1: 6.3f} +/- {2: 6.3f} mm\n Std  x = {3: 6.3f} +/- {4: 6.3f} mrad\n\n Mean y = {5: 6.3f} +/- {6: 6.3f} mm\n Std  y = {7: 6.3f} +/- {8: 6.3f} mrad\n".format( \
        self.__name, x_pro_mean*1e3, x_pro_mean_err*1e3, x_pro_std*1e3, x_pro_std_err*1e3, y_pro_mean*1e3, y_pro_mean_err*1e3, y_pro_std*1e3, y_pro_std_err*1e3 )

    
    return results_string


  def get_plots(self) :
    plot_dict = {}

    plot_dict["gradients"] = self.__gradient_reconstruction
    for i in range(5) :
      plot_dict["station_{0}".format(i)] = self.__station_distributions[i]

    plot_dict["x_trend"] = self.__x_graph
    plot_dict["y_trend"] = self.__y_graph

    return plot_dict






if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  parser = argparse.ArgumentParser( description='' )

  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS output root files' )

  parser.add_argument( '-O', '--output_filename', default='pulls_alignment_analysis', help='Set the output filename')

  parser.add_argument( '-D', '--output_directory', default='./', help='Set the output directory')

  parser.add_argument( '-T', '--TOF01_window', default=[0.0, 100.0], nargs=2, help='Specify the TOF0-TOF1 window.' )

  parser.add_argument( '-X', '--chi_squared_cut', default=5.0, help='Specify the max Chi-Squared per degree of Freedom.' )


  try :
    namespace = parser.parse_args()

    upstream_cut_control = Cuts("UP_")
    upstream_cut_control.set_tof_cut(namespace.TOF01_window[0], namespace.TOF01_window[1])
    upstream_cut_control.set_chi_squared_cut(namespace.chi_squared_cut)

    downstream_cut_control = Cuts("DOWN_")
    downstream_cut_control.set_tof_cut(namespace.TOF01_window[0], namespace.TOF01_window[1])
    downstream_cut_control.set_chi_squared_cut(namespace.chi_squared_cut)

    upstream_analyser = Analysis("UP_")
    downstream_analyser = Analysis("DOWN_")

  except BaseException as ex:
    raise
  else :

    files = []
    for filename in namespace.maus_root_files :
      files += glob.glob(filename)

    file_reader = event_loader.maus_reader(files)
    file_reader.set_print_progress("spill")
#    file_reader.set_max_num_events(1000)

    print "\n- Loading Spills...\n"
    try :
      while file_reader.next_selected_event() :

        try :
          scifi_event = file_reader.get_event( 'scifi' )
          tof_event = file_reader.get_event( 'tof' )

          if upstream_cut_control.is_cut_tof(tof_event) :
            continue

          upstream_track, downstream_track = get_tracks( scifi_event )

          if upstream_track is not None :
            the_up_track = fit_track( upstream_track )
            if not upstream_cut_control.is_cut_scifi( the_up_track ) :
              upstream_analyser.analyse_track(the_up_track)
              upstream_analyser.analyse_trackpoints(upstream_track)

          if downstream_track is not None :
            the_down_track = fit_track( downstream_track )
            if not downstream_cut_control.is_cut_scifi( the_down_track ) :
              downstream_analyser.analyse_track(the_down_track)
              downstream_analyser.analyse_trackpoints(downstream_track)


        except ValueError as ex :
          print "An Error Occured: " + str(ex)
          print "Skipping Event: " +\
                str(file_reader.get_current_event_number()) + " In Spill: " + \
                str(file_reader.get_current_spill_number()) + " In File: " + \
                str(file_reader.get_current_filenumber()) + "\n"
          continue

    except KeyboardInterrupt :
      print
      print " ###  Keyboard Interrupt  ###"
      print
    print "- {0:0.0f} Spills Loaded                                 ".format( \
                                            file_reader.get_total_num_spills())

    sys.stdout.write( "\n- Saving Plots and Data : Running\r" )
    sys.stdout.flush()

    details = upstream_analyser.conclude()
    details += downstream_analyser.conclude()

    filename = os.path.join(namespace.output_directory, namespace.output_filename)
    plots = { "upstream_cuts" : upstream_cut_control.get_plots(), "upstream_analysis" : upstream_analyser.get_plots(), \
              "downstream_cuts" : downstream_cut_control.get_plots(), "downstream_analysis" : downstream_analyser.get_plots() }
    tools.save_plots(plots, filename+'.root')

    text_file = os.path.join(namespace.output_directory, namespace.output_filename)
    with open(text_file+'.txt', 'w') as outfile :
      outfile.write(details)

    sys.stdout.write(   "- Saving Plots and Data : Done   \n" )

    print details

  print 
  print "Complete."
  print

