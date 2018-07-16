
from _selection_base import Selection_Base

from .. import LastAnalysis

import ROOT
import json
import analysis
import scipy.spatial
from analysis.beam_sampling import multivariate_gaussian
import numpy
import array


class VoronoiPhaseSpaceSelection(Selection_Base) :

  def __init__(self) :
    Selection_Base.__init__(self, "voronoi_phasespace_selection")

    self.__voronoi_data_file = LastAnalysis.LastData['analysis']['Voronoi_Tessellation']['voronoi_data']

    voronoi_data = None
    with open(self.__voronoi_data_file, 'r') as infile :
      voronoi_data = json.load(infile)

    self.__weights = voronoi_data['weights']
#    self.__densities = voronoi_data['densities']
#    self.__vertices = voronoi_data['vertices']
#    self.__regions = voronoi_data['regions']
#    self.__point_regions = voronoi_data['point_regions']
    self.__normalisation = voronoi_data['normalisation']
    
    self.__event_counter = 0


  def get_normalisation(self) :
    return self.__normalisation


  def weigh_event(self, event) :
    weight = self.__weights[self.__event_counter]
    self.__event_counter += 1

    if weight < 1.0e-3 : 
      return 0.0
    else :
      return weight


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass



#class VoronoiPhaseSpaceSelection(Selection_Base) :
#
#  def __init__(self, emittance, alpha, beta, momentum, mass=analysis.tools.MUON_MASS) :
#    Selection_Base.__init__(self, "voronoi_phasespace_selection")
#
#    self.__voronoi_data_file = LastAnalysis.LastData['analysis']['Voronoi_Tessellation']['voronoi_data']
#
#    voronoi_data = None
#    with open(self.__voronoi_data_file, 'r') as infile :
#      voronoi_data = json.load(infile)
#
#    self.__densities = voronoi_data['densities']
##    self.__vertices = voronoi_data['vertices']
##    self.__regions = voronoi_data['regions']
##    self.__point_regions = voronoi_data['point_regions']
#
#    cov_xx = beta*emittance*mass/momentum
#    cov_xp = -1.0*alpha*emittance*mass
#    cov_pp = ((1.0 + alpha**2) / beta) * emittance * momentum * mass
#
#    self.__means = numpy.array( [ 0.0, 0.0, 0.0, 0.0 ] )
#    self.__covariance = numpy.array( [[ cov_xx, cov_xp, 0.0, 0.0 ], \
#                                      [ cov_xp, cov_pp, 0.0, 0.0 ], \
#                                      [ 0.0, 0.0, cov_xx, cov_xp ], \
#                                      [ 0.0, 0.0, cov_xp, cov_pp ]] )
#
#
#    vector = numpy.array([0.0, 0.0, 0.0, 0.0])
#    max_val = multivariate_gaussian(vector, self.__means, self.__covariance)
#
#    vol = 0.0 
#    for den in self.__densities :
#      if den > 0.0 :
#        vol += 1.0/den
#    print "Densities Sum = ", sum(self.__densities)
#    print "Volume Sum = ", vol
##    self.__normalisation = sum(self.__densities)*(len(self.__densities)**2)*max_val/max(self.__densities)
#    self.__normalisation = len(self.__densities)**2/sum(self.__densities)
#
#    self.__weights_histogram = ROOT.TH1F("vornoi_weights", ";w;#", 1000, 0.0, 100.0)
#    self.__densities_histogram = ROOT.TH1F("vornoi_densities", ";#rho;#", 1000, 0.0, 1.0)
#
##    axis = self.__weights_histogram.GetXaxis();
##    bins = axis.GetNbins();
##    lower = axis.GetXmin();
##    upper = axis.GetXmax();
##    width = (lower - upper) / bins;
##    new_bins = array.array('d')
##    for i in range( bins, -1, -1 ) :
##      new_bins.append( ROOT.TMath.Power(10, lower + i * width) )
##    axis.Set(bins, new_bins)
#
#    self.__event_counter = 0
#
#
#  def weigh_event(self, event) :
#    weight = 0.0
#
#    hit = event.selection_trackpoint()
#    vector = numpy.array( [hit.get_x(), hit.get_px(), hit.get_y(), hit.get_py()] )
#
#    density = self.__densities[self.__event_counter]
#
#    expected = self.__normalisation*multivariate_gaussian(vector, self.__means, self.__covariance)
#
#    if density < 1.0e-12 :
#      weight = 0.0
#    else :
#      weight = expected / density
#
#    self.__event_counter += 1
#    self.__weights_histogram.Fill(weight)
#    self.__densities_histogram.Fill(density)
#    return weight
#
#
#  def _get_plots(self, plot_dict) :
#    plot_dict["voronoi_weights"] = self.__weights_histogram
#    plot_dict["voronoi_densities"] = self.__densities_histogram
#
#
#  def _get_data(self, data_dict) :
#    pass

