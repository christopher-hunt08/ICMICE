
from .. import LastAnalysis

from _analysis_base import Analysis_Base

import analysis
from analysis.beam_sampling import multivariate_gaussian
import os
import ROOT
import json
import math
import numpy
import scipy.spatial


class VoronoiTessellation(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "Voronoi_Tessellation")

    self.__output_filename = ""

    self.__ensemble_size = 10000
    self.__the_data = []
    self.__densities = []
    self.__regions = None
    self.__point_regions = None
    self.__vertices = None

    self.__voronoi_iteration = 0
    self.__current_ensemble = 0

    self.__weights_histogram = ROOT.TH1F("vornoi_weights", ";w;#", 1000, 0.0, 100.0)
    self.__densities_histogram = ROOT.TH1F("vornoi_densities", ";#rho;#", 1000, 0.0, 1.0)


  def analyse_event(self, analysis_event, weight=1.0) :
    hit = analysis_event.selection_trackpoint()

    vector = numpy.array( [hit.get_x(), hit.get_px(), hit.get_y(), hit.get_py()] )
    self.__the_data.append( vector )
    self.__current_ensemble += 1

    if self.__current_ensemble >= self.__ensemble_size :
      self.save_voronoi_file()

      self.__the_data = []
      self.__densities = []
      self.__regions = None
      self.__point_regions = None
      self.__vertices = None

      self.__voronoi_iteration += 1
      self.__current_ensemble = 0


  def _get_plots(self, plot_dict) :
    plot_dict["voronoi_weights"] = self.__weights_histogram
    plot_dict["voronoi_densities"] = self.__densities_histogram


  def _get_data(self, data_dict) :
    data_dict['voronoi_data'] = self.__output_filename
    data_dict['required_covariance_matrix'] = [ [ self.__covariance[i][j] for i in range(4) ] for j in range(4) ]


  def configure_arguments(self, parser) :
    parser.add_argument('--voronoi_selection', type=float, nargs=5, help='Specify the Emittance, Alpha, Beta, L and Momentum of the selection')
    parser.add_argument('--ensemble_size', default=self.__ensemble_size, type=int, help="Number of events to load for each ensemble.")


  def parse_arguments(self, namespace) :
    self.__output_filename = os.path.join(namespace.output_directory, namespace.output_filename+"-VoronoiData.json")
    self.__selection = namespace.voronoi_selection
    self.__ensemble_size = namespace.ensemble_size


  def save_voronoi_file(self) :
    V = scipy.spatial.Voronoi(self.__the_data)

    emittance, alpha, beta, L, momentum = self.__selection
    mass = analysis.tools.MUON_MASS

    cov_xx = beta*emittance*mass/momentum
    cov_xp = -1.0*alpha*emittance*mass
    cov_pp = ((1.0 + alpha**2 + L**2) / beta) * emittance * momentum * mass
    cov_px = -1.0*emittance*mass*L

    means = numpy.array( [ 0.0, 0.0, 0.0, 0.0 ] )
    self.__covariance = numpy.array( [[ cov_xx, cov_xp, 0.0, cov_px ], \
                                     [ cov_xp, cov_pp, -cov_px, 0.0 ], \
                                     [ 0.0, -cov_px, cov_xx, cov_xp ], \
                                     [ cov_px, 0.0, cov_xp, cov_pp ]] )

    self.__densities = []
    self.__weights = []

    self.__regions = V.regions
    self.__point_regions = V.point_region
    self.__vertices = V.vertices

    normalisation = 1.0/multivariate_gaussian(means, means, self.__covariance)

    for num, point in enumerate(self.__the_data) :
      region = self.__regions[self.__point_regions[num]]
      vector = numpy.array(point)
      expected = normalisation*multivariate_gaussian(vector, means, self.__covariance)
      density = 0.0
      weight = 0.0

      if len(region) == 0 or -1 in region : 
        density = 0.0
        weight = 0.0

      elif expected < 0.0001 :
        density = 0.0
        weight = 0.0

      else :
        coordinates = [ self.__vertices[i] for i in region ]
        hull = scipy.spatial.ConvexHull(coordinates)
        density = math.sqrt(1.0/hull.volume)
        weight = expected / density

      self.__densities.append( density )
      self.__weights.append( weight )

    print len(self.__weights)
    max_weight = max(self.__weights)
    sum_weight = sum(self.__weights)
#    norm = len(self.__weights)/sum_weight
    norm = len(self.__weights)/sum_weight
    self.__weights = [ weight*norm for weight in self.__weights ]
    mean_weight = numpy.mean(self.__weights)

    print len(self.__weights), sum_weight, max_weight, norm, mean_weight

    for density, weight in zip(self.__densities, self.__weights) :
      self.__densities_histogram.Fill( density )
      self.__weights_histogram.Fill( weight )


    data_dict = { "regions" : self.__regions, "point_regions": self.__point_regions.tolist(), \
                  "vertices" : self.__vertices.tolist(), "densities" : self.__densities, \
                  "weights" : self.__weights, "normalisation" : mean_weight/max_weight }

    with open( self.__output_filename+".{0:05d}".format(self.__voronoi_iteration), 'w' ) as outfile :
      json.dump(data_dict, outfile)


  def conclude(self) :
    if self.__current_ensemble >= (0.5*self.__ensemble_size) : # If we haven't got half the ensemble, throw it.
      self.save_voronoi_file()

      self.__the_data = []
      self.__densities = []
      self.__regions = None
      self.__point_regions = None
      self.__vertices = None

      self.__voronoi_iteration += 1
      self.__current_ensemble = 0



