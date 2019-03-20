
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


####################################################################################################
## Required to adjust the default KDE width
####################################################################################################
def silverman_3(kernel) :
  silverman = (kernel.n * (kernel.d + 2) / 4.)**(-1. / (kernel.d + 4))
  return silverman / 3.0

####################################################################################################



class KDECalculation(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "KDE_Calculation")

    self.__output_filename = ""
    self.__norm_samples = 10000

    self.__means = None
    self.__covariance = None

    self.__ensemble_size = 10000
    self.__the_data = []
    self.__densities = []

    self.__kde_iteration = 0
    self.__current_ensemble = 0
    self.__max_event_weight = 5.0

    self.__weights_histogram = ROOT.TH1F("kde_weights", ";w;#", 1000, 0.0, 100.0)
    self.__densities_histogram = ROOT.TH1F("kde_densities", ";#rho;#", 1000, 0.0, 1.0)


  def analyse_event(self, analysis_event, weight=1.0) :
    hit = analysis_event.selection_trackpoint()

    vector = numpy.array( [hit.get_x(), hit.get_px(), hit.get_y(), hit.get_py()] )
    self.__the_data.append( vector )
    self.__current_ensemble += 1

    if self.__current_ensemble >= self.__ensemble_size :
      self.save_kde_file()

      self.__the_data = []
      self.__densities = []

      self.__kde_iteration += 1
      self.__current_ensemble = 0


  def _get_plots(self, plot_dict) :
    plot_dict["kde_weights"] = self.__weights_histogram
    plot_dict["kde_densities"] = self.__densities_histogram


  def _get_data(self, data_dict) :
    data_dict['kde_data'] = self.__output_filename
    data_dict['required_covariance_matrix'] = [ [ self.__covariance[i][j] for i in range(4) ] for j in range(4) ]


  def configure_arguments(self, parser) :
    parser.add_argument('--kde_selection', type=float, nargs=5, help='Specify the Emittance, Alpha, Beta, L and Momentum of the selection')
    parser.add_argument('--ensemble_size', default=self.__ensemble_size, type=int, help="Number of events to load for each ensemble.")
#    parser.add_argument('--maximum_event_weight', type=float, default=self.__max_event_weight, help="Ensure no outlandish event weights are applied. Weights greater than the value specified will be cut.")
    parser.add_argument('--number_normalisation_samples', type=int, default=self.__norm_samples, help="Number of samples to take in order to estimate the normalalisation coefficient.")


  def parse_arguments(self, namespace) :
    self.__output_filename = os.path.join(namespace.output_directory, namespace.output_filename+"-KDEData.json")
    self.__selection = namespace.kde_selection
    self.__ensemble_size = namespace.ensemble_size
#    self.__max_event_weight = namespace.maximum_event_weight
    self.__norm_samples = namespace.number_normalisation_samples

    emittance, alpha, beta, L, momentum = self.__selection
    mass = analysis.tools.MUON_MASS

    cov_xx = beta*emittance*mass/momentum
    cov_xp = -1.0*alpha*emittance*mass
    cov_pp = ((1.0 + alpha**2 + L**2) / beta) * emittance * momentum * mass
    cov_px = -1.0*emittance*mass*L

    self.__means = numpy.array( [ 0.0, 0.0, 0.0, 0.0 ] )
    self.__covariance = numpy.array( [[ cov_xx, cov_xp, 0.0, cov_px ], \
                                     [ cov_xp, cov_pp, -cov_px, 0.0 ], \
                                     [ 0.0, -cov_px, cov_xx, cov_xp ], \
                                     [ cov_px, 0.0, cov_xp, cov_pp ]] )


  def save_kde_file(self) :
    kernel = scipy.stats.gaussian_kde(numpy.transpose(self.__the_data), bw_method='scott')

    self.__densities = []
    self.__weights = []

#    normalisation = 1.0/multivariate_gaussian(self.__means, self.__means, self.__covariance)
#    normalisation = 1.0

    print "Estimating Normalisation"
    norm = 1.0e20
    determinant = numpy.linalg.det(self.__covariance)
    covariance_inv = numpy.linalg.inv(self.__covariance)
    count = 0
    while count < self.__norm_samples :
      x = scipy.random.multivariate_normal(self.__means, self.__covariance)
      vector = numpy.array(numpy.array(x) - numpy.array(self.__means))
      amplitude = vector.transpose().dot(covariance_inv.dot(vector))
      if amplitude > 1.0 : continue
      count += 1

      test_norm = norm
      parent = kernel.pdf( x )
      daughter = multivariate_gaussian(x, self.__means, self.__covariance)

#      number_density = kernel.n * parent
#      uncertainty = 0.25*math.sqrt(daughter*(1.0-daughter) / number_density)
      uncertainty = 0.0
      if daughter > 1.0e-9 :
        test_norm = (parent+uncertainty) / daughter

      if test_norm < norm :
        norm = test_norm

    print "Normalisation Constant =", norm
    normalisation = norm
    print "Weighing Events"

    for num, point in enumerate(self.__the_data) :
      vector = numpy.array(point)
      expected = float(normalisation*multivariate_gaussian(vector, self.__means, self.__covariance))
      density = 1.0/kernel.pdf(point)[0]
#      weight = expected / density
      weight = float(expected * density)

      self.__densities.append( density )
      self.__weights.append( weight )


    print "Weighed", len(self.__weights), "Events"
#    norm = len(self.__weights)/sum(self.__weights)

#    self.__weights = [ weight*norm for weight in self.__weights ]

#    new_weights = []
#    for weight in self.__weights :
#      weight = weight*norm
#      if weight > self.__max_event_weight :
#        weight = self.__max_event_weight
#      new_weights.append(weight)
#
#    self.__weights = new_weights

    max_weight = max(self.__weights)
    max_density = max(self.__densities)
    sum_weight = sum(self.__weights)
    mean_weight = numpy.mean(self.__weights)

    print len(self.__weights), sum_weight, max_weight, norm, mean_weight

    for density, weight in zip(self.__densities, self.__weights) :
      self.__densities_histogram.Fill( density/max_density )
      self.__weights_histogram.Fill( weight )


    data_dict = { "densities" : self.__densities, "weights" : self.__weights, "normalisation" : mean_weight/max_weight }

    with open( self.__output_filename+".{0:05d}".format(self.__kde_iteration), 'w' ) as outfile :
      json.dump(data_dict, outfile)


  def conclude(self) :
    if self.__current_ensemble >= (0.5*self.__ensemble_size) : # If we haven't got half the ensemble, throw it.
      self.save_kde_file()

      self.__the_data = []
      self.__densities = []

      self.__kde_iteration += 1
      self.__current_ensemble = 0

