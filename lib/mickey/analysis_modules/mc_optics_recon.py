
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array


class MCTruthOpticsRecon(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "mc_truth_optics_reconstruction")

    self.__inspectors = []


  def analyse_event(self, analysis_event, weight=1.0) :
    pass


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
#    parser.add_argument( 
    pass


  def parse_arguments(self, namespace) :
    pass


  def conclude(self) :
    pass

