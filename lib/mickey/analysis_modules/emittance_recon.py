
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array


class EmittanceAnalysis(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "upstream_emittance_reconstruction")

    self.__inspectors = []
    self.__momentum_windows = [ (0, 0.0, 300.0) ]
    self.__emittance_graph = None


  def analyse_event(self, analysis_event, weight=1.0) :
    hit = analysis_event.upstream_reference_trackpoint()
    hit.set_weight(weight)

    p = hit.get_p()

    for num, low, high in self.__momentum_windows :
      if p >= low and p < high :
        self.__inspectors[num].add_hit(hit)


  def _get_plots(self, plot_dict) :
    for num, low, high in self.__momentum_windows :
      name = "bin_{0:d}".format(num)
      plot_dict[name] = self.__inspectors[num].get_plot_dictionary()

    plot_dict["emittance"] = self.__emittance_graph


  def _get_data(self, data_dict) :
    for num, low, high in self.__momentum_windows :
      name = "bin_{0:d}".format(num)
      data_dict[name] = self.__inspectors[num].get_data_dictionary()

      data_dict[name]['low_edge'] = low
      data_dict[name]['high_edge'] = high


  def configure_arguments(self, parser) :
    parser.add_argument( "--emittance_momentum_bins", nargs=3, type=float, default=[180.0, 260.0, 8], help="Specify start, end and number of bins" )


  def parse_arguments(self, namespace) :
    start = namespace.emittance_momentum_bins[0]
    end = namespace.emittance_momentum_bins[1]
    number = namespace.emittance_momentum_bins[2]

    width = (end-start) / number

    self.__momentum_windows = []

    for i in range( number ) :
      self.__momentum_windows.append(( i, (start+i*width), (start+(i+1)*width) ))
      self.__inspectors.append(inspectors.PhaseSpace2DInspector(i, 0))


  def conclude(self) :

    p = array.array("d")
    p_error = array.array("d")
    emittance = array.array("d")
    emittance_error = array.array("d")

    for num, low, high in self.__momentum_windows :
      p.append( 0.5*(high+low) )
      p_error.append( 0.5*(high-low) )
      em = self.__inspectors[num].covariance.get_emittance()
      emittance.append(em)
      emittance_error.append( em / math.sqrt(2*(self.__inspectors[num].covariance.length()-1)) )

    self.__emittance_graph = ROOT.TGraphErrors(len(p), p, emittance, p_error, emittance_error)
      

    

