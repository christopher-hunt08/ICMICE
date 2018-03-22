
from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array


class EmittanceSystematicErrors(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "mc_emittance_systematic_error_calculation", require_mc=True)

    self.__inspectors = []
    self.__momentum_windows = [ (0, 0.0, 300.0) ]
    self.__bias_graph = None
    self.__error_graph = None
    self.__bias_values = []
    self.__error_values = []


  def analyse_event(self, analysis_event, weight=1.0) :
    hit = analysis_event.upstream_reference_trackpoint()
    vhit = analysis_event.mc_upstream_reference_trackpoint()
    if vhit is None :
      print
      print "Expected Virtual Hit - Found Nothing!"
      print
      return

    hit.set_weight(weight)

    p = hit.get_p()

    for num, low, high in self.__momentum_windows :
      if p >= low and p < high :
        self.__inspectors[num].add_hit(hit, vhit)


  def _get_plots(self, plot_dict) :
    for num, low, high in self.__momentum_windows :
      name = "bin_{0:d}".format(num)
      plot_dict[name] = self.__inspectors[num].get_plot_dictionary()

    plot_dict["emittance_bias"] = self.__bias_graph
    plot_dict["emittance_systematic"] = self.__error_graph


  def _get_data(self, data_dict) :
    for num, low, high in self.__momentum_windows :
      name = "bin_{0:d}".format(num)
      data_dict[name] = self.__inspectors[num].get_data_dictionary()

      data_dict[name]['low_edge'] = low
      data_dict[name]['high_edge'] = high
      data_dict[name]['bias'] = self.__bias_values[num][0]
      data_dict[name]['bias_error'] = self.__bias_values[num][1]
      data_dict[name]['error'] = self.__error_values[num][0]
      data_dict[name]['error_low'] = self.__error_values[num][1]
      data_dict[name]['error_high'] = self.__error_values[num][2]


  def configure_arguments(self, parser) :
    parser.add_argument( "--emittance_momentum_bins", nargs=3, type=float, default=[175.0, 255.0, 8], help="Specify start, end and number of bins" )


  def parse_arguments(self, namespace) :
    start = namespace.emittance_momentum_bins[0]
    end = namespace.emittance_momentum_bins[1]
    number = namespace.emittance_momentum_bins[2]

    width = (end-start) / number

    self.__momentum_windows = []

    for i in range( number ) :
      self.__momentum_windows.append(( i, (start+i*width), (start+(i+1)*width) ))
      self.__inspectors.append(inspectors.EmittanceSystematicErrorInspector(i, 2000))


  def conclude(self) :

    p = array.array("d")
    p_error = array.array("d")
    bias = array.array("d")
    bias_err_low = array.array("d")
    bias_err_high = array.array("d")
    sys = array.array("d")
    sys_err_low = array.array("d")
    sys_err_high = array.array("d")
    zeros = array.array("d")

    self.__bias_values = []
    self.__error_values = []

    for num, low, high in self.__momentum_windows :
      p.append( 0.5*(high+low) )
      p_error.append( 0.5*(high-low) )

      b, bias_lower, bias_upper, s, sys_lower, sys_upper = self.__inspectors[num].get_systematic_error()

      bias.append(b)
      bias_err_low.append(b-bias_lower)
      bias_err_high.append(bias_upper-b)
      sys.append(s)
      sys_err_low.append(s-sys_lower)
      sys_err_high.append(sys_upper-s)
      zeros.append(0.0)

      self.__bias_values.append((b, 0.5*(bias_upper-bias_lower)))
      self.__error_values.append((s, sys_lower, sys_upper))

    self.__bias_graph = ROOT.TGraphAsymmErrors(len(p), p, bias, p_error, p_error, bias_err_low, bias_err_high)
    self.__error_graph = ROOT.TGraphAsymmErrors(len(p), p, sys, p_error, p_error, sys_err_low, sys_err_high)



