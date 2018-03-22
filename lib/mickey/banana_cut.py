

import ROOT
import math
from _cuts_base import Cut_Base

from analysis import tools

####################################################################################################
class Cut_banana_plot_mass(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "banana_plot_mass")
    self.__histogram = ROOT.TH2F( "banana_plot", ";p  [MeV/c];TOF  [ns]", 200, 100.0, 300.0, 200, 20.0, 40.0 )
    self.__momentum_histogram = ROOT.TH1F( "banana_momentum", "#Delta p  [MeV/c];# Muons", 400, -200.0, 200.0 )
    self.__mass_hypothesis = -1.0
    self.__momentum_range = None


  def _is_cut(self, analysis_event) :
    if self.__momentum_range is None :
      return False

    if analysis_event.num_upstream_tracks() == 0 :
      return True

    if analysis_event.num_tof1_spacepoints() == 0 or analysis_event.num_tof0_spacepoints() == 0 :
      return True

    p = analysis_event.upstream_reference_trackpoint().get_p()
    T = analysis_event.tof01()
    L = analysis_event.tof1_spacepoint().get_z() - analysis_event.tof0_spacepoint().get_z()

    p_min = p + self.__momentum_range[0]
    p_max = p + self.__momentum_range[1]

    t_max = L * math.sqrt(self.__mass_hypothesis**2 + p_min**2) / (300.0 * p_min)
    t_min = L * math.sqrt(self.__mass_hypothesis**2 + p_max**2) / (300.0 * p_max)

    if T > t_max or T < t_min :
      return True
    else :
      return False


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_upstream_tracks() >= 0 and analysis_event.num_tof1_spacepoints() >= 0 or analysis_event.num_tof0_spacepoints() >= 0 :

      p = analysis_event.upstream_reference_trackpoint().get_p()
      T = analysis_event.tof01()
      L = analysis_event.tof1_spacepoint().get_z() - analysis_event.tof0_spacepoint().get_z()
      b = L/(T*300.0)

      self.__histogram.Fill( p, T )

      if b > 1.0 : # Ruddy relativistic electrons!
        self.__momentum_histogram.Fill( 1.0e300 )
      else :
        self.__momentum_histogram.Fill( self.__mass_hypothesis / math.sqrt((1.0/(b*b)) - 1 ) - p )


  def _get_plots(self, plot_dict) :
    plot_dict["banana_plot"] = self.__histogram
    plot_dict["momentum_histogram"] = self.__momentum_histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--banana_plot_mass_cut", default=None, nargs=2, type=float, help="Cut on the banana plot using a momentum range." )


  def parse_arguments(self, namespace) :
      self.__mass_hypothesis = tools.MUON_MASS
      self.__momentum_range = namespace.banana_plot_mass_cut
      

