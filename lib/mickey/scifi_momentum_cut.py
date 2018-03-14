
import ROOT
from _cuts_base import Cut_Base

# chisq ndf
####################################################################################################
class Cut_scifi_upstream_momentum(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "upstream_scifi_momentum")
    self.__histogram = ROOT.TH1F( "upstream_scifi_momentum", ";#p  [MeV/c;# Muons", 500, 0.0, 500.0 )
    self.__p_cut = [0.0, 1000.0]


  def _is_cut(self, analysis_event) :
    if analysis_event.num_upstream_tracks() == 0 :
      return True

    if analysis_event.upstream_reference_trackpoint().get_p() < self.__p_cut[0] or analysis_event.upstream_reference_trackpoint().get_p() > self.__p_cut[1] :
      return True

    else :
      return False


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_upstream_tracks() > 0 :
      self.__histogram.Fill( analysis_event.upstream_reference_trackpoint().get_p() )


  def _get_plots(self, plot_dict) :
    plot_dict["momentum"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--upstream_momentum_cut", default=[0.0, 1000.0], type=float, nargs=2, help="Cut on upstream reference plane momentum." )


  def parse_arguments(self, namespace) :
    self.__p_cut = namespace.upstream_momentum_cut

