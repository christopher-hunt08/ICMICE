
import ROOT
import math

from _cuts_base import Cut_Base


####################################################################################################
class Cut_diffuser_aperture(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "diffuser_aperture")

    self.__histogram = ROOT.TH1F( "cut_diffuser_radius", ";r  [mm];# Events", 400, 0.0, 200.0 )
    self.__z_position = 13740.0
    self.__radius_cut = 90.0
    self.__tolerance = 1.0


  def _is_cut(self, analysis_event) :
    if analysis_event.num_global_tracks() == 0 :
      raise ValueError("No Global Track Was found")

    global_track = analysis_event.global_track()
    for i in range(len(global_track)) :
      tp = global_track[i]

      if math.fabs( tp.get_z() - self.__z_position ) < self.__tolerance :
        if tp.get_radius() < self.__radius_cut :
          return False
        else :
          return True

    raise ValueError("No virtual plane was located near the diffuser")


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_global_tracks() == 0 :
      return

    global_track = analysis_event.global_track()
    for i in range(len(global_track)) :
      tp = global_track[i]
      if math.fabs( tp.get_z() - self.__z_position ) < self.__tolerance :
        self.__histogram.Fill( tp.get_radius() )
        return

    self.__histogram.Fill(-1.0)


  def _get_plots(self, plot_dict) :
    plot_dict["tof01_number_spacepoints"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument('--diffuser_position', type=float, default=self.__z_position, help='z-Position of diffuser')
    parser.add_argument('--diffuser_radius_cut', type=float, default=self.__radius_cut, help='Radius to cut tracks')


  def parse_arguments(self, namespace) :
    self.__z_position = namespace.diffuser_position
    self.__radius_cut = namespace.diffuser_radius_cut
    


