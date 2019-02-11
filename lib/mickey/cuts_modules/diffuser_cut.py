
import ROOT
import math

from _cuts_base import Cut_Base


####################################################################################################
class DiffuserAperture(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "diffuser_aperture")

    self.__histogram = ROOT.TH1F( "cut_diffuser_radius", ";r  [mm];# Events", 400, 0.0, 200.0 )
    self.__residual_histogram = ROOT.TH1F( "cut_diffuser_position_residuals", ";z  [#{mu}m];# Events", 2000, -1000.0, 1000.0 )
    self.__z_position = 13740.0
    self.__radius_cut = 90.0
    self.__tolerance = 1.0
    self.__do_cut = False
    self.__missing_global_tracks = 0
    self.__missing_virtual_planes = 0


  def _is_cut(self, analysis_event) :
    if not self.__do_cut :
      return False

    if analysis_event.num_global_tracks() == 0 :
      self.__missing_global_tracks += 1
      return True
#      raise ValueError("No Global Track Was found")

    for track_i in range(analysis_event.num_global_tracks()) :
      global_track = analysis_event.global_track(track_i)
      if global_track.get_status() != 1 and global_track.get_status() != 3 :
        continue

      for i in range(len(global_track)) :
        tp = global_track[i]

        if math.fabs( tp.get_z() - self.__z_position ) < self.__tolerance :
          self.__residual_histogram.Fill( math.fabs( tp.get_z() - self.__z_position )*1000.0 ) # Just to check
          if tp.get_r() < self.__radius_cut :
            return False
          else :
            return True
#      break # Only try one track

    self.__missing_virtual_planes += 1
    return True
#    raise ValueError("No virtual plane was located near the diffuser")


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_global_tracks() == 0 :
      return

    global_track = analysis_event.global_track()
    for i in range(len(global_track)) :
      tp = global_track[i]
      if math.fabs( tp.get_z() - self.__z_position ) < self.__tolerance :
        self.__histogram.Fill( tp.get_r() )
        return

    self.__histogram.Fill(-1.0)


  def _get_plots(self, plot_dict) :
    if self.__radius_cut is not None :
      max_val = self.__histogram.GetMaximum()*1.05
      upper_line = ROOT.TLine(self.__radius_cut, 0.0, self.__radius_cut, max_val)
      self.__histogram.GetListOfFunctions().Add(upper_line)

    plot_dict["projected_diffuser_radius"] = self.__histogram
    plot_dict["position_residual"] = self.__residual_histogram


  def _get_data(self, data_dict) :
    data_dict['missing_virtual_planes'] = self.__missing_virtual_planes
    data_dict['missing_global_tracks'] = self.__missing_global_tracks


  def configure_arguments(self, parser) :
    parser.add_argument('--diffuser_position', type=float, default=self.__z_position, help='z-Position of diffuser')
    parser.add_argument('--diffuser_radius_cut', type=float, default=None, help='Radius to cut tracks')
    parser.add_argument('--position_tolerance', type=float, default=1.0, help='Tolerance when locating items in geometry/reconstruction')


  def parse_arguments(self, namespace) :
    self.__z_position = namespace.diffuser_position
    self.__radius_cut = namespace.diffuser_radius_cut
    self.__tolerance = namespace.position_tolerance

    if self.__radius_cut is None :
      self.__do_cut = False
    else :
      self.__do_cut = True
    


