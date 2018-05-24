
import ROOT
import math

from _cuts_base import Cut_Base


####################################################################################################
class VirtualCuts(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "virtual_cuts")

    self.__current_track = None
    self.__global_cuts = True

    self.__cut_transmision = [0, 0]
    self.__cut_pid = [0]

    self.__plane_cuts = {}

    self.__plane_list = []
    self.__plane_list_length = 0
    self.__virt_ensemble_size = 0


  def _is_cut(self, analysis_event) :


  def _get_plots(self, plot_dict) :


  def _get_data(self, data_dict) :


  def configure_arguments(self, parser) :
#    parser.add_argument('--virt_cut_pid', nargs='+', type=int, default=[13,-13], help='Cut on the PIDs of the particles in the beam')
    parser.add_argument('--virt_cut_transmission', type=float, nargs=2, default=None, help='Expect particles to pass through virtual planes')
    parser.add_argument('--virt_cut_radius', nargs='+', metavar='cut_info', action='append', default=None, help='Cut on the radius of beam at multiple planes.')
    parser.add_argument('--virt_pt_window', nargs='+', metavar='cut_info', action='append', default=None, help='Cut on a window of Pt of the particles in the beam at multple virtual planes')
    parser.add_argument('--virt_pz_window', nargs='+', action='append', metavar='cut_info', default=None, help='Cut on a window of Pz of the particles in the beam at multiple virtual planes')
    parser.add_argument('--virt_momentum_window', nargs='+', action='append', metavar='cut_info', default=None, help='Cut on a window of total momentum of the particles in the beam at multiple virtual planes')

    parser.add_argument('--not_global_cuts', action='store_true', help='Select whether a plane-based cut globally removes the track from the analysis.' )


  def parse_arguments(self, namespace) :
    self.__cut_pid = namespace.virt_cut_pid
    self.__cut_transmission = namespace.virt_cut_transmission
    
    if namespace.not_global_cuts :
      self.__global_cuts = False

    if namespace.virt_pt_window is not None :
      for option in namespace.virt_pt_window :
        cut_low = float(option[0])
        cut_high = float(option[1])
        for plane_id_str in option[2:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = _PlaneCut(plane_id)
          self.__plane_cuts[plane_id].set_pt_window(cut_low, cut_high)

    if namespace.virt_pz_window is not None :
      for option in namespace.virt_pz_window :
        cut_low = float(option[0])
        cut_high = float(option[1])
        for plane_id_str in option[2:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = _PlaneCut(plane_id)
          self.__plane_cuts[plane_id].set_pz_window(cut_low, cut_high)

    if namespace.virt_momentum_window is not None :
      for option in namespace.virt_momentum_window :
        cut_low = float(option[0])
        cut_high = float(option[1])
        for plane_id_str in option[2:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = _PlaneCut(plane_id)
          self.__plane_cuts[plane_id].set_momentum_window(cut_low, cut_high)

    if namespace.virt_cut_radius is not None :
      for option in namespace.virt_cut_radius :
        cut = float(option[0])
        for plane_id_str in option[1:] :
          plane_id = int(plane_id_str)
          if plane_id not in self.__plane_cuts :
            self.__plane_cuts[plane_id] = _PlaneCut(plane_id)
          self.__plane_cuts[plane_id].set_radius(cut)


####################################################################################################
class _PlaneCut() :

  def __init__(self, plane_id) :
    self.__plane_id = plane_id

    self.__radius = 1.0E+30
    self.__pt_window = (0.0, 1.0E+30)
    self.__pz_window = (0.0, 1.0E+30)
    self.__momentum_window = (0.0, 1.0E+30)


  def __repr__(self) :
    return "{0} : {7}, ({1} : {2}), ({3} : {4}), ({5} : {6})".format(self.__plane_id, self.__pt_window[0], self.__pt_window[1], self.__pz_window[0], self.__pz_window[1], self.__momentum_window[0], self.__momentum_window[1], self.__radius)

  def set_radius(self, radius) :
    self.__radius = radius


  def set_pt_window(self, pt_low, pt_high) :
    self.__pt_window = (pt_low, pt_high)


  def set_momentum_window(self, p_low, p_high) :
    self.__momentum_window = (p_low, p_high)


  def set_pz_window(self, pz_low, pz_high) :
    self.__pz_window = (pz_low, pz_high)


  def get_plane_id(self) :
    return self.__plane_id


  def is_cut(self, vhit) :
    rad = math.sqrt(vhit.GetPosition().x()**2 + vhit.GetPosition().y()**2)
    pt = math.sqrt(vhit.GetMomentum().x()**2 + vhit.GetMomentum().y()**2)
    p = math.sqrt(vhit.GetMomentum().x()**2 + vhit.GetMomentum().y()**2 + vhit.GetMomentum().z()**2)

    if rad > self.__radius :
      return True
    elif pt < self.__pt_window[0] :
      return True
    elif pt > self.__pt_window[1] :
      return True
    elif vhit.GetMomentum().z() < self.__pz_window[0] :
      return True
    elif vhit.GetMomentum().z() > self.__pz_window[1] :
      return True
    elif p < self.__momentum_window[0] :
      return True
    elif p > self.__momentum_window[1] :
      return True
  
    return False

