
from .. import LastAnalysis

from _analysis_base import Analysis_Base
from analysis import inspectors

import ROOT
import math
import array


class MCInspections(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "mc_inspections", require_mc=True)

    self.__inspectors = {}
    self.__inspect_primaries = False
    self.__PID = None


  def analyse_event(self, analysis_event, weight=1.0) :
    for plane_id, inspector in self.__inspectors.iteritems() :
      hit = analysis_event.mc_virtual_hit(plane_id)
      if hit is None :
        continue
      if self.__PID is not None :
        if hit.get_pid() != self.__PID :
          continue
      hit.set_weight(weight)
      inspector.add_hit(hit)

    if self.__inspect_primaries :
      hit = analysis_event.mc_primary()
      if hit is not None :
        if self.__PID is not None :
          if hit.get_pid() != self.__PID :
            return
        hit.set_weight(weight)
        self.__inspectors['primaries'].add_hit(hit)


  def _get_plots(self, plot_dict) :
    for plane_id, inspector in self.__inspectors.iteritems() :
      plot_dict[str(plane_id)] = inspector.get_plot_dictionary()


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--inspection_planes", nargs='+', type=int, default=[], help="Specify the virtual planes to inspect." )
    parser.add_argument( "--inspect_primaries", action='store_true', help="Specify whether in inspector the primary particles." )
    parser.add_argument( '--pid', type=int, default=None, help="Specify the Geant4 PID to analyse." )


  def parse_arguments(self, namespace) :
    for plane_id in namespace.inspection_planes :
      self.__inspectors[plane_id] = inspectors.PhaseSpace2DInspector(plane_id, 2000)

    if namespace.inspect_primaries :
      self.__inspect_primaries = True
      self.__inspectors['primaries'] = inspectors.PhaseSpace2DInspector(-1001, 2000)

    self.__PID = namespace.pid


  def conclude(self) :
    pass


