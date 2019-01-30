
import selection_modules
import _parsing
import LastAnalysis

import numpy

## Define the beam selection class.
# Controls the interface to the beam selection routines through a single argument
class BeamSelection(object) :

####################################################################################################
  def __init__(self) :
    self.__selectors = []
    self.__requires_parent = False
    self.__select_events = False
    self.__parent_analysis = None
    self.__selection_normalisation = 1.0


  def requires_parent(self) :
    return self.__requires_parent


  def get_plots(self) :
    plot_dict = {}

    plot_dict['parent_analysis'] = self.__parent_analysis.get_plots()

    for selector in self.__selectors :
      name, plots = selector.get_plots()
      plot_dict[name] = plots

    return plot_dict


  def get_data(self) :
    data_dict = {}

    data_dict['parent_analysis'] = self.__parent_analysis.get_data()

    for selector in self.__selectors :
      name, data = selector.get_data()
      data_dict[name] = data

    return data_dict


  def weigh_event(self, analysis_event, event_weight=1.0) :
    keep = True

    for selector in self.__selectors :
      event_weight *= selector.weigh_event(analysis_event)

    if self.__select_events :
      event_weight *= self.__selection_normalisation

      u = numpy.random.sample() # Random Numbers [0:1)
      if u >= event_weight :
        keep = False
        event_weight = 0.0
      else :
        keep = True
        event_weight = 1.0

    self.__parent_analysis.fill_plots(analysis_event, event_weight)

    return keep, event_weight


  def configure_arguments(self, parser) :
    parser.add_argument("--selection", action='append', nargs='+', help="Specify the selection algorithm to use, along with the required parameters. (See custom help command for more information.")
    parser.add_argument("--accept_reject", action="store_true", help="Use accept-reject method rather than event weights")
    parser.add_argument("--beam_selection_help", action="store_true", help="Display detailed beam selection help")


  def parse_arguments(self, namespace) :
    has_parent = not (LastAnalysis.LastData is None)
    if namespace.beam_selection_help :
      _parsing.print_beam_selection_help()
      raise SystemExit

    self.__select_events = namespace.accept_reject

    selection_routines = namespace.selection

    if selection_routines is not None :
      try :
        for values in selection_routines :
          if values[0] == "momentum" :
            if len(values) != 3 :
              raise ValueError("Momentum selection requires precisely 2 arguments")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.SelectMomentum(float(values[1]), float(values[2])))

          if values[0] == "longitudinal_momentum" :
            if len(values) != 3 :
              raise ValueError("Longitudinal momentum selection requires precisely 2 arguments")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.SelectLongMomentum(float(values[1]), float(values[2])))

          elif values[0] == "amplitude" :
            if len(values) != 2 :
              raise ValueError("Amplitude selection requires precisely 1 argument")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.SelectAmplitude(float(values[1])))

          elif values[0] == "amplitude_cut" :
            if len(values) != 2 :
              raise ValueError("Amplitude cut requires precisely 1 argument")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.CutAmplitude(float(values[1])))

          elif values[0] == "uncorrelated4d" :
            if len(values) != 3 :
              raise ValueError("Amplitude selection requires precisely 2 arguments")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.SelectUncorrelated4D(float(values[1]), float(values[2])))

          elif values[0] == "voronoi" :
            if len(values) != 1 :
              raise ValueError("Voronoi phasespace selection requires precisely 1 argument")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.VoronoiPhaseSpaceSelection())

          elif values[0] == "kde" :
            if len(values) != 1 :
              raise ValueError("KDE selection requires precisely 1 argument")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.KDESelection())

          elif values[0] == "analytic_beam" :
            if len(values) != 6 :
              raise ValueError("Analytic 4D phasespace selection requires precisely 4 argument")
            if not has_parent : raise RuntimeError(values[0])
            self.__selectors.append(selection_modules.SelectAnalyticBeam(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5])))

          else :
            raise ValueError("Unknown beam selection routine requested, '"+values[0]+"'.")
      except RuntimeError as ex :
        raise RuntimeError("Selection routine, '"+ex.args[0]+"', requires the parent analysis to be provided.")
 
    self.__parent_analysis = selection_modules.ParentAnalysis()

    self.__selection_normalisation = 1.0
    for selector in self.__selectors :
      self.__selection_normalisation *= selector.get_normalisation()


