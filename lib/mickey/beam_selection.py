
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


  def requires_parent(self) :
    return self.__requires_parent


  def get_plots(self) :
    plot_dict = {}

    for selector in self.__selectors :
      name, plots = selector.get_plots()
      plot_dict[name] = plots

    return plot_dict


  def get_data(self) :
    data_dict = {}

    for selector in self.__selectors :
      name, data = selector.get_data()
      data_dict[name] = data

    return data_dict


  def weigh_event(self, analysis_event) :
    event_weight = 1.0

    for selector in self.__selectors :
      event_weight *= selector.weigh_event(analysis_event)

    if self.__select_events :
      u = numpy.random.sample() # Random Numbers [0:1)
      if u >= event_weight :
        return False, event_weight

    return True, event_weight


  def configure_arguments(self, parser) :
    parser.add_argument("--selection", action='append', nargs='+', help="Specify the selection algorithm to use, along with the required parameters. (See custom help command for more information.")
    parser.add_argument("--accept_reject", action="store_true", help="Use accept-reject method rather than event weights")
    parser.add_argument("--beam_selection_help", action="store_true", help="Display extend beam selection help")


  def parse_arguments(self, namespace) :
    has_parent = not (LastAnalysis.LastData is None)
    if namespace.beam_selection_help :
      _parsing.print_beam_selection_help()
      raise SystemExit

    self.__select_events = namespace.accept_reject

    selection_routines = namespace.selection

    try :
      for values in selection_routines :
        print "Algorithm:", values

        if values[0] == "parent" :
          self.__selectors.append(selection_modules.ParentAnalysis())

        elif values[0] == "momentum" :
          if len(values) != 3 :
            raise ValueError("Momentum selection requires precisely 2 arguments")
          if not has_parent : raise RuntimeError(values[0])
          self.__selectors.append(selection_modules.SelectMomentum(float(values[1]), float(values[2])))

    except RuntimeError as ex :
      raise RuntimeError("Selection routine, '"+ex.args[0]+"', requires the parent analysis to be provided.")



