
import argparse

def get_engine_parser(parser) :
  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS output root files containing reconstructed data.')

  parser.add_argument( '-C', '--no_cuts', action='store_true', help="Don't run the cuts on the data - not always required if using a selection file")
  parser.add_argument( '-A', '--no_analysis', action='store_true', help="Don't run the analysis - e.g. just producing cut plots")

  parser.add_argument( '-q', '--quiet', action='store_true', help="Reduce the verbosity of the output")

  parser.add_argument( '-f', '--find_files', action="store_true", help='Assume the maus_root_files variables are directories and go looking inside for the individual root files')

  parser.add_argument( '-N', '--max_num_events', type=int, help='Maximum number of events to analyse.')

  parser.add_argument( '-O', '--output_filename', default=job_name, help='Set the output filename')

  parser.add_argument( '-D', '--output_directory', default='./', help='Set the output directory')
  parser.add_argument( '-P', '--print_plots', action='store_true', help="Flag to save the plots as individual pdf files" )
  parser.add_argument( '--selection_file', default=None, nargs='+', help='JSON files with a list of spill and event numbers to include in the analysis' )
  parser.add_argument( '-S', '--save_good_events', type=str, default=None, help='Save the good events from the specified analysis module to file' )
  parser.add_argument( '--list_modules', action='store_true', help='Prints out a list of all the modules used in this analysis' )
  parser.add_argument( '--mass_assumption', default=tools.MUON_MASS, type=float, help='Default mass to assume for all tracks' )
  parser.add_argument( '--last_analysis', default=None, help='Base name of the previous analysis run - some modules require this knowledge.' )
#    parser.add_argument( '-D', '--split_by_directory', action='store_true', \
#        help='Tells the script to split the supplied list of root files into "\
#              different data series, depending on the directory they are in.' )



