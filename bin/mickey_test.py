#!/usr/bin/env python

# This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
# MAUS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MAUS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MAUS.  If not, see <http://www.gnu.org/licenses/>.
#

import MAUS

import math
import sys
import os


import ROOT
import array

import mickey
from mickey import tof_cuts
from mickey import scifi_cuts
from mickey import diffuser_cut
from mickey import banana_cut
from mickey import analysis_modules



#def calculate_systematic_error(plot) :
#

def conclusion_analysis(plot_dict, data_dict) :
  outdir = data_dict['arguments']['output_directory']
  outfile = data_dict['arguments']['output_filename']

  p_z = array.array('d')
  p_z_error = array.array('d')
  residual = array.array('d')
  residual_error = array.array('d')

  for key, data in data_dict["analysis"]["upstream_emittance_reconstruction"].iteritems() :
    mc_data = data_dict["analysis"]["mc_upstream_emittance_reconstruction"][key]
    low = data['low_edge']
    high = data['high_edge']

    p_z.append( 0.5*(low + high) )
    p_z_error.append( 0.5*(high-low) )

    emittance = plot_dict["analysis"]["upstream_emittance_reconstruction"][key]["emittance"].GetMean()
    emittance_err = plot_dict["analysis"]["upstream_emittance_reconstruction"][key]["emittance"].GetRMS()
    
    residual.append( emittance - mc_data["emittance"] )
    residual_error.append( emittance_err )

  graph = ROOT.TGraphErrors( len(p_z), p_z, residual, p_z_error, residual_error )

  canvas = ROOT.TCanvas("temp_canvas")
  canvas.DrawFrame(170.0, -0.2, 260.0, 0.2)
  graph.SetTitle("")
  graph.GetXaxis().SetTitle("p_{z}   [MeV/c]")
  graph.GetYaxis().SetTitle("#Delta#epsilon  [mm]")
  graph.SetMarkerStyle(26)
  graph.SetMarkerColor(2)

  graph.Draw("p")
  canvas.SaveAs(os.path.join(outdir, outfile+"-residuals.pdf"))

  plot_dict["analysis"]['systematic_residuals'] = graph
  


    


if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  try :

    engine = mickey.Engine('testing_mickey')
#    engine.set_conclude_function(conclusion_analysis)

#    engine.add_cut( tof_cuts.Cut_tof01_spacepoints() )
#    engine.add_cut( tof_cuts.Cut_tof01_time() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_chisq_ndf() )
    engine.add_cut( scifi_cuts.Cut_scifi_upstream_pt() )
#    engine.add_cut( banana_cut.Cut_banana_plot_mass() )
#    engine.add_cut( diffuser_cut.Cut_diffuser_aperture() )

    engine.add_analysis( analysis_modules.EmittanceSystematicErrors() )

#    engine.add_analysis( analysis_modules.EmittanceAnalysis() )
#    engine.add_analysis( analysis_modules.MCEmittanceAnalysis() )

    namespace = engine.process_arguments()

    engine.go()

  except :
    raise

  print 
  print "Complete."
  print

