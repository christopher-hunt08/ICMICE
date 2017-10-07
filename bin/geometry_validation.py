#!/bin/env python

import MAUS
import argparse
import json
import ROOT
import Configuration  # MAUS configuration (datacards) module
import maus_cpp.globals as maus_globals # MAUS C++ calls
import maus_cpp.material # MAUS Geometry Navigator class
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


RANGE_Z = None
RANGE_R = None
GEOMETRY_FILE = None
STEP_SIZE = None


def main():
  """
  Spit out material information for the example geometry.
  """
  # set up datacards
  print "Welcome to MAUS geometry navigator example"
  configuration = Configuration.Configuration().getConfigJSON()
  conf_json = json.loads(configuration)
  conf_json["simulation_geometry_filename"] = GEOMETRY_FILE
  configuration = json.dumps(conf_json)

  # initialise field maps and geometry
  print "Building field maps and other initialisation (this can take a while)"
  maus_globals.birth(configuration)
  
  maus_cpp.material.set_position( 0.0, 0.0, RANGE_Z[0] )

  print "Analysing Geometry"
  print "Please Wait..."
  print

  z_steps = int((RANGE_Z[1] - RANGE_Z[0]) / STEP_SIZE)
  r_steps = int(RANGE_R / STEP_SIZE)

  material_plot = ROOT.TH2F("material_locations", "", z_steps, RANGE_Z[0], RANGE_Z[1], r_steps, 0.0, RANGE_R)
  material_list = ["Galactic"]

  pos_r = 0.0
  pos_z = RANGE_Z[0]

  for r in range( r_steps ) :
    pos_z = RANGE_Z[0]
    maus_cpp.material.set_position( pos_r, 0.0, pos_z )
    for z in range( z_steps ) :
      maus_cpp.material.step( 0.0, 0.0, STEP_SIZE )
      data = maus_cpp.material.get_material_data()
#      print pos_z, pos_r, data
#      print

      bin_x = material_plot.GetXaxis().FindBin(pos_z)
      bin_y = material_plot.GetYaxis().FindBin(pos_r)

      material = data['name']
      if material.startswith("G4_") :
        material = material[3:]

      if material in material_list :
#        material_plot.SetBinContent(pos_z, pos_r, float(material_list.index(material)))
        material_plot.SetBinContent(bin_x, bin_y, float(material_list.index(material)))
      else :
#        material_plot.SetBinContent(pos_z, pos_r, float(len(material_list)))
        material_plot.SetBinContent(bin_x, bin_y, float(len(material_list)))
        material_list.append(material)

      pos_z += STEP_SIZE

    pos_r += STEP_SIZE


  canvas = ROOT.TCanvas("temp", "temp")
  material_plot.Draw("COLZ")
  canvas.SaveAs("geometry_validation_out.pdf", "pdf")

  outfile = ROOT.TFile("geometry_validation_out.root", "RECREATE")
  material_plot.Write()
  outfile.Close()



  print 
  print "Found the following materials:"
  print
  for num in range(len(material_list)) :
    print num, "\t", material_list[num]

  print

  # Clean up
  maus_globals.death()
  # Finished
  print "Finished"



if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Examines MICE Geometry materials")

  parser.add_argument('--range_z', nargs=2, type=float, default=[-3000.0, 3000.0])
  parser.add_argument('--range_r', type=float, default=1000.0)
  parser.add_argument('--step_size', type=float, default=10.0)
  parser.add_argument('geometry')

  namespace = parser.parse_args()

  RANGE_Z = namespace.range_z
  RANGE_R = namespace.range_r
  STEP_SIZE = namespace.step_size
  GEOMETRY_FILE = namespace.geometry

  main()
