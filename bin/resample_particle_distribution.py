
import MAUS
import sys
import os
import numpy
import scipy.stats
import argparse





PARTICLE_ID_MAP = { -13 : 2, 13 : -2 }




if __name__ == "__main__" :
  parser = argparse.ArgumentParser("KDE Sampler")

  parser.add_argument( "data_file", help="An icool03 type beam file of particles to perform the resampling")

  parser.add_argument( "-N", "--number_particles", default=100, help="Number of particles to generate" )

  parser.add_argument( "-O", "--output_filename", default="resampled_beam", help="Output File Name" )

  parser.add_argument( "-D", "--output_directory", default="./", help="Output Directory" )

  parser.add_argument( "-P", "--particle_id", default=-13, type=int, help="Specify the PID of the beam to produce (Geant4 notation)" )


  namespace = parser.parse_args()


  vectors = None
  weights = None

  vectors = numpy.loadtxt(namespace.data_file, delimiter=' ', usecols=(6, 7, 8, 9, 10, 11), skiprows=2)
  weights = numpy.loadtxt(namespace.data_file, delimiter=' ', usecols=(5), skiprows=2)


#  kernel = scipy.stats.gaussian_kde(numpy.transpose(vectors), weights=weights)
  kernel = scipy.stats.gaussian_kde(numpy.transpose(vectors))

  outfilename = os.path.join(namespace.output_directory, namespace.output_filename)+".txt"
  with open( outfilename, "w" ) as outfile :

    outfile.write( "Custom Input file\n" )
    outfile.write( "0. 0. 0. 0. 0. 0. 0. 0.\n" )

    for count in range( namespace.number_particles ) :
      particle = kernel.resample(1)
#      print particle

      # 0  : Event Number
      value_1 = 1 # Particle Number
      value_2 = PARTICLE_ID_MAP[namespace.particle_id] # Particle ID
      value_3 = 0 # Status
      value_4 = 0.0 # Time
      value_5 = 0.0 # Local Weight
      value_6 = float(particle[0]) # x [m]
      value_7 = float(particle[1]) # y [m]
      value_8 = float(particle[2]) # z [m]
      value_9 = float(particle[3]) # px [GeV/c]
      value_10 = float(particle[4]) # py [GeV/c]
      value_11 = float(particle[5]) # pz [GeV/c]
      value_12 = 0 # sx
      value_13 = 0 # sy
      value_14 = 0 # sz

      outfile.write( str( count ) + ' ' +\
                     str( value_1 ) + ' ' + \
                     str( value_2 ) + ' ' + \
                     str( value_3 ) + ' ' + \
                     str( value_4 ) + ' ' +\
                     str( value_5 ) + ' ' +\
                     str( value_6 ) + ' ' +\
                     str( value_7 ) + ' ' +\
                     str( value_8 ) + ' ' +\
                     str( value_9 ) + ' ' +\
                     str( value_10 ) + ' ' +\
                     str( value_11 ) + ' ' +\
                     str( value_12 ) + ' ' +\
                     str( value_13 ) + ' ' +\
                     str( value_14 ) + '\n' )



