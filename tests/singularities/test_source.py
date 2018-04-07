#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def source_vs_sourcearr(source, sourcearr):

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.plot(source[:,1], source[:,4], '-x')
  ax.plot(sourcearr[:,1], sourcearr[:,4])
  ax.grid()
  ax.legend(['Quad source', 'Source array'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.plot(source[:,1], source[:,5], '-x')
  ax.plot(sourcearr[:,1], sourcearr[:,5])
  ax.grid()
  ax.legend(['Quad source', 'Source array'])
  plt.show()

if __name__ == "__main__":

  source = ld('source1.dat')
  sourcearr = ld('sourcearr1.dat')

  source_vs_sourcearr(source, sourcearr)
