#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def source_vs_sourcearr(source, source2, sourcearr):

  # Plot x velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('Vx')
  ax.plot(source[:,0], source[:,3], '-x')
  ax.plot(source2[:,0], source2[:,3], '-+', markersize=10)
  ax.plot(sourcearr[:,0], sourcearr[:,3])
  ax.grid()
  ax.legend(['Quad source', 'Tri sources', 'Source array'])
  plt.show()

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('Vy')
  ax.plot(source[:,0], source[:,4], '-x')
  ax.plot(source2[:,0], source2[:,4], '-+', markersize=10)
  ax.plot(sourcearr[:,0], sourcearr[:,4])
  ax.grid()
  ax.legend(['Quad source', 'Tri sources', 'Source array'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('Vz')
  ax.plot(source[:,0], source[:,5], '-x')
  ax.plot(source2[:,0], source2[:,5], '-+', markersize=10)
  ax.plot(sourcearr[:,0], sourcearr[:,5])
  ax.grid()
  ax.legend(['Quad source', 'Tri sources', 'Source array'])
  plt.show()

if __name__ == "__main__":

  source = ld('source3.dat')
  source2 = ld('source4.dat')
  sourcearr = ld('sourcearr2.dat')

  source_vs_sourcearr(source, source2, sourcearr)
