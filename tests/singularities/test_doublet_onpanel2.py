#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def below_and_above(below, above):

  # Plot Phi
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('Phi')
  ax.set_ylabel('z')
  ax.plot(below[:,3], below[:,2], '-x')
  ax.plot(above[:,3], above[:,2], '-o')
  ax.grid()
  ax.legend(['Below panel', 'Above panel'])
  plt.show()

if __name__ == "__main__":

  below = ld('quad_doublet2_below.dat')
  above = ld('quad_doublet2_above.dat')

  below_and_above(below, above)
