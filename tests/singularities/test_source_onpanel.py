#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def below_and_above(below, above, title):

  # Plot Phi
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('Phi')
  ax.set_ylabel('z')
  ax.set_title(title)
  ax.plot(below[:,3], below[:,2], '-x')
  ax.plot(above[:,3], above[:,2], '-o')
  ax.grid()
  ax.legend(['Below panel', 'Above panel'])
  plt.show()

if __name__ == "__main__":

  below1 = ld('quad_source1_below.dat')
  above1 = ld('quad_source1_above.dat')
  below2 = ld('quad_source2_below.dat')
  above2 = ld('quad_source2_above.dat')

  below_and_above(below1, above1, "Inside panel border")
  below_and_above(below2, above2, "Outside panel border")
