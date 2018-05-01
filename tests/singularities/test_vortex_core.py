#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def plot_velocity(vortex1, vortex2):

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.plot(vortex1[:,1], vortex1[:,5], '-o')
  ax.plot(vortex2[:,1], vortex2[:,5], '-x')
  ax.grid()
  ax.legend(['rcore = 1e-3', 'rcore = 1e-12'])
  plt.show()

if __name__ == "__main__":

  vortex1 = ld('vortex_core1.dat')
  vortex2 = ld('vortex_core2.dat')
  plot_velocity(vortex1, vortex2)
