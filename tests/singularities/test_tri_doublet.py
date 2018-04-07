#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def tri_doublet_vs_vring(tri_doublet, vring):

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.plot(tri_doublet[:,1], tri_doublet[:,4], '-x')
  ax.plot(vring[:,1], vring[:,4])
  ax.grid()
  ax.legend(['Tri doublets', 'Vortex ring'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.plot(tri_doublet[:,1], tri_doublet[:,5], '-x')
  ax.plot(vring[:,1], vring[:,5])
  ax.grid()
  ax.legend(['Tri doublets', 'Vortex ring'])
  plt.show()

if __name__ == "__main__":

  tri_doublet = ld('tri_doublets1.dat')
  vring1 = ld('vring1.dat')

  tri_doublet_vs_vring(tri_doublet, vring1)
