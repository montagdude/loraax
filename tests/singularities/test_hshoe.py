#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def hshoe_vs_vring(hshoe, vring, title):

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.set_title(title)
  ax.plot(hshoe[:,1], hshoe[:,4])
  ax.plot(vring[:,1], vring[:,4])
  ax.grid()
  ax.legend(['Horseshoe vertex', 'Vortex ring'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.set_title(title)
  ax.plot(hshoe[:,1], hshoe[:,5])
  ax.plot(vring[:,1], vring[:,5])
  ax.grid()
  ax.legend(['Horseshoe vertex', 'Vortex ring'])
  plt.show()

if __name__ == "__main__":

  hshoe = ld('hshoe1.dat')
  vring1 = ld('vring1.dat')
  vring2 = ld('vring2.dat')
  vring3 = ld('vring3.dat')

  hshoe_vs_vring(hshoe, vring1, 'Vortex ring length = 1')
  hshoe_vs_vring(hshoe, vring2, 'Vortex ring length = 2')
  hshoe_vs_vring(hshoe, vring3, 'Vortex ring length = 10')
