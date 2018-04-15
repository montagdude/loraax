#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def hshoe_vs_vring_vs_doublet(hshoe, vring, doublet, title):

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.set_title(title)
  ax.plot(hshoe[:,1], hshoe[:,4])
  ax.plot(vring[:,1], vring[:,4], '-o')
  ax.plot(doublet[:,1], doublet[:,4], '-x')
  ax.grid()
  ax.legend(['Horseshoe vertex', 'Vortex ring', 'Quad doublet'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.set_title(title)
  ax.plot(hshoe[:,1], hshoe[:,5])
  ax.plot(vring[:,1], vring[:,5], '-o')
  ax.plot(doublet[:,1], doublet[:,5], '-x')
  ax.grid()
  ax.legend(['Horseshoe vertex', 'Vortex ring', 'Quad doublet'])
  plt.show()

if __name__ == "__main__":

  hshoe = ld('hshoe1.dat')
  vring1 = ld('vring1.dat')
  vring2 = ld('vring2.dat')
  vring3 = ld('vring3.dat')
  doublet1 = ld('quad_doublet3.dat')
  doublet2 = ld('quad_doublet4.dat')
  doublet3 = ld('quad_doublet5.dat')

  hshoe_vs_vring_vs_doublet(hshoe, vring1, doublet1, 'Ring/doublet length = 1')
  hshoe_vs_vring_vs_doublet(hshoe, vring2, doublet2, 'Ring/doublet length = 2')
  hshoe_vs_vring_vs_doublet(hshoe, vring3, doublet3, 'Ring/doublet length = 10')
