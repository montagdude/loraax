#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def doublet_vs_vring(doublet, vring):

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.plot(doublet[:,1], doublet[:,4], '-x')
  ax.plot(vring[:,1], vring[:,4])
  ax.grid()
  ax.legend(['Quad doublet', 'Vortex ring'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.plot(doublet[:,1], doublet[:,5], '-x')
  ax.plot(vring[:,1], vring[:,5])
  ax.grid()
  ax.legend(['Quad doublet', 'Vortex ring'])
  plt.show()

if __name__ == "__main__":

  doublet = ld('quad_doublet1.dat')
  vring2 = ld('vring3.dat')

  doublet_vs_vring(doublet, vring2)
