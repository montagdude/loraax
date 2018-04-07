#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def doublet_vs_doubletarr(doublet, doubletarr):

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.plot(doublet[:,1], doublet[:,4], '-x')
  ax.plot(doubletarr[:,1], doubletarr[:,4])
  ax.grid()
  ax.legend(['Quad doublet', 'Doublet array'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.plot(doublet[:,1], doublet[:,5], '-x')
  ax.plot(doubletarr[:,1], doubletarr[:,5])
  ax.grid()
  ax.legend(['Quad doublet', 'Doublet array'])
  plt.show()

if __name__ == "__main__":

  doublet = ld('doublet1.dat')
  doubletarr = ld('doubletarr1.dat')

  doublet_vs_doubletarr(doublet, doubletarr)
