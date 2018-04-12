#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def mirror_vs_twopanels(mirror1, mirror2):

  # Plot x velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('Vx')
  ax.plot(mirror1[:,0], mirror1[:,3], '-x')
  ax.plot(mirror2[:,0], mirror2[:,3])
  ax.grid()
  ax.legend(['Mirror', 'Two panels'])
  plt.show()

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('Vy')
  ax.plot(mirror1[:,0], mirror1[:,4], '-x')
  ax.plot(mirror2[:,0], mirror2[:,4])
  ax.grid()
  ax.legend(['Mirror', 'Two panels'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('Vz')
  ax.plot(mirror1[:,0], mirror2[:,5], '-x')
  ax.plot(mirror1[:,0], mirror2[:,5])
  ax.grid()
  ax.legend(['Mirror', 'Two panels'])
  plt.show()

if __name__ == "__main__":

  mirror1 = ld('mirror1.dat')
  mirror2 = ld('mirror2.dat')

  mirror_vs_twopanels(mirror1, mirror2)
