#!/usr/bin/env python

from data_loader import load_data as ld
from matplotlib import pyplot as plt

def exact_vs_fd(panel1, panel2):

  # Plot x velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vx')
  ax.plot(panel1[:,1], panel1[:,3], '-o')
  ax.plot(panel2[:,1], panel2[:,3], '-x')
  ax.grid()
  ax.legend(['Exact', 'Finite difference potential'])
  plt.show()

  # Plot y velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vy')
  ax.plot(panel1[:,1], panel1[:,4], '-o')
  ax.plot(panel2[:,1], panel2[:,4], '-x')
  ax.grid()
  ax.legend(['Exact', 'Finite difference potential'])
  plt.show()

  # Plot z velocity
  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('y')
  ax.set_ylabel('Vz')
  ax.plot(panel1[:,1], panel1[:,5], '-o')
  ax.plot(panel2[:,1], panel2[:,5], '-x')
  ax.grid()
  ax.legend(['Exact', 'Finite difference potential'])
  plt.show()

if __name__ == "__main__":

  panel1 = ld('tri_panel1.dat')
  panel2 = ld('tri_panel2.dat')

  exact_vs_fd(panel1, panel2)
