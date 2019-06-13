#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 02:15:56 2019

@author: ChrisLin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

t = 0
particle = 2
r = np.zeros((3, particle))

df = pd.read_csv('test_data.csv', sep=',',header=0)
rx = np.array(df['rx'].values)
ry = np.array(df['ry'].values)
rz = np.array(df['rz'].values)

rx = np.reshape(rx,(rx.shape[0]/particle,particle))
ry =  np.reshape(ry,(ry.shape[0]/particle,particle))
rz = np.reshape(rz,(rz.shape[0]/particle,particle))

def data():
    global t
    for i in range(particle):
        r[0][i] = rx[t][i]
        r[1][i] = ry[t][i]
        r[2][i] = rz[t][i]
    t = t + 1
    



Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

def update_r(num):
    global t, r
    data()
    graph._offsets3d = (r[0], r[1], r[2])
    return graph,


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection="3d")
graph = ax.scatter(r[0], r[1], r[2], color='darkblue')

ax.set_xlim3d(0, 1)
ax.set_ylim3d(0, 1)
ax.set_zlim3d(0, 1)

# Creating the Animation object
ani = animation.FuncAnimation(fig, update_r, frames=98, interval=1, blit=False)
ani.save('testABC.mp4', writer=writer)
plt.show()
