#!/usr/bin/env python
# coding: utf-8

# In[50]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
from matplotlib import cm


# In[51]:


df = pd.read_csv('z1.csv')
df = df.iloc[:, 1:]
z1 = df.values.T
z1.shape

df = pd.read_csv('z2.csv')
df = df.iloc[:, 1:]
z2 = df.values.T

df = pd.read_csv('z3.csv')
df = df.iloc[:, 1:]
z3 = df.values.T

df = pd.read_csv('z22.csv')
df = df.iloc[:, 1:]
z22 = df.values.T

df = pd.read_csv('z33.csv')
df = df.iloc[:, 1:]
z33 = df.values.T

df = pd.read_csv('z4.csv')
df = df.iloc[:, 1:]
z4 = df.values.T


# In[52]:


x = np.linspace(1959, 2008, num=50, dtype=int)
y = np.linspace(1, 365, num=365, dtype=int)
x, y = np.meshgrid(x, y)
x.shape


# In[53]:


fig = plt.figure(figsize=(10,7), facecolor='white')
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, y, z3, color='deepskyblue', alpha=0.2)
ax.plot_surface(x, y, z2,  color='deepskyblue', alpha=0.2)

ax.plot_surface(x, y, z1,  color='red', alpha=0.2)

ax.grid(False)



ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))


ax.set_xlabel('Year', labelpad=10)
ax.set_ylabel('Day', labelpad=10)
ax.set_zlabel('Temperature', labelpad=10)


ax.set_zlim([12,30])
ax.view_init(elev=15, azim=-135)
#ax.view_init(elev=15, azim=45)
plt.tight_layout()
plt.savefig('sydney_temperature1.pdf', format='pdf')
plt.show()


# In[54]:


signed_dev = np.zeros_like(z4)

below = z4 > z22
above = z4 < z33

signed_dev[above] = z4[above] - z22[above]
signed_dev[below] = z4[below] - z33[below]

max_abs = np.max(np.abs(signed_dev))
if max_abs > 0:
    signed_dev_norm = signed_dev / max_abs
else:
    signed_dev_norm = signed_dev

fig = plt.figure(figsize=(10,7), facecolor='white')
ax = fig.add_subplot(111, projection='3d')

ax.set_facecolor('white')
fig.patch.set_facecolor('white')

ax.xaxis.pane.set_facecolor((1,1,1,1))
ax.yaxis.pane.set_facecolor((1,1,1,1))
ax.zaxis.pane.set_facecolor((1,1,1,1))

ax.xaxis.pane.set_edgecolor('white')
ax.yaxis.pane.set_edgecolor('white')
ax.zaxis.pane.set_edgecolor('white')

ax.grid(False)


ax.plot_surface(
    x, y, z22,
    color='#6BAED6',
    alpha=0.20,
    linewidth=0,
    shade=False
)

ax.plot_surface(
    x, y, z33,
    color='#6BAED6',
    alpha=0.20,
    linewidth=0,
    shade=False
)


ax.plot_surface(
    x, y, z4,
    facecolors=cm.coolwarm((signed_dev_norm + 1)/2),
    linewidth=0,
    antialiased=True,
    shade=False
)


ax.set_xlabel('Year', labelpad=10)
ax.set_ylabel('Day', labelpad=10)
ax.set_zlabel('Temperature', labelpad=10)

ax.set_zlim([12,30])
ax.view_init(elev=15, azim=-135)

#mappable = cm.ScalarMappable(cmap=cm.seismic)
mappable = cm.ScalarMappable(cmap=cm.coolwarm)
mappable.set_array(signed_dev)
mappable.set_clim(-max_abs, max_abs)

cbar = fig.colorbar(mappable, shrink=0.6, pad=0.1)
#cbar.set_label('Signed Deviation')

ax.set_zlim([12,30])
ax.view_init(elev=15, azim=-135)
ax.grid(False)


plt.tight_layout()
plt.savefig('sydney_temperature2.pdf', bbox_inches='tight')
plt.show()


# In[55]:


df = pd.read_csv('wu1.csv')
df = df.iloc[:, 1:]
z1 = df.values.T
z1.shape

df = pd.read_csv('wu2.csv')
df = df.iloc[:, 1:]
z2 = df.values.T

df = pd.read_csv('wu3.csv')
df = df.iloc[:, 1:]
z3 = df.values.T

fig = plt.figure(figsize=(10,7), facecolor='white')
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, y, z3, color='deepskyblue', alpha=0.1)
ax.plot_surface(x, y, z2,  color='deepskyblue', alpha=0.1)

ax.plot_surface(x, y, z1,  color='red', alpha=0.1)
#ax.plot_surface(x, y, z3,edgecolors='black', alpha=0.5)

ax.grid(False)


ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))


ax.set_xlabel('Year', labelpad=10)
ax.set_ylabel('Day', labelpad=10)
ax.set_zlabel('Temperature', labelpad=10)


ax.set_zlim([12,30])
ax.view_init(elev=15, azim=-135)
#ax.view_init(elev=15, azim=45)
plt.tight_layout()
plt.savefig('sydney_temperature_wu.pdf', format='pdf')
plt.show()


# In[ ]:




