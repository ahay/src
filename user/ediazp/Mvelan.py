#!/usr/bin/env python
##   Copyright (C) 2014 Colorado School of Mines
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import rsf.api as rsf
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import sys
plt.rcParams.update({'font.size': 20,'legend.fontsize': 20})


def sffile1d(array,o,d,n):
  a = rsf.File(array) # here i am creating an rsf file from numpy array
  return rsf.put(n1=n,o1=o,d1=d,n2=1,n3=1)[a]

def sffile2d(array,o,d,n,o2,d2,n2):
  a = rsf.File(array)
  return rsf.put(n1=n,o1=o,d1=d,n2=n2,o2=o2,d2=d2,n3=1)[a]


class interp:
  def __init__(self,o,d,n):
    self.o = o
    self.d = d
    self.n = n
    self.xout = np.linspace(o,o+(n-1)*d,n)

  def int1d(self,xin,yin,extend=True):
    if len(yin) <2:
      x,y= self.xout,yin[0]*np.ones(len(self.xout),'f')
    else:
      x,y = self.scipy(xin,yin,extend=True)
    return x,y

  def scipy(self,xin,yin,extend=True):
    self.f = interpolate.interp1d(xin,yin,kind='linear',bounds_error=False)
    yout = self.f(self.xout)
    xout = self.xout
    if extend:
      for i in range(len(yout)):
        if np.isnan(yout[i]):
          if i <len(yout)/2:
            yout[i] = yin[0]
          else:
            yout[i] = yin[len(yin)-1]
    else: 
      pops = []
      for i in range(len(yout)):
        if np.isnan(yout[i]):
           pops += [i]
      yout = np.delete(yout,pops)
      xout = np.delete(xout,pops)
    return xout,yout

class interact:
  def __init__(self,ov,dv,nv,
                    ot,dt,nt,     
                    oo,do,no,
                    half):
    self.ov = ov
    self.dv = dv
    self.nv = nv
    self.mv = ov + (nv-1)*dv

    self.ot = ot
    self.dt = dt
    self.nt = nt
    self.mt = ot + (nt-1)*dt

    self.oo = oo
    self.do = do
    self.no = no
    self.mo = oo + (no-1)*do

    self.interp = interp(ot,dt,nt)
    self.fign = 0
    self.half = half

  def plot(self,array,y,x,gather):
    px = []
    py = []
    extend = [self.ox,self.mx,self.my,self.oy]
    
    fig = plt.figure(1)
    ax = fig.add_subplot(121)
    im = ax.imshow(array.T,aspect='auto',extent=extend)
    ax.set_ylabel("t(s)")
    ax.set_xlabel("vel(km/s)")
    ax.set_xlim(self.ox, self.mx)
    ax.set_ylim(self.oy, self.my)
    ax.invert_yaxis() 
    ax.plot(x,y,'yo')
    
    ax2 = fig.add_subplot(122, sharey=ax)
    extend2 = [oo+(no-1)*do,oo,self.my,self.oy]
    ax2.set_xlabel("offset(km)")
    ax2.imshow(gather.T,aspect='auto',extent=extend2,cmap='gray')
    for label in ax2.get_yticklabels():
      label.set_visible(False)

    plt.show()
    
  def one(self,semblance,gather,offset=None,xold=None,yold=None):  
    global px, py
    px = []
    py = []
    vmin = gather.min()*0.1
    vmax = gather.max()*0.1

    extend2 = [self.oo,self.mo,self.mt,self.ot]
 
    fig,ax,ax2 = self.getplots(semblance,gather)
    pickable_artists = []
    if xold==None:
      pt, = ax.plot(self.ov +self.nv*0.5*self.dv, self.ot + self.nt*0.5*self.dt, 'yo')  # 5 points tolerance
    else:
      ptold, = ax.plot(xold,yold, 'ro')  # 5 points tolerance
      pt, = ax.plot(self.ov +self.nv*0.5*self.dv, self.ot + self.nt*0.5*self.dt, 'yo')  # 5 points tolerance
    pickable_artists.append(pt)
    
    def onclick(event):
      global px,py
      if event.inaxes is not None and not hasattr(event, 'already_picked'):
        ax = event.inaxes
        remove = [artist for artist in pickable_artists if artist.contains(event)[0]]
        x, y = ax.transData.inverted().transform_point([event.x, event.y])
        if not remove:
          px += [x]
          py += [y]
          # add a pt
          pt, = ax.plot(x, y, 'yo', picker=5)
          pickable_artists.append(pt)
        else:
          for artist in remove:
            try:
              artist.remove()
            except:
              pass
          hit=False
          for i in range(len(px)):
            if abs(py[i] - y) <5*self.dt:
              pyy = py[i] 
              pxx = px[i]
              hit=True
              break
          if hit:
            py.remove(pyy)
            px.remove(pxx)
        t,v = self.interp.int1d(py,px)
        b = self.nmo(v,gather,offset)
        ax2.imshow(b.T,aspect='auto',extent=extend2,cmap='gray',vmin=vmin,vmax=vmax)
        plt.draw()

    fig.canvas.mpl_connect('button_release_event', onclick)
    plt.show()
    t,v = self.interp.int1d(py,px)
    return t,v,px,py


  def getplots(self,semblance,gather):
    vmin = gather.min()*0.3
    vmax = gather.max()*0.3
    extend = [self.ov,self.mv,self.mt,self.ot]
    self.fign+=1
    fig = plt.figure(self.fign)
    ax = fig.add_subplot(121)
    ax.imshow(semblance.T,aspect='auto',extent=extend)
    ax.set_xlim(self.ov, self.mv)
    ax.set_ylim(self.ov, self.mv)
    ax.set_ylabel("t(s)")
    ax.set_xlabel("vel(km/s)")
    ax.invert_yaxis() 

    ax2 = fig.add_subplot(122, sharey=ax)
    extend2 = [self.oo,self.mo,self.mt,self.ot]
    ax2.set_xlabel("offset(km)")
    ax2.imshow(gather.T,aspect='auto',extent=extend2,cmap='gray',vmin=vmin,vmax=vmax)
    for label in ax2.get_yticklabels():
      label.set_visible(False)
    return fig,ax,ax2    

  def nmo(self,v,gather,offsets=None):
    vel = sffile1d(v,self.ot,self.dt,self.nt)
    c = sffile2d(gather,self.ot,self.dt,self.nt,self.oo,self.do,self.no)
    if not offsets==None:
      of = sffile1d(offsets,self.oo,self.do,self.no)
      d = rsf.nmo(velocity=vel,half=self.half,offset=of)[c]
    else:
      d = rsf.nmo(velocity=vel,half=self.half)[c]
    b = np.array(d)
    # close temp rsf files
    d.close()
    vel.close()
    c.close()
    return b    
    

def get_axis(File,axis):
  o = File.float("o%d"%axis)
  d = File.float("d%d"%axis)
  n = File.int("n%d"%axis)    
  return o,d,n

def put_axis(File,axis,o,d,n):
  File.put("o%d"%axis,o)
  File.put("d%d"%axis,d)
  File.put("n%d"%axis,n)    



par = rsf.Par()
Fsemb = rsf.Input()
Fcmp  = rsf.Input("cmp")
off = par.bool("useoffset",True) # if irregular offset, pass it
half = par.bool("half",True)  # half or full offset?
if off:
  Foff = rsf.Input("offset")
  n1 = Foff.int("n1")
  offsets = np.zeros(n1,'f')
  Foff.read(offsets)
  
Fout = rsf.Output()



o1,d1,n1 = get_axis(Fsemb,1)
o2,d2,n2 = get_axis(Fsemb,2)
o3,d3,n3 = get_axis(Fsemb,3)

oo,do,no = get_axis(Fcmp,2)

put_axis(Fout,2,o3,d3,n3)
put_axis(Fout,3,o3,d3,1)



interactivePlot = interact(o2,d2,n2,o1,d1,n1,oo,do,no,half)
panel = np.zeros((n2,n1),'f')
gather = np.zeros((no,n1),'f')


for i3 in range(n3):
  print('CMP analysis %02d of %02d'%(i3+1,n3),file=sys.stderr)

  Fsemb.read(panel)
  Fcmp.read(gather)
  
  if off and i3==0:
    x,y,xold,yold = interactivePlot.one(panel,gather,offsets)
  if off and i3>0:
    x,y,xold,yold = interactivePlot.one(panel,gather,offsets,xold,yold)
  if off==False and i3==0:
    x,y,xold,yold = interactivePlot.one(panel,gather)
  if off==False and i3>0:
    x,y,xnew,ynew = interactivePlot.one(panel,gather,xold=xold,yold=yold)
    xold=xnew
    yold=ynew
  Fout.write(y)

Fsemb.close()
Fout.close()
