#!/usr/bin/env python
'Two layer NN training'

##   Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
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

from __future__ import print_function
import sys
import numpy
try:
        from numpy import *
        import rsf.api as rsf
except Exception as e:
        print('ERROR : need numpy')
        sys.exit(1)

import time

starttime = time.time()

par=rsf.Par()
dat=rsf.Input()
lbl=rsf.Input("label")
valdat=rsf.Input("valdata")
vallbl=rsf.Input("vallabel")
wt1=rsf.Input("weight1")
wt2=rsf.Input("weight2")
bs1=rsf.Input("bias1")
bs2=rsf.Input("bias2")
lr=par.float("lr")
nepoch=par.float("niter")
act=par.float("act") # Activation function - 0:sigmoid 1:tanh 2:relu 3:identity
opt=par.float("opt") # Optimization method - 0:SGD 1:momentum 2:Adam
seed=par.float("seed")
stop=par.float("stop")
lossfunc=par.float("lossfunc") # Loss function - 0:MSE 1:L1
reg=par.float("reg") # Regularization - 0:L2 1:L1
alpha=par.float("alpha") # Regularization coeff. If not, set alpha=0

numpy.random.seed(int(seed))

n1=dat.int('n1')
n2=dat.int('n2')
npatch=dat.int('n3')
u1=zeros((npatch,n2,n1),'f')
dat.read(u1)
u1=numpy.transpose(u1,(0,2,1))
u1.astype('float64')
u2=zeros((npatch,1,n1),'f')
lbl.read(u2)
u2=numpy.transpose(u2,(0,2,1))
u2.astype('float64')

n1val=valdat.int('n1')
u1val=zeros((n2,n1val),'f')
valdat.read(u1val)
u1val=numpy.transpose(u1val,(1,0))
u1val.astype('float64')
u2val=zeros((n1val,1),'f')
vallbl.read(u2val)
u2val.astype('float64')

n3=wt1.int('n1')
u3=zeros((n2,n3),'f')
wt1.read(u3)
u3=numpy.transpose(u3,(1,0))
u3.astype('float64')
u4=zeros((n3),'f')
bs1.read(u4)
u4.astype('float64')

u5=zeros((n3),'f')
wt2.read(u5)
u5.astype('float64')
u6=zeros((1),'f')
bs2.read(u6)
u6.astype('float64')

loss=rsf.Output()
wt1out=rsf.Output("weight1out")
wt2out=rsf.Output("weight2out")
bs1out=rsf.Output("bias1out")
bs2out=rsf.Output("bias2out")
valloss=rsf.Output("valloss")

data = list(zip(u1,u2))
data_val = list(zip(u1val,u2val))

loss_history, loss_val_history = [], []

vw1=numpy.zeros_like(u3)
vb1=numpy.zeros_like(u4)
vw2=numpy.zeros_like(u5)
vb2=numpy.zeros_like(u6)

sw1=numpy.zeros_like(u3)
sb1=numpy.zeros_like(u4)
sw2=numpy.zeros_like(u5)
sb2=numpy.zeros_like(u6)

i = 1
prev_loss = 10000
prev_loss_val = 10000
count = 0
beta=0.9
beta2=0.999
epsilon=1e-8
t=0

while i<= nepoch:
   # Validation
   numpy.random.shuffle(data_val)
   vloss = 0

   for xi, yi in data_val:
      # Forward
      z1 = numpy.dot(u3,xi) + u4
      if act==0:
         a1 = 1 / (1 + numpy.exp(-z1))
      elif act==1:
         a1 = (numpy.exp(z1) - numpy.exp(-z1)) / (numpy.exp(z1) + numpy.exp(-z1))
      elif act==2:
         a1 = z1 * (z1 > 0)
      else:
         a1 = z1       

      z2 = numpy.dot(u5,a1) + u6
      vloss += numpy.square(z2-yi)

   loss_val_history.append(vloss/u2val.shape[0])

   # Training
   numpy.random.shuffle(data)
   tloss = 0

   for xi, yi in data:
      t+=1
      # Forward
      xi1 = numpy.transpose(xi,(1,0))
      z1 = numpy.dot(u3,xi1) + numpy.expand_dims(u4,-1)
      if act==0:
         a1 = 1 / (1 + numpy.exp(-z1))
      elif act==1:
         a1 = (numpy.exp(z1) - numpy.exp(-z1)) / (numpy.exp(z1) + numpy.exp(-z1))
      elif act==2:
         a1 = z1 * (z1 > 0)
      else:
         a1 = z1
         
      z2 = numpy.dot(numpy.reshape(u5,(1,n3)),numpy.reshape(a1,(n3,n1))) + numpy.reshape(u6,(1,1))
      
      # Backward
      if lossfunc==0:
         err_output = numpy.reshape(z2,(n1))-numpy.reshape(yi,(n1))
      else:
         var = numpy.reshape(z2,(n1))-numpy.reshape(yi,(n1))
         err_output = (var>=0)*1.0+(var<0)*(-1.0)
      if reg==0:
         grad_W2 = numpy.dot(numpy.expand_dims(err_output,0),numpy.transpose(a1))/n1 + 2*alpha*u5
      else:
         grad_W2 = numpy.dot(numpy.expand_dims(err_output,0),numpy.transpose(a1))/n1 + alpha*((numpy.expand_dims(u5,0)>0)*1.0+(numpy.expand_dims(u5,0)<0)*(-1.0)) 
           
      if opt==0:
         u5 -= lr*numpy.reshape(grad_W2,(n3))
      elif opt==1:
         vw2 = beta*vw2+(1-beta)*numpy.reshape(grad_W2,(n3))
         u5 -= lr*vw2
      elif opt==2:
         vw2 = beta*vw2+(1-beta)*numpy.reshape(grad_W2,(n3))
         vcorrw2 = vw2/(1-numpy.power(beta,t))
         sw2 = beta2*sw2+(1-beta2)*numpy.power(numpy.reshape(grad_W2,(n3)),2)      
         scorrw2 = sw2/(1-numpy.power(beta2,t))
         u5 -= lr*vcorrw2/numpy.sqrt(scorrw2+epsilon)         

      grad_b2 = numpy.squeeze(numpy.sum(numpy.expand_dims(err_output,0),axis=1,keepdims=True))/n1
      
      if opt==0:
         u6 -= lr*numpy.reshape(grad_b2,(1))
      elif opt==1:
         vb2 = beta*vb2+(1-beta)*numpy.reshape(grad_b2,(1))
         u6 -= lr*vb2
      elif opt==2:
         vb2 = beta*vb2+(1-beta)*numpy.reshape(grad_b2,(1))
         vcorrb2 = vb2/(1-numpy.power(beta,t))
         sb2 = beta2*sb2+(1-beta2)*numpy.power(numpy.reshape(grad_b2,(1)),2)
         scorrb2 = sb2/(1-numpy.power(beta2,t))
         u6 -= lr*vcorrb2/numpy.sqrt(scorrb2+epsilon)
      
      if act==0:
         derivative = a1 * (1 - a1)
      elif act==1:
         derivative = 1-tanh(a1)**2
      elif act==2:
         derivative = 1 * (a1 > 0)
      else:
         derivative = numpy.ones_like(a1)
         
      err_hidden = numpy.expand_dims(err_output,0) * derivative * numpy.expand_dims(u5,-1)
      
      if reg==0:
         grad_W1 = numpy.dot(err_hidden,xi)/n1 + 2*alpha*u3
      else:
         grad_W1 = numpy.dot(err_hidden,xi)/n1 + alpha*((u3>0)*1.0+(u3<0)*(-1.0))

      if opt==0:
         u3 -= lr * grad_W1
      elif opt==1:
         vw1 = beta*vw1+(1-beta)*grad_W1
         u3 -= lr*vw1
      elif opt==2:
         vw1 = beta*vw1+(1-beta)*grad_W1
         vcorrw1 = vw1/(1-numpy.power(beta,t))
         sw1 = beta2*sw1+(1-beta2)*numpy.power(grad_W1,2)
         scorrw1 = sw1/(1-numpy.power(beta2,t))
         u3 -= lr*vcorrw1/numpy.sqrt(scorrw1+epsilon)
      
      grad_b1 = numpy.squeeze(numpy.sum(err_hidden,axis=1,keepdims=True))/n1
      
      if opt==0:
         u4 -= lr*numpy.reshape(grad_b1,(n3))
      elif opt==1:
         vb1 = beta*vb1+(1-beta)*numpy.reshape(grad_b1,(n3))
         u4 -= lr*vb1
      elif opt==2:
         vb1 = beta*vb1+(1-beta)*numpy.reshape(grad_b1,(n3))
         vcorrb1 = vb1/(1-numpy.power(beta,t))
         sb1 = beta2*sb1+(1-beta2)*numpy.power(numpy.reshape(grad_b1,(n3)),2)
         scorrb1 = sb1/(1-numpy.power(beta2,t))
         u4 -= lr*vcorrb1/numpy.sqrt(scorrb1+epsilon)
      
      tloss += numpy.mean(numpy.square(z2-numpy.reshape(yi,(n1))),keepdims=True)
   
   loss_history.append(tloss/u2.shape[0])

   i+=1
 
   # Early stopping
   if vloss/u2val.shape[0] < prev_loss_val and tloss/u2.shape[0] < prev_loss:
     prev_loss_val = vloss/u2val.shape[0]
     prev_loss = tloss/u2.shape[0]
     count = 0
   else:
     count += 1
   if count==0: # Save the best parameters
     u3fin=u3
     u4fin=u4
     u5fin=u5
     u6fin=u6

   # Break the loop when no improvement
   if count == int(stop):
     break

sys.stderr.write('Training epoch: '+str(i-1)+'\n')

loss_history = numpy.array(loss_history)

loss_val_history = numpy.array(loss_val_history)

loss.put('n1',int(i-1))
loss.put('n2',1)
loss.put('n3',1)
loss.put('n4',1)
loss.write(loss_history)

valloss.put('n1',int(i-1))
valloss.put('n2',1)
valloss.put('n3',1)
valloss.put('n4',1)
valloss.write(loss_val_history)

wt1out.put('n1',n2)
wt1out.put('n2',n3)
wt1out.put('n3',1)
wt1out.put('n4',1)
wt1out.write(u3fin)

bs1out.put('n1',n3)
bs1out.put('n2',1)
bs1out.put('n3',1)
bs1out.put('n4',1)
bs1out.write(u4fin)

wt2out.put('n1',n3)
wt2out.put('n2',1)
wt2out.put('n3',1)
wt2out.put('n4',1)
wt2out.write(u5fin)

bs2out.put('n1',1)
bs2out.put('n2',1)
bs2out.put('n3',1)
bs2out.put('n4',1)
bs2out.write(u6fin)

sys.stderr.write('Total time: '+str(time.time()-starttime)+' seconds\n')

