#!/usr/bin/env python
'Generate geometry parameters for 2d/3d RTM'

##   Copyright (C) 2016 University of Texas at Austin
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
import os, sys, numpy
import rsf.prog
import rsf.api as rsf
from subprocess import call

def write2mat(mat,nz,nx,ny,sou_z,sou_ox,sou_oy,sou_jx,sou_jy,sou_nx,sou_ny,rec_z,rec_nx,rec_ny,npad,noff,roll):
    if roll==0: # fixed-spread acquisitiohn
        if ny > 1:
            mod_oz = 0
            mod_nz = nz
            rec_ny_new = rec_ny
            rec_nx_new = rec_nx
            for iy in range(0,sou_ny):
                sou_y = sou_oy + iy*sou_jy
                if sou_y < rec_ny/2:
                    rec_oy = 0
                elif sou_y > ny-rec_ny/2:
                    rec_oy = ny-rec_ny 
                else:
                    rec_oy = sou_y-rec_ny/2
                #mod_oy = rec_oy
                #mod_ny = rec_ny_new
                if rec_oy < npad:
                    mod_oy = 0
                else:
                    mod_oy = rec_oy-npad
                if rec_oy+rec_ny_new > ny-npad:
                    mod_ny = ny-mod_oy;
                else:
                    mod_ny = rec_oy+rec_ny_new+npad-mod_oy;
                for ix in range(0,sou_nx):
                    sou_x = sou_ox + ix*sou_jx
                    if sou_x < rec_nx/2:
                        rec_ox = 0
                    elif sou_x > nx-rec_nx/2:
                        rec_ox = nx-rec_nx 
                    else:
                        rec_ox = sou_x-rec_nx/2
                    #mod_ox = rec_ox
                    #mod_nx = rec_nx_new
                    if rec_ox < npad:
                        mod_ox = 0
                    else:
                        mod_ox = rec_ox-npad
                    if rec_ox+rec_nx_new > nx-npad:
                        mod_nx = nx-mod_ox;
                    else:
                        mod_nx = rec_ox+rec_nx_new+npad-mod_ox;
                    mat.append([mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx_new,rec_ny_new]) 
        else:
            mod_oz = 0
            mod_nz = nz
            mod_oy = 0
            mod_ny = 1
            sou_y  = 0
            rec_oy = 0
            rec_nx_new = rec_nx
            rec_ny_new = 1
            for ix in range(0,sou_nx):
                sou_x = sou_ox + ix*sou_jx
                if sou_x < rec_nx/2:
                    rec_ox = 0
                elif sou_x > nx-rec_nx/2:
                    rec_ox = nx-rec_nx 
                else:
                    rec_ox = sou_x-rec_nx/2
                #mod_ox = rec_ox
                #mod_nx = rec_nx_new
                if rec_ox < npad:
                    mod_ox = 0
                else:
                    mod_ox = rec_ox-npad
                if rec_ox+rec_nx_new > nx-npad:
                    mod_nx = nx-mod_ox;
                else:
                    mod_nx = rec_ox+rec_nx_new+npad-mod_ox;
                mat.append([mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx_new,rec_ny_new]) 

    elif roll==1: # towed marine streamer to negative
        if ny > 1:
            mod_oz = 0
            mod_nz = nz
            for iy in range(0,sou_ny):
                sou_y = sou_oy + iy*sou_jy
                if sou_y < rec_ny/2:
                    rec_oy = 0
                elif sou_y > ny-rec_ny/2:
                    rec_oy = ny-rec_ny 
                else:
                    rec_oy = sou_y-rec_ny/2
                #mod_oy = rec_oy
                #mod_ny = rec_ny_new
                if rec_oy < npad:
                    mod_oy = 0
                else:
                    mod_oy = rec_oy-npad
                if rec_oy+rec_ny_new > ny-npad:
                    mod_ny = ny-mod_oy;
                else:
                    mod_ny = rec_oy+rec_ny_new+npad-mod_oy;
                for ix in range(0,sou_nx):
                    sou_x = sou_ox + ix*sou_jx
                    if sou_x < rec_nx+noff:
                        rec_ox = 0
                        rec_nx_new = sou_x - noff
                    else:
                        rec_ox = sou_x - rec_nx - noff
                        rec_nx_new = rec_nx
                    #mod_ox = rec_ox
                    #mod_nx = rec_nx_new
                    if rec_ox < npad:
                        mod_ox = 0
                    else:
                        mod_ox = rec_ox-npad
                    if rec_ox+rec_nx_new > nx-npad:
                        mod_nx = nx-mod_ox;
                    else:
                        mod_nx = rec_ox+rec_nx_new+npad-mod_ox;
                    mat.append([mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx_new,rec_ny_new]) 
        else:
            mod_oz = 0
            mod_nz = nz
            mod_oy = 0
            mod_ny = 1
            sou_y  = 0
            rec_oy = 0
            rec_ny_new = 1
            for ix in range(0,sou_nx):
                sou_x = sou_ox + ix*sou_jx
                if sou_x < rec_nx+noff:
                    rec_ox = 0
                    rec_nx_new = sou_x - noff
                else:
                    rec_ox = sou_x - rec_nx - noff
                    rec_nx_new = rec_nx
                #mod_ox = rec_ox
                #mod_nx = rec_nx_new
                if rec_ox < npad:
                    mod_ox = 0
                else:
                    mod_ox = rec_ox-npad
                if rec_ox+rec_nx_new > nx-npad:
                    mod_nx = nx-mod_ox;
                else:
                    mod_nx = rec_ox+rec_nx_new+npad-mod_ox;
                mat.append([mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx_new,rec_ny_new]) 

    else: # towed marine streamer to positive
        if ny > 1:
            mod_oz = 0
            mod_nz = nz
            for iy in range(0,sou_ny):
                sou_y = sou_oy + iy*sou_jy
                if sou_y < rec_ny/2:
                    rec_oy = 0
                elif sou_y > ny-rec_ny/2:
                    rec_oy = ny-rec_ny 
                else:
                    rec_oy = sou_y-rec_ny/2
                #mod_oy = rec_oy
                #mod_ny = rec_ny_new
                if rec_oy < npad:
                    mod_oy = 0
                else:
                    mod_oy = rec_oy-npad
                if rec_oy+rec_ny_new > ny-npad:
                    mod_ny = ny-mod_oy;
                else:
                    mod_ny = rec_oy+rec_ny_new+npad-mod_oy;
                for ix in range(0,sou_nx):
                    sou_x = sou_ox + ix*sou_jx
                    rec_ox = sou_x + noff
                    if sou_x > nx-rec_nx-noff:
                        rec_nx_new = nx-sou_x-noff 
                    else:
                        rec_nx_new = rec_nx 
                    #mod_ox = rec_ox
                    #mod_nx = rec_nx_new
                    if rec_ox < npad:
                        mod_ox = 0
                    else:
                        mod_ox = rec_ox-npad
                    if rec_ox+rec_nx_new > nx-npad:
                        mod_nx = nx-mod_ox;
                    else:
                        mod_nx = rec_ox+rec_nx_new+npad-mod_ox;
                    mat.append([mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx_new,rec_ny_new]) 
        else:
            mod_oz = 0
            mod_nz = nz
            mod_oy = 0
            mod_ny = 1
            sou_y  = 0
            rec_oy = 0
            rec_ny_new = 1
            for ix in range(0,sou_nx):
                sou_x = sou_ox + ix*sou_jx
                rec_ox = sou_x + noff
                if sou_x > nx-rec_nx-noff:
                    rec_nx_new = nx-sou_x-noff
                else:
                    rec_nx_new = rec_nx 
                #mod_ox = rec_ox
                #mod_nx = rec_nx_new
                if rec_ox < npad:
                    mod_ox = 0
                else:
                    mod_ox = rec_ox-npad
                if rec_ox+rec_nx_new > nx-npad:
                    mod_nx = nx-mod_ox;
                else:
                    mod_nx = rec_ox+rec_nx_new+npad-mod_ox;
                mat.append([mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx_new,rec_ny_new]) 
 
if __name__ == "__main__":
    try:
        nz     = int(sys.argv[1])  # dimension in z
        nx     = int(sys.argv[2])  # dimension in x
        ny     = int(sys.argv[3])  # dimension in y
        sou_z  = int(sys.argv[4])  # source position in depth
        sou_ox = int(sys.argv[5])  # source starting location in x
        sou_oy = int(sys.argv[6])  # source starting location in y
        sou_jx = int(sys.argv[7])  # source interval in x
        sou_jy = int(sys.argv[8])  # source interval in y
        sou_nx = int(sys.argv[9])  # number of sources in x
        sou_ny = int(sys.argv[10]) # number of sources in y
        rec_z  = int(sys.argv[11]) # receiver position in depth
        rec_nx = int(sys.argv[12]) # number of receivers in x
        rec_ny = int(sys.argv[13]) # number of receivers in y
        npad   = int(sys.argv[14]) # computational domain padding
        noff   = int(sys.argv[15]) # near offset
        roll   = int(sys.argv[16]) # acquisition pattern: 0-> fixed-spread, 1-> towed-streamer to left, 2-> twoed streamer to right
        print "nz=",nz,", nx=",nx,", ny=",ny,", sou_z=",sou_z,", sou_ox=",sou_ox,", sou_oy=",sou_oy,", sou_jx=",sou_jx,", sou_jy=",sou_jy,", sou_nx=",sou_nx,", sou_ny=",sou_ny,", rec_z=",rec_z,", rec_nx=",rec_nx,", rec_ny=",rec_ny,", npad=",npad,", noff=",noff,", roll=",roll
    except:
        print 'Usage:',sys.argv[0],'nz nx ny sou_z sou_ox sou_oy sou_jx sou_jy sou_nx sou_ny rec_z rec_nx rec_ny npad noff(if roll>0) roll'
        print 'Output format:','mod_oz mod_ox mod_oy mod_nz mod_nx mod_ny sou_z sou_x sou_y rec_z rec_ox rec_oy rec_nx rec_ny'
        sys.exit("Execution failed.")

    # double check dimension
    if sou_nx > (nx-sou_ox)/sou_jx:
        sou_nx = (nx-sou_ox)/sou_jx
    if sou_ny > 1 & sou_ny > (ny-sou_oy)/sou_jy:
        sou_ny = (ny-sou_oy)/sou_jy

    # do the work
    dim1=14
    dim2=sou_nx*sou_ny
    mat = []
    write2mat(mat,nz,nx,ny,sou_z,sou_ox,sou_oy,sou_jx,sou_jy,sou_nx,sou_ny,rec_z,rec_nx,rec_ny,npad,noff,roll)
    mat2 = numpy.array(mat)
    output = rsf.Output()
    output.put("n1",dim1)
    output.put("n2",dim2)
    output.settype('int')
    output.write(mat2)

    sys.exit(0)
