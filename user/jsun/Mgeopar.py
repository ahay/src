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
import os, sys
import rsf.prog
from subprocess import call

def write2txt(filename,nz,nx,ny,sou_z,sou_nx,sou_n6,rec_z,rec_nx,rec_ny,nbell):
    f = open(filename,'w')
    
    if ny>1:
        jx = (nx-1-2*nbell)/(sou_nx-1)
        jy = (ny-1-2*nbell)/(sou_ny-1)
        for iy in range(0,sou_ny):
            sou_y = nbell + iy*jy
            if sou_y<rec_ny/2:
                rec_oy = 0
            elif sou_y>ny-rec_ny/2:
                rec_oy = ny-rec_ny 
            else:
                rec_oy = sou_y-rec_ny/2
            for ix in range(0,sou_nx):
                sou_x = nbell + ix*jx
                if sou_x<rec_nx/2:
                    rec_ox = 0
                elif sou_x>nx-rec_nx/2:
                    rec_ox = nx-rec_nx 
                else:
                    rec_ox = sou_x-rec_nx/2
                f.write('%8d %8d %8d %8d %8d %8d %8d %8d\n' %(sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx,rec_ny))
    else:
        sou_y  = 0
        rec_oy = 0
        rec_ny = 1
        jx = (nx-1-2*nbell)/(sou_nx-1)
        for ix in range(0,sou_nx):
            sou_x = nbell + ix*jx
            if sou_x<rec_nx/2:
                rec_ox = 0
            elif sou_x>nx-rec_nx/2:
                rec_ox = nx-rec_nx 
            else:
                rec_ox = sou_x-rec_nx/2
            f.write('%8d %8d %8d %8d %8d %8d %8d %8d\n' %(sou_z,sou_x,sou_y,rec_z,rec_ox,rec_oy,rec_nx,rec_ny))
    
    f.close()

def write2rsf(file_in,file_out,dim1,dim2):

    dd = os.path.join(rsf.prog.RSFROOT,'bin','sfdd')
    command="echo in=%s n1=%d n2=%d data_format=ascii_int | %s type=int form=native > %s" %(file_in,dim1,dim2,dd,file_out)
    print command

    call(command, shell=True)

if __name__ == "__main__":
    try:
        nz     = int(sys.argv[1])  # dimension in z
        nx     = int(sys.argv[2])  # dimension in x
        ny     = int(sys.argv[3])  # dimension in y
        sou_z  = int(sys.argv[4])  # source position in depth
        sou_nx = int(sys.argv[5])  # number of sources in x
        sou_ny = int(sys.argv[6])  # number of source in y
        rec_z  = int(sys.argv[7])  # receiver position in depth
        rec_nx = int(sys.argv[8])  # number of receivers in x
        rec_ny = int(sys.argv[9])  # number of receivers in y
        nbell  = int(sys.argv[10]) # bell width
        print "nz=",nz,", nx=",nx,", ny=",ny,", sou_z=",sou_z,", sou_nx=",sou_nx,", sou_ny=",sou_ny,", rec_z=",rec_z,", rec_nx=",rec_nx,", rec_ny=",rec_ny,", nbell=",nbell
    except:
        print 'Usage:',sys.argv[0],'nz nx ny sou_z sou_nx sou_ny rec_z rec_nx rec_ny nbell'
        sys.exit("Execution failed.")

    # do the work
    dim1=8
    dim2=sou_nx*sou_ny
    write2txt('geo.txt',nz,nx,ny,sou_z,sou_nx,sou_ny,rec_z,rec_nx,rec_ny,nbell)
    write2rsf('geo.txt','geo.rsf',dim1,dim2)

    sys.exit(0)
