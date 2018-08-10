#!/usr/bin/env python
'''Plot wavefields over a background image
'''

# Copyright (C) 2018 Carlos Alberto da Costa Filho
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.

# Guides to m8r contained in
# http://www.ahay.org/rsflog/index.php?/archives/
#    173-Extending-Python-interface.html
# http://www.ahay.org/rsflog/index.php?/archives/
#    264-Running-Madagascar-in-an-interactive-console.html

from __future__ import print_function
import sys
import timeit
import m8r
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import animation
from matplotlib.colors import colorConverter
from mpl_toolkits.axes_grid1 import make_axes_locatable

def sf_warning(message):
    "Prints warning."
    print("sfpgreywfl: %s" % message, file=sys.stderr)

def sf_error(message):
    "Prints fatal warning and quits."
    sf_warning(message)
    sys.exit(1)

def read_axis(infile, axis=1):
    "Returns n, d, o, u, l of specified axis a from infile of type m8r.Input"
    return (infile.int("n"+str(axis)), infile.float("d"+str(axis)),
            infile.float("o"+str(axis)), infile.string("unit"+str(axis)),
            infile.string("label"+str(axis)))

def create_tranparent_cmap(cmap, treshold=0.01):
    '''
    Creates a linear colormap from `color1` to `color2` in which everything
    below a certain `treshold` is set to transparent.

    Parameters
    ----------
    cmap : string
        For valid cmaps see https://matplotlib.org/users/colormaps.html
    treshold : float, optional
        Treshold below which colors will be transparent. Must be a number
        between 0 and 1.

    Returns
    -------
    cmap : matplotlib.colors.LinearSegmentedColormap
        Linear colormap from `color1` to `color2` with values lower then
        `treshold` set to transparent.
    '''

    cmap = cm.get_cmap(cmap, lut=256)
    cmap._init()
    N = cmap.N+3
    alphas = [a if a > treshold else 0 for a in abs(np.linspace(-1.0, 1.0, N))]
    N2 = len([a for a in alphas if a < treshold]) # No. of elements with 0 alpha
    N3 = (N-N2)/2 # No. of elements with nonzero alpha
    N4 = (N+N2)/2

    # These next commands will taper alphas, so that no abrupt transition
    # between zero and nonzero elements occur
    window = np.kaiser(N-N2, 2)
    alphas[:N3] *= window[N3:]
    alphas[N4:] *= window[:N3]
    cmap._lut[:, -1] = alphas
    return cmap

if __name__ == "__main__":
    par = m8r.Par()
    inp = m8r.Input()
    shape = inp.shape()
    if len(shape) != 3:
        sf_error("Must have 3 axes")
    if inp.type != 'float':
        sf_error("Only supports float")
    bgf = par.string("bg", None) # Background for animation
    if bgf:
        bg  = m8r.Input(bgf)
        if bg.type != 'float':
            sf_error("Only supports float")
        bgshape = bg.shape()
        if len(bgshape) != 2 or bgshape != shape[1:]:
            sf_error("Background must have the same two last dimensions of input")
    else:
        sf_warning("No background")
    savefile = par.string("savefile") # Save animation to file. If not present, display animation.

    # Parameters
    title = par.string("title") # Plot title
    pclip = par.float("pclip", 100)  # Clip amplitude percentage from (0-100)
    aclip = par.bool("absclip", False)  # Clipping is done for all gathers rather than per frame (y/n)
    cbar  = par.bool("scalebar", True) # Colorbar (y/n)
    cbarl = par.string("barlabel") # Colorbar label
    ttxt = par.bool("timetext", True) # Time text (y/n)
    bmap = par.string("bgcmap", 'viridis') # Background colormap. See https://matplotlib.org/users/colormaps.html
    wmap = par.string("wflcmap", 'gray') # Wavefield colormap (should be sequential)
    jsnap = par.int("jsnap", 1)  # Number of timesteps at which to plot wavefield
    aspect = par.int("aspect", 1)  # Aspect ratio
    tmin  = par.float("tmin", None)  # Minimum time
    tmax  = par.float("tmax", None)  # Maximum time
    figx  = par.float("figx", 10) # Figure x size in inches
    figy  = par.float("figy", 8) # Figure y size in inches
    xints = par.int("xints", None) # Max number of x intervals
    yints = par.int("yints", None) # Max number of y intervals
    fps   = par.float("fps", None) # Frames per second (when saving file)
    verb  = par.bool("verb", False) # Verbosity flag

    # Error/bounds checking
    if pclip < 0 or pclip > 100:
        sf_warning("pclip must be between 0 and 100, using 100")
        pclip = 100.0
    pclip /= 100.0
    if cbarl:
        cbar = True
    if jsnap < 1:
        sf_warning("setting jsnap=1")
        jsnap = 1
    if figx < 0:
        sf_warning("setting figx=6")
        figx = 6
    if figy < 0:
        sf_warning("setting figy=8")
        figy = 6
    if xints < 0 and xints is not None:
        sf_warning("setting xints=None")
        xints = None
    if yints < 0 and yints is not None:
        sf_warning("setting yints=None")
        yints = None

    inp_d = inp[:]
    if bgf:
        bg_d = bg[:].T
    else:
        bg_d = np.zeros((inp_d.shape[2], inp_d.shape[1]))

    # Read spacial axes
    (nz, dz, oz, uz, lz) = read_axis(inp, 1)
    (nx, dx, ox, ux, lx) = read_axis(inp, 2)
    zmax = oz + dz*(nz-1)
    xmax = ox + dx*(nx-1)
    inp_d = np.swapaxes(inp_d, 1, 2)

    # Read time axes
    (nt, dt, ot, ut, lt) = read_axis(inp, 3)
    t = np.arange(nt)*dt + ot
    if tmin is None:
        tmin = ot
    if tmax is None:
        tmax = t[-1]
    inp_d = inp_d[np.logical_and(tmin <= t, t <= tmax),:,:]
    inp_d = inp_d[::jsnap,:,:]
    ot = tmin
    dt *= jsnap
    nt = inp_d.shape[0]

    if aclip:
        inp_d /= np.max(np.abs(inp_d))
    else:
        for i in range(nt):
            inp_d[i,:,:] /= np.max(np.abs(inp_d[i,:,:]))

    if bgf:
        (uz2, lz2) = read_axis(bg, 1)[3:]
        (ux2, lx2) = read_axis(bg, 2)[3:]
    else:
        uz2 = lz2 = ux2 = lx2 = ""
    ux = ux if ux else ux2
    lx = lx if lx else lx2
    uz = uz if uz else uz2
    lz = lz if lz else lz2

    mpl.rcParams.update({'font.size': 16})
    dpi = 90

    # Initialize figure and axes
    fig = plt.figure(figsize=(figx, figy), dpi=dpi)
    axi = fig.add_subplot(111, xlim=(ox, xmax), ylim=(zmax, oz))
    if cbar:
        div = make_axes_locatable(axi)
        axc = div.new_horizontal(size="5%", pad=0.15)
        fig.add_axes(axc)
    img_bg = axi.imshow(bg_d, aspect=aspect, extent=(ox, xmax, zmax, oz), cmap=bmap)
    if cbar:
        cb_bg = plt.colorbar(img_bg, cax=axc)
        cb_bg.ax.set_ylabel(cbarl) if cbarl else None

    # Transparent colormap
    # cmap_t = create_tranparent_cmap('white', 'black')
    cmap_t = create_tranparent_cmap(wmap)
    img = axi.imshow(np.zeros((nz, nx), 'f'), aspect=aspect, interpolation='lanczos',
                     extent=(ox, xmax, zmax, oz), vmin=-pclip, vmax=pclip, cmap=cmap_t)

    axi.set_xlabel(r'%s (%s)' % (lx, ux)) if ux else axi.set_xlabel(r'%s' % lx)
    axi.set_ylabel(r'%s (%s)' % (lz, uz)) if uz else axi.set_ylabel(r'%s' % lz)
    if title:
        axi.set_title(title, fontsize=22)
    if xints:
        axi.locator_params(axis='x', nbins=xints)
    if yints:
        axi.locator_params(axis='y', nbins=yints)

    # Time text box
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
    if ttxt:
        time_text = axi.text(0.05, 0.05, '', bbox=bbox_props, transform=axi.transAxes)
    fig.tight_layout()

    if verb:
        msg = ''
        timings = []
        tic = timeit.default_timer()
        tic_total = tic

    def init():
        if ttxt:
            return img, time_text
        else:
            return img,

    # Overlay this at each animation step
    def animate(i):
        if verb:
            global tic
            global msg
            toc = timeit.default_timer()
            sys.stderr.write('\b'*(len(msg)))
            sys.stderr.flush()
        img.set_data(inp_d[i,:,:])
        if ttxt:
            time_text.set_text('t = %.2fs' % (ot+i*dt,) )
        
        if verb:
            timings.append(toc-tic)
            elapsed = sum(timings)/len(timings)
            m, s = divmod(elapsed*(nt-i-1), 60)
            h, m = divmod(m, 60)
            msg = "ETA: %02d:%02d:%02d" % (h, m, s)
            sys.stderr.write(msg)
            sys.stderr.flush()

            tic = timeit.default_timer()
        if ttxt:
            return img, time_text
        else:
            return img,

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=nt, interval=10*dt,
                                   blit=True)

    # Plot or save animations
    if savefile:
        if savefile.endswith('html'):
            from IPython.display import HTML
            html = HTML(anim.to_jshtml(fps=30)).data
            with open(savefile, 'w') as f:
                f.write(html)
        elif savefile.endswith('gif'):
            anim.save(savefile, writer='imagemagick')
        else:
            anim.save(savefile, fps=fps)
        sys.stderr.write('\b'*(len(msg)))
        sys.stderr.flush()
        toc = timeit.default_timer()
        if verb:
            m, s = divmod(toc-tic_total, 60)
            h, m = divmod(m, 60)
            msg = "%02d:%02d:%02d" % (h, m, s)
            sf_warning("Elapsed time - %s" % msg)
    else:
        plt.show()
        sys.stderr.write('\b'*(len(msg)))
        sys.stderr.flush()
    sys.stderr.write('')
