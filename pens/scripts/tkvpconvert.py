#!/usr/bin/env python
'''GUI for vpconvert'''

##   Copyright (C) 2012 University of Texas at Austin
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

import sys, os

try:
    import Tkinter as tk
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror, showinfo
except:
    sys.stderr.write('Please install Tkinter!\n\n')
    sys.exit(1)
      
import vpconvert

def select_file(entry):
    '''Select a file into entry'''
    filename = askopenfilename()
    if not filename:
        return

    entry.delete(0, tk.END)
    entry.insert(0, filename)

def convert(vpl,format,opts,pars):
    '''Convert a VPL file'''
    if not os.path.isfile(vpl):
        showerror("ERROR", "Can't find " + vpl)
        return

    new = '.'.join([os.path.splitext(vpl)[0],format.lower()])
    opts = ' '.join(map(lambda x: '='.join([x, str(pars[x].get())]), 
                        pars.keys())) + ' ' + opts
    
    fail = vpconvert.convert(vpl,new,format,None,opts,False)
    
    run = "%s to %s using \"%s\"" % (vpl,new,opts)
    if fail:
        showerror("Could not convert",run)
    else:
        showinfo("Converted",run)

def main():
    root = tk.Tk()
    root.title("Vplot Converter")

    # Vplot File: _________________ [...]
    frame = tk.Frame(root)
    tk.Label(frame, text="Vplot File:").pack(side=tk.LEFT)
    vpl = tk.Entry(frame, width=60)
    vpl.pack(side=tk.LEFT)
    tk.Button(frame, text="...",
              command=lambda: select_file(vpl)).pack(side=tk.LEFT)
    frame.pack()

    pars = {}

    # Format
    fmt = tk.StringVar()
    fmt.set("png") # default value

    frame = tk.Frame(root)
    tk.Label(frame, text="Format:").pack(side=tk.LEFT)
    subframe = tk.Frame(frame)
    formats = vpconvert.pens.keys()
    formats.sort()
    nf = len(formats)
    for i in range(nf):
        format = formats[i]
        rb = tk.Radiobutton(subframe,text=format.upper(),value=format,variable=fmt)
        rb.grid(row=i%2,column=i/2,sticky=tk.W)
    subframe.pack(side=tk.LEFT)
    frame.pack(fill=tk.X,pady=10)

    pars['format'] = fmt

    # Fat
    fat = tk.IntVar()
    fat.set(1)


    frame = tk.Frame(root)
    
    tk.Label(frame, text="Fat:").pack(side=tk.LEFT)
    scale = tk.Scale(frame,orient=tk.HORIZONTAL,variable=fat,from_=1,to=10,length=330)
    scale.pack(side=tk.LEFT)

    frame.pack(fill=tk.X,pady=10)

    pars['fat'] = fat

    # bgcolor

    bgcolor = tk.StringVar()
    bgcolor.set("light")

    frame = tk.Frame(root)
    tk.Label(frame, text="Background Color:").pack(side=tk.LEFT)
    tk.Radiobutton(frame,text="Light",value="light",variable=bgcolor).pack(side=tk.LEFT)
    tk.Radiobutton(frame,text="Dark",value="dark",bg="gray50",fg="white",variable=bgcolor).pack(side=tk.LEFT)
    tk.Radiobutton(frame,text="Black",value="black",bg="black",fg="white",variable=bgcolor).pack(side=tk.LEFT)
    tk.Radiobutton(frame,text="White",value="white",bg="white",variable=bgcolor).pack(side=tk.LEFT)

    pars['bgcolor'] = bgcolor

    # Serifs
    serifs = tk.IntVar()
    serifs.set(1)

    tk.Checkbutton(frame,text="Serifs",variable=serifs,padx=10).pack(side=tk.RIGHT)

    frame.pack(fill=tk.X,pady=10)

    pars['serifs'] = serifs

    # Other options
    options = tk.StringVar()
    options.set("")

    frame = tk.Frame(root)
    tk.Label(frame, text="Other Options:").pack(side=tk.LEFT)
    tk.Entry(frame,textvariable=options,width=60).pack(side=tk.LEFT)
    frame.pack(fill=tk.X,pady=10)

    # [Convert] [Quit]
    frame = tk.Frame(root)
    run = tk.Button(frame,text="Convert",background="yellow",
                    command=lambda: convert(vpl.get(),
                                            fmt.get(),
                                            options.get(),
                                            pars))
    run.pack(side=tk.LEFT)
    tk.Button(frame,text="Quit",command=root.quit).pack(side=tk.RIGHT)
    frame.pack(fill=tk.X,pady=10)

    root.mainloop()

if __name__ == "__main__":
    main()

