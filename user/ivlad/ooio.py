#! /usr/bin/env python
"""
NAME
	ooio
DESCRIPTION
	Object-Oriented I/O
SOURCE
	user/ivlad/ooio.py
"""
# Copyright (C) 2010 Ioan Vlad
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os, sys, socket, datetime, time, getpass, struct
import rsf.path
import rsf.version as version

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsf.user.ivlad as ivlad

################################################################################

stdout = 'stdout'

RSF_hdr_stream_end = 2*'\014' + '\004'

dtype_avl = ('int','float') # Acceptable Value List
################################################################################

class BaseFile:
    'Base class for both standalone files and metafiles'

    def set_intent(self, intent):
        # Can be used both in init, and to change intent later
        assert intent in ('in', 'out')
        self.intent = intent

    def __init__(self, intent):
        self.bny_attribs = [] # Unprintable object attribs (i.e. SEG-Y reel hdr)
        self.set_intent(intent)

    def print_self(self, varname):
        attribs = self.__dict__
        akeys = list(attribs.keys())
        akeys.sort()
        max_nm_len = max(list(map(len,akeys)))
        ivlad.msg(ivlad.hr)
        ivlad.msg(varname + ': instance of ' + self.__class__.__name__)
        indent = '  '
        spc = ' '
        for key in akeys:
            dots_to_print = max_nm_len - len(key) + 2
            line_to_print = indent + key + spc + dots_to_print * '.' + spc
            if key in self.bny_attribs:
                line_to_print += '<' + str(len(attribs[key])) + \
        	                    ' binary bytes, not printed>'
            else:
                val = attribs[key]
                if type(val) == list and len(val) == 1:
                    val=val[0] # Better readability without brackets
            line_to_print += str(val)
            ivlad.msg(line_to_print)

################################################################################

class File(BaseFile):
    'Single file, header, binary or together, i.e. segy, .rsf, .rsf@, etc'

    def path_analysis(self,path):
        self.full_nm = os.path.abspath(path) # /abs/path/file.rsf
        self.parent_dir = os.path.abspath(os.path.dirname(path)) # /abs/path
        self.nm = os.path.basename(path) # file.rsf
        [self.nm_noext, self.ext] = os.path.splitext(self.nm) # [file, .rsf]
        self.writing_to_stdout = (self.nm == stdout and self.ext == '')

    def chk_ok(self, expected_sz=None):
        if self.intent == 'in':
            # Check file exists and is readable
            pass
            if expected_sz != None:
                # check if file has expected size
                pass
        elif self.intent == 'out':
            # Check if parent directory is writable
            ivlad.chk_dir_wx(self.parent_dir)
            if expected_sz != None:
                # Check there is space on disk
                pass

    def __init__(self, path, intent, expected_sz=None):
        BaseFile.__init__(self, intent)
        if path != None:
            self.path_analysis(path)
            self.chk_ok(expected_sz)
        if intent == 'in':
            self.sz = os.path.getsize(self.full_nm)
        elif intent == 'out':
            self.sz = None

    def open_new(self):
        # Opening is not done automatically upon init because we need to
        # handle stdout (first hdr, then dat)
        if self.writing_to_stdout:
            self.handle = sys.stdout
        else:
            self.handle = open(self.full_nm, 'w', encoding='utf-8')

    def close(self):
        if not self.writing_to_stdout: # stdout does not need to be closed
            self.handle.close()
        self.handle = None

################################################################################

class MetaFile(BaseFile):
    'A collection of files with logical relationship between them'

    def __init__(self, intent, fmt, valid_ext):
        BaseFile.__init__(self, intent)
        self.total_sz = None # Total size in bytes for all files
        self.all_exist = None # Whether all files listed exist
        self.open_files  = {}
        assert type(fmt) == str
        self.fmt = fmt # File format name: RSF, ISF, etc
        if type(valid_ext) == str:
            self.valid_ext = [valid_ext]
        elif type(valid_ext) == list:
            self.valid_ext = valid_ext

    def close_all(self):
        for filenm in list(self.open_files.keys()):
            open_files[filenm].close()

################################################################################

def first_line_rsf_hdr_entry():
    'Returns the first line of a program entry in the history file'
    # Defined outside the class because called directly from outside

    vers = version.label
    prog = sys.argv[0]
    cwd  = os.getcwd()
    user = getpass.getuser()
    host = socket.gethostname()
    date = datetime.date.today()
    now = time.strftime('%a %b %d %X %Y')
    return '%s %s %s: %s@%s %s\n\n' % (vers, prog, cwd, user, host, now)

################################################################################

class RSFheader(File):

    def __init__(self, path, intent, par=None):
        File.__init__(self, path, intent)
        # For reproducibility, it would be nice to know the command-line
        # arguments, and this works with rsf.apibak, but I need to know
        # how to make it work with rsf.api as well.
        # if intent == 'out':
            # self.cli_par = par._Par__args
        # To write the CLI args to header, have
        # self.write_dict(self.cli_par, 'From CLI')
        # in RSFheader.open()

    def write_1st_line(self):
        'Writes to history file the first line of a program entry'

        self.handle.write(first_line_rsf_hdr_entry())

    def add2hist(self, item, val, comment=None):
        'Write variable to a RSF history file'
        newline = '\n'
        indent = 4 * ' '
        to_print = indent
        if comment != None:
            to_print +=  '# ' + comment + newline + indent
        to_print += item + '='
        # To do: properly handle comments containing newline characters
        # i.e. replace '\n' with '\n# '
        if type(val) in (str, str):
            if newline in val: # multiline string
                quot="'''"
                if quot in val:
                    quot='"""'
                    if quot in val:
                        pass # replace it in val
                        # using escaped characters, i.e. \'\'\'
                to_print += quot + newline
                if val[-1] != newline:
                    to_print += newline
                to_print += quot
            else:
                quot='"'
                if quot in val:
                    quot="'"
                    if quot in val:
                        pass # Use \'
                to_print += quot+val+quot
        elif type(val) == datetime.date:
            quot='"'
            to_print += quot + str(val) + quot
        else:
            to_print += str(val)
        self.handle.write( to_print + newline )

    def write_dict(self, adict, common_comment=None):
        # adict must be the dictionary of arguments passed to the program
        keys = list(adict.keys())
        keys.sort()
        for key in keys:
            keycomment = common_comment
            if type(adict[key]) == dict:
                val = adict[key]['val']
                if common_comment == None:
                    keycomment = adict[key]['com']
                else:
                    keycomment += adict[key]['com']
            else:
                val = adict[key]
            self.add2hist( key, val, keycomment )

    def open(self):
        self.open_new()
        self.write_1st_line()
        # self.write_dict(self.cli_par, 'From CLI')
        self.add2hist('in',self.dat)

################################################################################

class RSFfile(MetaFile):
    '''Regularly Sampled Format metafile (header+data).
    See http://ahay.org/wiki/RSF_Comprehensive_Description'''

    def set_defaults_and_constants(self):
        self.ndim_max = 9
        self.defaults = {'n':1,'o':0.0,'d':'1.0','unit':'unknown',
        'label':'unknown'}
        self.hdr_flushed = False
        # To do: implement XDR and ASCII encodings
        self.encoding = 'native'
        self.esize = 4

    def chk_io_type(self):
        self.writing_to_stdout = self.hdr.writing_to_stdout
        if self.intent == 'out':
            if self.writing_to_stdout:
                self.writing_to_stdout_pipe = ivlad.stdout_is_pipe()
                self.writing_to_stdout_file = not self.writing_to_stdout_pipe
            else:
                self.writing_to_stdout_pipe = False
                self.writing_to_stdout_file = False

    def init_dat(self,dpath_suffix=''):
        assert self.intent == 'out'
        if self.writing_to_stdout_pipe:
            self.hdr.dat = 'stdin'
            self.dat = File(None, intent='out')
            self.dat.writing_to_stdout = True
        else:
            # We will write a data file to disk. Get or create its name:
            if self.writing_to_stdout_file:
                dat_nm = ivlad.data_file_nm()
            elif not self.writing_to_stdout:
                dat_nm = self.hdr.nm
            dat_nm = dpath_suffix + dat_nm + '@'
            self.dpath = os.path.abspath(rsf.path.datapath())
            ivlad.chk_dir_wx(self.dpath)
            self.hdr.dat = os.path.join(self.dpath, dat_nm)
            self.dat = File(self.hdr.dat, intent='out')

    def open_hdr(self):
        self.hdr.open()
        self.hdr.add2hist('esize',self.esize)
        self.hdr.add2hist('data_format',self.encoding + '_' + self.dtype)
        if not self.hdr.writing_to_stdout:
            self.open_files[self.hdr.full_nm] = self.hdr.handle

    def flush_hdr(self):
        'Writes everything in the header, opens data file for writing'
        all_ax_dict = {}
        for i in range(self.ndim):
            for k in list(self.ax[i].keys()):
                all_ax_dict[k+str(i+1)] = self.ax[i][k]
        self.hdr.write_dict(all_ax_dict)
        if self.writing_to_stdout_pipe:
            self.hdr.handle.write(RSF_hdr_stream_end)
        elif not self.writing_to_stdout:
            self.hdr.close()
            del self.open_files[self.hdr.full_nm]
        self.dat.open_new()
        self.hdr_flushed = True

    def set_ndim(self, ndim):
        if ndim > 0 and ndim < self.ndim_max:
            self.ndim = ndim

    def __init__(self, hdr_nm, par, ndim, intent, dtype='float',
                                                            dpath_suffix=''):
        MetaFile.__init__(self, intent, fmt='RSF', valid_ext = '.rsf')
        self.set_defaults_and_constants()
        assert dtype in dtype_avl
        self.dtype = dtype
        self.ax = []
        for i in range(ndim):
            self.ax.append(
                {'n':None,'o':None,'d':None, 'unit':None, 'label':None})
        self.set_ndim(ndim)
        self.set_intent(intent)
        self.hdr = RSFheader(hdr_nm, intent, par)
        self.chk_io_type()
        if intent == 'out':
            self.init_dat(dpath_suffix)
            self.open_hdr()

        # Note to self: add a parent_file argument
        # When parent_file is a RSF file, history should be copied
	    # When parent file is a non-RSF file, just the name should be kept

    def __list2dict(self, ilist, dictkey):
        if ilist != None:
            ilist_tmp = ivlad.trunc_or_append(self.ndim, ilist,
                                              self.defaults[dictkey])
            for i in range(self.ndim):
                self.ax[i][dictkey] = ilist_tmp[i]

    def set_hdr_info(self, n=None, o=None, d=None, unit=None, lbl=None):
        assert self.intent == 'out'
        self.__list2dict(n, 'n')
        self.__list2dict(o, 'o')
        self.__list2dict(d, 'd')
        self.__list2dict(unit, 'unit')
        self.__list2dict(lbl, 'label')

    def write(self,val):
        if not self.hdr_flushed:
            self.flush_hdr()
        if sys.version_info[0] >= 3:
            self.dat.handle.write(struct.pack(ivlad.fmt[self.dtype],val).decode('latin1'))
        else:
            self.dat.handle.write(struct.pack(ivlad.fmt[self.dtype],val))

