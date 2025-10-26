"""
  RsfPy - Python tools for Madagascar RSF data file reading/writing and scientific array handling.

  Copyright (C) 2025 Jilin University

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""



import io, re, warnings, os
from subprocess import Popen, PIPE, SubprocessError, run as Run

def _check_input_source(src, mode='rb'):
    """
    Check if src is a valid readable/writable file path or IOBase object.
    """
    if isinstance(src, str):
        try:
            fp = open(src, mode)
            return fp
        except Exception as e:
            warnings.warn(f"File not accessible with {mode}: {src}, {e}")
            return None
    elif isinstance(src, io.IOBase):
        return src
    else:
        warnings.warn(f"Invalid input type: {type(src)}")
        return None


def _str_match_re(str_in: str, pattern: str = r'\s+(?=(?:[^"]*"[^"]*")*[^"]*$)', strip: str | None = None) -> dict:
    """
    Match strings in the input using the specified regex pattern.
    """
    out_dict = {}
    tokens = re.split(pattern, str_in.strip(strip))
    for token in tokens:
        if "=" not in token:
            continue
        k, v = token.split("=", 1)
        out_dict[k] = v.strip('"').strip("'")
    return out_dict


def _get_datapath(cwd=os.getcwd()):

    top = _datapath()
    path = os.path.dirname(top)
    if top[:2] != './':
        # create a hierarchical structure
        tree = _dirtree(cwd)
        for level in tree:
            if level:
                path = os.path.join(path, level)
        _mkdir(path)
    path = os.path.join(path, os.path.basename(top))
    return path

def _datapath():
    '''
         Path for binary and temporary files.

           Datapath rules:
           1. check datapath= on the command line (Not applicable)
           2. check .datapath file in the current directory
           3. check DATAPATH environmental variable
           4. check .datapath in the home directory
           5. use '.' (not a SEPlib behavior)

       Note that option 5 results in:
       if you move the header to another directory then you must also move
       the data to the same directory. Confusing consequence.

        :return: str
           datapath
    '''
    path = os.environ.get('DATAPATH')
    if not path:
        try:
            pathfile = open('.datapath','r')
        except:
            try:
                pathfile = open(os.path.join(os.environ.get('HOME'),
                                             '.datapath'),'r')
            except:
                pathfile = None
        if pathfile:
            for line in pathfile.readlines():
                check = re.match(r"(?:%s\s+)?datapath=(\S+)" % os.uname()[1],
                                 line)
                if check:
                    path = check.group(1)
            pathfile.close()
    if not path:
        path = './' # the ultimate fallback
    return path

def _dirtree(cwd=None):
    'hierarcical directory structure'
    if not cwd:
        cwd = os.getcwd()
    tree = (os.path.basename(os.path.dirname(os.path.dirname(cwd))),
            os.path.basename(os.path.dirname(cwd)),
            os.path.basename(cwd))
    return tree


def _mkdir(dir):
    'Recursive directory making'
    while os.path.basename(dir) == '.':
        dir = os.path.dirname(dir)
    if dir and not os.path.isdir(dir):
        _mkdir(os.path.dirname(dir))
        os.mkdir(dir)
    return dir


def flow(cmd, source=None, verb=False):
    '''
    RSF Flow
    Usage:
      > sfin = flow("sfmath n1=1 output='1' | sfin")\n
      > print(sfin.read().decode())\n
        in:
            in="stdin"\n
            esize=4 type=float form=native\n
            n1=1           d1=1           o1=0\n
            1 elements 4 bytes\n
    :param source: str | BaseIO | None
        Input data file
    :param cmd: str
    :return: BytesIO
    '''
    out = io.BytesIO()

    if source is None:
        fsrc = None
    elif type(source) == str:
        fsrc = open(source,'rb')
    elif isinstance(source,io.IOBase):
        fsrc = source
    else: raise TypeError("Wrong type of source: %s"%type(source))

    if fsrc is not None:
        fsrc.seek(0)
        subprc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
        [sout, serr] = subprc.communicate(fsrc.read())
        if subprc.returncode != 0:
            raise SubprocessError("In Command: '%s':\n%s"%(cmd,serr.decode()))
    else :
        subprc = Run(cmd, stdin=None, stdout=PIPE, stderr=PIPE, shell=True, check=True)
        sout = subprc.stdout
        serr = subprc.stderr
        if subprc.returncode != 0:
            raise SubprocessError("In Command: '%s':\n"%(cmd,serr.read().decode()))

    out.write(sout)
    out.seek(0)
    return out