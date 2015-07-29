## Copyright (C) 2007 Ioan Vlad
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{par}, @var{stat}]=} rsf_par (@var{file}, @var{name}, @var{default})
##
## Reads parameter @var{name} from RSF header @var{file}
## Returns @var{default} if @var{name} is not present in @var{file}
## If no error was encountered, @var{stat}.err = false.
## Else, @var{stat}.err = true and @var{stat}.msg contains the error message.
##
## @example
## @code{[n1, stat] = rsf_par('junk.rsf','n1',1);}
## or simply
## @code{n1 = rsf_par('junk.rsf','n1',1);}
## @end example
## @seealso{rsf_create, rsf_dim, rsf_read, rsf_write}
## @end deftypefn

function [par, stat] = rsf_par(file, name, default)

    par = [];
    stat.err = true; % assume we will have an error
    stat.msg = [];

    docstr = '\nUsage: [par, stat] = rsf_par(file, name, default)';

    if nargin<3
        stat.msg = ['Too few arguments.' docstr];
        return
    endif

    if nargin>3
        stat.msg = ['Too many arguments.' docstr];
        return
    endif

    if !isstr(file)
        stat.msg = ['Argument 1 of rsf_dim must be a string.' docstr];
        return
    endif

    if !isstr(name)
        stat.msg = ['Argument 2 of rsf_dim must be a string.' docstr];
        return
    endif

    [sys_stat, stdout] = system(['sfget <' file ' ' name ' parform=n']);

    if sys_stat!=0
        stat.msg = ['sfget failed with code ' num2str(sys_stat)];
        return
    else
        stat.err = false; % completed without error
        if strncmp(stdout, 'sfget: No key ', 14)
            par = default;
        else
            par = str2num(stdout);
        endif
    endif

endfunction
