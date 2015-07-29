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
## @deftypefn {Function File} {@var{stat}=} rsf_create (@var{file}, @var{arg2})
##
## Writes RSF header with desired info to disk.
## Returns @var{default} if @var{name} is not present in @var{file}
## If @var{arg2} is a filename, it copies that header
## If @var{arg2} is a vector with dimensions, creates a header with those 
## dimensions. If success, @var{stat}.err = false.
## Else, @var{stat}.err = true and @var{stat}.msg contains the error message.
##
## @example
## @code{[n1, stat] = rsf_par('junk.rsf','n1',1);}
## or simply
## @code{n1 = rsf_par('junk.rsf','n1',1);}
## @end example
## @seealso{rsf_dim, rsf_par, rsf_read, rsf_write}
## @end deftypefn


function stat = rsf_create( out_filename, arg2 )
    % Writes RSF header with desired info to disk
    %
    % Usage: stat = rsf_create( out_filename, inp_filename | dims )

    docstr = '\nUsage: stat = rsf_create( out_filename, inp_filename | dims )';

    stat.err = true;  % assume we will have an error
    stat.msg = [];

    if nargin<2
        stat.msg = ['Too few arguments.' docstr];
        return
    endif

    if nargin>2
        stat.msg = ['Too many arguments.' docstr];
        return
    endif

    if !isstr(out_filename)
        stat.msg = ['First argument must be a filename string'];
        return
    endif

    if isstr(arg2) % copy the header

        command = ['sffileflush <' arg2 ' >' ];

        [dummy_var, sys_stat] = system( command );

        if sys_stat!=0
            stat.msg = ['sffileflush failed with code ' num2str(sys_stat)];
        else
            stat.err = false; % completed without error
        endif

    elseif size(size(arg2),2)==2 & (size(arg2,1)==1 | size(arg2,2)==1)

        % create header from arg2, which was just tested to be a vector
        % (2-D array with one of the sizes of length 1)

        command = 'sfcreate';
        for i=1:length(arg2)
            command = [command ' n' num2str(i) '=' arg2(i)];
        endfor
        command = [command '>' out_filename];

        [dummy_var, sys_stat] = system( command );

        if sys_stat!=0
            stat.msg = ['sfcreate failed with code ' num2str(sys_stat)];
        else
            stat.err = false; % completed without error
        endif

    else
        stat.msg = ['Second argument must be either string or vector'];
    endif

endfunction
