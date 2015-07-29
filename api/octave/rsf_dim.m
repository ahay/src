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
## @deftypefn {Function File} {[@var{dims}, @var{stat}]=} rsf_dim (@var{in})
##
## Returns a vector @var{dims} with the dimensions of RSF header @var{in}
## without trailing dimensions of length 1.
## If no error was encountered, @var{stat}.err = false.
## Else, @var{stat}.err = true and @var{stat}.msg contains the error message
##
## @example
## @code{[n, stat] = rsf_par('junk.rsf');}
## or simply
## @code{n = rsf_par('junk.rsf');}
## @end example
## @seealso{rsf_create, rsf_par, rsf_read, rsf_write}
## @end deftypefn

function [dims, stat] = rsf_dim(in)

    dims = [];
    stat.err = true; % assume we will have an error
    stat.msg = [];

    docstr = '\nUsage: [dims, stat] = rsf_dim(filename)';

    if nargin==0
        stat.msg = ['Please specify filename.' docstr];
        return
    endif

    if nargin>1
        stat.msg = ['Too many arguments.' docstr];
        return
    endif

    if !isstr(in)
        stat.msg = 'Input to rsf_dim must be a string';
        return
    endif

    command = ['sffiledims <' in ' parform=n'];

    [sys_stat, stdout] = system( command );

    if sys_stat!=0
        stat.msg = ['sffiledims failed with code ' num2str(sys_stat)];
        return
    endif

    % sffiledims output looks like: "ndims:n1,n2,n3"
    % We want an array: [n1;n2;n3]
    % could have easily used awk to process it, but who knows
    % if awk is present on all systems

    % put the "n1,n2,n3" part into dstring
    dstring = split( stdout, ':' )(2,:);

    % convert dstring to an array
    dims = str2num( split( dstring, ',' ));

    stat.err = false; % completed without error

endfunction
