% Used by the Octave API

function [stat, dims] = rsf_dim(in)
    % Returns a vector "dims" with the dimensions of a rsf file
    % without trailing dimensions of length 1
    % 
    % Usage: [dims stat] = rsf_dim(filename)

    dims = [];
    stat.err = true; % assume we will have an error
    stat.msg = [];

    docstr = 'Usage: [dims stat] = rsf_dim(filename)';

    if nargin==0
        stat.msg = ['Please specify filename.\n' docstr];
        return
    endif

    if nargin>1
        stat.msg = ['Too many arguments.\n' docstr];
        return
    endif

    if !isstr(in)
        stat.msg = 'Input to rsf_dim must be a string';
        return
    endif

    command = ['sffiledims <' in ' parform=n'];

    [stdout, sys_stat] = system( command );

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
