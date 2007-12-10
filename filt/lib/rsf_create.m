% Used by the Octave API

function stat = rsf_create( out_filename, arg2 )
    % Writes RSF header with desired info to disk
    %
    % Usage: stat = rsf_create( out_filename, inp_filename | dims )

    docstr = 'Usage: error = rsf_create( out_filename, inp_filename | dims )';

    stat.err = true;  % assume we will have an error
    stat.msg = [];

    if nargin<2
        stat.msg = ['Too few arguments.\n' docstr];
        return
    endif

    if nargin>2
        stat.msg = ['Too many arguments.\n' docstr];
        return
    endif

    if !isstr(out_filename)
        stat.msg = ['First argument must be a filename string'];
        return
    endif

    if isstr(arg2) % copy the header

        command = ['sffileflush <' arg2 ' >' ];

        [dummy_var sys_stat] = system( command );

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

        [dummy_var sys_stat] = system( command );

        if sys_stat!=0
            stat.msg = ['sfcreate failed with code ' num2str(sys_stat)];
        else
            stat.err = false; % completed without error
        endif

    else
        stat.msg = ['Second argument must be either string or vector'];
    endif

endfunction
