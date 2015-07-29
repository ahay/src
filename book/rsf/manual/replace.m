function replace(oldtxt,newtxt,file)

% REPLACE will lookinto a file and change all occurrences of a string with another string.
% This function is extremely fragile and untested (and dangerous). Use with extreme caution.
%
% Use \\ to produce a backslash character and %% to produce the percent
% character.
%
% SYNTAX
%   status = replaceText('oldtxt','newtxt','file');
%   status = replaceText('oldtxt','newtxt','folder');
% 
% INPUTS 
%   'oldtxt': string to be replaced
%   'newtxt': string to be changed to. 
%   'folder':   full path to folder, whose entire files will be subject to change. 
% 
% 
% 
% 
% 
% 
% 

[F,loc] = getfilelist(file);


initial_dir = pwd;
eval(sprintf('cd %s',loc));

Lo = length(oldtxt);
% Ln = length(newtxt);


for i = 1:length(F)
    
    ischanged = 0;
    
    % read file
    lines = readfile(F(i).name);
   
    % find and replace
    L = length(lines);
    for l = 1:L   
        idx = strfind(lines{l},oldtxt); % find occurrences
        if ~isempty(idx)
            fprintf('Found one.\n')
            ischanged = 1;
            for ii = length(idx):-1:1
                lines{l}(idx(ii):idx(ii)+Lo-1) = []; % delete oldtxt
                lines{l} = squeezein(lines{l},idx(ii),newtxt);
            end
            tmp = lines{l};
        end
    end
    
    if ischanged, 
        % Make backup if required
        makebackup(F(i).name); 
        % Write to file
        writefile(F(i).name,lines);
    end
    

end   
   
% Go back to oldtxt location
eval(sprintf('cd %s;',initial_dir));
end


%% #### Local functions ####

%% getfilelist
function [F,loc] = getfilelist(file)
%
% OUTPUTS
%   F:   structure array with the file list.
%   loc: absolute path to files.
    
    F = dir(file);
    if isempty(F), 
        error('%s: No such file of directory.',upper(mfilename)); 
    end
    
    
    for i = length(F):-1:1  % get rid of ., .., and other folders
        if F(i).isdir, F(i) = []; end
    end

    if isempty(F)
        fprintf('%s: No such file or directory.\n',upper(mfilename));
        return;
    end
    
    % Find loc
    switch isdir(file)  % the user specified folder
        case 1
            idx = strfind(file,'/');
            if isempty(idx) || ~strcmp(idx(end),'/')
                file = [file,'/'];
%                 idx = [idx length(file)];
            end            
            loc = getabspath(file);
            
        case 0          % the user specified file
           idx = strfind(file,'/');
           if isempty(idx)
               loc = pwd;
           else
               folder = file(1:idx(end));
               loc = getabspath(folder);
           end
    end
    
    
end

%% makebackup 
function makebackup(fname)
    success = copyfile(fname, sprintf('%s.bak',fname));
    if ~success, 
        error('%s: error occured while copying. Exiting.',upper(mfilename));
    end
end

%% readfile   
function lines = readfile(fname)

    fid = fopen(fname,'r');
    
    l = 0;
    while ~feof(fid)
        % read line
        l = l+1;
        lines{l} = fgetl(fid);        %#ok<AGROW>
%         fprintf('l=%d: %s\n',l,lines{l});
        % Care for special characters
        special = ['\' '%'];
        for s = 1:length(special)
            idx = strfind(lines{l},special(s));
            if ~isempty(idx)
                for ii = length(idx):-1:1
%                     fprintf('length of line: %d, position of %s\:%d\n', length(lines{l}),special(s),idx(ii))
                     lines{l} = squeezein(lines{l},idx(ii),special(s));
                end
%                 fprintf('l=%d: %s\n',l,lines{l});
            end
        end
    end
    fclose(fid);
end

%% writefile  
function writefile(fname,lines)
    
    fid = fopen(fname,'w');
    L = length(lines);
    
    for l = 1:L
        fprintf(fid,[lines{l} '\n']);
    end
    fclose(fid);
end

%% getabspath 
function loc = getabspath(folder)
    if strcmp(folder(1),filesep)
        loc = folder;
    else
        loc = [pwd,filesep,folder];
    end
end

%% squeezein
function line = squeezein(line,pos,str)
    if pos==1
        line = [str line];
    else
        line = [line(1:pos-1) str line(pos:end)];
    end
end