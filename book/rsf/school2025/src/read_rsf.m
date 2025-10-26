function [data, header] = read_rsf(file)
%READ_RSF  Read RSF header and data (binary or ASCII).
%
%   [DATA, HEADER] = READ_RSF(FILE) parses FILE up to the RSFHSPLITER
%   (0x0C0C04) as text header, collects key=val pairs into HEADER
%   (last occurrence wins), converts certain HEADER fields to numeric,
%   then reads the data portion according to HEADER.in and HEADER.data_format.
%
%   HEADER fields:
%     .in           – 'stdin' or external filename
%     .n1, .n2, …   – integer dimensions
%     .o1, .o2, …   – float origins
%     .d1, .d2, …   – float deltas
%     .esize        – integer element size
%     .data_format  – e.g. 'native_float', 'xdr_complex', 'ascii_int'
%
%   DATA is returned as an N-D array sized [n1,n2,…].

    % 1) Define the RSF header/data splitter
    splitter = uint8([12;12;4]);  % RSFHSPLITER = 0x0C0C04

    % 2) Read entire file as bytes
    fid = fopen(file,'rb');
    if fid < 0
        error('Cannot open file: %s', file);
    end
    allBytes = fread(fid, Inf, 'uint8=>uint8');
    fclose(fid);

    % 3) Split header vs. raw bytes
    hdrStr    = char(allBytes');
    splitStr  = char(splitter');
    idx       = strfind(hdrStr, splitStr);
    % idx = strfind(allBytes', splitter');
    if ~isempty(idx)
        hdrBytes  = allBytes(1:idx(1)-1);
        dataBytes = allBytes(idx(1)+numel(splitter):end);
    else
        hdrBytes  = allBytes;
        dataBytes = [];
    end

    % 4) Parse header text into key/value pairs
    hdrText = char(hdrBytes');
    tokens  = regexp(hdrText,'\s+','split');
    tokens(cellfun('isempty',tokens)) = [];
    header = struct();
    for k = 1:numel(tokens)
        t  = tokens{k};
        eq = strfind(t,'=');
        if isempty(eq), continue; end
        key = t(1:eq(1)-1);
        val = t(eq(1)+1:end);
        % strip surrounding quotes
        if numel(val)>=2 && ((val(1)=='''' && val(end)=='''') || ...
                             (val(1)=='"'  && val(end)=='"' ))
            val = val(2:end-1);
        end
        header.(key) = val;  % last one wins
    end

    % 5) Convert certain header fields to numeric
    fnames = fieldnames(header);
    for i = 1:numel(fnames)
        f = fnames{i};
        v = header.(f);
        % n#: dimensions → int32
        if regexp(f,'^n[1-9]$','once')
            header.(f) = int32(str2double(v));
        % o#, d#: origins/deltas → double
        elseif regexp(f,'^[od][1-9]$','once')
            header.(f) = str2double(v);
        % esize → int32
        elseif strcmp(f,'esize')
            header.(f) = int32(str2double(v));
        % other keys remain as strings
        end
    end

    % 6) Determine data source: 'stdin' or external file
    if ~isfield(header,'in')
        error('Header field "in" not found.');
    end
    src = header.in;
    if strcmp(src,'stdin')
        rawBytes = dataBytes;
    else
        fid2 = fopen(src,'rb');
        if fid2 < 0
            error('Cannot open binary file: %s', src);
        end
        rawBytes = fread(fid2, Inf, 'uint8=>uint8');
        fclose(fid2);
    end

    % 7) Gather dimensions n1,n2,…
    dims = [];
    idx = 1;
    while true
        fld = sprintf('n%d',idx);
        if isfield(header,fld)
            dims(end+1) = double(header.(fld)); %#ok<AGROW>
            idx = idx + 1;
        else
            break;
        end
    end
    if isempty(dims)
        error('No dimension fields (n1,n2,…) found.');
    end
    nElem = prod(dims);

    % 8) Parse data_format
    if ~isfield(header,'data_format')
        error('Header field "data_format" not found.');
    end
    % parts = split(header.data_format,'_');
    parts = strsplit(header.data_format, '_');

    if numel(parts)~=2
        error('Invalid data_format: %s', header.data_format);
    end
    order = parts{1};   % native, xdr, ascii
    btype = parts{2};   % int, float, complex

    % 9) Handle ASCII formats
    if strcmp(order,'ascii')
        txt    = char(rawBytes');
        toks   = regexp(txt,'\s+','split');
        toks(cellfun('isempty',toks)) = [];
        switch btype
          case 'int'
            if numel(toks)<nElem
                error('Not enough ASCII tokens for int data.');
            end
            vec = int32(str2double(toks(1:nElem)));
          case 'float'
            if numel(toks)<nElem
                error('Not enough ASCII tokens for float data.');
            end
            vec = single(str2double(toks(1:nElem)));
          case 'complex'
            if numel(toks)<2*nElem
                error('Not enough ASCII tokens for complex data.');
            end
            vec = complex(zeros(nElem,1,'single'));
            for k = 1:nElem
                re = str2double(toks{2*k-1});
                imt = toks{2*k};  % e.g. '3.14i'
                if imt(end)=='i'
                    imt = imt(1:end-1);
                end
                im = str2double(imt);
                vec(k) = complex(single(re), single(im));
            end
          otherwise
            error('Unknown ASCII type: %s', btype);
        end
        data = reshape(vec, dims);
        return
    end

    % 10) Handle binary formats
    switch btype
      case 'int'
        prec     = 'int32'; bytesPer = 4;
      case 'float'
        prec     = 'single'; bytesPer = 4;
      case 'complex'
        prec     = 'single'; bytesPer = 8;
      otherwise
        error('Unknown binary type: %s', btype);
    end

    if numel(rawBytes) < nElem*bytesPer
        error('Binary data too short: need %d bytes, got %d.', ...
              nElem*bytesPer, numel(rawBytes));
    end
    rawBytes = rawBytes(1:nElem*bytesPer);

    if strcmp(btype,'complex')
        tmp = typecast(rawBytes,'single');
        if strcmp(order,'xdr')
            tmp = swapbytes(tmp);
        end
        re  = tmp(1:2:end);
        im  = tmp(2:2:end);
        vec = complex(re,im);
    else
        vec = typecast(rawBytes, prec);
        if strcmp(order,'xdr')
            vec = swapbytes(vec);
        end
    end

    % 11) Reshape to [n1,n2,…]
    data = reshape(vec, dims);
end
