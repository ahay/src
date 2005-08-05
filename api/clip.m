function clip(in,out,clip)
%CLIP Clip the data

dims = rsf_dim(in);
n1 = dims(1);           % trace length
n2 = prod(dims(2:end)); % number of traces
trace = 1:n1;           % allocate trace
rsf_create(out,in)      % create an output file

for i2 = 1:n2           % loop over traces
    rsf_read(trace,in,'same');
    trace(trace >   clip) =  clip;
    trace(trace < - clip) = -clip;
    rsf_write(trace,out,'same');
end

