function clip(in,out,clip)
%CLIP Clip the data
% CLIP(in,out,clip) reads data from RSF file "inp", clips it, and writes it to RSF file "out"
% Processing is done trace by trace to save memory 

path(path,'/home/fomels/RSF/filt/lib/mex');

dims = rsf_dim(in);
n1 = dims(1);           % trace length
n2 = prod(dims(2:end)); % number of traces
trace = 1:n1;           % allocate trace
rsf_create(out,in)      % create an output file

for i2 = 1:n2           % loop over traces
    i2
    rsf_read(trace,in,'same');
    trace
    trace(trace >   clip) =  clip;
    trace(trace < - clip) = -clip;
    trace
    rsf_write(trace,out,'same');
end
