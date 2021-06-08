// Chapel-Madagascar API
use m8r;

// Clipping value as config.
// This would allow the cli syntax --clip=
//config var clip: real(32) = 100.0;

proc main(args: [] string)
{
    // Initialize Madagascar
    sf_init(args);

    // Get input and output files
    var fin: sf_file = sf_input("in");
    var fout: sf_file = sf_output("out");

    // Input file must be float
    if SF_FLOAT != sf_gettype(fin) then
        sf_error("Need float input");
    
    // Number of samples
    var n1: int(32);
    if !sf_histint(fin,"n1",n1) then
        sf_error("No n1= in input");    

    // Number of traces
    var n2 = sf_leftsize(fin,1);

    // Get params
    var clip: real(32);
    if !sf_getfloat("clip", clip) then
        sf_error("Need clip=");

    // Trace storage
    var trace: [0..n1-1] real(32);

    // For every trace
    for i2 in 0..n2-1 {
        // Read trace
        sf_floatread(trace, n1, fin);

        for val in trace {
            if val >  clip then val =  clip;
            if val < -clip then val = -clip; 
        }

        // Write trace
        sf_floatwrite(trace, n1, fout);
    }   

    // Close files
    sf_fileclose(fin);
    sf_fileclose(fout);

    // End Madagascar
    sf_close();
} 
