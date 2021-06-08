// Time-domain acoustic FD modeling

// Chapel-Madagascar API
use m8r;

proc main(args: [] string)
{
    // Laplacian coefficients
    const c0 = -30.0/12: real(32);
    const c1 =  16.0/12: real(32);
    const c2 = -1.0/12:  real(32);

    sf_init(args); // init RSF
    var verb: bool; // verbose flag
    if( !sf_getbool("verb", verb) ) then verb = false;

    // Setup I/O files
    var Fw: sf_file = sf_input("in");
    var Fv: sf_file = sf_input("vel");
    var Fr: sf_file = sf_input("ref");
    var Fo: sf_file = sf_output("out");

    // Read/Write axes
    var nt:  int(32); sf_histint  ( Fw, "n1", nt ); 
    var dt: real(32); sf_histfloat( Fw, "d1", dt );
    var nz:  int(32); sf_histint  ( Fv, "n1", nz );
    var dz: real(32); sf_histfloat( Fv, "d1", dz );
    var nx:  int(32); sf_histint  ( Fv, "n2", nx );
    var dx: real(32); sf_histfloat( Fv, "d2", dx );

    sf_putint  ( Fo, "n1", nz );
    sf_putfloat( Fo, "d1", dz );
    sf_putint  ( Fo, "n2", nx );
    sf_putfloat( Fo, "d2", dx );
    sf_putint  ( Fo, "n3", nt );
    sf_putfloat( Fo, "d3", dt );

    var dt2: real(32) = dt**2;
    var idz: real(32) = 1.0/(dz**2);
    var idx: real(32) = 1.0/(dx**2);

    // Read wavelet, velocity and reflectivity
    var ww: [0..nt-1]       real(32); sf_floatread( ww, nt, Fw );
    var vv: [0..<nx,0..<nz] real(32); sf_floatread( vv, nz*nx, Fv );
    var rr: [0..<nx,0..<nz] real(32); sf_floatread( rr, nz*nx, Fr );

    // Allocate temporary arrays
    var um: [0..<nz,0..<nx] real(32);
    var uo: [0..<nz,0..<nx] real(32);
    var up: [0..<nz,0..<nx] real(32);
    var ud: [0..<nz,0..<nx] real(32);

    // Main loop
    for it in 0..nt-1 {

        // 4th order laplacian
        forall (ix,iz) in {2..nx-3,2..nz-3} {
            ud[ix,iz] = c0* uo[ix,  iz  ]                 *(idx + idz) +
                        c1*(uo[ix,  iz-1] + uo[ix,  iz+1])*idx +
                        c2*(uo[ix,  iz-2] + uo[ix,  iz+2])*idx +
                        c1*(uo[ix-1,iz  ] + uo[ix+1,iz  ])*idz +
                        c2*(uo[ix-2,iz  ] + uo[ix+2,iz  ])*idz;  
        }

	// Inject wavelet
        forall (ix,iz) in {0..<nx,0..<nz} {
            ud[ix,iz] -= ww[it] * rr[ix,iz];
        }

        // Scale by velocity
        forall (ix,iz) in {0..<nx,0..<nz} {
            ud[ix,iz] *= vv[ix,iz]**2;
        }

        // Time step
        forall (ix,iz) in {0..<nx,0..<nz} {
            up[ix,iz] = 2.0*uo[ix,iz] -
                            um[ix,iz] +
                            ud[ix,iz] * dt2;                        
        }

	// Interchange wavefiels
	um = uo;
        uo = up;

	// Write wavefield to output
        sf_floatwrite( uo, nz*nx, Fo );
    }

    // Close files
    sf_fileclose( Fw );
    sf_fileclose( Fv );
    sf_fileclose( Fr );
    sf_fileclose( Fo );

    // End Madagascar
    sf_close();
} 
