struct device dev = {
    0,0,0,0,
    0.0,
    0.0,
    -1,  /* num_color */
    1,   /* lost */
    false, /* need_end_erase */
    false, /* smart_clip */
    false, /* smart_raster */
    false, /* smart_background */
    false, /* cachepipe */

    DEFAULT_FONT, /* txfont */
    DEFAULT_PREC, /* txprec */
    OVLY_NORMAL,  /* txovly */

    0,0,   /* xorigin, yorigin */

    BREAK_BREAK, /* brake */

    /* control routines */
    opendev,		/* open */
    nulldev,		/* reset */
    genmessage,		/* message */
    nullclose,		/* erase */
    nullclose,		/* close */
    
    /* high level output */
    genvector,		/* vector */
    genmarker,		/* marker */
    gentext,		/* text */
    genarea,		/* area */
    genraster,		/* raster */
    genpoint,		/* point */
    nullattributes,	/* attributes */
    
    /* input */
    gen_dovplot,            /* reader */
    nullgetpoint,	    /* getpoint */
    nullinteract,	    /* interact */
    
    /* low level output */
    nullplot,		/* plot */
    nullclose,		/* startpoly */
    nullmidpoly,		/* midpoly */
    nullclose		/* endpoly */
};

struct s_txalign txalign;
