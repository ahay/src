/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/*-----------------------------------------------------------------------
 * segd.h - include file for SEGD traces
 *          Version 2.1, 10/10/94, Marc Schaming
 *          Version 2.2, 17/08/95, Celine Girard et Marc Schaming
 *          Version 2.3,  04/1997, Marc Schaming
 *          Version 2.4,  02/1998, Stewart A. Levin (SEP)
 *                        changed general_header_1.f to unsigned char[2]
 *                        changed non-struct typedefs to struct typedefs
 *                        to allow for uniform future expansion
 *          Version 2.5,  02/2001, Stewart A. Levin (SEP)
 *                        Little endian fixups for channel_set_header.cs
 *                        and dem_trace_header,tn
 *                        

 * declarations for:
 *		typedef struct {} general_header_1 :    general header
 *		typedef struct {} general_header_2 :    general header
 *		typedef struct {} general_header_n :    general header
 *		typedef struct {} gen_head_sn358 :      general header for sn358 
 *		typedef struct {} channel_set_header :  channel set descriptor
 *		typedef struct {} sample_skew :         sample skew block
 *		typedef struct {} extended_header   :   extended header  
 *		typedef struct {} external_header   :   external header  
 *		typedef struct {} general_trailer  :    general trailer  

 *		typedef struct {} dem_trace_header :    demultiplexed trace header
 *		typedef struct {} trace_header_ext :    trace header extension
 *		typedef struct {} general_trailer :     general trailer
 *
 *
 * Reference:
 *	SEG Subcommitee on Digital Tape Formats, "Digital field
 *		tape format standards - SEG-D"
 * Digital field tape format standards - SEG-D, REVISION 1 (1994)
 *    Geophysics, vol. 59, p.668-684
 *	
 *---------------------------------------------------------------------*/

/*---------------- General header block #1 -------------------------*/

typedef struct {
	unsigned char f[2];   /* 01-02 File number (0-9999) */
	unsigned short y;     /* 03-04 Format code */
	char k1_k2;           /* 05    General constants */
	char k3_k4;           /* 06    General constants */
	char k5_k6;           /* 07    General constants */
	char k7_k8;           /* 08    General constants */
	char k9_k10;          /* 09    General constants */
	char k11_k12;         /* 10    General constants */
	unsigned char yr;     /* 11    Year (0-99) */
	unsigned char gh_dy1; /* 12    Number blocks in general header */
	                      /* 12     - day of year (x--) */
	unsigned char dy;     /* 13    Day of year (xx) */
	unsigned char h;      /* 14    Hour of day */
	unsigned char mi;     /* 15    Minute of hour */
	unsigned char se;     /* 16    Second of minute */
	unsigned char m[3];   /* 17    Manufacturer's code */
	                      /* 18-19 and serial number */
	unsigned char b[3];   /* 20-22 Bytes per scan (multiplexed only) */
	unsigned char i;      /* 23    Base scan interval */
	unsigned char p_sbx;  /* 24    Polarity */
	                      /* 24    - Number of scans per block */
	unsigned char sb;     /* 25    Number of scans per block */
	unsigned char z_r1;   /* 26    Record type */
	                      /* 26    - Record length */
	unsigned char r;      /* 27    Record length */
	unsigned char str;    /* 28    Scan types per record  */
	unsigned char cs;     /* 29    Channels sets per scan type */ 
	unsigned char sk;     /* 30    Skew blocks */
	unsigned char ec;     /* 31    Extended header length */
	unsigned char ex;     /* 32    External header length */
} general_header_1;


/*---------------- General header block #2 -------------------------*/

typedef struct {
	unsigned char ef[3]; /* 01-03 Extended file number */
	unsigned char en[2];  /* 04-05 Extended channel sets and scan types */
	unsigned char ecx[2]; /* 06-07 Extended header blocks */
	unsigned char eh[2];  /* 08-09 external header blocks */
	char x1;              /* 10    undefined */
	unsigned char rev[2]; /* 11-12 SEG-D revision number */
	unsigned short gt;    /* 13-14 General trailer number */
	unsigned char erl[3]; /* 15-17 Extended record length */
	char x2;              /* 18    undefined */
	unsigned char bn;     /* 19    General header block number */
	char x3[13];          /* 20-32 undefined */
} general_header_2;


/*---------------- General header block #n -------------------------*/

typedef struct {
	char x1[3];           /* 01-03 undefined */
	unsigned char sln[5]; /* 04-08 Source line number */
	unsigned char spn[5]; /* 09-13 Source point number */
	unsigned char spi;    /* 14    Source point index */
	unsigned char pc;     /* 15    Phase control */
	unsigned char v;      /* 16    Type vibrator */
	short pa;             /* 17-18 Phase angle */
	unsigned char bn;     /* 19    General header block number */
	unsigned char ss;     /* 20    Source set number */
	char x2[12];          /* 21-32 undefined */
} general_header_n;


/*---------------- General header extension (Sercel SN358) ---------*/

typedef struct {
	unsigned char fc1;    /* 01    first and last channel of seismic param 1 */
	unsigned char lc1;    /* 01-02 last channel of seismic param 1 */
	unsigned char fc2;    /* 03-04 first channel of seismic param 2 */
	unsigned char f_lc2;  /* 03-04 first and last channel of seismic param 2 */
	unsigned char lc2;    /* 04-05 last channel of seismic param 2 */
	unsigned char fc3;    /* 06-07 first channel of seismic param 3 */
	unsigned char f_lc3;  /* 07    fist and last channel of seismic param 3 */
	unsigned char lc3;    /* 07-08 last channel of seismic param 3 */
	unsigned char fc4;    /* 09-10 first channel of seismic param 4 */
	unsigned char f_lc4;  /* 10    first and last channel of seismic param 4 */
	unsigned char lc4;    /* 10-11 fist channel of seismic param 4 */
	unsigned char f_lac1; /* 12    first and last aux channel of scan type 1 */
	unsigned char fsc1;   /* 13-14 first seismic channel of scan type 1 */
	unsigned char f_lsc1; /* 14    first and last seismic channel of scan type 1 */
	unsigned char lsc1;   /* 14-15 last seismic channel of scan type 1 */
	unsigned char sam_int1;/* 16    sample interval of scan type 1 */
	unsigned char fac2;   /* 17    first and last aux channel of scan type 2 */
	unsigned char fsc2;   /* 18-19 first seismic channel of scan type 2 */
	unsigned char f_lsc2; /* 19    first and last seismic channel of scan type 2 */
	unsigned char lsc2;   /* 19-20 last seismic channel of scan type 2 */
	unsigned char sam_int2;/* 21    sample interval of scan type 2 */
	unsigned char bl_sig_le;/* 22    block signature length (n*0.1 s) */
	unsigned char rec_length[2];  /* 23-24 record length n*0.1 s (00.0-99.9) s */	
	unsigned char dyn_swit_del[2];   /* 25-26 dynamically switching delay (00.0-99.9 s) */
	unsigned char rec_del[2];       /* 27-28 recording delay (00.0-99.9s) */
	unsigned char ty_a_cha12;/* 29    type of auxiliary channel 1, channel 2 */
	unsigned char ty_a_cha34;/* 30    type of auxiliary channel 3, channel 4 */
	unsigned char ty_a_cha56;/* 31    type of auxiliary channel 5, channel 6 */
	unsigned char ty_a_cha78;/* 32    type of auxiliary channel 7, channel 8 */
	unsigned char mode_num;/* 33    mode number */ 
	unsigned char an_sys_co;/* 34    analog sys count (1-2), tape transport num */
	unsigned char reel_num[2];       /* 35-36 reel number (0-9999) */
	unsigned char file_num[2];       /* 37-38 file logical number */
	unsigned char sp_num[2];         /* 39-40 shot point number */ 
	unsigned char lc_nf;  /* 41    lc=0 lowcut fi. out, nf=0 notch fi. out */ 
	unsigned char fg_ic;  /* 42    fg=1 fist fixed gain, ic=1 internal osc. */
	unsigned char of1;    /* 43-45 oscillator frequency (0-9999.9 Hz) */ 
	unsigned char of2;    /* 43-45 oscillator frequency (0-9999.9 Hz) */ 
	unsigned char of3;    /* 43-45 oscillator frequency (0-9999.9 Hz) */ 
	unsigned char osc_att;/* 46    oscillator attenuator */ 
	unsigned char t_sig_ph;/* 47    test signal phase (1:+,0:-) */
	unsigned char m_gain; /* 48    main gain amplifier (0-15) 16: IFP */
	unsigned char dum[16];/* 49-64 additional information */ 
} gen_head_sn358;


/*---------------- Scan type header (Channel set descriptor ) ------*/

typedef struct {
	unsigned char st;     /* 01    Scan type */
	unsigned char cn;     /* 02    Channel set number */
	unsigned short tf;    /* 03-04 Channel set start time */
	unsigned short te;    /* 05-06 Channel set end time */
	unsigned char mp[2];  /* 07-08 Descaling exponent */
	unsigned char cs[2];  /* 09-10 Channels in this channel set */
	unsigned char c;      /* 11    Channel type identification */
	unsigned char sc_j;   /* 12    Sample/channel gain */
	                      /* 12    - Gain control method */
	unsigned short af;    /* 13-14 Alias filter frequency */
	unsigned short as;    /* 15-16 Alias filter slope */
	unsigned short lc;    /* 17-18 Low cut filter frequency */
	unsigned short ls;    /* 19-20 Low cut filter slope */
	unsigned short nt[3]; /* 21-26 Notch filter frequency */
	unsigned short ecs;   /* 27-28 Extended channel set number */
	unsigned char efh;    /* 29    Extended header flag */
	unsigned char vs;     /* 30    Vertical stack */
	unsigned char cab;    /* 31    Streamer number */
	unsigned char ary;    /* 32    Array forming */
} channel_set_header;


/*---------------- Sample skew -------------------------------------*/
          
typedef struct { 
	unsigned char skew[32]; 
} sample_skew;


/*---------------- Extended header ---------------------------------*/

typedef struct {
	unsigned char dummy[32];
} extended_header;


/*---------------- External header ---------------------------------*/

typedef struct {
	unsigned char dummy[32];
} external_header;


/*---------------- General trailer ---------------------------------*/

typedef struct {
	unsigned short gt;    /* 01-02 General trailer number */
	char x1[8];           /* 03-10 undefined */
	unsigned char c;      /* 11    Channel type identification */
	char x2[21];          /* 12-32 undefined */
} general_trailer;


/*--------------- Demultiplexed trace header -----------------------*/

typedef struct {
	unsigned short f;     /* 01-02 file number */
	unsigned char st;     /* 03    scan type (1-2) */
	unsigned char cn;     /* 04    channel set */
	unsigned char tn[2];    /* 05-06 trace number */
	unsigned char t[3];   /* 07-09 timing word of the first sample if
                                the data were written in the multiplexed
                                format */
	unsigned char the;    /* 10    trace header extensions */
	unsigned char ss;     /* 11    sample skew of the first sample
                                of the trace. It is a part of
                                the fractional part of the base scan interval. Res :
                                1/256 scan interval */
	unsigned char tr;     /* 12    trace edit */
	unsigned char tw[3];  /* 13-15 time from time break to the end of 
                                the internal time break window (binary 
                                number, inc 1ms) */
	unsigned char en[2];  /* 16-17 extended channel set number */
	unsigned char efn[3]; /* 18-20 extended file number */
} dem_trace_header;


/*---------------- Trace header extension --------------------------*/

typedef struct {
	unsigned char rln[3]; /* 01-03 Receiver line number */
	unsigned char rpn[3]; /* 04-06 Receiver point number */
	unsigned char rpi;    /* 07    Receiver point index */
	unsigned char nbs[3]; /* 08-10 Number of samples per traces */
	char x[22];           /* 11-32 undefined */
} trace_header_ext;
