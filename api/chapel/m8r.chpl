/* RSF-Chapel interface */
/*
  Copyright (C) 2021 Jorge Monsegny
 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

module m8r {

    require "rsf.h";

    use SysBasic;
    use CTypes;

    //Extern declaration of C Madagascar procedures
    //and datatypes w=wrapper nw=no wrapper

    //Functions for Madagascar initialization and closure
    extern proc sf_init(argc: c_int, argv: [] c_ptr(c_char)): void; //w
    extern proc sf_close() : void; //nw

    //Functions for opening and closing files
    extern type sf_file;
    extern proc sf_input(tag : c_string) : sf_file; //w
    extern proc sf_output(tag : c_string) : sf_file; //w
    extern proc sf_fileclose(file : sf_file) : void; //nw

    //Functions to get/set file type
    extern type sf_datatype = c_int;
    extern const SF_UCHAR :sf_datatype;
    extern const SF_CHAR :sf_datatype;
    extern const SF_INT :sf_datatype;
    extern const SF_FLOAT :sf_datatype;
    extern const SF_COMPLEX :sf_datatype;
    extern const SF_SHORT :sf_datatype;
    extern const SF_DOUBLE :sf_datatype;
    extern const SF_LONG :sf_datatype;
    extern proc sf_gettype(file : sf_file) : sf_datatype; //nw
    extern proc sf_settype(file : sf_file, type_arg : sf_datatype) : void; //nw

    //Functions to get/set file form
    extern type sf_dataform = c_int;
    extern const SF_ASCII :sf_dataform;
    extern const SF_XDR :sf_dataform;
    extern const SF_NATIVE :sf_dataform;
    extern proc sf_getform(file : sf_file) : sf_dataform; //nw
    extern proc sf_setform(file : sf_file, form : sf_dataform) : void; //nw
    extern proc sf_setformat(file : sf_file, format : c_string) : void; //nw

    //Terminal functions
    extern proc sf_warning(format : c_string) : void; //w
    extern proc sf_error(format : c_string) : void; //w

    //Functions to get keys from files 
    extern proc sf_histint(file : sf_file, 
                           key : c_string, 
                           ref par : int(32)) : bool; //nw
    extern proc sf_histints(file : sf_file,
                            key : c_string,
                            par : c_ptr(c_int),
                            n : int(32)) : bool; //w
    extern proc sf_histlargeint(file : sf_file,
                                key : c_string,
                                ref par : int(64)) : bool; //nw
    extern proc sf_histfloat(file : sf_file, 
                             key : c_string, 
                             ref par : real(32)) : bool; //nw
    extern proc sf_histfloats(file : sf_file,
                              key : c_string,
                              par : c_ptr(c_float),
                              n : int(32)) : bool; //w
    extern proc sf_histbool(file : sf_file, 
                            key : c_string, 
                            ref par : bool) : bool; //nw
    extern proc sf_histbools(file : sf_file,
                             key : c_string,
                             par : c_ptr(bool),
                             n : int(32)) : bool; //w
    extern proc sf_histdouble(file : sf_file, 
                              key : c_string, 
                              ref par : real(64)) : bool; //nw
    extern proc sf_histstring(file : sf_file,
                              key : c_string) : c_string; //nw

    //Functions to put keys in files
    extern proc sf_putint(file : sf_file, 
                          key : c_string, 
                          par : c_int) : void; //nw
    extern proc sf_putlargeint(file : sf_file, 
                               key : c_string, 
                               par : c_long) : void; //nw
    extern proc sf_putints(file : sf_file, key : c_string, 
                           par : [] c_int, n : c_int) : void; //w
    extern proc sf_putfloat(file : sf_file, 
                            key : c_string, 
                            par : c_float) : void; //nw
    extern proc sf_putstring(file : sf_file, 
                             key : c_string, 
                             par : c_string) : void; //nw
    extern proc sf_putline(file : sf_file, line : c_string) : void; //nw

    //Size functions
    extern const SF_MAX_DIM : int;
    extern proc sf_leftsize(file : sf_file, dim : int(32)) : int(32); //nw
    extern proc sf_filedims(file : sf_file, n : c_ptr(c_int)) : int(32); //w
    extern proc sf_largefiledims(file : sf_file, 
                                 n : c_ptr(c_long)) : int(32); //w 
    extern proc sf_filesize(file : sf_file) : c_int; //nw
    extern proc sf_esize(file : sf_file) : c_int; //nw

    //Functions that read command line parameters
    extern proc sf_getint(key : c_string, ref par : int(32)) : bool; //nw
    extern proc sf_getlargeint(key : c_string, ref par : int(64)) : bool; //nw
    extern proc sf_getints(key : c_string, par : c_ptr(c_int), n : int(32)) : bool; //w
    extern proc sf_getfloat(key : c_string, ref par : real(32)) : bool; //nw
    extern proc sf_getdouble(key : c_string, ref par : real(64)) : bool; //nw
    extern proc sf_getfloats(key : c_string, par : c_ptr(c_float), n : int(32)) : bool; //w
    extern proc sf_getstring(key : c_string) : c_string; //nw
    extern proc sf_getstrings(key : c_string, par : [] c_ptr(c_char), n : int(32)) : bool; //w looks like the strings are colon separated: par=str1:str2:...
    extern proc sf_getbool(key : c_string, ref par : bool) : bool; //nw
    extern proc sf_getbools(key : c_string, par : c_ptr(bool), n : int(32)) : bool; //w

    //Complex external type
    extern type sf_complex;

    //Read data routines.  
    extern proc sf_complexread(arr : c_ptr(sf_complex), 
                               size : c_int, 
                               file : sf_file) : void; //w  
    extern proc sf_floatread(arr : c_ptr(c_float), 
                             size : c_int, 
                             file : sf_file) : void; //w
    extern proc sf_intread(arr : c_ptr(c_int), 
                           size : c_int, 
                           file : sf_file) : void; //w
    extern proc sf_shortread(arr : c_ptr(c_short),
                             size : c_int,
                             file : sf_file) : void; //w
    extern proc sf_charread(arr : c_ptr(c_char),
                            size : c_int,
                            file : sf_file) : void; //w
    extern proc sf_uncharread(arr : c_ptr(c_uchar),
                              size : c_int,
                              file : sf_file) : void; //w

    //Data writing routines
    //sf_complexwrite is not implemented
    extern proc sf_complexwrite(arr : c_ptr(sf_complex),
                                size : c_int,
                                file : sf_file) : void; //w
    extern proc sf_floatwrite(arr : c_ptr(c_float), 
                              size : c_int, 
                              file : sf_file) : void; //w
    extern proc sf_intwrite(arr : c_ptr(c_int), 
                            size : c_int, 
                            file : sf_file) : void; //w
    extern proc sf_shortwrite(arr : c_ptr(c_short),
                              size : c_int,
                              file : sf_file) : void; //w
    extern proc sf_charwrite(arr : c_ptr(c_char),
                             size : c_int,
                             file : sf_file) : void; //w
    extern proc sf_uncharwrite(arr : c_ptr(c_uchar),
                               size : c_int,
                               file : sf_file) : void; //w

    //Some procedures wrappers

    //Functions for Madagascar initialization and closure
    proc sf_init(args: [] string) : void
    {
        var argv: [0..args.size-1] c_ptr(c_char);        
        for i in (0..args.size-1) {
            var tmp: c_string = args[i].localize().c_str();            
            argv[i] = tmp:c_ptr(c_char);          
        }
        sf_init(args.size:c_int, argv);
    }

    //Functions for opening and closing files
    proc sf_input(tag: string): sf_file
    {
        return sf_input(tag.c_str());
    }

    proc sf_output(tag: string): sf_file
    {   
        return sf_output(tag.c_str());
    }

    //Terminal functions
    proc sf_warning(args...?n) : void
    {   
        var msg: string;
        for param i in 0..n-1 {//param must be here because args is heterogenous
            msg += args(i):string;
        }
        sf_warning(msg.c_str());
    }

    proc sf_error(args...?n) : void
    {
        var msg: string;
        for param i in 0..n-1 {//param must be here because args is heterogenous
            msg += args(i):string;
        }
        sf_error(msg.c_str());
    }

    //Size functions
    proc sf_filedims(file : sf_file, n : [] int(32)) : int
    {
        var cn = c_ptrTo(n);
        return sf_filedims(file, cn);        
    }

    proc sf_largefiledims(file : sf_file, n : [] int(64)) : int(32)
    {
        var cn = c_ptrTo(n);
        return sf_largefiledims(file, cn);
    }

    //Functions to get keys from files
    proc sf_histints(file : sf_file, key : string, 
                     par : [] int(32), n : int(32)) : bool
    {
        var cpar = c_ptrTo(par);
        return sf_histints(file, key.c_str(), cpar, n);
    }

    proc sf_histfloats(file : sf_file, key : string, 
                       par : [] real(32), n : int(32)) : bool
    {
        var cpar = c_ptrTo(par);
        return sf_histfloats(file, key.c_str(), cpar, n);
    }

    proc sf_histbools(file : sf_file, key : string, 
                      par : [] bool, n : int(32)) : bool
    {
        var cpar = c_ptrTo(par);
        return sf_histbools(file, key.c_str(), cpar, n);
    }


    //Functions that read command line parameters
    proc sf_getints(key : string, par : [] int(32), n : int(32)) : bool
    {
        var cpar = c_ptrTo(par);
        return sf_getints(key.c_str(), cpar, n);
    }

    proc sf_getfloats(key : string, par : [] real(32), n : int(32)) : bool
    {
        var cpar = c_ptrTo(par);
        return sf_getfloats(key.c_str(), cpar, n);
    }

    proc sf_getstrings(key : string, par : [] string, n : int(32)) : bool
    {
        var cpar: [0..n-1] c_ptr(c_char);
        var ret = sf_getstrings(key.c_str(), cpar, n);
        for i in 0..n-1 {
            par[i] = cpar[i]:c_string;
        }
        return ret;
    }

    proc sf_getbools(key : string, par : [] bool, n : int(32)) : bool
    {
        var cpar = c_ptrTo(par);
        return sf_getbools(key.c_str(), cpar, n);
    }

    //Data reading functions
    proc sf_complexread(arr : [] complex(64), size : int, 
                        file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_complexread(carr:c_ptr(sf_complex), csize, file);
    }

    proc sf_floatread(arr : [] real(32), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_floatread(carr, csize, file);
    }

    proc sf_intread(arr : [] int(32), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_intread(carr, csize, file);
    }

    proc sf_shortread(arr : [] int(16), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_shortread(carr, csize, file);
    }

    proc sf_charread(arr : [] int(8), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_charread(carr, csize, file);
    }

    proc sf_uncharread(arr : [] uint(8), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_ucharread(carr, csize, file);
    }

    //Data writing functions
    proc sf_complexwrite(arr : [] complex(64), size : int, 
                         file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_complexwrite(carr:c_ptr(sf_complex), csize, file);
    }

    proc sf_floatwrite(arr : [] real(32), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_floatwrite(carr, csize, file);
    }
  
    proc sf_intwrite(arr : [] int(32), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_intwrite(carr, csize, file);
    }

    proc sf_shortwrite(arr : [] int(16), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_shortwrite(carr, csize, file);
    }
   
    proc sf_charwrite(arr : [] int(8), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_charwrite(carr, csize, file);
    }

    proc sf_uncharwrite(arr : [] uint(8), size : int, file : sf_file) : void
    {
        var carr = c_ptrTo(arr);
        var csize: c_int = size: int(32);
        sf_uncharwrite(carr, csize, file);
    }
}
