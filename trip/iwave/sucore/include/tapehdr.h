/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/* tapehdr.h - include file for SEGY traces as bytes (only for segyread,write)
 *
 * Reference:
 *      K. M. Barry, D. A. Cavers and C. W. Kneale, "Special Report:
 *              Recommended Standards for Digital Tape Formats",
 *              Geophysics, vol. 40, no. 2 (April 1975), P. 344-352.
 *      
 * $Author: john $
 * $Source: /usr/local/cwp/src/su/include/RCS/tapehdr.h,v $
 * $Revision: 1.1 $ ; $Date: 1996/09/09 16:18:41 $
 */ 

#ifndef TAPEHDR_H
#define TAPEHDR_H

static struct {
        char *key;      char *type;     int offs;
} tapehdr[] = {
           {"tracl",             "P",            0},
           {"tracr",             "P",            4},
            {"fldr",             "P",            8},
           {"tracf",             "P",            12},
              {"ep",             "P",            16},
             {"cdp",             "P",            20},
            {"cdpt",             "P",            24},
            {"trid",             "U",            28},
             {"nvs",             "U",            30},
             {"nhs",             "U",            32},
            {"duse",             "U",            34},
          {"offset",             "P",            36},
           {"gelev",             "P",            40},
           {"selev",             "P",            44},
          {"sdepth",             "P",            48},
            {"gdel",             "P",            52},
            {"sdel",             "P",            56},
           {"swdep",             "P",            60},
           {"gwdep",             "P",            64},
          {"scalel",             "U",            68},
          {"scalco",             "U",            70},
              {"sx",             "P",            72},
              {"sy",             "P",            76},
              {"gx",             "P",            80},
              {"gy",             "P",            84},
          {"counit",             "U",            88},
           {"wevel",             "U",            90},
          {"swevel",             "U",            92},
             {"sut",             "U",            94},
             {"gut",             "U",            96},
           {"sstat",             "U",            98},
           {"gstat",             "U",            100},
           {"tstat",             "U",            102},
            {"laga",             "U",            104},
            {"lagb",             "U",            106},
           {"delrt",             "U",            108},
            {"muts",             "U",            110},
            {"mute",             "U",            112},
              {"ns",             "U",            114},
              {"dt",             "U",            116},
            {"gain",             "U",            118},
             {"igc",             "U",            120},
             {"igi",             "U",            122},
            {"corr",             "U",            124},
             {"sfs",             "U",            126},
             {"sfe",             "U",            128},
            {"slen",             "U",            130},
            {"styp",             "U",            132},
            {"stas",             "U",            134},
            {"stae",             "U",            136},
           {"tatyp",             "U",            138},
           {"afilf",             "U",            140},
           {"afils",             "U",            142},
          {"nofilf",             "U",            144},
          {"nofils",             "U",            146},
             {"lcf",             "U",            148},
             {"hcf",             "U",            150},
             {"lcs",             "U",            152},
             {"hcs",             "U",            154},
            {"year",             "U",            156},
             {"day",             "U",            158},
            {"hour",             "U",            160},
          {"minute",             "U",            162},
             {"sec",             "U",            164},
          {"timbas",             "U",            166},
            {"trwf",             "U",            168},
          {"grnors",             "U",            170},
          {"grnofr",             "U",            172},
          {"grnlof",             "U",            174},
            {"gaps",             "U",            176},
           {"otrav",             "U",            178},
};

static struct {
        char *key;      char *type;     int offs;
} tapebhdr[] = {
           {"jobid",             "P",            0},
           {"lino",              "P",            4},
           {"reno",              "P",            8},
           {"ntrpr",             "U",            12},
           {"nart",              "U",            14},
           {"hdt",               "U",            16},
           {"dto",               "U",            18},
           {"hns",               "U",            20},
           {"nso",               "U",            22},
           {"format",            "U",            24},
           {"fold",              "U",            26},
           {"tsort",             "U",            28},
           {"vscode",            "U",            30},
           {"hsfs",              "U",            32},
           {"hsfe",              "U",            34},
           {"hslen",             "U",            36},
           {"hstyp",             "U",            38},
           {"schn",              "U",            40},
           {"hstas",             "U",            42},
           {"hstae",             "U",            44},
           {"htatyp",            "U",            46},
           {"hcorr",             "U",            48},
           {"bgrcv",             "U",            50},
           {"rcvm",              "U",            52},
           {"mfeet",             "U",            54},
           {"polyt",             "U",            56},
           {"vpol",              "U",            58},
};

#endif /* end TAPEHDR_H */
