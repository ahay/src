/* vplot filter for LibGD. */
/*
  Copyright (C) 2009 The University of Texas at Austin
  
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

 #include <gd.h>

#include "../include/attrcom.h"
#include "../include/extern.h"
#include "../include/params.h"
#include "../include/err.h"
#include "../include/erasecom.h"
#include "../include/closestat.h"
#include "../include/enum.h"

#include "../genlib/genpen.h" 
#include "../utilities/util.h"

#include "dovplot.h"
#include "init_vplot.h"

#ifdef FFMPEG
#include <avcodec.h>

#ifdef FFMS_USE_FFMPEG_COMPAT
#  define VERSION_CHECK(LIB, cmp, u1, u2, u3, major, minor, micro) ((LIB) cmp (AV_VERSION_INT(major, minor, micro)))
#else
#  define VERSION_CHECK(LIB, cmp, major, minor, micro, u1, u2, u3) ((LIB) cmp (AV_VERSION_INT(major, minor, micro)))
#endif

#endif

#include "gdpen.h"

char            name[] = "gdpen";
#include "gddoc.h"

#include "_device.h"

#define PPI 100.0  /* pixels per inch */
#define NCOLOR 256 /* number of colors */
#define MAXVERT 1000

static gdImagePtr image;
static bool light = false;
static int color_table[NCOLOR], gdcolor, bgcolor, delay, nx, ny;
static const char *image_type;

#ifdef FFMPEG
static AVCodec *codec = NULL;
static AVCodecContext *codec_ctx = NULL;
static AVFrame *mpeg_frame = NULL;
#if LIBAVCODEC_VERSION_MAJOR >= 54
static AVPacket mpeg_pkt;
static int mpeg_gout;
#else
static int frame_out_size;
#endif
static int frame_size, frame_outbuf_size;
static uint8_t *frame_outbuf, *frame_buf;
static int bitrate;
static int frame_num = 0;
static int frame_mark = 0; /* Current time mark */
static int frame_delay = 4; /* 25 fps */

static void ffmpeg_init (void);
static void ffmpeg_finish (void);
static void ffmpeg_write (void);
#endif

void opendev (int argc, char* argv[])
/*< open >*/
{
    char newpath[60];
    const char *color;
    int value;

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    /* control routines */
    dev.reset = gdreset;
    dev.erase = gderase;
    dev.close = gdclose;

    dev.area = gdarea;
    dev.attributes = gdattr;
/*    dev.vector = gdvector; */
    dev.plot = gdplot;
    dev.point = gdpoint;

    dev.pixels_per_inch = PPI;
    nx = VP_STANDARD_HEIGHT * dev.pixels_per_inch;
    ny = VP_SCREEN_RATIO * nx;
    dev.aspect_ratio = 1.;

    sf_getfloat ("aspect", &dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &dev.pixels_per_inch);
    /* pixels per inch */
    nx *= dev.pixels_per_inch / PPI;
    ny *= (1.0 / dev.aspect_ratio) *
	(dev.pixels_per_inch / PPI);
    sf_getint ("n1", &nx);
    sf_getint ("n2", &ny);
    /* image size */

    dev.need_end_erase = true;
    dev.smart_clip = true; 
    dev.num_col = NCOLOR;

    dev.xmax = nx-1;
    dev.ymax = ny-1;

    image = gdImageCreate(nx,ny);

    if (NULL == (color = sf_getstring("bgcolor"))) color="black";
    /* background color (black,white,light,dark) 
       'light' and 'dark' cause the background to be transparent (in PNG and GIF) */

    switch (color[0]) {
	case 'b': /* black */
	case 'd': /* dark */
	default:
	    light = false;
	    bgcolor = gdImageColorAllocate(image, 0, 0, 0);
	    break;
	case 'w': /* white */
	case 'l': /* light */
	    light = true;
	    bgcolor = gdImageColorAllocate(image, 255, 255, 255);
	    break;
    }

    if (color[0]=='l' || color[0]=='d')
	gdImageColorTransparent(image,bgcolor);

    if (isatty(fileno(pltout)))
    {
	sprintf (newpath, "%s", "image_file");
	pltout = fopen (newpath, "wb");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }

    for (value = 0; value < NCOLOR; value++)
    {
	color_table[value] = -1;
    }

    if (NULL == (image_type = sf_getstring("type"))) 
	image_type = "png";
    /* image type (png, jpeg, gif, mpeg) */

    if (!sf_getint("delay",&delay)) delay=10;
    /* animation delay (if type=="gif" or "mpeg") */

#ifndef FFMPEG
    if (image_type[0] == 'm' || image_type[0] == 'M') {
        ERR (FATAL, name, "This pen has not been compiled with MPEG support.\n");
    }
#else
    if (!sf_getint ("bitrate", &bitrate)) bitrate=400000;
    /* MPEG bitrate */
    if (image_type[0] == 'm' || image_type[0] == 'M')
        ffmpeg_init ();
    if (delay < frame_delay)
        delay = frame_delay; /* Can't have more than 25 fps */
#endif
}

void gdreset (void)
/*< reset >*/
{
    int value, color, r, g, b;

    /* reset color table */
    color = color_table[0];
    if (-1 != color) gdImageColorDeallocate(image, color);
    color_table[0] = bgcolor;

    for (value = 1; value < 8; value++)
    {
	color = color_table[value];
	if (-1 != color) gdImageColorDeallocate(image, color);
	r = MAX_GUN * ((value & 2) / 2);
	g = MAX_GUN * ((value & 4) / 4);
	b = MAX_GUN * ((value & 1) / 1);

	if (light) {
	    color_table[value] = gdImageColorAllocate(image,255-r,255-g,255-b);
	} else {
	    color_table[value] = gdImageColorAllocate(image,r,g,b);
	}
    }
 
    for (value = 8; value < NCOLOR; value++)
    {
	color = color_table[value];
	if (-1 != color) {
	    gdImageColorDeallocate(image, color);
	    color_table[value] = -1;
	}
    }
}

static void gd_write (void)
{
    static bool called=false;

#ifdef GIFANIM
    static gdImagePtr oldimage;
#endif

    /* Reset clipping, otherwise GD will write a partial frame */
    gdImageSetClip(image, dev.xmin, dev.ymin, dev.xmax - 1, dev.ymax - 1);

    switch (image_type[0]) {
	case 'j':
	case 'J':
	    if(called) return;
	    gdImageJpeg(image, pltout, -1);
	    break;
#ifdef GIFANIM
	case 'g':
	case 'G':
	    if (called) {
		gdImageGifAnimAdd(image, pltout, 0, 0, 0, delay, 1, oldimage);
                gdImageDestroy (oldimage);
	    } else {
		gdImageGifAnimBegin(image, pltout, 1, 0);
		gdImageGifAnimAdd(image, pltout, 0, 0, 0, delay, 1, NULL);
	    }
	    oldimage = image;
	    image = gdImageCreate(nx,ny);
	    gdImagePaletteCopy(image, oldimage);
	    break;
#endif
#ifdef FFMPEG
	case 'm':
	case 'M':
            frame_num++;
            /* Since MPEG has a very limited set of frame rates,
               we create a lot of dummy frames to "slow" it down */
            while (frame_mark < frame_num * delay) {
                ffmpeg_write ();
                frame_mark += frame_delay;
                if ((frame_mark % 500) == 0)
                    ERR (FATAL, name, "Encoding video: finished %d seconds.\n", frame_mark / 100);
            }
	    break;
#endif
	case 'p':
	case 'P':
	default:
	    if(called) return;
	    gdImagePng(image, pltout);
	    break;
    }

    called=true;
}

void gderase (int command)
/*< erase >*/
{
    switch (command)
    {
	case ERASE_START:
	    gdImageFilledRectangle(image, 0, 0, nx, ny, color_table[0]);
	    break;
	case ERASE_MIDDLE:
	    gd_write();
	    gdImageFilledRectangle(image, 0, 0, nx, ny, color_table[0]);
	    break;
	case ERASE_END:
	    gd_write();
#ifdef GIFANIM
	    if (image_type[0]=='g' || image_type[0]=='G') 
		gdImageGifAnimEnd(pltout);
#endif
#ifdef FFMPEG
	    if (image_type[0] == 'm' || image_type[0] == 'M')
                ffmpeg_finish ();
#endif
	    fclose (pltout);
	    break;
	case ERASE_BREAK:
	    gd_write();
	    break;
	default:
	    break;
    }
}

void gdclose (int status)
/*< close >*/
{
    switch (status)
    {
	case CLOSE_NORMAL:
	    break;
	case CLOSE_FLUSH:
	    fflush (pltout);
	    break;
	default:
	    break;
    }
}

void gdattr (int command, int value, int v1, int v2, int v3)
/*< attr >*/
{
    int xmin, ymin, xmax, ymax;

    switch (command)
    {
	case SET_COLOR:
	    gdcolor = color_table[value];
	    break;
	case SET_COLOR_TABLE:
	    color_table[value] = gdImageColorAllocate(image, v1, v2, v3);
	    break;
	case SET_WINDOW:
	    xmin = value;
	    ymin = dev.ymax - v3;
	    xmax = v2;
	    ymax = dev.ymax - v1;
 	    gdImageSetClip(image, xmin, ymin, xmax, ymax);
	    break;
	default:
	    break;
    }
}

void gdplot(int x, int y, int draw)
/*< plot >*/
{
    static int oldx = 0, oldy = 0;

    if (draw)
    {
	gdImageLine(image, oldx, dev.ymax - oldy, x, dev.ymax - y, gdcolor);
    } else {
	dev.lost = 0;
    }

    oldx = x;
    oldy = y;
}

void gdpoint(int x, int y)
/*< point >*/
{
    gdImageSetPixel(image, x, dev.ymax-y, gdcolor);
}

void gdarea (int npts, struct vertex *head)
/*< area >*/
{
    int i;
    gdPoint vlist[MAXVERT];

    /* translate data structures */
    for (i = 0; i < npts && i < MAXVERT; i++)
    {
	vlist[i].x = head->x;
	vlist[i].y = dev.ymax - head->y;
	head = head->next;
    }

    gdImageFilledPolygon(image, vlist, npts, gdcolor);
}

#ifdef FFMPEG
/* Open and initialize codec + buffers */
static void ffmpeg_init (void) {
    /* must be called before using avcodec lib */
#if LIBAVCODEC_VERSION_MAJOR < 54
    avcodec_init();
#endif
#if LIBAVCODEC_VERSION_INT < 0x363b64 /* 54.59.100 @ ffmpeg-1.0 */
#define AV_CODEC_ID_MPEG1VIDEO CODEC_ID_MPEG1VIDEO 
#endif
#if LIBAVCODEC_VERSION_MAJOR < 58
    /*register all the codecs */
    avcodec_register_all ();
#endif

    /* find the mpeg1 video encoder */
    codec = avcodec_find_encoder (AV_CODEC_ID_MPEG1VIDEO);
    if (!codec) {
        ERR (FATAL, name, "Could not initialize MPEG1 codec\n");
    }

#if LIBAVCODEC_VERSION_MAJOR < 54
    codec_ctx = avcodec_alloc_context ();
#else
    codec_ctx = avcodec_alloc_context3 (codec);
#endif
#if LIBAVCODEC_VERSION_INT < AV_VERSION_INT(55,28,1)
#define av_frame_alloc  avcodec_alloc_frame
#endif

#if LIBAVUTIL_VERSION_MAJOR < 52
#define AV_PIX_FMT_YUV420P PIX_FMT_YUV420P
#endif

    mpeg_frame = av_frame_alloc(); 

    codec_ctx->bit_rate = bitrate;
    /* resolution must be a multiple of two */
    codec_ctx->width = (nx / 2) * 2;
    codec_ctx->height = (ny / 2) * 2;
    /* frames per second */
    codec_ctx->time_base= (AVRational){1,100/frame_delay};
    codec_ctx->gop_size = 10; /* number of frames in a group */
    codec_ctx->has_b_frames = 0;
    codec_ctx->max_b_frames = 0;
    codec_ctx->pix_fmt = AV_PIX_FMT_YUV420P;

    /* open it */
#if LIBAVCODEC_VERSION_MAJOR < 54
    if (avcodec_open (codec_ctx, codec) < 0) {
#else
    if (avcodec_open2 (codec_ctx, codec, NULL) < 0) {
#endif
        ERR (FATAL, name, "Could not open MPEG1 codec\n");
    }

    /* alloc image and output buffer */
    frame_outbuf_size = 1000000;
    frame_outbuf = malloc (frame_outbuf_size);
    frame_size = codec_ctx->width * codec_ctx->height;
    frame_buf = malloc ((frame_size * 3) / 2); /* size for YUV 420 */

    mpeg_frame->data[0] = frame_buf;
    mpeg_frame->data[1] = mpeg_frame->data[0] + frame_size;
    mpeg_frame->data[2] = mpeg_frame->data[1] + frame_size / 4;
    mpeg_frame->linesize[0] = codec_ctx->width;
    mpeg_frame->linesize[1] = codec_ctx->width / 2;
    mpeg_frame->linesize[2] = codec_ctx->width / 2;
}

#define SCALEBITS 10
#define ONE_HALF  (1 << (SCALEBITS - 1))
#define FIX(x)    ((int) ((x) * (1<<SCALEBITS) + 0.5))
/* Convert RGB to YUP420P and do encoding */
static void ffmpeg_write (void) {
    int wrap, x, y, index;
    int r, g, b, r1, g1, b1;
    uint8_t *lum, *cb, *cr;

    lum = mpeg_frame->data[0];
    cb = mpeg_frame->data[1];
    cr = mpeg_frame->data[2];

    wrap = codec_ctx->width;
    for (y = 0; y < codec_ctx->height; y += 2) {
        for (x = 0; x < codec_ctx->width; x += 2) {
            index = gdImageGetPixel (image, x, y);
            r = gdImageRed (image, index);
            g = gdImageGreen (image, index);
            b = gdImageBlue (image, index);
            r1 = r;
            g1 = g;
            b1 = b;
            lum[0] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            index = gdImageGetPixel (image, x + 1, y);
            r = gdImageRed (image, index);
            g = gdImageGreen (image, index);
            b = gdImageBlue (image, index);
            r1 += r;
            g1 += g;
            b1 += b;
            lum[1] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            lum += wrap;

            index = gdImageGetPixel (image, x, y + 1);
            r = gdImageRed (image, index);
            g = gdImageGreen (image, index);
            b = gdImageBlue (image, index);
            r1 += r;
            g1 += g;
            b1 += b;
            lum[0] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            index = gdImageGetPixel (image, x + 1, y + 1);
            r = gdImageRed (image, index);
            g = gdImageGreen (image, index);
            b = gdImageBlue (image, index);
            r1 += r;
            g1 += g;
            b1 += b;
            lum[1] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            cb[0] = ((- FIX(0.16874) * r1 - FIX(0.33126) * g1 +
                      FIX(0.50000) * b1 + 4 * ONE_HALF - 1) >> (SCALEBITS + 2)) + 128;
            cr[0] = ((FIX(0.50000) * r1 - FIX(0.41869) * g1 -
                     FIX(0.08131) * b1 + 4 * ONE_HALF - 1) >> (SCALEBITS + 2)) + 128;
            cb++;
            cr++;
            lum += -wrap + 2;
        }
        lum += wrap;
    }

   

    /* encode the image */
#if LIBAVCODEC_VERSION_MAJOR < 54
    frame_out_size = avcodec_encode_video (codec_ctx, frame_outbuf,
                                           frame_outbuf_size, mpeg_frame);
    fwrite (frame_outbuf, 1, frame_out_size, pltout);
#else
    av_init_packet (&mpeg_pkt);
    mpeg_pkt.data = NULL;
    mpeg_pkt.size = 0;
#if LIBAVCODEC_VERSION_MAJOR < 58
    if (avcodec_encode_video2 (codec_ctx, &mpeg_pkt, mpeg_frame, &mpeg_gout) < 0)
        ERR (FATAL, name, "MPEG encoding error\n");
#else
    if (avcodec_send_frame(codec_ctx, mpeg_frame) < 0)
	ERR (FATAL, name, "MPEG send frame error\n");
    mpeg_gout = avcodec_receive_packet(codec_ctx, &mpeg_pkt);
    if (mpeg_gout == AVERROR(EAGAIN) || mpeg_gout == AVERROR_EOF ||
	mpeg_gout < 0)
	ERR (FATAL, name, "MPEG receive packet error\n");
#endif    
    if (mpeg_gout) {
        fwrite (mpeg_pkt.data, 1, mpeg_pkt.size, pltout);

#if VERSION_CHECK(LIBAVCODEC_VERSION_INT, <, 57, 8, 0, 57, 12, 100)
	av_free_packet (&mpeg_pkt);
#else
	av_packet_unref (&mpeg_pkt);
#endif
    }
#endif
}

static void ffmpeg_finish (void) {
    int i = 0;
    /* get the delayed frames */
#if LIBAVCODEC_VERSION_MAJOR < 54
    for (; frame_out_size; i++) {
        fflush (pltout);
        frame_out_size = avcodec_encode_video (codec_ctx,
                                               frame_outbuf,
                                               frame_outbuf_size,
                                               NULL);
        fwrite (frame_outbuf, 1, frame_out_size, pltout);
    }
#else
    for (mpeg_gout = 1; mpeg_gout; i++) {
#if LIBAVCODEC_VERSION_MAJOR < 58
        if (avcodec_encode_video2 (codec_ctx, &mpeg_pkt, NULL, &mpeg_gout) < 0)
            ERR (FATAL, name, "MPEG encoding error\n");
#else
	if (avcodec_send_frame(codec_ctx, NULL) < 0)
	    ERR (FATAL, name, "MPEG send frame error\n");
	mpeg_gout = avcodec_receive_packet(codec_ctx, &mpeg_pkt);
	if (mpeg_gout == AVERROR(EAGAIN) || mpeg_gout == AVERROR_EOF ||
	    mpeg_gout < 0)
	    ERR (FATAL, name, "MPEG receive packet error\n");
#endif
        if (mpeg_gout) {
            fwrite (mpeg_pkt.data, 1, mpeg_pkt.size, pltout);
#if VERSION_CHECK(LIBAVCODEC_VERSION_INT, <, 57, 8, 0, 57, 12, 100)
	av_free_packet (&mpeg_pkt);
#else
	av_packet_unref (&mpeg_pkt);
#endif
        }
    }
#endif
    /* add sequence end code to have a real mpeg file */
    frame_outbuf[0] = 0x00;
    frame_outbuf[1] = 0x00;
    frame_outbuf[2] = 0x01;
    frame_outbuf[3] = 0xb7;
    fwrite (frame_outbuf, 1, 4, pltout);
}
#endif
