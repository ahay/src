
/* functions for dealing with the list of frames and getting old images */

struct cmap_ {
unsigned char   red[65536] ;
unsigned char   green[65536];
unsigned char   blue[65536];
unsigned long  map[65536];
};

typedef struct xtf_ {
    int file_num;		/* file to read (in list of files) */
    long  file_position;        /* starting poition in the file */
    int  end_file;              /* file in which the frame ends */
    long  end_pos;              /* position in that file */
    long total_len;             /* total length in bytes */
    int   frame_num;            /* number of this frame */
    int has_image;              /* A stored image of this frame exists */
    int break_end;              /* This frame ended with a break (not erase)*/
		char filename[256];         /* File name */
    XImage* image;		/* An XIMage of this frame */
    Pixmap  pixmap;		/* A pixmap of this frame */
    struct cmap_ cmap;		/* the colormap for this frame */
    struct xtf_ *next;		/* The next frame */
    struct xtf_ *prev;		/* The previous frame */
  
    
} xtFrame ;


/* 
 * A list of frames that can be plotted, the parent list of which this list
 * is an orderd subset is the parent member.
 */

typedef struct xtfl_ {
    xtFrame* start;
    xtFrame* end;
    int num_frame;
    struct xtfl_ *parent;
} xtFrameList;

#ifdef __STDC__    

extern xtFrameList* xt_new_list( xtFrameList*, int*, int );
extern xtFrame* xt_frame( xtFrameList*, int, int );
extern xtFrame* xt_first_frame( xtFrameList* );
extern xtFrame* xt_last_frame( xtFrameList* );
extern xtFrame* xt_frame_num( xtFrameList*, int  );
extern xt_store_image( xtFrame* );
extern xt_put_image( xtFrame* );
extern xt_clear_images( xtFrameList* );    

#else
extern xtFrameList* xt_new_list();
extern xtFrame* xt_frame();
extern xtFrame* xt_first_frame();
extern xtFrame* xt_last_frame();
extern xtFrame* xt_frame_num();
extern xt_store_image();
extern xt_put_image();
extern xt_clear_images();    

#endif
