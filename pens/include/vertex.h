#ifndef _pens_vertex_h
#define _pens_vertex_h

struct vertex
{
	int x;
	int y;
	struct vertex *next;		/* pointer to next vertex */
	struct vertex *last;		/* pointer to last vertex */
	struct vertex *soft;		/* pointer to some other vertex */
};

#endif
