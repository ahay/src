#include <rsf.h>

enum{
    EMPTY=-100,
	NOLINK=-10,
	ENDOFLIST=-1,
	UNCLAIMED=0,
	CLAIMED=1,
	INTERIOR=0
	};

#define POLYS  5000	/* Memory alloted for polygon */
#define MAXPOL 500	/* Maximum number of polygons made */

static int poly[POLYS][5];
static int pols[MAXPOL];	/* Point to the start of each polygon */
static int polsc[MAXPOL];	/* Cycle length */
static int npols;	/* How many polygons we have */

static int pedge[POLYS];	/* Edge points of the polygon */
static int nedge, nedge2;	/* number of edge points */

static int point, endlist;

static int edge (int x, int y, int xwmax, int xwmin, int ywmax, int ywmin);
static void insert (int where, int x, int y, int z);
static void delete (int where);

/* Read in the data for polystart */
void vp_polyfix (int x, int y, bool *first, bool *allgone)  
{
    static int oldx, oldy;

    if (*first) {
	*first = false;
	point = 1;
	*allgone = false;
	oldx = x + 1;
	oldy = y + 1;		
    }

    /* ignore repeated points */
    if ((x == oldx) && (y == oldy)) return;
    
    poly[point][0] = x;		/* X of vertex */
    poly[point][1] = y;		/* Y of vertex */
    oldx = x;
    oldy = y;
    point++;
    if (point >= POLYS)
	sf_error ("%s: Not enough memory for polygon.",__FILE__);
}

/* Start working on the polygons */
void vp_polystart (int xwmax, int xwmin, int ywmax, int ywmin)		
{
    int i, j, k, l, ii;
    int firstpoint;
    bool flag, double_check;
    int temp1, temp2, temp;

    endlist = point;		/* Last element in use */
    /* initialize array */
    for (i = 0; i <= point; i++) {
	poly[i][4] = INTERIOR;
	poly[i][2] = i + 1;
	poly[i][3] = i - 1;
    }
    poly[0][0] = poly[0][1] = EMPTY;
    poly[point][0] = poly[point][1] = EMPTY;

    /* fix links of free list and last point */
    poly[0][3] = NOLINK;
    poly[0][2] = point;
    poly[point][2] = ENDOFLIST;
    poly[point][3] = 0;		/* point to start of free list */

    poly[point - 1][2] = 1;	/* point to start of polygon */
    poly[1][3] = point - 1;

    /* Now we have our polygon stored as a circular double-linked list! */

    /*
    if (smart_clip)
	goto skip_it;
    */

    nedge = 0;
    i = firstpoint = 1;

    do {
	poly[i][4] = edge (poly[i][0], poly[i][1], xwmax, xwmin, ywmax, ywmin);
	if (poly[i][4])	{ /* It's on an edge */
	    nedge++;
	    pedge[nedge] = i;
	}
	i = poly[i][2];
    } while (firstpoint != i);


    /* Look at each edge segment. See if any edge points occur inside
     * segments joining two other edge points. */
    for (ii = 1; ii <= nedge; ii++) {
	i = pedge[ii];
	k = (poly[i][4] & poly[poly[i][2]][4]);
	if (k) {
	    /* We have an edge link. That is, it connects to edge vertices. */
	    /* See if any other edge points lie inbetween these two. */
	    for (j = 1; j <= nedge; j++) {
		if (k & poly[pedge[j]][4]) {
		    /* This point is on the correct edge */

		    l = ((k & (8 + 2)) > 0);
		    
		    /* l = 0 if top or bottom, 1 if left or right */

		    if (poly[i][l] > poly[poly[i][2]][l]) {
			if ((poly[pedge[j]][l] > poly[poly[i][2]][l])
			    && (poly[pedge[j]][l] < poly[i][l])) {
			    nedge++;
			    insert (i, poly[pedge[j]][0],
				    poly[pedge[j]][1],
				    poly[pedge[j]][4]);
			    pedge[nedge] = poly[i][2];
			}
		    } else {
			if ((poly[pedge[j]][l] < poly[poly[i][2]][l])
			    && (poly[pedge[j]][l] > poly[i][l]))
			{
			    nedge++;
			    insert (i, poly[pedge[j]][0],
				    poly[pedge[j]][1],
				    poly[pedge[j]][4]);
			    pedge[nedge] = poly[i][2];
			}
		    }
		}
	    }
	}
    }

    double_check = true;
    while (double_check) {
	double_check = false;

	for (j = 1; j <= nedge; j++) {
	    flag = true;
	    while (flag) {
		flag = false;

		if ( (poly[pedge[j]][4] & poly[poly[pedge[j]][2]][4]) && 
		    ((poly[pedge[j]][0] != poly[poly[pedge[j]][2]][0]) || 
		     (poly[pedge[j]][1] != poly[poly[pedge[j]][2]][1]))) {
		    for (k = j; k <= nedge; k++) {
			if ((poly[pedge[j]][0] == poly[pedge[k]][0]) && 
			    (poly[pedge[j]][1] == poly[pedge[k]][1])) {
			    /* OK, it is in the same position. See if either
			     * its succesor or predecessor is the same as the
			     * other end of the link. */
			    if ((poly[poly[pedge[j]][2]][0] ==
				 poly[poly[pedge[k]][3]][0]) && 
				(poly[poly[pedge[j]][2]][1] ==
				 poly[poly[pedge[k]][3]][1])) {
				/* links go in opposite directions. (The easy
				 * case) Break the links and instead link
				 * points at the same location. Delete
				 * duplicated point. */
				temp1 = poly[pedge[j]][2];
				/* remember pedge[j]'s succesor */
				poly[pedge[j]][2] = pedge[k];
				/* pedge[j]'s succesor is now pedge[k] */
				temp2 = poly[pedge[k]][3];
				/* remember pedge[k]'s predecessor */
				poly[pedge[k]][3] = pedge[j];
				poly[temp1][3] = temp2;
				poly[temp2][2] = temp1;
				/* Done. We have just split one polygon into
				 * two! */
				/* Clean up by removing repeated vertices */

				delete (pedge[k]);
				pedge[k] = 0;
				/* Point this to a place we know is marked as
				 * an interior point. Always fails checks to
				 * see if it is on the edge we want. */
				delete (temp2);
				flag = true;
				double_check = true;
			    } else {
				if ((j != k) &&
				    (poly[poly[pedge[j]][2]][0] ==
				     poly[poly[pedge[k]][2]][0]) && 
				    (poly[poly[pedge[j]][2]][1] ==
				     poly[poly[pedge[k]][2]][1])) {
				    /* the hard case. Both links go the same
				     * direction. Do as before, but Reverse
				     * the direction of one of the two
				     * pieces. */
				    temp1 = poly[pedge[j]][2];
				    temp2 = poly[pedge[k]][2];
				    i = temp1;
				    do {  /* Reverse one piece first */
					temp = poly[i][2];
					poly[i][2] = poly[i][3];
					poly[i][3] = temp;
					i = temp;
				    } while ((temp2 != i) && (temp1 != i));
				    poly[temp1][2] = temp2;
				    poly[temp2][3] = temp1;
				    if (i == temp2) {
					poly[pedge[j]][2] = pedge[k];
					poly[pedge[k]][3] = pedge[j];
					/* We have not created a new polygon,
					 * merely re-ordered one! */
				    } else {
					poly[pedge[j]][3] = pedge[k];
					poly[pedge[k]][2] = pedge[j];
					/* We have merged two polygons back
					 * into one! */
				    }
				    delete (pedge[k]);
				    pedge[k] = 0;
				    delete (temp2);
				    flag = true;
				    double_check = true;
				}
			    }
			}
			if (flag) break;
		    }
		}
	    }
	}
    }
}

#ifdef kjhkjh


/* Our polygon has been fragmented into multiple smaller polygons as
 * necessary! Output the resulting polygons */

    scan ();

    for (i = 1; i <= npols; i++)
    {
	dev.startpoly (polsc[i]);
	j = pols[i];
	do
	{
	    dev.midpoly (poly[j][0], poly[j][1]);
	    j = poly[j][2];
	}
	while (j != pols[i]);
	dev.endpoly (i == npols);
    }
}


/*
 * See how many polygons we have, and where they start.
 * Put this information in pols.
 */
scan ()
{
int             polyon[POLYS];	/* Array to remember whom this point belongs
				 * to */
int             joe;
int             cycle;
int             firstpoint;
int             where;

    for (joe = 1; joe < endlist; joe++)
    {
	if (poly[joe][0] == EMPTY)
	{
	    polyon[joe] = EMPTY;
	}
	else
	{
	    polyon[joe] = UNCLAIMED;
	}
    }

    npols = 0;
    for (joe = 1; joe < endlist; joe++)
    {
	if (polyon[joe] == UNCLAIMED)
	{
	    /* Found the start of a polygon */
	    firstpoint = joe;
	    polyon[firstpoint] = CLAIMED;
	    cycle = 1;
	    where = poly[firstpoint][2];
	    while (firstpoint != where)
	    {
		cycle++;
		polyon[where] = CLAIMED;
		where = poly[where][2];
		if (where == NOLINK)
		    ERR (FATAL, name, "Oh, No! Polygon list damaged! (This shouldn't be able to happen.)");
	    }
	    if (cycle < 3)
	    {			/* Not really a polygon, remove it */
		where = poly[firstpoint][2];
		while (firstpoint != where)
		{
		    delete (where);
		    where = poly[firstpoint][2];
		}
		delete (firstpoint);
	    }
	    else
	    {
		/* We found another polygon */
		npols++;
		if (npols > MAXPOL)
		    ERR (FATAL, name, "Too many polygons!");
		pols[npols] = firstpoint;
		polsc[npols] = cycle;
	    }
	}
    }
}

#endif

/* Find out which edges this point is on */
static int edge (int x, int y, int xwmax, int xwmin, int ywmax, int ywmin)
{
    int bottom, left;

    bottom = left = 0;
    if      (x == xwmin) left = 2;
    else if (x == xwmax) left = 8;

    if      (y == ywmin) bottom = 1;
    else if (y == ywmax) bottom = 4;

    return (bottom + left);
}

/* insert the given vertex between the element where and its forward link
 * remove an element from the free list */
static void insert (int where, int x, int y, int z)
{
    int temp;

    
    temp = poly[0][2];		/* free element */
    if (temp == endlist) {
	endlist++; /* Need to make our list one longer. */
	if (endlist >= POLYS)
	    sf_error ("%s: Ran out of memory on the polygon",__FILE__);
	poly[endlist][3] = temp;
	poly[endlist][2] = ENDOFLIST;
	poly[temp][2] = endlist;
    }
    /* OK, Now you can remove it, it isn't at the end anymore. */
    poly[0][2] = poly[temp][2];
    poly[poly[temp][2]][3] = 0;

    poly[temp][0] = x;		/* Put vertex in for element */
    poly[temp][1] = y;
    poly[temp][4] = z;

    /* update links */

    poly[temp][2] = poly[where][2];
    poly[where][2] = temp;
    
    poly[temp][3] = poly[poly[temp][2]][3];
    poly[poly[temp][2]][3] = temp;
}

/* Get rid of element number where; put it on the free list
 * This must work even if element where points to itself! */
static void delete (int where)
{
    int temp;

    poly[where][0] = poly[where][1] = EMPTY;
    poly[where][4] = INTERIOR;
    temp = poly[where][3];
    poly[where][3] = 0;
    poly[temp][2] = poly[where][2];
    poly[where][2] = poly[0][2];
    poly[0][2] = where;
    poly[poly[where][2]][3] = where;
    poly[poly[temp][2]][3] = temp;
}
