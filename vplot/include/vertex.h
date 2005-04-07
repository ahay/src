struct vertex
{
	int x;
	int y;
	struct vertex *next;		/* pointer to next vertex */
	struct vertex *last;		/* pointer to last vertex */
	struct vertex *soft;		/* pointer to some other vertex */
};
