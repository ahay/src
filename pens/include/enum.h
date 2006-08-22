
/*
 * enumeration values
 */

/*
 * flags
 */
#define NEVER_SET	-1
#ifndef NO
#define NO  		0
#endif
#ifndef YES
#define YES 		1
#endif

/*
 * size
 */
#define RELATIVE 	0	/* ABSOLUTE is also a scale */

/*
 * breaks
 */
#define	BREAK_BREAK	0
#define	BREAK_ERASE	1
#define	BREAK_IGNORE	2

/*
 * erases
 */
#define	DO_LITERALS		1
#define FORCE_INITIAL	2

/*
 * Color table entry states
 */
#define UNSET 	0		/* Must be zero */
#define SET		1
#define MAPPED	2
#define MAP_SET 3		/* MAPPED & SET */

/*
 * Color table array elements
 */
#define STATUS	0
#define MAP		1
#define _RED		2
#define _GREEN	3
#define _BLUE	4
#define _GREY	5

/*
 * Number of Color table array elements
 */
#define _NUM_PRIM	6
