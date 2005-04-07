/*
 * device routine table structure
 */

struct device{

	/* control routines */
	int (*open)();
	int (*reset)();
	int (*message)();
	int (*erase)();
	int (*close)();

	/* high level output */
	int (*vector)();
	int (*marker)();
	int (*text)();
	int (*area)();
	int (*raster)();
	int (*point)();
	int (*attributes)();

	/* input */
	int (*getpoint)();
	int (*interact)();

	/* low level output
	 * - these are called by the generic high level output routines
	 */
	int (*plot)();
	int (*startpoly)();
	int (*midpoly)();
	int (*endpoly)();
};
