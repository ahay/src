#ifndef PROGRESS_MONITOR_H
#define PROGRESS_MONITOR_H 1

#include <stdio.h>

class ProgressMonitor
{
public:
  ProgressMonitor ():_object (0)
  {
  }

  /** Call with "this" as first argument.
      Call before you do any work, with 0. as fraction.
      First object to call controls progress.
      Nested calls from other objects will be ignored.
      @in object pointer to your object "this"
      @in fraction Ranges from 0 to 1. When value is 1 or greater,
      then this object will be released for use by others.
      @return Returns 1 if user has canceled the activity;
      otherwise returns 0;
  */

  virtual int set (const void *object, double fraction)
  {
    int result = 1;
    if (0 == _object)
      _object = object;
    if (object == _object) {
      if (fraction >= 1.) {
        _object = 0;
      }
      if (fraction < 0.)
        fraction = 0.;
      if (fraction > 1.)
        fraction = 1.;
      result = update (fraction);
    }
    return result;
  }

  virtual ~ ProgressMonitor () {
  }

protected:
  /** Derived classes should override this update method
      Return 1 if the user wants to termininate processing.
      Otherwise, return 0;
   */
  virtual int update (double fraction)
  {
    fprintf (stdout, "Completed %g\n", fraction);
    fflush (stdout);
    return 0;
  }
private:
  const void *_object;
};

#endif
