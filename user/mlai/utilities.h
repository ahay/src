#ifndef UTILITIES_H
#define UTILITIES_H

#include <sys/stat.h>
#include <stdio.h>
#define SECONDS_TO_SLEEP 30
#define DO_DEBUG_FILE_NAME "/tmp/DEBUG"

static inline int doesFileExist(const char* filename)
/*< Checks if a file exist using stat() function. Returns 1 
 * if the file exist otherwise returns 0.
 * >*/
{
    struct stat buffer;
    int exist = stat(filename,&buffer);
    if(exist == 0)
        return 1;
    else // -1
        return 0;
}

static inline void sleepIfFileExist(int seconds_to_sleep, const char* file_name)
/*<Sleeps if a file exists>*/
{	
    if (doesFileExist(file_name)){
    	
		printf("The file %s exists.  Waiting %d seconds for someone to attach ...\n", 
				file_name, 
				seconds_to_sleep);
		sleep(seconds_to_sleep);
		printf("Done waiting!\n");
    }
}

#endif /* UTILITIES_H */
