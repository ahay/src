/* Simplified system command */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "system.h"
#include "error.h"

void sf_system(const char *command)
/*< System command >*/
{
    pid_t pid;
    int status;

    pid = fork();

    if (pid < 0) {
	sf_error("Failed to fork");
    } else if (pid==0) { /* child */
	execl("/bin/sh","sh","-c",command,(char*) NULL);
	_exit(127);
    } else { /* parent */
	waitpid(pid,&status,0);
	if (status < 0) sf_error("Failed to run \"%s\"",command);
    }
}

