#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <pthread.h>

/* Time delay (secs) between checks for job requests */
#define SLEEP_TIME 1

/* Number of processors, ID of this processor */
static int NTask;
static int ThisTask;

/* Flag to signal job completion */
volatile int job_running;

/* Maximum length of command to execute */
#define COMMAND_LENGTH 4096
static char command[COMMAND_LENGTH];

/* MPI tags */
#define JOB_REQUEST_TAG 1
#define TERMINATION_TAG 2

/*
  Check to see what format specifier a string contains,
  if any. Returns index of the character at the end of the
  format string or 0 if there isn't one.
*/
size_t identify_format(char *str, size_t len)
{
  size_t i = 0;

  /* Advance to first '%' sign */
  while((i < len-1) && (str[i] != '%'))
    i += 1;
  
  /* Return 0 if we didn't find a % or we're at the end */
  if((str[i] != '%') || (i >= len-1))
    return 0;
  
  /* Advance to type character/flags etc */
  i += 1;
  if(i >= len)return 0;

  /* Skip over any flags */
  while((i < len-1) && 
	((str[i]=='+') ||
	 (str[i]=='-') ||
	 (str[i]==' ') ||
	 (str[i]=='#') ||
	 (str[i]=='0')))
    i += 1;
  
  /* Skip width */
  while((i < len-1) && 
	(str[i] >= '0') && 
	(str[i] <= '9'))
    i += 1;

  /* Skip precision */
  if(str[i] == '.')
    {
      /* Skip over dot */
      i+= 1;
      if(i>=len-1)
	return 0;
      
      /* Skip any digits */
      while((i < len-1) && 
	    (str[i] >= '0') && 
	    (str[i] <= '9'))
	i += 1;
    }

  /* Skip modifier */
  while((i < len-1) && 
	((str[i]=='h') ||
	 (str[i]=='l') ||
	 (str[i]=='L')))
    i += 1;

  /* Should now have type character */
  switch (str[i])
    {
    case 'd':
    case 'i':
    case 'f':
    case 'e':
    case 'E':
    case 'g':
    case 'G':
      return i;
    break;
    default:
      return 0;
      break;
    }
}



void *run_job(void *ptr)
{
  char cmd_exec[COMMAND_LENGTH];
  char tmp[COMMAND_LENGTH];
  int ijob = *((int *) ptr);
  size_t exec_offset = 0;
  size_t len = strlen(command);
  size_t offset = 0;
  size_t fpos;

  /* 
     Construct command by substituting in job index
     wherever we find an int or double format specifier
  */
  while(offset < len)
    {
      fpos = identify_format(command+offset, len-offset);
      if(fpos == 0)
	{
	  /* No format strings left, so just copy */
	  strncpy(cmd_exec+exec_offset, command+offset, 
		 COMMAND_LENGTH-exec_offset);
	  offset = len;
	}
      else
	{
	  /* Have to sub in the job number */
	  strncpy(tmp, command+offset, fpos+1);
	  tmp[fpos+1] = (char) 0;

	  switch (command[offset+fpos])
	    {
	    case 'd':
	    case 'i':
	      {
		sprintf(cmd_exec+exec_offset, tmp, ijob);
		exec_offset = strlen(cmd_exec);
		offset += fpos + 1;
	      }
	    break;
	    case 'f':
	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
	      {
		sprintf(cmd_exec+exec_offset, tmp, (double) ijob);
		exec_offset = strlen(cmd_exec);
		offset += fpos + 1;
	      }
	    break;
	    default:
	      /* Can't handle this format, so just copy */
	      strncpy(cmd_exec+exec_offset, command+offset, 
		     COMMAND_LENGTH-exec_offset);
	      offset = len;
	    }
	}
    }
  
  /* Run the command */
  printf("Running job %d on process %d\n", ijob, ThisTask);
  system(cmd_exec);
  printf("Job %d on process %d finished\n", ijob, ThisTask);
  job_running = 0;
  return NULL;
}



int main(int argc, char *argv[])
{
  int ifirst, ilast;
  int njobs_tot;
  int next_to_assign;
  int iproc;
  int job_received;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  /* Read command line args */
  if(ThisTask==0)
    {
      if(argc != 4)
	{
	  printf("Usage: parallel_tasks ifirst ilast command\n");
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
      sscanf(argv[1], "%d", &ifirst);
      sscanf(argv[2], "%d", &ilast);
      strncpy(command, argv[3], COMMAND_LENGTH);
      command[COMMAND_LENGTH-1] = (char) 0;

      if((ifirst < 0) || (ilast < 0))
	{
	  printf("Job index must be non-negative!\n");
	  MPI_Abort(MPI_COMM_WORLD, 1);  
	}
    }

  if(ThisTask==0)
    printf("Parallel tasks - command is: %s\n", command);

  /* Broadcast arguments */
  MPI_Bcast(&ifirst, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ilast,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(command, COMMAND_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);

  /* Initial number of jobs to assign */
  njobs_tot      = ilast-ifirst+1;
  next_to_assign = ifirst;

  if(ThisTask==0)
    {
      /*
	Processor 0 needs to handle job requests and 
	run its own jobs at the same time. Jobs on this
	process are run in a separate thread so we can continue
	responding to job requests. The main thread spends most
	of its time sleeping to avoid taking cpu time from the jobs.
      */
      int nfinished  = 0;
      int proc0_done = 0;
      int local_job  = -1;
      pthread_t job_thread;
      job_running    = 0;
      while((nfinished < NTask-1) || (proc0_done==0))
	{
	  int flag;
	  MPI_Status probe_status;
	  MPI_Status recv_status;
	  /* If no local job is running, try to start one */
	  if(job_running == 0)
	    {
	      /* Wait for thread which ran the previous job to finish */
	      if(local_job >= 0)
		{
		  pthread_join(job_thread, NULL);
		  local_job = -1;
		}

	      /* Check if we have jobs left to assign */
	      if(next_to_assign <= ilast)
		{
		  /* Launch the next job in a new thread */
		  job_running = 1;
		  local_job = next_to_assign;
		  pthread_create(&job_thread, NULL, 
				 &run_job, 
				 (void *) &local_job);
		  next_to_assign += 1;
		}
	      else
		{
		  /* No more jobs left */
		  proc0_done = 1;
		}
	    }

	  /* 
	     Check for job requests from other processes.
	     May be more than one message waiting
	  */
	  flag = 1;
	  while(flag)
	    {
	      MPI_Iprobe(MPI_ANY_SOURCE, JOB_REQUEST_TAG, MPI_COMM_WORLD, 
			 &flag, &probe_status);
	      if(flag)
		{
		  /* We have a job request to deal with  */
		  int ireq;
		  int ijob;
		  MPI_Recv(&ireq, 1, MPI_INT, 
			   probe_status.MPI_SOURCE, probe_status.MPI_TAG, 
			   MPI_COMM_WORLD, &recv_status);
		  /* Check if we have jobs to hand out */
		  if(next_to_assign <= ilast)
		    {
		      ijob = next_to_assign;
		      next_to_assign += 1;
		    }
		  else
		    {
		      /* No more jobs for this processor */
		      nfinished += 1;
		      ijob = -1;
		    }
		  /* Send the job index back */
		  MPI_Send(&ijob, 1, MPI_INT, probe_status.MPI_SOURCE, 
			   JOB_REQUEST_TAG, MPI_COMM_WORLD);
		}
	    }

	  /* Go back to sleep */
	  sleep(SLEEP_TIME);
	}
      /* 
	 At this point, all jobs are complete so signal other tasks.
	 This is just to avoid having them sit at 100% cpu in MPI_Barrier.
      */
      int i;
      int dummy = 1;
      for(i=1;i<NTask;i+=1)
	MPI_Send(&dummy, 1, MPI_INT, i, TERMINATION_TAG, MPI_COMM_WORLD);
    }
  else
    {
      /*
	Processors other than 0 request and execute jobs
	until there are none left (signalled by -ve job ID).
      */
      while(1)
	{
	  int ireq = 1;
	  int ijob;
	  MPI_Status status;
	  /* Ask for a job */
	  MPI_Sendrecv(&ireq, 1, MPI_INT, 0, JOB_REQUEST_TAG,
		       &ijob, 1, MPI_INT, 0, JOB_REQUEST_TAG,
		       MPI_COMM_WORLD, &status);
	  if(ijob >= 0)
	    {
	      /* Run the job if we got one */
	      run_job(&ijob);
	    }
	  else
	    {
	      /* If there are no jobs left, we're done */
	      break;
	    }
	}
      /*
	Now wait for termination signal from task 0
      */
      int flag = 0;
      while(!flag)
	{
	  MPI_Status probe_status;
	  MPI_Status recv_status;	  
	  MPI_Iprobe(0, TERMINATION_TAG, MPI_COMM_WORLD, 
		     &flag, &probe_status);
	  if(flag)
	    {
	      /* Termination message is waiting */
	      int dummy;
	      MPI_Recv(&dummy, 1, MPI_INT, 
		       probe_status.MPI_SOURCE, probe_status.MPI_TAG, 
		       MPI_COMM_WORLD, &recv_status);
	    }
	  else
	    {
	      /* Wait a bit before trying again */
	      sleep(SLEEP_TIME);
	    }
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask==0)
    printf("All jobs complete.\n");

  MPI_Finalize();
  
  return 0;
}
