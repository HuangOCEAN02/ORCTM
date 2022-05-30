#include <sys/time.h>

/* Replament if MPI_Wtime is not present,
   UNIX systems only */

#ifdef NOMPI

double mpi_wtime_()
{
   struct timeval tv;
   double d;

   gettimeofday(&tv,0);
   d = tv.tv_sec + tv.tv_usec*1.e-6;

  return d;
}

#endif
