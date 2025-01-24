/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: msieve.h 1008 2016-11-11 17:48:13Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _MSIEVE_H_
#define _MSIEVE_H_

#ifdef __cplusplus
extern "C" {
#endif

	/* Lightweight factoring API */

#include <util.h>
#include "savefile.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* version info */

#define MSIEVE_MAJOR_VERSION 1
#define MSIEVE_MINOR_VERSION 54

#ifndef MSIEVE_SVN_VERSION
#define MSIEVE_SVN_VERSION "unknown"
#endif

/* The final output from the factorization is a linked
   list of msieve_factor structures, one for each factor
   found. */

enum msieve_factor_type {
	MSIEVE_COMPOSITE,
	MSIEVE_PRIME,
	MSIEVE_PROBABLE_PRIME
};

typedef struct msieve_factor {
	enum msieve_factor_type factor_type;
	char *number;
	struct msieve_factor *next;
} msieve_factor;

/* These flags are used to 'communicate' status and
   configuration info back and forth to a factorization
   in progress */

enum msieve_flags {
	MSIEVE_DEFAULT_FLAGS = 0,		/* just a placeholder */
	MSIEVE_FLAG_USE_LOGFILE = 0x01,	    /* append log info to a logfile */
	MSIEVE_FLAG_LOG_TO_STDOUT = 0x02,   /* print log info to the screen */
	MSIEVE_FLAG_STOP_SIEVING = 0x04,    /* tell library to stop sieving
					       when it is safe to do so */
	MSIEVE_FLAG_FACTORIZATION_DONE = 0x08,  /* set by the library if a
						   factorization completed */
	MSIEVE_FLAG_SIEVING_IN_PROGRESS = 0x10, /* set by the library when 
						   any sieving operations are 
						   in progress */
	MSIEVE_FLAG_SKIP_QS_CYCLES = 0x20,  /* do not perform exact tracking of
	                                    the number of cycles while sieving
					    is in progress; for distributed
					    sieving where exact progress info
					    is not needed, sieving clients can
					    save a lot of memory with this */
	MSIEVE_FLAG_NFS_POLY1 = 0x40,     /* if input is large enough, perform
	                                    stage 1 polynomial selection for NFS */
	MSIEVE_FLAG_NFS_POLYSIZE = 0x80,  /* if input is large enough, perform
	                                    NFS polynomial size optimization */
	MSIEVE_FLAG_NFS_POLYROOT = 0x100, /* if input is large enough, perform
	                                    NFS polynomial root optimization */
	MSIEVE_FLAG_NFS_SIEVE = 0x200,   /* if input is large enough, perform
	                                    sieving for NFS */
	MSIEVE_FLAG_NFS_FILTER = 0x400,  /* if input is large enough, perform
	                                    filtering phase for NFS */
	MSIEVE_FLAG_NFS_LA = 0x800,      /* if input is large enough, perform
	                                    linear algebra phase for NFS */
	MSIEVE_FLAG_NFS_SQRT = 0x1000,    /* if input is large enough, perform
	                                    square root phase for NFS */
	MSIEVE_FLAG_NFS_LA_RESTART = 0x2000,/* restart the NFS linear algbra */
	MSIEVE_FLAG_DEEP_ECM = 0x4000,   /* perform nontrivial-size ECM */
	MSIEVE_FLAG_NFS_ONLY = 0x8000    /* go straight to NFS */
};


/* One factorization is represented by a msieve_obj
   structure. This contains all the static information
   that gets passed from one stage of the factorization
   to another. If this was C++ it would be a simple object */

typedef struct {
	char *input;		  /* pointer to string version of the 
				     integer to be factored */
	msieve_factor *factors;   /* linked list of factors found (in
				     ascending order */
	volatile uint32 flags;	  /* input/output flags */
	savefile_t savefile;      /* data for savefile */
	char *logfile_name;       /* name of the logfile that will be
				     used for this factorization */
	uint32 seed1, seed2;      /* current state of random number generator
				     (updated as random numbers are created) */
	char *nfs_fbfile_name;    /* name of factor base file */
	uint32 max_relations;      /* the number of relations that the sieving
	                              stage will try to find. The default (0)
				      is to keep sieving until all necessary 
				      relations are found. */
	uint32 which_gpu;         /* ordinal ID of GPU to use */


	uint32 cache_size1;       /* bytes in level 1 cache */
	uint32 cache_size2;       /* bytes in level 2 cache */
	enum cpu_type cpu;

	uint32 num_threads;

#ifdef HAVE_MPI
	uint32 mpi_size;          /* number of MPI processes, each with
                                     num_threads threads */
	uint32 mpi_rank;          /* from 0 to mpi_size - 1 */

	uint32 mpi_nrows;         /* a 2-D MPI lanczos grid */
	uint32 mpi_ncols;
	MPI_Comm mpi_la_grid;
	MPI_Comm mpi_la_row_grid; /* communicator for the current MPI row */
	MPI_Comm mpi_la_col_grid; /* communicator for the current MPI col */
	uint32 mpi_la_row_rank;
	uint32 mpi_la_col_rank;
#endif

	char *mp_sprintf_buf;    /* scratch space for printing big integers */

	const char *nfs_args;   /* arguments for NFS */
} msieve_obj;

msieve_obj * msieve_obj_new(char *input_integer,
			    uint32 flags,
			    char *savefile_name,
			    char *logfile_name,
			    char *nfs_fbfile_name,
			    uint32 seed1,
			    uint32 seed2,
			    uint32 max_relations,
			    enum cpu_type cpu,
			    uint32 cache_size1,
			    uint32 cache_size2,
			    uint32 num_threads,
			    uint32 which_gpu,
			    const char *nfs_args);

msieve_obj * msieve_obj_free(msieve_obj *obj);

void msieve_run(msieve_obj *obj);
				
#define MSIEVE_DEFAULT_LOGFILE "msieve.log"
#define MSIEVE_DEFAULT_SAVEFILE "msieve.dat"
#define MSIEVE_DEFAULT_NFS_FBFILE "msieve.fb"

#ifdef __cplusplus
}
#endif

#endif /* _MSIEVE_H_ */
