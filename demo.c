/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: demo.c 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_msieve.h>
#include <signal.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

msieve_obj *g_curr_factorization = NULL;

/*--------------------------------------------------------------------*/
void handle_signal(int sig) {

	msieve_obj *obj = g_curr_factorization;

	printf("\nreceived signal %d; shutting down\n", sig);
	
	if (obj && (obj->flags & MSIEVE_FLAG_SIEVING_IN_PROGRESS))
		obj->flags |= MSIEVE_FLAG_STOP_SIEVING;
	else
		_exit(0);
}

/*--------------------------------------------------------------------*/
void get_random_seeds(uint32 *seed1, uint32 *seed2) {

	uint32 tmp_seed1, tmp_seed2;

	/* In a multithreaded program, every msieve object
	   should have two unique, non-correlated seeds
	   chosen for it */

#if !defined(WIN32) && !defined(_WIN64)

	FILE *rand_device = fopen("/dev/urandom", "r");

	if (rand_device != NULL) {

		/* Yay! Cryptographic-quality nondeterministic randomness! */

		fread(&tmp_seed1, sizeof(uint32), (size_t)1, rand_device);
		fread(&tmp_seed2, sizeof(uint32), (size_t)1, rand_device);
		fclose(rand_device);
	}
	else

#endif
	{
		/* <Shrug> For everyone else, sample the current time,
		   the high-res timer (hopefully not correlated to the
		   current time), and the process ID. Multithreaded
		   applications should fold in the thread ID too */

		uint64 high_res_time = read_clock();
		tmp_seed1 = ((uint32)(high_res_time >> 32) ^
			     (uint32)time(NULL)) * 
			    (uint32)getpid();
		tmp_seed2 = (uint32)high_res_time;
	}

	/* The final seeds are the result of a multiplicative
	   hash of the initial seeds */

	(*seed1) = tmp_seed1 * ((uint32)40499 * 65543);
	(*seed2) = tmp_seed2 * ((uint32)40499 * 65543);
}

/*--------------------------------------------------------------------*/
void print_usage(char *progname) {

	printf("\nMsieve v. %d.%02d (SVN %s)\n", MSIEVE_MAJOR_VERSION, 
					MSIEVE_MINOR_VERSION,
					MSIEVE_SVN_VERSION);

	printf("\nusage: %s [options] [one_number]\n", progname);
	printf("\nnumbers starting with '0' are treated as octal,\n"
		"numbers starting with '0x' are treated as hexadecimal\n");
	printf("\noptions:\n"
	         "   -s <name> save intermediate results to <name>\n"
		 "             instead of the default %s\n"
	         "   -l <name> append log information to <name>\n"
		 "             instead of the default %s\n"
	         "   -i <name> read one or more integers to factor from\n"
		 "             <name> (default worktodo.ini) instead of\n"
		 "             from the command line\n"
		 "   -m        manual mode: enter numbers via standard input\n"
	         "   -q        quiet: do not generate any log information,\n"
		 "             only print any factors found\n"
	         "   -d <min>  deadline: if still sieving after <min>\n"
		 "             minutes, shut down gracefully (default off)\n"
		 "   -r <num>  stop sieving after finding <num> relations\n"
		 "   -p        run at idle priority\n"
	         "   -v        verbose: write log information to screen\n"
		 "             as well as to logfile\n"
		 "   -z        you are Paul Zimmermann\n"
#ifdef HAVE_CUDA
		 "   -g <num>  use GPU <num>, 0 <= num < (# graphics cards)>\n"
#endif
	         "   -t <num>  use at most <num> threads\n"
		 "\n"
		 " elliptic curve options:\n"
		 "   -e        perform 'deep' ECM, seek factors > 15 digits\n\n"
		 " quadratic sieve options:\n"
		 "   -c        client: only perform sieving\n\n"
		 " number field sieve options:\n\n"
		 "           [nfs_phase] \"arguments\"\n\n"
		 " where the first part is one or more of:\n"
		 "   -n        use the number field sieve (80+ digits only;\n"
		 "             performs all NFS tasks in order)\n"
	         "   -nf <name> read from / write to NFS factor base file\n"
		 "             <name> instead of the default %s\n"
		 "   -np       perform only NFS polynomial selection\n"
		 "   -np1      perform stage 1 of NFS polynomial selection\n"
		 "   -nps      perform NFS polynomial size optimization\n"
		 "   -npr      perform NFS polynomial root optimization\n"
		 "   -ns       perform only NFS sieving\n"
		 "   -nc       perform only NFS combining (all phases)\n"
		 "   -nc1      perform only NFS filtering\n"
		 "   -nc2      perform only NFS linear algebra\n"
		 "   -ncr      perform only NFS linear algebra, restarting\n"
		 "             from a previous checkpoint\n"
		 "   -nc3      perform only NFS square root\n\n"
		 " the arguments are a space-delimited list of:\n"
		 " polynomial selection options:\n"
#ifdef HAVE_CUDA
		 "   sortlib=X       use GPU sorting library X\n"
		 "   gpu_mem_mb=X    use X megabytes of GPU memory\n"
#endif
		 "   polydegree=X    select polynomials with degree X\n"
		 "   min_coeff=X     minimum leading coefficient to search\n"
		 "                   in stage 1\n"
		 "   max_coeff=X     maximum leading coefficient to search\n"
		 "                   in stage 1\n"
		 "   stage1_norm=X   the maximum norm value for stage 1\n"
		 "   stage2_norm=X   the maximum norm value for stage 2\n"
		 "   min_evalue=X    the minimum score of saved polyomials\n"
		 "   poly_deadline=X stop searching after X seconds (0 means\n"
		 "                   search forever)\n"
		 "   X,Y             same as 'min_coeff=X max_coeff=Y'\n"
		 " line sieving options:\n"
		 "   X,Y             handle sieve lines X to Y inclusive\n"
		 " filtering options:\n"
		 "   filter_mem_mb=X  try to limit filtering memory use to\n"
		 "                    X megabytes\n"
		 "   filter_maxrels=X limit the filtering to using the first\n"
		 "                    X relations in the data file\n"
		 "   filter_lpbound=X have filtering start by only looking\n"
		 "                    at ideals of size X or larger\n"
		 "   target_density=X attempt to produce a matrix with X\n"
		 "                    entries per column\n"
		 "   max_weight=X     have filtering start by looking at ideals\n"
		 "                    of max weight >= X\n"
		 "   X,Y              same as 'filter_lpbound=X filter_maxrels=Y'\n"
		 " linear algebra options:\n"
		 "   skip_matbuild=1  start the linear algebra but skip building\n"
		 "                    the matrix (assumes it is built already)\n"
		 "   la_block=X       use a block size of X (512<=X<=65536)\n"
		 "   la_superblock=X  use a superblock size of X\n"
		 "   cado_filter=1    assume filtering used the CADO-NFS suite\n"
#ifdef HAVE_MPI
		 "   mpi_nrows=X      use a grid with X rows\n"
		 "   mpi_ncols=X      use a grid with X columns\n"
		 "   X,Y              same as 'mpi_nrows=X mpi_ncols=Y'\n"
		 "                    (if unspecified, default grid is\n"
		 "                    1 x [argument to mpirun])\n"
#endif
		 " square root options:\n"
		 "   dep_first=X start with dependency X, 1<=X<=64\n"
		 "   dep_last=Y  end with dependency Y, 1<=Y<=64\n"
		 "   X,Y         same as 'dep_first=X dep_last=Y'\n"
		 ,
		 MSIEVE_DEFAULT_SAVEFILE, 
		 MSIEVE_DEFAULT_LOGFILE,
		 MSIEVE_DEFAULT_NFS_FBFILE);
}

/*--------------------------------------------------------------------*/
void factor_integer(char *buf, uint32 flags,
		    char *savefile_name,
		    char *logfile_name,
		    char *nfs_fbfile_name,
		    uint32 *seed1, uint32 *seed2,
		    uint32 max_relations,
		    enum cpu_type cpu,
		    uint32 cache_size1,
		    uint32 cache_size2,
		    uint32 num_threads,
		    uint32 which_gpu,
		    const char *nfs_args) {
	
	char *int_start, *last;
	msieve_obj *obj;
	msieve_factor *factor;

	/* point to the start of the integer or expression;
	   if the start point indicates no integer is present,
	   don't try to factor it :) */

	last = strchr(buf, '\n');
	if (last)
		*last = 0;
	int_start = buf;
	while (*int_start && !isdigit(*int_start) &&
			*int_start != '(' ) {
		int_start++;
	}
	if (*int_start == 0)
		return;

	g_curr_factorization = msieve_obj_new(int_start, flags,
					savefile_name, logfile_name,
					nfs_fbfile_name,
					*seed1, *seed2, max_relations,
					cpu, cache_size1, cache_size2,
					num_threads, which_gpu,
					nfs_args);
	if (g_curr_factorization == NULL) {
		printf("factoring initialization failed\n");
		return;
	}

	msieve_run(g_curr_factorization);

	if (!(g_curr_factorization->flags & MSIEVE_FLAG_FACTORIZATION_DONE)) {
		printf("\ncurrent factorization was interrupted\n");
		exit(0);
	}

	/* If no logging is specified, at least print out the
	   factors that were found */

	if (!(g_curr_factorization->flags & (MSIEVE_FLAG_USE_LOGFILE |
					MSIEVE_FLAG_LOG_TO_STDOUT))) {
		factor = g_curr_factorization->factors;

		printf("\n");
		printf("%s\n", buf);
		while (factor != NULL) {
			char *factor_type;

			if (factor->factor_type == MSIEVE_PRIME)
				factor_type = "p";
			else if (factor->factor_type == MSIEVE_COMPOSITE)
				factor_type = "c";
			else
				factor_type = "prp";

			printf("%s%d: %s\n", factor_type, 
					(int32)strlen(factor->number), 
					factor->number);
			factor = factor->next;
		}
		printf("\n");
	}

	/* save the current value of the random seeds, so that
	   the next factorization will pick up the pseudorandom
	   sequence where this factorization left off */

	*seed1 = g_curr_factorization->seed1;
	*seed2 = g_curr_factorization->seed2;

	/* free the current factorization struct. The following
	   avoids a race condition in the signal handler */

	obj = g_curr_factorization;
	g_curr_factorization = NULL;
	if (obj)
		msieve_obj_free(obj);
}

#ifdef WIN32
DWORD WINAPI countdown_thread(LPVOID pminutes) {
	DWORD minutes = *(DWORD *)pminutes;

	if (minutes > 0x7fffffff / 60000)
		minutes = 0;            /* infinite */

	Sleep(minutes * 60000);
	raise(SIGINT);
	return 0;
}

#else
void *countdown_thread(void *pminutes) {
	uint32 minutes = *(uint32 *)pminutes;

	if (minutes > 0xffffffff / 60)
		minutes = 0xffffffff / 60;   /* infinite */

	sleep(minutes * 60);
	raise(SIGINT);
	return NULL;
}
#endif

/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {

	char buf[500];
	uint32 seed1, seed2;
	char *savefile_name = NULL;
	char *logfile_name = NULL;
	char *infile_name = "worktodo.ini";
	char *nfs_fbfile_name = NULL;
	uint32 flags;
	char manual_mode = 0;
	int i;
	int32 deadline = 0;
	uint32 max_relations = 0;
	enum cpu_type cpu;
	uint32 cache_size1; 
	uint32 cache_size2; 
	uint32 num_threads = 0;
	uint32 which_gpu = 0;
	const char *nfs_args = NULL;
		
	get_cache_sizes(&cache_size1, &cache_size2);
	cpu = get_cpu_type();

	if (signal(SIGINT, handle_signal) == SIG_ERR) {
	        printf("could not install handler on SIGINT\n");
	        return -1;
	}
	if (signal(SIGTERM, handle_signal) == SIG_ERR) {
	        printf("could not install handler on SIGTERM\n");
	        return -1;
	}     
#ifdef HAVE_MPI
	{
		int32 level;
		if ((i = MPI_Init_thread(&argc, &argv,
				MPI_THREAD_FUNNELED, &level)) != MPI_SUCCESS) {
			printf("error %d initializing MPI, aborting\n", i);
			MPI_Abort(MPI_COMM_WORLD, i);
		}
	}
#endif

	flags = MSIEVE_FLAG_USE_LOGFILE;

	i = 1;
	buf[0] = 0;
	while (i < argc) {
		if (argv[i][0] == (char)('-')) {
			switch(argv[i][1]) {
			case 'h':
			case '?':
				print_usage(argv[0]);
				return 0;

			case 'i':
			case 's':
			case 'l':
				if (i + 1 < argc && argv[i+1][0] != '-') {
					if (tolower(argv[i][1]) == 'i')
						infile_name = argv[i+1];
					else if (tolower(argv[i][1]) == 's') {
						char *p;
						savefile_name = argv[i+1];
						if((p=strstr(savefile_name, ".gz"))) {
							savefile_name = strdup(argv[i+1]);
							savefile_name[p-argv[i+1]] = 0;
						}
					} 
					else {
						logfile_name = argv[i+1];
					}
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
					
			case 'm':
				manual_mode = 1;
				i++;
				break;

			case 'e':
				flags |= MSIEVE_FLAG_DEEP_ECM;
				i++;
				break;

			case 'z':
				flags |= MSIEVE_FLAG_NFS_ONLY;
				i++;
				break;

			case 'n': 
				switch (argv[i][2]) {
				case 0:
					flags |= MSIEVE_FLAG_NFS_POLY1 |
						 MSIEVE_FLAG_NFS_POLYSIZE |
						 MSIEVE_FLAG_NFS_POLYROOT |
						 MSIEVE_FLAG_NFS_SIEVE |
						 MSIEVE_FLAG_NFS_FILTER |
						 MSIEVE_FLAG_NFS_LA |
						 MSIEVE_FLAG_NFS_SQRT;
					break;

				case 'f':
					break;

				case 'p':
					if (argv[i][3] == '1')
						flags |= MSIEVE_FLAG_NFS_POLY1;
					else if (argv[i][3] == 's')
						flags |= MSIEVE_FLAG_NFS_POLYSIZE;
					else if (argv[i][3] == 'r')
						flags |= MSIEVE_FLAG_NFS_POLYROOT;
					else
						flags |= MSIEVE_FLAG_NFS_POLY1 |
							 MSIEVE_FLAG_NFS_POLYSIZE |
							 MSIEVE_FLAG_NFS_POLYROOT;
					break;
				
				case 's':
					flags |= MSIEVE_FLAG_NFS_SIEVE;
					break;

				case 'c':
					switch (argv[i][3]) {
					case 0:
						flags |= MSIEVE_FLAG_NFS_FILTER |
						         MSIEVE_FLAG_NFS_LA |
						         MSIEVE_FLAG_NFS_SQRT;
						break;

					case '1':
						flags |= MSIEVE_FLAG_NFS_FILTER;
						break;

					case 'r':
						flags |= MSIEVE_FLAG_NFS_LA_RESTART;
					case '2': /* fall through */
						flags |= MSIEVE_FLAG_NFS_LA;
						break;

					case '3':
						flags |= MSIEVE_FLAG_NFS_SQRT;
						break;
					}
					break;

				default:
					print_usage(argv[0]);
					return -1;
				}

				if (i + 1 < argc && argv[i+1][0] != '-') {
					if (argv[i][2] == 'f')
						nfs_fbfile_name = argv[++i];
					else
						nfs_args = argv[++i];
				}
				i++;
				break;
					
			case 'q':
				flags &= ~(MSIEVE_FLAG_USE_LOGFILE |
					   MSIEVE_FLAG_LOG_TO_STDOUT);
				i++;
				break;
					
			case 'd':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					deadline = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
					
			case 'r':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					max_relations = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
					
			case 't':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					num_threads = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
#ifdef HAVE_CUDA					
			case 'g':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					which_gpu = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
#endif					
			case 'c':
				flags |= MSIEVE_FLAG_SKIP_QS_CYCLES;
				i++;
				break;

			case 'v':
				flags |= MSIEVE_FLAG_LOG_TO_STDOUT;
				i++;
				break;

			case 'p':
				set_idle_priority();
				i++;
				break;

			default:
				print_usage(argv[0]);
				return -1;
			}
		}
		else {
			if (isdigit(argv[i][0]) || argv[i][0] == '(' )
				strncpy(buf, argv[i], sizeof(buf));
			i++;
		}
	}

	get_random_seeds(&seed1, &seed2);

	if (deadline) {
#if defined(WIN32) || defined(_WIN64)
		DWORD thread_id;
		CreateThread(NULL, 0, countdown_thread, 
				&deadline, 0, &thread_id);
#else
		pthread_t thread_id;
		pthread_create(&thread_id, NULL, 
				countdown_thread, &deadline);
#endif
	}

	if (isdigit(buf[0]) || buf[0] == '(' ) {
		factor_integer(buf, flags, savefile_name, 
				logfile_name, nfs_fbfile_name,
				&seed1, &seed2,
				max_relations, 
				cpu, cache_size1, cache_size2,
				num_threads, which_gpu,
				nfs_args);
	}
	else if (manual_mode) {
		while (1) {
			printf("\n\nnext number: ");
			fflush(stdout);
			buf[0] = 0;
			fgets(buf, (int)sizeof(buf), stdin);
			factor_integer(buf, flags, savefile_name, 
					logfile_name, nfs_fbfile_name,
					&seed1, &seed2,
					max_relations, 
					cpu, cache_size1, cache_size2,
					num_threads, which_gpu, nfs_args);
			if (feof(stdin))
				break;
		}
	}
	else {
		FILE *infile = fopen(infile_name, "r");
		if (infile == NULL) {
			printf("cannot open input file '%s'\n", infile_name);
			return 0;
		}

		while (1) {
			buf[0] = 0;
			fgets(buf, (int)sizeof(buf), infile);
			factor_integer(buf, flags, savefile_name, 
					logfile_name, nfs_fbfile_name,
					&seed1, &seed2,
					max_relations, 
					cpu, cache_size1, cache_size2,
					num_threads, which_gpu, nfs_args);
			if (feof(infile))
				break;
		}
		fclose(infile);
	}

#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	return 0;
}
