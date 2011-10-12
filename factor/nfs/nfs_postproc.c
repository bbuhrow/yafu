#include "nfs.h"
#include "gmp_xface.h"

uint32 do_msieve_filtering(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN)
{
	FILE *tmp, *logfile;
	uint32 relations_needed;
	uint32 flags = 0;

	flags = flags | MSIEVE_FLAG_USE_LOGFILE;
	if (VFLAG > 0)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_FILTER;
	obj->flags = flags;

	//msieve: filter
	if (VFLAG >= 0)
		printf("nfs: commencing msieve filtering\n");

	logfile = fopen(fobj->flogname, "a");
	if (logfile == NULL)
		printf("could not open yafu logfile for appending\n");
	else
	{
		logprint(logfile, "nfs: commencing msieve filtering\n");
		fclose(logfile);
	}

	tmp = fopen(fobj->nfs_obj.fbfile,"r");
	if (tmp == NULL)
		ggnfs_to_msieve(fobj, job);
	else
		fclose(tmp);

	printf("%s\n",mp_print(mpN, 10, NULL, gstr1.s));
	relations_needed = nfs_filter_relations(obj, mpN);

	return relations_needed;
}


void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj)
{
	int i;
	FILE *logfile;
	char c[3];

	// extract the factors
	for (i=0;i<factor_list->num_factors;i++)
	{
		z tmpz;
		mpz_t tmp;

		//init locals
		zInit(&tmpz);
		mpz_init(tmp);
		
		//convert the factor
		mp_t2gmp(&factor_list->final_factors[i]->factor,tmp);

		//divide it out
		mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, tmp);

		//check if its prime and log accordingly
		if (mpz_probab_prime_p(tmp, NUM_WITNESSES))
		{
			//need to convert to yafu bigint to store
			gmp2mp(tmp, &tmpz);
			tmpz.type = PRP;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"prp");
		}
		else
		{
			gmp2mp(tmp, &tmpz);
			tmpz.type = COMPOSITE;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"C");
		}

		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "%s%d = %s\n",c,
				mpz_sizeinbase(tmp, 10), mpz_get_str(gstr1.s, 10, tmp));
			fclose(logfile);
		}		

		//free locals
		zFree(&tmpz);
		mpz_clear(tmp);
	}

	//log anything left over
	if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) > 0) 
	{
		z tmpz;
		char c[3];

		zInit(&tmpz);
		if (mpz_probab_prime_p(fobj->nfs_obj.gmp_n, NUM_WITNESSES))
		{
			gmp2mp(fobj->nfs_obj.gmp_n, &tmpz);
			tmpz.type = PRP;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"prp");			
		}
		else
		{
			gmp2mp(fobj->nfs_obj.gmp_n, &tmpz);
			tmpz.type = COMPOSITE;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"C");
		}
		
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "%s%d = %s\n",c,
				mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10), 
				mpz_get_str(gstr1.s, 10, fobj->nfs_obj.gmp_n));
			fclose(logfile);
		}		

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		zFree(&tmpz);
	}

	return;
}

