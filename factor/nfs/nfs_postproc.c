/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 12/6/2012
----------------------------------------------------------------------*/

#include "nfs_impl.h"
#include "gmp_xface.h"
#include "arith.h"

#ifdef USE_NFS

uint32_t do_msieve_filtering(fact_obj_t *fobj, msieve_obj *obj, nfs_job_t *job)
{
	FILE *tmp, *logfile;
	uint32_t relations_needed;
	uint32_t flags = 0;
	char nfs_args[100];

	flags = flags | MSIEVE_FLAG_USE_LOGFILE;
	if (fobj->VFLAG > 0)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_FILTER;
	obj->flags = flags;

	if (job->use_max_rels > 0)
	{
		sprintf(nfs_args, "filter_maxrels=%u", job->use_max_rels);
		obj->nfs_args = nfs_args;
		if (fobj->VFLAG >= 0)
			printf("nfs: commencing msieve filtering with reduced relation set\n");
	}
	else
	{
		if (fobj->VFLAG >= 0)
			printf("nfs: commencing msieve filtering\n");
	}
	
	if(fobj->nfs_obj.td > 0){
	    if(strlen(nfs_args)>0)
            sprintf(nfs_args + strlen(nfs_args), " target_density=%d", fobj->nfs_obj.td);
        else
            sprintf(nfs_args, "target_density=%d",fobj->nfs_obj.td);

	    obj->nfs_args = nfs_args;
    }

    if (fobj->LOGFLAG)
    {
        logfile = fopen(fobj->flogname, "a");
        if (logfile == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open yafu logfile for appending\n");
        }
        else
        {
            if (job->use_max_rels > 0)
                logprint(logfile, "nfs: commencing msieve filtering with reduced relation set\n");
            else
                logprint(logfile, "nfs: commencing msieve filtering\n");
            fclose(logfile);
        }
    }

	tmp = fopen(fobj->nfs_obj.fbfile,"r");
	if (tmp == NULL)
		ggnfs_to_msieve(fobj, job);
	//else
	//	fclose(tmp);
	else // test if the fb file is for this job
	{
		char line[GSTR_MAXSIZE];
		mpz_t num, r;
		mpz_init(num);
		mpz_init(r);
		

		while (fgets(line,GSTR_MAXSIZE,tmp))
		{
			if (line[0] == 'N')
			{
				mpz_set_str(num, line + 2, 0);
				mpz_tdiv_r(r, fobj->nfs_obj.gmp_n, num);
				if (mpz_cmp_ui(r, 0) == 0) // match, do nothing
				{	
					fclose(tmp);
					break;
				}
				else
				{
					if (fobj->VFLAG > 0)
						printf("nfs: warning: .fb file didn't match current job, overwriting\n");
					fclose(tmp);
					ggnfs_to_msieve(fobj, job);
					break;
				}
			}
		}
		mpz_clear(r);
		mpz_clear(num);
	}

	gmp_sprintf(obj->input, "%Zd", fobj->nfs_obj.gmp_n);
	gmp_printf("%s\n", obj->input); //mp_print(mpN, 10, NULL, gstr1.s));
	relations_needed = nfs_filter_relations(obj, fobj->nfs_obj.gmp_n);

	// reset the args list
	obj->nfs_args = NULL;

	return relations_needed;
}


void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj)
{
	int i;
	FILE *logfile;
	char c[4];

	// extract the factors
	for (i=0;i<factor_list->num_factors;i++)
	{
		mpz_t tmp;

		//init locals
		mpz_init(tmp);
		
		//convert the factor
		mp_t2gmp(&factor_list->final_factors[i]->factor,tmp);

		//divide it out
		mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, tmp);

		//check if its prime and log accordingly
		if (is_mpz_prp(tmp, fobj->NUM_WITNESSES))
		{
			//need to convert to yafu bigint to store
			add_to_factor_list(fobj->factors, tmp, fobj->VFLAG, fobj->NUM_WITNESSES, 0);
			strncpy(c,"prp",4);
		}
		else
		{
			add_to_factor_list(fobj->factors, tmp, fobj->VFLAG, fobj->NUM_WITNESSES, 0);
			strncpy(c,"C",2);
		}

        if (fobj->LOGFLAG)
        {
            logfile = fopen(fobj->flogname, "a");
            if (logfile == NULL)
            {
                printf("fopen error: %s\n", strerror(errno));
                printf("could not open yafu logfile for appending\n");
            }
            else
            {
                char* s = mpz_get_str(NULL, 10, tmp);
                logprint(logfile, "%s%d = %s\n", c, gmp_base10(tmp), s);
                fclose(logfile);
                free(s);
            }
        }

		//free locals
		mpz_clear(tmp);
	}

	//log anything left over
	if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) > 0) 
	{
		char c[4];

		if (is_mpz_prp(fobj->nfs_obj.gmp_n, fobj->NUM_WITNESSES))
		{
			add_to_factor_list(fobj->factors, fobj->nfs_obj.gmp_n, fobj->VFLAG, fobj->NUM_WITNESSES, 0);
			strncpy(c,"prp",4);			
		}
		else
		{
			add_to_factor_list(fobj->factors, fobj->nfs_obj.gmp_n, fobj->VFLAG, fobj->NUM_WITNESSES, 0);
			strncpy(c,"C",2);
		}
		
        if (fobj->LOGFLAG)
        {
            logfile = fopen(fobj->flogname, "a");
            if (logfile == NULL)
            {
                printf("fopen error: %s\n", strerror(errno));
                printf("could not open yafu logfile for appending\n");
            }
            else
            {
                char* s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
                logprint(logfile, "%s%d = %s\n", c,
                    gmp_base10(fobj->nfs_obj.gmp_n), s);
                fclose(logfile);
                free(s);
            }
        }

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
	}

	return;
}

#endif
