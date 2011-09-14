#include "nfs.h"


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
	z tmp1, tmp2, *N = &fobj->nfs_obj.n;
	int i;
	FILE *logfile;
	char c[3];

	// extract the factors
	zInit(&tmp1);
	zInit(&tmp2);
	zCopy(N,&tmp1);
	for (i=0;i<factor_list->num_factors;i++)
	{
		z tmpz;
		zInit(&tmpz);				
		
		mp_t2z(&factor_list->final_factors[i]->factor,&tmpz);

		zDiv(&tmp1,&tmpz,N,&tmp2);
		if (isPrime(&tmpz))
		{
			tmpz.type = PRP;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"prp");
		}
		else
		{
			tmpz.type = COMPOSITE;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"C");
		}

		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "%s%d = %s\n",c,ndigits(&tmpz),z2decstr(&tmpz,&gstr1));
			fclose(logfile);
		}		

		zFree(&tmpz);		
		zCopy(N,&tmp1);
	}

	if (zCompare(&tmp1,&zOne) > 0)
	{
		char c[3];
		if (isPrime(&tmp1))
		{
			tmp1.type = PRP;
			add_to_factor_list(fobj, &tmp1);
			strcpy(c,"prp");
			zCopy(&zOne, N);
		}
		else
		{
			tmp1.type = COMPOSITE;
			add_to_factor_list(fobj, &tmp1);
			strcpy(c,"C");
			zCopy(&zOne, N);
		}

		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "%s%d = %s\n",c,ndigits(&tmp1),z2decstr(&tmp1,&gstr1));
			fclose(logfile);
		}		
	}

	zFree(&tmp1);
	zFree(&tmp2);

	return;
}

