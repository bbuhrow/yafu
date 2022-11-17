@
Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@<Candidate search@>=

#if defined( ASM_SEARCH0) && !defined(AVX512_SIEVE_SEARCH)
        {
            unsigned char* srbs;
            u32_t i;
            srbs = sieve_report_bounds[s][j_offset / CANDIDATE_SEARCH_STEPS];
            ncand = lasieve_search0(sieve_interval, horizontal_sievesums,
                horizontal_sievesums + j_per_strip,
                srbs, srbs + n_i / CANDIDATE_SEARCH_STEPS,
                cand, fss_sv);
            for (i = 0; i < ncand; i++)fss_sv[i] += horizontal_sievesums[cand[i] >> i_bits];
    }
#else
        {
            unsigned char* srbs;
            u32_t i;

            __m512i vincr = _mm512_set_epi16(
                31, 30, 29, 28, 27, 26, 25, 24,
                23, 22, 21, 20, 19, 18, 17, 16,
                15, 14, 13, 12, 11, 10, 9, 8,
                7, 6, 5, 4, 3, 2, 1, 0);
            __m512i vincr32 = _mm512_set1_epi16(32);

            srbs = sieve_report_bounds[s][j_offset / CANDIDATE_SEARCH_STEPS];
            ncand = 0;
            for (i = 0; i < n_i; i += CANDIDATE_SEARCH_STEPS) {
                unsigned char st;
                u32_t j;

                st = *(srbs++);
                for (j = 0; j < j_per_strip; j++) {
                    unsigned char* i_o, * i_max, st1;

                    i_o = sieve_interval + (j << i_bits) + i;
                    i_max = i_o + CANDIDATE_SEARCH_STEPS;
                    if (st <= horizontal_sievesums[j]) {
                        // in this condition we load all
                        // CANDIDATE_SEARCH_STEPS into cand/fss_sv.

#if defined(AVX512_SIEVE_SEARCH) && (CANDIDATE_SEARCH_STEPS == 128)
                        __m512i vidx = _mm512_add_epi16(vincr,
                            _mm512_set1_epi16((j << i_bits) + i));
                        _mm512_store_si512(&cand[ncand], vidx);
                        vidx = _mm512_add_epi16(vidx, vincr32);
                        _mm512_store_si512(&cand[ncand + 32], vidx);
                        vidx = _mm512_add_epi16(vidx, vincr32);
                        _mm512_store_si512(&cand[ncand + 64], vidx);
                        vidx = _mm512_add_epi16(vidx, vincr32);
                        _mm512_store_si512(&cand[ncand + 96], vidx);
                        __m512i v_io = _mm512_load_si512(i_o);
                        _mm512_store_si512(&fss_sv[ncand],
                            _mm512_add_epi8(v_io,
                                _mm512_set1_epi8(horizontal_sievesums[j])));
                        v_io = _mm512_load_si512(i_o + 64);
                        _mm512_store_si512(&fss_sv[ncand + 64],
                            _mm512_add_epi8(v_io,
                                _mm512_set1_epi8(horizontal_sievesums[j])));
                        ncand += 128;
#else
                        while (i_o < i_max) {
                            cand[ncand] = i_o - sieve_interval;
                            fss_sv[ncand++] = *(i_o++) + horizontal_sievesums[j];
                        }
#endif
                        continue;
                    }
                    st1 = st - horizontal_sievesums[j];

#if defined(AVX512_SIEVE_SEARCH) && (CANDIDATE_SEARCH_STEPS == 128)
#define HAVE_SSIMD
#endif

#ifndef HAVE_SSIMD
#ifdef GNFS_CS32 
#define bc_t unsigned long
#define BC_MASK 0x80808080
#else
#define bc_t unsigned long long
#define BC_MASK 0x8080808080808080
#endif
                    {
                        if (st1 < 0x80) {
                            bc_t bc, * i_oo;

                            bc = st1;
                            bc = (bc << 8) | bc;
                            bc = (bc << 16) | bc;
#ifndef GNFS_CS32
                            bc = (bc << 32) | bc;
#endif
                            bc = BC_MASK - bc;
                            for (i_oo = (bc_t*)i_o; i_oo < (bc_t*)i_max; i_oo++) {
                                bc_t v = *i_oo;
                                if (((v & BC_MASK) | ((v + bc) & BC_MASK)) == 0)continue;
                                for (i_o = (unsigned char*)i_oo; i_o < (unsigned char*)(i_oo + 1); i_o++) {
                                    if (*i_o >= st1) {
                                        cand[ncand] = i_o - sieve_interval;
                                        fss_sv[ncand++] = *i_o + horizontal_sievesums[j];

                                    }
                                }
                            }
                        }
                        else {
                            bc_t* i_oo;

                            for (i_oo = (bc_t*)i_o; i_oo < (bc_t*)i_max; i_oo++) {
                                if ((*i_oo & BC_MASK) == 0)continue;
                                for (i_o = (unsigned char*)i_oo; i_o < (unsigned char*)(i_oo + 1); i_o++) {
                                    if (*i_o >= st1) {
                                        cand[ncand] = i_o - sieve_interval;
                                        fss_sv[ncand++] = *i_o + horizontal_sievesums[j];

                                    }
                                }
                            }
                        }
                    }
#else

#if defined(AVX512_SIEVE_SEARCH) && (CANDIDATE_SEARCH_STEPS == 128)

                    {
                        // removed the looping as only 2 iterations are
                        // needed as long as CANDIDATE_SEARCH_STEPS = 128.

                        {
                            // we don't actually need any of the
                            // max functions when we can directly
                            // compare bytes to st1 (with AVX512_BW)
                            __m512i x;
                            __m512i vi_o1 = _mm512_load_si512(i_o);
                            __m512i vi_o2 = _mm512_load_si512(i_o + 64);

                            x = _mm512_set1_epi8(st1);

                            __mmask64 m1 = _mm512_cmpge_epu8_mask(vi_o1, x);
                            __mmask64 m2 = _mm512_cmpge_epu8_mask(vi_o2, x);

                            // search the sparse set bits for 
                            // candidate locations.
                            u64_t ma = m1;
                            while (ma > 0)
                            {
                                int id = _tzcnt_u64(ma);
                                cand[ncand] = i_o + id - sieve_interval;
                                fss_sv[ncand++] = i_o[id] + horizontal_sievesums[j];

                                ma = _blsr_u64(ma);
                            }

                            // search the sparse set bits for 
                            // candidate locations.
                            u64_t mb = m2;
                            while (mb > 0)
                            {
                                int id = _tzcnt_u64(mb);
                                cand[ncand] = i_o + 64 + id - sieve_interval;
                                fss_sv[ncand++] = i_o[64 + id] + horizontal_sievesums[j];

                                mb = _blsr_u64(mb);
                            }

                            i_o += 128;
                        }
                    }
#else

                    {
                        unsigned long long x;

                        x = st1 - 1;
                        x |= x << 8;
                        x |= x << 16;
                        x |= x << 32;
                        while (i_o < i_max) {
                            asm volatile("movq (%%eax),%%mm7\n"
                                "1:\n"
                                "movq (%%esi),%%mm1\n"
                                "movq 8(%%esi),%%mm0\n"
                                "pmaxub 16(%%esi),%%mm1\n"
                                "pmaxub 24(%%esi),%%mm0\n"
                                "pmaxub %%mm7,%%mm1\n"
                                "pmaxub %%mm1,%%mm0\n"
                                "pcmpeqb %%mm7,%%mm0\n"
                                "pmovmskb %%mm0,%%eax\n"
                                "cmpl $255,%%eax\n"
                                "jnz 2f\n"
                                "leal 32(%%esi),%%esi\n"
                                "cmpl %%esi,%%edi\n"
                                "ja 1b\n"
                                "2:\n"
                                "emms":"=S"(i_o) : "a"(&x), "S"(i_o),
                                "D"(i_max));
                            if (i_o < i_max) {
                                unsigned char* i_max2 = i_o + 32;
                                while (i_o < i_max2) {
                                    if (*i_o >= st1) {
                                        cand[ncand] = i_o - sieve_interval;
                                        fss_sv[ncand++] = *i_o + horizontal_sievesums[j];

                                    }
                                    i_o++;
                                }
                            }
                        }
                    }
#endif


#endif

                    }
                }
            }
#endif

#if 0
{
  char *ofn;
  FILE  *of;
  asprintf(&ofn,"cdump.%u.%u.j%u.ot%u",special_q,r[root_no],
	   j_offset,oddness_type);
  if((of=fopen(ofn,"w")) != NULL) {
    u32_t i;
    fprintf(of,"%u candidates\n",ncand);
    for(i=0;i<ncand;i++)
      fprintf(of,"%u %u\n",cand[i],fss_sv[i]);
    fclose(of);
  } else errprintf("Cannot open debug file %s: %m\n",ofn);
  free(ofn);
}
#endif

@
@<MMX Candidate searcher@>=
#ifndef HAVE_SSIMD
#ifdef GNFS_CS32 /* Use 32 bit registers for candidate search. */
#define bc_t unsigned long
#define BC_MASK 0x80808080
#else
#define bc_t unsigned long long
#define BC_MASK 0x8080808080808080
#endif
{
  if(st1<0x80) {
    bc_t bc,*i_oo;

    bc=st1;
    bc=(bc<<8)|bc;
    bc=(bc<<16)|bc;
#ifndef GNFS_CS32
    bc=(bc<<32)|bc;
#endif
    bc=BC_MASK-bc;
    for(i_oo=(bc_t*)i_o;i_oo<(bc_t*)i_max;i_oo++) {
      bc_t v=*i_oo;
      if( ( (v & BC_MASK) | ((v+bc) & BC_MASK)) == 0) continue;
      for(i_o=(unsigned char *)i_oo;i_o<(unsigned char *)(i_oo+1);i_o++) {
	if(*i_o>=st1) {
	  @<Store survivor@>@;
	}
      }
    }
  } else {
    bc_t *i_oo;

    for(i_oo=(bc_t*)i_o;i_oo<(bc_t*)i_max;i_oo++) {
      if( (*i_oo & BC_MASK) == 0) continue;
      for(i_o=(unsigned char *)i_oo;i_o<(unsigned char *)(i_oo+1);i_o++) {
	if(*i_o>=st1) {
	  @<Store survivor@>@;
	}
      }
    }
  }
}
#else
{
  unsigned long long x;

  x=st1-1;
  x|=x<<8;
  x|=x<<16;
  x|=x<<32;
  while(i_o<i_max) {
    asm volatile ("movq (%%eax),%%mm7\n"
		  "1:\n"
		  "movq (%%esi),%%mm1\n"
		  "movq 8(%%esi),%%mm0\n"
		  "pmaxub 16(%%esi),%%mm1\n"
		  "pmaxub 24(%%esi),%%mm0\n"
		  "pmaxub %%mm7,%%mm1\n"
		  "pmaxub %%mm1,%%mm0\n"
		  "pcmpeqb %%mm7,%%mm0\n"
		  "pmovmskb %%mm0,%%eax\n"
		  "cmpl $255,%%eax\n"
		  "jnz 2f\n"
		  "leal 32(%%esi),%%esi\n"
		  "cmpl %%esi,%%edi\n"
		  "ja 1b\n"
		  "2:\n"
		  "emms": "=S" (i_o) : "a" (&x), "S" (i_o),
		  "D" (i_max));
    if(i_o<i_max) {
      unsigned char *i_max2=i_o+32;
      while(i_o<i_max2) {
	if(*i_o>=st1) {
	  @<Store survivor@>@;
	}
	i_o++;
      }
    }
  }
}
#endif

@
@<Store survivor@>=
cand[ncand]=i_o-sieve_interval;
fss_sv[ncand++]=*i_o+horizontal_sievesums[j];

@
@<Final candidate search@>=
{
  u32_t i,nc1;
  unsigned char *srbs;
  static u32_t bad_pvl=0;
  double sr_inv;

  srbs=sieve_report_bounds[s][j_offset/CANDIDATE_SEARCH_STEPS];
  n_prereports+=ncand;
  if (ncand) sr_inv=1./(M_LN2*sieve_multiplier[s]);
  for(i=0,nc1=0;i<ncand;i++) {
    u16_t st_i,t_j,ii,jj,j;
    double pvl;

    j=cand[i]>>i_bits;
#ifndef DEBUG_SIEVE_REPORT_BOUNDS
    if(sieve_interval[cand[i]]+horizontal_sievesums[j]<
       srbs[(cand[i]&(n_i-1))/CANDIDATE_SEARCH_STEPS])
      continue;
#endif
    jj=j_offset+j;
    ii=cand[i]&(n_i-1);
    st_i=2*ii+ ( oddness_type==2 ? 0 : 1 );
    t_j=2*jj+ ( oddness_type==1 ? 0 : 1 );
#if 1
    pvl=log(fabs(rpol_eval(tpoly_f[s],poldeg[s],
			   (double)st_i-(double)i_shift,(double)t_j)));
#else
    pvl=log(fabs(rpol_eval0(tpoly_f[s],poldeg[s],
                            (i32_t)st_i-(i32_t)i_shift,t_j)));
#endif
    if(special_q_side == s)
      pvl-=special_q_log;
    pvl*=sieve_multiplier[s];
    pvl-=sieve_report_multiplier[s]*FB_maxlog[s];
    if((double)(sieve_interval[cand[i]]+horizontal_sievesums[j]) >=pvl) {
/*
In fss_sv2 we save an approximation of the number of bits of the cofactor:
*/
      pvl+=sieve_report_multiplier[s]*FB_maxlog[s];
      pvl-=(double)(sieve_interval[cand[i]]+horizontal_sievesums[j]);
      if (pvl<0.) pvl=0.;
      pvl*=sr_inv;
//      pvl/=(M_LN2*sieve_multiplier[s]);
//      pvl/=M_LN2;
      fss_sv2[nc1]=(unsigned char)(pvl);

#ifdef DEBUG_SIEVE_REPORT_BOUNDS
      @<Test correctness of sieve report bounds@>@;
#endif
      fss_sv[nc1]=fss_sv[i];
      cand[nc1++]=cand[i];
    }
  }
  rpol_eval_clear();
  ncand=nc1;
}

@
@<Test correctness of sieve report bounds@>=
if(sieve_interval[cand[i]]+horizontal_sievesums[j]<
   srbs[(cand[i]&(n_i-1))/CANDIDATE_SEARCH_STEPS]) {
     double pvl1;

     pvl=fabs(rpol_eval(tpoly_f[s],poldeg[s],
			(double)st_i-(double)i_shift,(double)t_j));
     fprintf(stderr,"Bad pvl min %u at (%f,%f),spq=%u\npvl: %.5g->",
	     bad_pvl++,(double)st_i-(double)i_shift,(double)t_j,
	     special_q,pvl);
     pvl=log(pvl);
     fprintf(stderr,"%.3f->",pvl);
     pvl=sieve_multiplier[s]*pvl;
     fprintf(stderr,"%.3f->",pvl);
     if(special_q_side == s) pvl-=sieve_multiplier[s]*special_q_log;
     fprintf(stderr,"%.3f->",pvl);
     pvl-=sieve_report_multiplier[s]*FB_maxlog[s];
     fprintf(stderr,"%.3f\nLower bound was %u sv was %u=%u+%u\n",pvl,
	     (u32_t)srbs[(cand[i]&(n_i-1))/CANDIDATE_SEARCH_STEPS],
	     (u32_t)sieve_interval[cand[i]]+(u32_t)horizontal_sievesums[j],
	     (u32_t)sieve_interval[cand[i]],(u32_t)horizontal_sievesums[j]);
   }
