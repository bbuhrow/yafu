/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/
#ifndef _YECM_H_
#define _YECM_H_

#include "factor.h"
#include "ytools.h"
#include "arith.h"
#include <stdio.h>
#include <gmp_u64_xface.h>


/* produced using ecm -v -v -v for the various B1 bounds (default B2).
/	Thanks A. Schindel !
/
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */

static int ecm_levels[NUM_ECM_LEVELS] = {
    15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 };

static uint64_t ecm_std_b1[NUM_ECM_LEVELS] = { 2000, 11000, 50000, 250000, 1000000,
    3000000, 11000000, 43000000, 110000000, 260000000, 850000000, 2900000000 };

static double ecm_data[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
    /*t15, 2000,	*/	{30,		12,			7,			5,			3,			2,			2,			2,			2,		1,		1},
    /*t20, 11000,	*/	{844,		74,			21,			8,			5,			3,			2,			2,			2,		2,		1},
    /*t25, 50000,	*/	{58129,		1539,		214,		50,			20,			11,			7,			5,			4,		3,		3},
    /*t30, 250000,	*/	{6711967,	49962,		3288,		430,		118,		54,			26,			14,			10,		8,		6},
    /*t35, 1E+06,	*/	{1.20E+09,	2292278,	68422,		4914,		904,		322,		122,		54,			34,		23,		15},
    /*t40, 3E+06,	*/	{2.90E+12,	1.40E+08,	1849287,	70293,		8613,		2350,		681,		242,		135,	82,		47},
    /*t45, 11E+06,	*/	{9.00E+99,	1.10E+10,	6.10E+07,	1214949,	97057,		20265,		4480,		1263,		613,	333,	168},
    /*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	2.50E+09,	2.50E+07,	1270662,	199745,		33652,		7404,		3133,	1512,	661},
    /*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	1.30E+11,	5.90E+08,	1.90E+07,	2246256,	283939,		48714,		17769,	7643,	2865},
    /*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	5.80E+16,	1.60E+10,	3.20E+08,	2.80E+07,	2655154,	350439,		111196,	42017,	13611},
    /*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	8.20E+21,	2.70E+13,	6.10E+09,	4.00E+08,	2.70E+07,	2768535,	751771,	250214,	69408} };

/* produced using ecm -v -v -v for the various B1 bounds using param 0 (default B2).
/ version 7.0.5
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
static double ecm_data_param0[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
    /*t15, 2000,	*/	{34,		7,			3,			2,			2,			1,			1,			1,			1,		1,		1},
    /*t20, 11000,	*/	{1282,		86,		    21,			8,			5,			3,			2,			2,			2,		2,		1},
    /*t25, 50000,	*/	{95034,	    1864,		214,		50,			20,			11,			7,			5,			4,		3,		3},
    /*t30, 250000,	*/	{1.2e+07,	62295,		3283,		430,		118,		54,			26,			14,			10,		8,		6},
    /*t35, 1E+06,	*/	{2.20E+09,	2924742,	69076,		4911,		910,		324,		122,		55,			34,		23,		15},
    /*t40, 3E+06,	*/	{9.00E+99,	1.80E+08,	1847472,	70940,		8615,		2351,		686,		246,		135,	82,		47},
    /*t45, 11E+06,	*/	{9.00E+99,	1.50E+10,	6.20E+07,	1226976,	97096,		20272,		4482,		1286,		614,	335,	168},
    /*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	2.50E+09,	2.50E+07,	1281819,	201449,		33676,		7557,		3135,	1521,	661},
    /*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	1.30E+11,	5.80E+08,	1.90E+07,	2247436,	284176,		49831,		17884,	7650,	2867},
    /*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	5.90E+16,	1.60E+10,	3.10E+08,	2.80E+07,	2657998,	361851,		111314,	42057,	13623},
    /*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	8.40E+21,	2.70E+13,	6.20E+09,	3.90E+08,	2.70E+07,	2844041,	752662,	250476,	69471} };


/* produced using ecm -v -v -v for the various B1 bounds using param 1 (default B2).
/ version 7.0.5
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
static double ecm_data_param1[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
    /*t15, 2000,	*/	{43,		8,			4,			2,			2,			1,			1,			1,			1,		1,		1},
    /*t20, 11000,	*/	{1743,		107,		24,			9,			5,			4,			3,			2,			2,		2,		1},
    /*t25, 50000,	*/	{134769,	2402,		261,		58,			22,			13,			8,			5,			4,		3,		3},
    /*t30, 250000,	*/	{1.7e+07,	82576,		4108,		513,		137,		61,			29,			16,			11,		8,		6},
    /*t35, 1E+06,	*/	{3.30E+09,	3967858,	88265,		6022,		1071,		374,		138,		61,			38,		25,		17},
    /*t40, 3E+06,	*/	{9.00E+99,	2.50E+08,	2402639,	87544,		10283,		2753,		788,		278,		151,	91,		52},
    /*t45, 11E+06,	*/	{9.00E+99,	2.0E+10,	8.10E+07,	1534319,	118226,		24017,		5208,		1459,		692,	373,	185},
    /*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	3.30E+09,	3.10E+07,	1565171,	241048,		39497,		8704,		3583,	1709,	737},
    /*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	2.30E+11,	7.50E+08,	2.40E+07,	2713723,	336066,		57844,		20479,	8656,	3220},
    /*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	1.50E+17,	2.00E+10,	3.90E+08,	3.40E+07,	3167410,	419970,		128305,	47888,	15391},
    /*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	2.10E+22,	6.70E+13,	7.70E+09,	4.80E+08,	3.20E+07,	3346252,	872747,	288516,	78923} };


/* produced using ecm -v -v -v for the various B1 bounds using param 3 (default B2).
/ version 7.0.5
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
static double ecm_data_param3[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
    /*t15, 2000,	*/	{43,		8,			4,			2,			2,			1,			1,			1,			1,		1,		1},
    /*t20, 11000,	*/	{1743,		107,		24,			9,			5,			4,			3,			2,			2,		2,		1},
    /*t25, 50000,	*/	{134769,	2402,		261,		58,			22,			13,			8,			5,			4,		3,		3},
    /*t30, 250000,	*/	{1.7e+07,	82576,		4108,		513,		137,		61,			29,			16,			11,		8,		6},
    /*t35, 1E+06,	*/	{3.30E+09,	3967858,	88265,		6022,		1071,		374,		138,		61,			38,		25,		17},
    /*t40, 3E+06,	*/	{9.00E+99,	2.50E+08,	2402639,	87544,		10283,		2753,		788,		278,		151,	91,		52},
    /*t45, 11E+06,	*/	{9.00E+99,	2.0E+10,	8.10E+07,	1534319,	118226,		24017,		5208,		1459,		692,	373,	185},
    /*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	3.30E+09,	3.10E+07,	1565171,	241048,		39497,		8704,		3583,	1709,	737},
    /*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	2.30E+11,	7.50E+08,	2.40E+07,	2713723,	336066,		57844,		20479,	8656,	3220},
    /*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	1.50E+17,	2.00E+10,	3.90E+08,	3.40E+07,	3167410,	419970,		128305,	47888,	15391},
    /*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	2.10E+22,	6.70E+13,	7.70E+09,	4.80E+08,	3.20E+07,	3346252,	872747,	288516,	78923} };


/* produced using ecm -v -power 1 -param 0 for the various B1 bounds (B2=B1*100).
/  Thanks tbusby!
/
/                          2k       11k      50k      250k     1M       3M       11M      43M      110M     260M     850M     2900M */
static double avx_ecm_data[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
    /*t15, 2000,    */    {32,      8,       4,       3,       2,       2,       2,       1,       1,       1,       1,       1     },
    /*t20, 11000,   */    {1203,    98,      26,      11,      7,       5,       4,       3,       3,       2,       2,       2     },
    /*t25, 50000,   */    {88933,   2174,    281,     74,      32,      20,      12,      8,       7,       5,       4,       4     },
    /*t30, 250000,  */    {1.1E+07, 73653,   4455,    671,     206,     100,     50,      28,      20,      15,      11,      8     },
    /*t35, 1E+06,   */    {2.0E+09, 3495192, 95330,   7965,    1676,    635,     250,     115,     73,      51,      33,      23    },
    /*t40, 3E+06,   */    {9.0E+99, 2.2E+08, 2610023, 117616,  16589,   4858,    1492,    550,     308,     194,     111,     68    },
    /*t45, 11E+06,  */    {9.0E+99, 1.8E+10, 8.8e+07, 2088438, 194098,  43492,   10273,   3014,    1481,    835,     419,     228   },
    /*t50, 44E+06,  */    {9.0E+99, 9.0E+99, 3.7E+09, 4.4E+07, 2626303, 445887,  80191,   18579,   7942,    3995,    1748,    838   },
    /*t55, 110E+06, */    {9.0E+99, 9.0E+99, 2.7E+11, 1.0E+09, 4.1E+07, 5150239, 699132,  126904,  46946,   20992,   7961,    3352  },
    /*t60, 260E+06, */    {9.0E+99, 9.0E+99, 9.0E+99, 2.9E+10, 6.7E+08, 6.6e+07, 6727121, 949988,  302925,  119976,  39223,   14455 },
    /*t65, 850E+06, */    {9.0E+99, 9.0E+99, 9.0E+99, 9.0E+99, 1.4E+10, 9.7e+08, 7.1E+07, 7721832, 2115419, 739454,  207648,  66688 } };
/*t70, 260E+06,     {9.0E+99, 9.0E+99, 9.0E+99, 9.0E+99, 9.0E+99, 1.5e+10, 7.6E+08, 6.7E+07, 1.6E+07, 4880638, 1173260, 327240} */



/* ============================ interface to gmpecm/pm1/pp1 ============================ */
extern void pollard_loop(fact_obj_t* fobj);
extern void williams_loop(fact_obj_t* fobj);
extern int ecm_loop(fact_obj_t* fobj);


/* ============================ interface to avxecm ============================ */
extern void vec_ecm_main(fact_obj_t* fobj, uint32_t numcurves, uint64_t B1,
    uint64_t B2, int threads, int* numfactors, int verbose,
    int save_b1, uint32_t* curves_run);
extern void vecPP1(fact_obj_t* fobj);
extern void vecPM1(fact_obj_t* fobj);

/* ============================ interface to utilities and table data ============================ */
extern uint32_t get_curves_for_tlevel(int tlevel, int b1_method, int b2_method);
extern void get_ecm_method(fact_obj_t* fobj, uint64_t b1, int* b1_method, int* b2_method);
extern uint32_t get_curves_required(ecm_obj_t* ecm_obj, double target_tlevel, uint64_t b1, int b1_method, int b2_method);
extern void print_std_ecm_work_done(ecm_obj_t* ecm_obj, int disp_levels, FILE* log, int VFLAG, int LOGFLAG);
extern void record_curves_completed(ecm_obj_t* ecm_obj, int curves, uint64_t b1, int b1_method, int b2_method);
int print_B1B2(fact_obj_t* fobj, FILE* fid);

#endif // #ifndef _YECM_H_
