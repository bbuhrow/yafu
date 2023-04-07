#include <arith.h>
#include <ytools.h>
#include <math.h>

#define pi 3.1415926535897932384626433832795

double * parsewords(mpz_t a, int b, int fftlen)
{
	// given mpz input a, make a double-float array with
	// b bits of the input per double-float.
	mpz_t aa;
	int s = mpz_sizeinbase(a, 2);
	int i = 0;
	int bits = 0;
	double* d = (double *)xmalloc(fftlen * sizeof(double));
	uint64_t mask = (1ULL << b) - 1;

	mpz_init(aa);
	mpz_set(aa, a);

	while (mpz_cmp_ui(aa, 0) > 0)
	{
		d[i++] = (double)(mpz_get_64(aa) & mask);
		mpz_tdiv_q_2exp(aa, aa, b);

		if (i > fftlen)
		{
			printf("input a is too long for input fftlen (%d words) @ %d bits/word\n",
				fftlen, b);
			exit(1);
		}
	}

	for (; i < fftlen; i++)
	{
		d[i] = 0.0;
	}

	mpz_clear(aa);
	return d;
}

double* parsewords_balanced(mpz_t a, int b, int fftlen)
{
	// given mpz input a, make a double-float array with
	// b bits of the input per double-float.
	mpz_t aa;
	int s = mpz_sizeinbase(a, 2);
	int i = 0;
	int bits = 0;
	double* d = (double*)xmalloc(fftlen * sizeof(double));
	int64_t W = (1ULL << b);
	int64_t mask = W - 1;
	int64_t b2 = W / 2;
	uint64_t* n = a->_mp_d;

	mpz_init(aa);
	mpz_set(aa, a);
	int64_t carry = 0;
	while ((mpz_cmp_ui(aa, 0) > 0) && (i < fftlen))
	{
		int64_t q = (mpz_get_64(aa) + carry) & mask;
		if ((q == 0) && carry)
		{
			d[i++] = 0.0;
			carry = 1;
		}
		else if (q >= b2)
		{
			d[i++] = (double)(q - W);
			carry = 1;
		}
		else
		{
			d[i++] = (double)q;
			carry = 0;
		}
		mpz_tdiv_q_2exp(aa, aa, b);

		if (i > fftlen)
		{
			printf("input a is too long for input fftlen (%d words) @ %d bits/word\n",
				fftlen, b);
			exit(1);
		}
	}
	if ((i >= fftlen) && carry)
	{
		printf("input a is too long for input fftlen (%d words) @ %d bits/word\n",
			fftlen, b);
		exit(1);
	}
	d[i++] = (double)carry;

	for (; i < fftlen; i++)
	{
		d[i] = 0.0;
	}

	mpz_clear(aa);
	return d;
}

double* fastparse_balanced(mpz_t a, int b, int fftlen)
{
	// given mpz input a, make a double-float array with
	// b bits of the input per double-float.
	int s = mpz_sizeinbase(a, 2);
	int i = 0;
	int bits = 0;
	double* d = (double*)xmalloc(fftlen * sizeof(double));
	int64_t W = (1ULL << b);
	int64_t mask = W - 1;
	int64_t b2 = W / 2;
	uint64_t* n = a->_mp_d;
	int64_t carry = 0;
	int j = 0;

	//printf("W,m = %u,%x\n", W, mask);
	//gmp_printf("a: %Zx\n", a);

	while (j < (s-b))
	{
		int x = j / 64;
		int y = j & 63;
		int64_t q = (n[x] >> y) & mask;

		//printf("in loop: j,x,y,s = %d,%d,%d,%d, q = %lx\n", j, x, y, s, q);

		if ((y + b) > 64)
			q += n[x + 1] << (64 - y);
		q = (q + carry) & mask;

		//int64_t q = (mpz_get_64(aa) + carry) & mask;
		if ((q == 0) && carry)
		{
			d[i++] = 0.0;
			carry = 1;
		}
		else if (q >= b2)
		{
			d[i++] = (double)(q - W);
			carry = 1;
		}
		else
		{
			d[i++] = (double)q;
			carry = 0;
		}
		//mpz_tdiv_q_2exp(aa, aa, b);

		if (i > fftlen)
		{
			printf("input a is too long for input fftlen (%d words) @ %d bits/word\n",
				fftlen, b);
			exit(1);
		}
		j += b;
	}

	int x = j / 64;
	int y = j & 63;

	//printf("last: j,x,y,s = %d,%d,%d,%d\n", j, x, y, s);

	int64_t q = ((n[x] >> y) + carry) & mask;

	//printf("q = %lx\n", q);

	if ((q == 0) && carry)
	{
		d[i++] = 0.0;
		carry = 1;
	}
	else if (q >= b2)
	{
		d[i++] = (double)(q - W);
		carry = 1;
	}
	else
	{
		d[i++] = (double)q;
		carry = 0;
	}

	if ((i >= fftlen) && carry)
	{
		printf("input a is too long for input fftlen (%d words) @ %d bits/word\n",
			fftlen, b);
		exit(1);
	}
	d[i++] = (double)carry;

	for (; i < fftlen; i++)
	{
		d[i] = 0.0;
	}

	return d;
}

__inline void cmplx_mul(double* rc, double* ic, double ra, double ia, double rb, double ib)
{
	// (ac - bd) + i(ad + bc)
	*rc = ra * rb - ia * ib;
	*ic = ra * ib + ia * rb;
	return;
}


void fftr(double* dar, double* dai, int fftlen, int invert)
{
	/*
	https://cp-algorithms.com/algebra/fft.html#implementation
	recursive implementation

	using cd = complex<double>;
	const double PI = acos(-1);

	void fft(vector<cd> &a, bool invert) {
		int n = a.size();
		if (n == 1)
			return;

		vector<cd> a0(n / 2), a1(n / 2);
		for (int i = 0; 2 * i < n; i++) {
			a0[i] = a[2 * i];
			a1[i] = a[2 * i + 1];
		}
		fft(a0, invert);
		fft(a1, invert);

		double ang = 2 * PI / n * (invert ? -1 : 1);
		cd w(1), wn(cos(ang), sin(ang));
		for (int i = 0; 2 * i < n; i++) {
			a[i] = a0[i] + w * a1[i];
			a[i + n / 2] = a0[i] - w * a1[i];
			if (invert) {
				a[i] /= 2;
				a[i + n / 2] /= 2;
			}
			w *= wn;
		}
	}

	*/
	int debug = 0;
	int i;

	if (fftlen == 1)
		return;

	double* a0r = xmalloc(fftlen / 2 * sizeof(double));
	double* a0i = xmalloc(fftlen / 2 * sizeof(double));
	double* a1r = xmalloc(fftlen / 2 * sizeof(double));
	double* a1i = xmalloc(fftlen / 2 * sizeof(double));

	for (i = 0; 2 * i < fftlen; i++) {
		a0r[i] = dar[2 * i];
		a1r[i] = dar[2 * i + 1];
		a0i[i] = dai[2 * i];
		a1i[i] = dai[2 * i + 1];
	}
	fftr(a0r, a0i, fftlen / 2, invert);
	fftr(a1r, a1i, fftlen / 2, invert);

	double ang = 2 * pi / fftlen * (invert ? -1 : 1);
	//cd w(1), wn(cos(ang), sin(ang));
	double wr = 1.0, wi = 0.0;
	double wnr = cos(ang);
	double wni = sin(ang);

	for (int i = 0; 2 * i < fftlen; i++) {
		//a[i] = a0[i] + w * a1[i];
		double tr, ti;
		cmplx_mul(&tr, &ti, wr, wi, a1r[i], a1i[i]);
		dar[i] = a0r[i] + tr;
		dai[i] = a0i[i] + ti;

		//a[i + n / 2] = a0[i] - w * a1[i];
		dar[i + fftlen / 2] = a0r[i] - tr;
		dai[i + fftlen / 2] = a0i[i] - ti;

		if (invert) {
			dar[i] /= 2;
			dar[i + fftlen / 2] /= 2;
			dai[i] /= 2;
			dai[i + fftlen / 2] /= 2;
		}
		//w *= wn;
		cmplx_mul(&tr, &ti, wr, wi, wnr, wni);
		wr = tr;
		wi = ti;
	}

	if (debug)
	{
		printf("transform: ");
		for (i = 0; i < fftlen; i++)
		{
			printf("%lf,%lf ", dar[i], dai[i]);
		}
		printf("\n");
	}

	free(a0r);
	free(a0i);
	free(a1r);
	free(a1i);


	/*
	// in-place iterative implementation
	using cd = complex<double>;
	const double PI = acos(-1);

	int reverse(int num, int lg_n) {
		int res = 0;
		for (int i = 0; i < lg_n; i++) {
			if (num & (1 << i))
				res |= 1 << (lg_n - 1 - i);
		}
		return res;
	}

	void fft(vector<cd> & a, bool invert) {
		int n = a.size();
		int lg_n = 0;
		while ((1 << lg_n) < n)
			lg_n++;

		for (int i = 0; i < n; i++) {
			if (i < reverse(i, lg_n))
				swap(a[i], a[reverse(i, lg_n)]);
		}

		for (int len = 2; len <= n; len <<= 1) {
			double ang = 2 * PI / len * (invert ? -1 : 1);
			cd wlen(cos(ang), sin(ang));
			for (int i = 0; i < n; i += len) {
				cd w(1);
				for (int j = 0; j < len / 2; j++) {
					cd u = a[i+j], v = a[i+j+len/2] * w;
					a[i+j] = u + v;
					a[i+j+len/2] = u - v;
					w *= wlen;
				}
			}
		}

		if (invert) {
			for (cd & x : a)
				x /= n;
		}
	}

	// without explicit reversing
	using cd = complex<double>;
	const double PI = acos(-1);

	void fft(vector<cd> & a, bool invert) {
		int n = a.size();

		for (int i = 1, j = 0; i < n; i++) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1)
				j ^= bit;
			j ^= bit;

			if (i < j)
				swap(a[i], a[j]);
		}

		for (int len = 2; len <= n; len <<= 1) {
			double ang = 2 * pi / len * (invert ? -1 : 1);
			cd wlen(cos(ang), sin(ang));
			for (int i = 0; i < n; i += len) {
				cd w(1);
				for (int j = 0; j < len / 2; j++) {
					cd u = a[i+j], v = a[i+j+len/2] * w;
					a[i+j] = u + v;
					a[i+j+len/2] = u - v;
					w *= wlen;
				}
			}
		}

		if (invert) {
			for (cd & x : a)
				x /= n;
		}
	}


	Additionally we can precompute the bit-reversal permutation beforehand.
	This is especially useful when the size n is the same for all calls.
	But even when we only have three calls (which are necessary for multiplying
	two polynomials), the effect is noticeable. Also we can precompute all
	roots of unity and their powers.



	*/



}

int reverse(int num, int lg_n) {
	int res = 0;
	for (int i = 0; i < lg_n; i++) {
		if (num & (1 << i))
			res |= 1 << (lg_n - 1 - i);
	}
	return res;
}

void ffti(double* dar, double* dai, int fftlen, int invert)
{
	/*
	https://cp-algorithms.com/algebra/fft.html#implementation

	// in-place iterative implementation
	using cd = complex<double>;
	const double PI = acos(-1);

	int reverse(int num, int lg_n) {
		int res = 0;
		for (int i = 0; i < lg_n; i++) {
			if (num & (1 << i))
				res |= 1 << (lg_n - 1 - i);
		}
		return res;
	}

	void fft(vector<cd> & a, bool invert) {
		int n = a.size();
		int lg_n = 0;
		while ((1 << lg_n) < n)
			lg_n++;

		for (int i = 0; i < n; i++) {
			if (i < reverse(i, lg_n))
				swap(a[i], a[reverse(i, lg_n)]);
		}

		for (int len = 2; len <= n; len <<= 1) {
			double ang = 2 * PI / len * (invert ? -1 : 1);
			cd wlen(cos(ang), sin(ang));
			for (int i = 0; i < n; i += len) {
				cd w(1);
				for (int j = 0; j < len / 2; j++) {
					cd u = a[i+j], v = a[i+j+len/2] * w;
					a[i+j] = u + v;
					a[i+j+len/2] = u - v;
					w *= wlen;
				}
			}
		}

		if (invert) {
			for (cd & x : a)
				x /= n;
		}
	}

	*/
	int debug = 0;
	int i;
	int lg_n = 0;

	while ((1 << lg_n) < fftlen)
		lg_n++;

	for (i = 0; i < fftlen; i++) {
		int ri = reverse(i, lg_n);
		if (i < ri)
		{
			//swap(a[i], a[reverse(i, lg_n)]);
			double tmpr = dar[i];
			double tmpi = dai[i];
			dar[i] = dar[ri];
			dai[i] = dai[ri];
			dar[ri] = tmpr;
			dai[ri] = tmpi;
		}
	}

	int len;
	for (len = 2; len <= fftlen; len <<= 1) {
		double ang = 2 * pi / len * (invert ? -1 : 1);
		//cd wlen(cos(ang), sin(ang));
		double wlenr = cos(ang);
		double wleni = sin(ang);

		for (i = 0; i < fftlen; i += len) {
			//cd w(1);
			double wr = 1.0, wi = 0.0;
			int j;
			for (j = 0; j < len / 2; j++) {
				//cd u = a[i + j], v = a[i + j + len / 2] * w;
				double ur = dar[i + j];
				double ui = dai[i + j];
				double vr;
				double vi;
				cmplx_mul(&vr, &vi, dar[i + j + len / 2], dai[i + j + len / 2], wr, wi);

				//a[i + j] = u + v;
				//a[i + j + len / 2] = u - v;
				dar[i + j] = ur + vr;
				dai[i + j] = ui + vi;
				dar[i + j + len / 2] = ur - vr;
				dai[i + j + len / 2] = ui - vi;

				//w *= wlen;
				cmplx_mul(&vr, &vi, wr, wi, wlenr, wleni);
				wr = vr;
				wi = vi;
			}
		}
	}

	if (invert) {
		//for (cd& x : a)
		//	x /= n;
		for (i = 0; i < fftlen; i++) {
			dar[i] /= fftlen;
			dai[i] /= fftlen;
		}
	}

	if (debug)
	{
		printf("transform: ");
		for (i = 0; i < fftlen; i++)
		{
			printf("%lf,%lf ", dar[i], dai[i]);
		}
		printf("\n");
	}


	/*
	

	// without explicit reversing
	using cd = complex<double>;
	const double PI = acos(-1);

	void fft(vector<cd> & a, bool invert) {
		int n = a.size();

		for (int i = 1, j = 0; i < n; i++) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1)
				j ^= bit;
			j ^= bit;

			if (i < j)
				swap(a[i], a[j]);
		}

		for (int len = 2; len <= n; len <<= 1) {
			double ang = 2 * pi / len * (invert ? -1 : 1);
			cd wlen(cos(ang), sin(ang));
			for (int i = 0; i < n; i += len) {
				cd w(1);
				for (int j = 0; j < len / 2; j++) {
					cd u = a[i+j], v = a[i+j+len/2] * w;
					a[i+j] = u + v;
					a[i+j+len/2] = u - v;
					w *= wlen;
				}
			}
		}

		if (invert) {
			for (cd & x : a)
				x /= n;
		}
	}


	Additionally we can precompute the bit-reversal permutation beforehand.
	This is especially useful when the size n is the same for all calls.
	But even when we only have three calls (which are necessary for multiplying
	two polynomials), the effect is noticeable. Also we can precompute all
	roots of unity and their powers.



	*/



}


const double tw_r[21] = { 
	-1.0, 									   // c(3.1415926535897932384626433832795/2^0)
	0.0, 									   // c(3.1415926535897932384626433832795/2^1)
	0.70710678118654752440084436210485,		   // c(3.1415926535897932384626433832795/2^2)
	0.92387953251128675612818318939679,		   // c(3.1415926535897932384626433832795/2^3)
	0.98078528040323044912618223613424,		   // c(3.1415926535897932384626433832795/2^4)
	0.99518472667219688624483695310948,		   // c(3.1415926535897932384626433832795/2^5)
	0.9987954562051723927147716047591,		   // c(3.1415926535897932384626433832795/2^6)
	0.99969881869620422011576564966617,		   // c(3.1415926535897932384626433832795/2^7)
	0.99992470183914454092164649119638,		   // c(3.1415926535897932384626433832795/2^8)
	0.99998117528260114265699043772857,		   // c(3.1415926535897932384626433832795/2^9)
	0.99999529380957617151158012570012,		   // c(3.1415926535897932384626433832795/2^10)
	0.99999882345170190992902571017153,		   // c(3.1415926535897932384626433832795/2^11)
	0.99999970586288221916022821773877,		   // c(3.1415926535897932384626433832795/2^12)
	0.99999992646571785114473148070739,		   // c(3.1415926535897932384626433832795/2^13)
	0.99999998161642929380834691540291,		   // c(3.1415926535897932384626433832795/2^14)
	0.99999999540410731289097193313961,		   // c(3.1415926535897932384626433832795/2^15)
	0.99999999885102682756267330779455,		   // c(3.1415926535897932384626433832795/2^16)
	0.99999999971275670684941397221864,		   // c(3.1415926535897932384626433832795/2^17)
	0.99999999992818917670977509588385,		   // c(3.1415926535897932384626433832795/2^18)
	0.99999999998204729417728262414778,		   // c(3.1415926535897932384626433832795/2^19)
	0.999999999995511823544310584173		   // c(3.1415926535897932384626433832795/2^20)
};

const double tw_i[21] = {
    0.0, 									   // s(3.1415926535897932384626433832795/2^0)
    1.0, 									   // s(3.1415926535897932384626433832795/2^1)
    0.707106781186547524400844362104,		   // s(3.1415926535897932384626433832795/2^2)
    0.382683432365089771728459984029,	       // s(3.1415926535897932384626433832795/2^3)
    0.195090322016128267848284868476,	       // s(3.1415926535897932384626433832795/2^4)
    0.098017140329560601994195563888,	       // s(3.1415926535897932384626433832795/2^5)
    0.049067674327418014254954976941,	       // s(3.1415926535897932384626433832795/2^6)
    0.024541228522912288031734529458,	       // s(3.1415926535897932384626433832795/2^7)
    0.012271538285719926079408261950,	       // s(3.1415926535897932384626433832795/2^8)
    0.006135884649154475359640234589,	       // s(3.1415926535897932384626433832795/2^9)
    0.003067956762965976270145365489,	       // s(3.1415926535897932384626433832795/2^10)
    0.001533980186284765612303697149,	       // s(3.1415926535897932384626433832795/2^11)
    0.000766990318742704526938568357,	       // s(3.1415926535897932384626433832795/2^12)
    0.000383495187571395589072461680,	       // s(3.1415926535897932384626433832795/2^13)
    0.000191747597310703307439909561,	       // s(3.1415926535897932384626433832795/2^14)
    0.000095873799095977345870517210,	       // s(3.1415926535897932384626433832795/2^15)
    0.000047936899603066884549003990,	       // s(3.1415926535897932384626433832795/2^16)
    0.000023968449808418218729186576,	       // s(3.1415926535897932384626433832795/2^17)
    0.000011984224905069706421521561,	       // s(3.1415926535897932384626433832795/2^18)
    0.000005992112452642427842879711,	       // s(3.1415926535897932384626433832795/2^19)
    0.000002996056226334660750454811	       // s(3.1415926535897932384626433832795/2^20)
};



void fft(double* dar, double* dai, int fftlen, int invert)
{
	/*
	https://cp-algorithms.com/algebra/fft.html#implementation

	// in-place iterative implementation
	// without explicit reversing
	using cd = complex<double>;
	const double PI = acos(-1);

	void fft(vector<cd> & a, bool invert) {
		int n = a.size();

		for (int i = 1, j = 0; i < n; i++) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1)
				j ^= bit;
			j ^= bit;

			if (i < j)
				swap(a[i], a[j]);
		}

		for (int len = 2; len <= n; len <<= 1) {
			double ang = 2 * pi / len * (invert ? -1 : 1);
			cd wlen(cos(ang), sin(ang));
			for (int i = 0; i < n; i += len) {
				cd w(1);
				for (int j = 0; j < len / 2; j++) {
					cd u = a[i+j], v = a[i+j+len/2] * w;
					a[i+j] = u + v;
					a[i+j+len/2] = u - v;
					w *= wlen;
				}
			}
		}

		if (invert) {
			for (cd & x : a)
				x /= n;
		}
	}

	*/
	int debug = 0;
	int i, j;
	
	for (i = 1, j = 0; i < fftlen; i++) {
		int bit = fftlen >> 1;
		for (; j & bit; bit >>= 1)
			j ^= bit;
		j ^= bit;

		if (i < j)
		{
			//swap(a[i], a[j]);
			double tmpr = dar[i];
			double tmpi = dai[i];
			dar[i] = dar[j];
			dai[i] = dai[j];
			dar[j] = tmpr;
			dai[j] = tmpi;
		}
	}

	int len;
	int l;
	for (len = 2, l = 0; len <= fftlen; len <<= 1, l++) {
		//double ang = 2 * pi / len * (invert ? -1 : 1);
		//double wlenr = cos(ang);
		//double wleni = sin(ang);

		//cd wlen(cos(ang), sin(ang));
		double wlenr = tw_r[l];
		double wleni = (invert ? -1.0 : 1.0) * tw_i[l];

		//printf("len = %d, l = %d, wlenr = %lf, wleni = %lf\n",
		//	len, l, wlenr, wleni);

		for (i = 0; i < fftlen; i += len) {
			//cd w(1);
			double wr = 1.0, wi = 0.0;
			for (j = 0; j < len / 2; j++) {
				//cd u = a[i + j], v = a[i + j + len / 2] * w;
				double ur = dar[i + j];
				double ui = dai[i + j];
				double vr;
				double vi;
				cmplx_mul(&vr, &vi, dar[i + j + len / 2], dai[i + j + len / 2], wr, wi);

				//a[i + j] = u + v;
				//a[i + j + len / 2] = u - v;
				dar[i + j] = ur + vr;
				dai[i + j] = ui + vi;
				dar[i + j + len / 2] = ur - vr;
				dai[i + j + len / 2] = ui - vi;

				//w *= wlen;
				cmplx_mul(&vr, &vi, wr, wi, wlenr, wleni);
				wr = vr;
				wi = vi;
			}
		}
	}

	if (invert) {
		//for (cd& x : a)
		//	x /= n;
		for (i = 0; i < fftlen; i++) {
			dar[i] /= fftlen;
			dai[i] /= fftlen;
		}
	}

	if (debug)
	{
		printf("transform: ");
		for (i = 0; i < fftlen; i++)
		{
			printf("%lf,%lf ", dar[i], dai[i]);
		}
		printf("\n");
	}


	/*


	


	Additionally we can precompute the bit-reversal permutation beforehand.
	This is especially useful when the size n is the same for all calls.
	But even when we only have three calls (which are necessary for multiplying
	two polynomials), the effect is noticeable. Also we can precompute all
	roots of unity and their powers.



	*/



}



//int64_t round(double d)
//{
//	if (d < 0)
//		return (int64_t)(d - 0.49) + 1;
//	else
//		return (int64_t)(d + 0.51);
//}

void reconstruct_int(mpz_t c, int64_t* ci, int bits_per_word, int fftlen)
{
	mpz_t base, t;
	int i;
	mpz_init(base);
	mpz_init(t);
	mpz_set_ui(c, 0);
	mpz_set_ui(base, 1);

	for (i = 0; i < fftlen; i++)
	{
		mpz_set_si(t, ci[i]);
		mpz_mul_2exp(t, t, bits_per_word * i);
		mpz_add(c, c, t);
	}

	if (mpz_cmp_ui(c, 0) < 0)
	{
		mpz_set_si(t, 1.0);
		mpz_mul_2exp(t, t, bits_per_word * i);
		mpz_add(c, c, t);
		mpz_sub_ui(c, c, 1);
	}
	return;
}

void fast_reconstruct_int(mpz_t c, int64_t* ci, int bits_per_word, int fftlen)
{
	int i = 0;
	int s = bits_per_word * fftlen;
	int j = 0;
	int b = bits_per_word;
	mpz_set_ui(c, 1);
	mpz_mul_2exp(c, c, bits_per_word * fftlen);
	uint64_t* n = c->_mp_d;

	while (j < (s - b))
	{
		int x = j / 64;
		int y = j & 63;

		n[x] += (ci[i] << j);

		if ((y + b) > 64)
			n[x + 1] = ci[i] >> (64 - y);
		i++;
	}

	//if (mpz_cmp_ui(c, 0) < 0)
	//{
	//	mpz_set_si(t, 1);
	//	mpz_mul_2exp(t, t, bits_per_word * i);
	//	mpz_add(c, c, t);
	//	mpz_sub_ui(c, c, 1);
	//}
	return;
}

void fftmul(mpz_t c, mpz_t a, mpz_t b, int bits_per_word, int fftlen)
{
	double* dar, * dbr, * dbi, * dai;
	int64_t *ci;
	int sza, szb;
	mpz_t cc;
	int i, q;
	int debug = 0;
	int test = 1;
	struct timeval start, stop;
	double ttest, tsetup;
	int iterations = 1024;
	gettimeofday(&start, NULL);

	mpz_init(cc);

	if (test)
	{
		mpz_t t, m;
		mpz_init(t);
		mpz_init(m);

		mpz_set(cc, a);
		mpz_set_ui(t, 1);
		mpz_mul_2exp(m, t, bits_per_word * fftlen);
		mpz_sub_ui(m, m, 1);
		if (debug)
		{
			gmp_printf("mod = %Zd\n", m);
		}

		for (i = 0; i < iterations; i++)
		{
			mpz_mul(cc, cc, b);
			mpz_mod(cc, cc, m);
			if (debug)
			{
				gmp_printf("iteration %d: %Zd\n", i, cc);
			}
		}

		mpz_clear(t);

		gettimeofday(&stop, NULL);
		ttest = ytools_difftime(&start, &stop);
		printf("test took %lf seconds\n", ttest);
	}

	ci = (int64_t*)xmalloc(fftlen * sizeof(int64_t));

	if (debug)
	{
		gmp_printf("a: %Zd\n", a);
	}

	dar = fastparse_balanced(a, bits_per_word, fftlen);
	
	if (debug)
	{
		for (i = 0; i < fftlen; i++)
		{
			printf("%lf ", dar[i]);
		}
		printf("\n");

		gmp_printf("b: %Zd\n", b);
	}

	dbr = fastparse_balanced(b, bits_per_word, fftlen);

	if (debug)
	{
		for (i = 0; i < fftlen; i++)
		{
			printf("%lf ", dbr[i]);
		}
		printf("\n");
	}

	dai = (double*)xmalloc(fftlen * sizeof(double));
	dbi = (double*)xmalloc(fftlen * sizeof(double));

	for (i = 0; i < fftlen; i++)
	{
		dai[i] = dbi[i] = 0.0;
	}

	gettimeofday(&stop, NULL);

	tsetup = ytools_difftime(&start, &stop) - ttest;

	printf("setup took %lf seconds\n", tsetup);

	fft(dbr, dbi, fftlen, 0);

	for (q = 0; q < iterations; q++)
	{
		fft(dar, dai, fftlen, 0);

		// multiply
		for (i = 0; i < fftlen; i++)
		{
			// (ac - bd) + i(ad + bc)
			double tr, ti;
			cmplx_mul(&tr, &ti, dar[i], dai[i], dbr[i], dbi[i]);
			dar[i] = tr;
			dai[i] = ti;
		}

		if (debug)
		{
			printf("multiply ta*tb: ");
			for (i = 0; i < fftlen; i++)
			{
				printf("%lf,%lf ", dar[i], dai[i]);
			}
			printf("\n");
		}

		fft(dar, dai, fftlen, 1);

		if (debug)
		{
			printf("inverse transform: ");
			for (i = 0; i < fftlen; i++)
			{
				printf("%lf,%lf ", dar[i], dai[i]);
			}
			printf("\n");
		}

		for (i = 0; i < fftlen; i++)
		{
			ci[i] = round(dar[i]);
		}

		if (debug)
		{
			printf("rounded: ");
			for (i = 0; i < fftlen; i++)
			{
				printf("%ld ", ci[i]);
			}
			printf("\n");
		}

		// carry prop balanced
		int64_t carry = 0;
		int64_t W = (1ULL << bits_per_word);
		int64_t mask = W - 1;
		int done = 0;
		int nc = 0;
		while (!done)
		{
			done = 1;
			for (i = 0; i < fftlen; i++)
			{
				int64_t q = ci[i];
				ci[i] = (ci[i] + carry) % W;
				carry = (q + carry) / W;
			}

			if (debug)
			{
				printf("carry: ");
				for (i = 0; i < fftlen; i++)
				{
					printf("%ld ", ci[i]);
				}
				printf("%ld\n", carry);
			}

			if (carry != 0)
			{
				done = 0;
			}
			//else if (ci[fftlen - 1] < ( -W / 2) )
			//{
			//	ci[fftlen - 1] += W;
			//	carry = -1;
			//	done = 0;
			//}
			//else if (ci[fftlen - 1] >= (W/2))
			//{
			//	ci[fftlen - 1] -= W;
			//	carry = 1;
			//	done = 0;
			//}

			nc++;
		}

		int64_t p = 0;
		carry = 0;
		int64_t b2 = (1ULL << (bits_per_word - 1));

		for (i = 0; i < fftlen; i++)
		{
			dar[i] = (double)ci[i];

			//p = (ci[i] + carry);
			//
			//if ((p == 0) && carry)
			//{
			//	dar[i] = 0.0;
			//	//carry = 1;
			//}
			//else if (p >= b2)
			//{
			//	dar[i] = (double)(p - W);
			//	carry = 1;
			//}
			//else if (p < -b2)
			//{
			//	dar[i] = (double)(p + W);
			//	carry = -1;
			//}
			//else
			//{
			//	dar[i] = (double)p;
			//	carry = 0;
			//}
		}
		//if (carry)
		//	dar[i-1] += carry*W;

		if (debug)
		{
			printf("convert: ");
			for (i = 0; i < fftlen; i++)
			{
				printf("%lf ", dar[i]);
			}
			printf("\n");
			//exit(1);

			reconstruct_int(c, ci, bits_per_word, fftlen);
			gmp_printf("c = %Zd\n", c);
		}
	}

	gettimeofday(&stop, NULL);

	printf("fftmul took %lf seconds\n", ytools_difftime(&start, &stop) - ttest - tsetup);

	reconstruct_int(c, ci, bits_per_word, fftlen);

	gettimeofday(&stop, NULL);

	double tt = ytools_difftime(&start, &stop);

	printf("total time: %lf seconds\n", tt);

	if (test)
	{
		mpz_t t;
		mpz_init(t);
		mpz_sub(t, c, cc);
		if (mpz_cmp_ui(t, 0) != 0)
		{
			gmp_printf("ans is incorrect\n"); // should be : \nref = % Zd\n", cc);

			if (debug)
			{
				gmp_printf("ref = % Zd\n", cc);
			}
		}
		mpz_clear(t);
	}

	free(dar);
	free(dai);
	free(dbr);
	free(dbi);
	free(ci);
	mpz_clear(cc);
	return;
}























