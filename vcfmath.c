#include <math.h>
#include <assert.h>
#include "vcfmath.h"

static double g_q2p[256];

double *bcf_m_get_pdg3(const bcf_hdr_t *h, bcf1_t *b, const char *fmt)
{
	double *pdg;
	int i, k, id_PL;
	bcf_fmt_t *ptr = 0;
	if (g_q2p[0] == 0.)
		for (i = 0; i < 256; ++i)
			g_q2p[i] = pow(10., -0.1 * i);
	id_PL = bcf_id2int(h, BCF_DT_ID, fmt);
	if (id_PL < 0) return 0;
	bcf_unpack(b, BCF_UN_FMT);
	pdg = malloc(b->n_sample * 3 * sizeof(double));
	for (i = 0; i < b->n_fmt; ++i) {
		if (b->d.fmt[i].id == id_PL) {
			ptr = &b->d.fmt[i];
			break;
		}
	}
	if (ptr->type == BCF_BT_INT8) {
		int8_t *p = (int8_t*)ptr->p;
		for (i = k = 0; i < b->n_sample; ++i, p += ptr->size) // FIXME: check if *p < 0
			pdg[k++] = g_q2p[p[2]], pdg[k++] = g_q2p[p[1]], pdg[k++] = g_q2p[p[0]];
	} else if (ptr->type == BCF_BT_INT16) {
		int16_t *p = (int16_t*)ptr->p;
		for (i = k = 0; i < b->n_sample; ++i, p += ptr->size) {
			pdg[k++] = g_q2p[p[2]<255?p[2]:255];
			pdg[k++] = g_q2p[p[1]<255?p[1]:255];
			pdg[k++] = g_q2p[p[0]<255?p[0]:255];
		}
	} else if (ptr->type == BCF_BT_FLOAT) {
		float *p = (float*)ptr->p;
		for (i = k = 0; i < b->n_sample; ++i, p += ptr->size)
			pdg[k++] = p[2], pdg[k++] = p[1], pdg[k++] = p[1];
	}
	return pdg;
}

#define TINY 1e-20

int bcf_m_lk2(int n, const double *pdg3, double *y, const uint8_t *_ploidy)
{
	double *z[2], *tmp, *yswap;
	int j, M, last_min, last_max;
	const uint8_t *ploidy = 0;
	const double *pdg;

	if (_ploidy) {
		for (j = M = 0; j < n; ++j) {
			assert(_ploidy[j] > 0 && _ploidy[j] <= 2);
			M += _ploidy[j];
		}
		if (M != n<<1) ploidy = _ploidy;
	} else M = n<<1;
	z[0] = y;
	memset(y, 0, sizeof(double) * (M + 1));
	z[1] = yswap = calloc(M + 1, sizeof(double));
	z[0][0] = 1.;
	last_min = last_max = 0;
	for (j = M = 0, pdg = pdg3; j < n; ++j, pdg += 3) {
		int _min = last_min, _max = last_max;
		long k, M0;
		for (; _min < _max && z[0][_min] < TINY; ++_min) z[0][_min] = z[1][_min] = 0.;
		for (; _max > _min && z[0][_max] < TINY; --_max) z[0][_max] = z[1][_max] = 0.;
		M0 = M; M += ploidy? ploidy[j] : 2;
		if (ploidy && ploidy[j] == 1) { // haploid
			double p[2], sum;
			p[0] = pdg[0]; p[1] = pdg[2];
			++_max;
			if (_min == 0) k = 0, z[1][k] = (M0+1-k) * p[0] * z[0][k];
			for (k = _min < 1? 1 : _min; k <= _max; ++k)
				z[1][k] = (M0+1-k) * p[0] * z[0][k] + k * p[1] * z[0][k-1];
			for (k = _min, sum = 0.; k <= _max; ++k) sum += z[1][k];
			for (k = _min; k <= _max; ++k) z[1][k] /= sum;
			if (_min >= 1) z[1][_min-1] = 0.;
			if (j < n - 1) z[1][_max+1] = 0.;
		} else { // diploid
			double p[3], sum;
			p[0] = pdg[0]; p[1] = 2 * pdg[1]; p[2] = pdg[2];
			_max += 2;
			if (_min == 0) k = 0, z[1][k] = (M0-k+1) * (M0-k+2) * p[0] * z[0][k];
			if (_min <= 1) k = 1, z[1][k] = (M0-k+1) * (M0-k+2) * p[0] * z[0][k] + k*(M0-k+2) * p[1] * z[0][k-1];
			for (k = _min < 2? 2 : _min; k <= _max; ++k)
				z[1][k] = (M0-k+1)*(M0-k+2) * p[0] * z[0][k] + k*(M0-k+2) * p[1] * z[0][k-1] + k*(k-1) * p[2] * z[0][k-2];
			for (k = _min, sum = 0.; k <= _max; ++k) sum += z[1][k];
			for (k = _min; k <= _max; ++k) z[1][k] /= sum;
			if (_min >= 1) z[1][_min-1] = 0.;
			if (_min >= 2) z[1][_min-2] = 0.;
			if (j < n - 1) z[1][_max+1] = z[1][_max+2] = 0.;
		}
		tmp = z[0]; z[0] = z[1]; z[1] = tmp;
		last_min = _min; last_max = _max;
	}
	if (z[0] != y) memcpy(y, z[0], sizeof(double) * (M + 1));
	free(yswap);
	return M;
}

void bcf_m_prior(int type, double theta, int M, double *phi)
{
	int i;
	if (type == BCF_MP_COND2) {
		for (i = 0; i <= M; ++i)
			phi[i] = 2. * (i + 1) / (M + 1) / (M + 2);
	} else if (type == BCF_MP_FLAT) {
		for (i = 0; i <= M; ++i)
			phi[i] = 1. / (M + 1);
	} else {
		double sum;
		for (i = 0, sum = 0.; i < M; ++i)
			sum += (phi[i] = theta / (M - i));
		phi[M] = 1. - sum;
	}
}

