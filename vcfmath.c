#include <math.h>
#include "vcfmath.h"

static double g_q2p[256];

double *bcf_m_get_pdg3(const bcf_hdr_t *h, bcf1_t *b)
{
	double *pdg;
	int i, k, id_PL, id_GL;
	bcf_fmt_t *ptr;
	if (g_q2p[0] == 0.)
		for (i = 0; i < 256; ++i)
			g_q2p[i] = pow(10., -0.1 * i);
	id_PL = bcf_id2int(h, BCF_DT_ID, "PL");
	id_GL = bcf_id2int(h, BCF_DT_ID, "GL");
	if (id_PL < 0 && id_GL < 0) return 0;
	bcf_unpack(b, BCF_UN_FMT);
	pdg = malloc(b->n_sample * 3 * sizeof(double));
	for (i = 0; i < b->n_fmt; ++i) {
		if (b->d.fmt[i].id == id_PL || b->d.fmt[i].id == id_GL) {
			ptr = &b->d.fmt[i];
			break;
		}
	}
	if (ptr->type == BCF_BT_INT8) {
		int8_t *p = (int8_t*)ptr->p;
		for (i = k = 0; i < b->n_sample; ++i, p += ptr->size) // FIXME: check if *p < 0
			pdg[k++] = g_q2p[p[0]], pdg[k++] = g_q2p[p[1]], pdg[k++] = g_q2p[p[2]];
	} else if (ptr->type == BCF_BT_INT16) {
		int16_t *p = (int16_t*)ptr->p;
		for (i = k = 0; i < b->n_sample; ++i, p += ptr->size) {
			pdg[k++] = g_q2p[p[0]<255?p[0]:255];
			pdg[k++] = g_q2p[p[1]<255?p[1]:255];
			pdg[k++] = g_q2p[p[2]<255?p[2]:255];
		}
	} else if (ptr->type == BCF_BT_FLOAT) {
		float *p = (float*)ptr->p;
		for (i = k = 0; i < b->n_sample; ++i, p += ptr->size)
			pdg[k++] = p[0], pdg[k++] = p[1], pdg[k++] = p[2];
	}
	return pdg;
}
