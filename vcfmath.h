#ifndef VCFMATH_H
#define VCFMATH_H

#include "vcf.h"

#define BCF_MP_FILE  0
#define BCF_MP_FLAT  1
#define BCF_MP_COND2 2
#define BCF_MP_WF    3

double *bcf_m_get_pdg3(const bcf_hdr_t *h, bcf1_t *b);
int bcf_m_lk2(int n, const double *pdg3, double *y, const uint8_t *_ploidy);
void bcf_m_prior(int type, double theta, int M, double *phi);

#endif
