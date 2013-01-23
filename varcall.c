#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include "vcfmath.h"
#include "vcf.h"

#define FLAG_VCFIN   1
#define FLAG_BCFOUT  2
#define FLAG_CALL    4
#define FLAG_NOINDEL 8
#define FLAG_VARONLY 16

int main(int argc, char *argv[])
{
	int i, c, clevel = -1, flag = 0, n_samples = -1, *imap = 0, prior_type = BCF_MP_WF, M;
	char *fn_ref = 0, *fn_out = 0, moder[8], modew[8], **samples = 0;
	bcf_hdr_t *h, *hsub = 0;
	htsFile *in, *out;
	hts_idx_t *idx = 0;
	bcf1_t *b;

	double *phi = 0, theta = 0.001, p_var_thr = 0.5;

	while ((c = getopt(argc, argv, "l:bSt:o:T:s:GvcIP:u:p:")) >= 0) {
		if (c == 'l') clevel = atoi(optarg), flag |= 2;
		else if (c == 'S') flag |= FLAG_VCFIN;
		else if (c == 'b') flag |= FLAG_BCFOUT;
		else if (c == 'c') flag |= FLAG_CALL;
		else if (c == 'G') n_samples = 0;
		else if (c == 't') fn_ref = optarg, flag |= 1;
		else if (c == 'o') fn_out = optarg;
		else if (c == 's') samples = hts_readlines(optarg, &n_samples);
		else if (c == 'I') flag |= FLAG_NOINDEL;
		else if (c == 'u') theta = atof(optarg);
		else if (c == 'p') p_var_thr = atof(optarg);
		else if (c == 'v') flag |= FLAG_VARONLY | FLAG_CALL;
		else if (c == 'P') {
			if (strcmp(optarg, "wf") == 0) prior_type = BCF_MP_WF;
			else if (strcmp(optarg, "cond2") == 0) prior_type = BCF_MP_COND2;
			else if (strcmp(optarg, "flat") == 0) prior_type = BCF_MP_FLAT;
			else abort(); // not implemented
			break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\nUsage:   vcfview [options] <in.bcf>|<in.vcf>|<in.vcf.gz>\n\n");
		fprintf(stderr, "Options: -b           BCF output\n");
		fprintf(stderr, "         -S           VCF input\n");
		fprintf(stderr, "         -o FILE      output file name [stdout]\n");
		fprintf(stderr, "         -l INT       compression level [%d]\n", clevel);
		fprintf(stderr, "         -t FILE      list of reference names and lengths [null]\n");
		fprintf(stderr, "         -s FILE/STR  list of samples (STR if started with ':'; FILE otherwise) [null]\n");
		fprintf(stderr, "         -G           drop individual genotype information\n");
		fprintf(stderr, "         -I           exclude INDELs\n\n");
		fprintf(stderr, "         -c           variant calling\n");
		fprintf(stderr, "         -v           only output variants (force -c)\n");
		fprintf(stderr, "         -t FLOAT     theta [%.g]\n", theta);
		fprintf(stderr, "         -p FLOAT     variant if P(var)>FLOAT [%g]\n", p_var_thr);
		fprintf(stderr, "         -P FILE/STR  prior: wf, cond2, flat or FILE [wf]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	h = vcf_hdr_read(in);
	if (h == 0) {
		fprintf(stderr, "[E::%s] fail to read the VCF/BCF2 header\n", __func__);
		hts_close(in);
		return 1;
	}
	if (n_samples >= 0) {
		if (n_samples) imap = (int*)malloc(n_samples * sizeof(int));
		hsub = bcf_hdr_subset(h, n_samples, samples, imap);
	}
	b = bcf_init1();

	strcpy(modew, "w");
	if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
	if (flag&2) strcat(modew, "b");
	out = hts_open(fn_out? fn_out : "-", modew, 0);
	vcf_hdr_write(out, hsub? hsub : h);

	if (optind + 1 < argc && !(flag&FLAG_VCFIN)) { // BAM input and has a region
		if ((idx = bcf_index_load(argv[optind])) == 0) {
			fprintf(stderr, "[E::%s] fail to load the BCF index\n", __func__);
			return 1;
		}
	}
	for (i = optind; i < argc; ++i) {
		hts_itr_t *iter;
		if (idx && i == optind) continue; // when there is an index, don't do linear reading
		if (idx == 0 && i != optind) break; // when there is no index, don't do indexed reading
		if (i == optind) { // linear reading
			iter = hts_itr_query(0, HTS_IDX_REST, 0, 0);
		} else { // indexed reading
			if ((iter = bcf_itr_querys(idx, h, argv[i])) == 0) {
				fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
				continue;
			}
		}
		while (bcf_itr_next((BGZF*)in->fp, iter, b) >= 0) {
			bcf_hdr_t *h1 = hsub? hsub : h;
			if ((flag&FLAG_NOINDEL) && !bcf_is_snp(b)) continue;
			if (imap) bcf_subset(h, b, n_samples, imap);
			if (flag&FLAG_CALL) {
				double *pdg, *y, p_var, p_ref;
				long double sum = 0;
				int k;
				if (phi == 0) {
					M = b->n_sample * 2; // FIXME: read ploidy from input!
					phi = malloc((M + 1) * sizeof(double));
					bcf_m_prior(prior_type, theta, M, phi);
				}
				pdg = bcf_m_get_pdg3(h1, b, "PL");
				y = malloc((M + 1) * sizeof(double));
				if (pdg == 0) {
					double *pdg2[2], *y2[2];
					y2[0] = malloc((M + 1) * sizeof(double));
					y2[1] = malloc((M + 1) * sizeof(double));
					pdg2[0] = bcf_m_get_pdg3(h1, b, "FL");
					pdg2[1] = bcf_m_get_pdg3(h1, b, "RL");
					bcf_m_lk2(b->n_sample, pdg2[0], y2[0], 0);
					bcf_m_lk2(b->n_sample, pdg2[1], y2[1], 0);
					for (k = 0; k <= M; ++k) y[k] = y2[0][k] * y2[1][k];
					free(y2[0]); free(y2[1]); free(pdg2[0]); free(pdg2[1]);
				} else bcf_m_lk2(b->n_sample, pdg, y, 0);
				for (k = 0, sum = 0.; k <= M; ++k) sum += (long double)phi[k] * y[k];
				for (k = 0, p_var = 0.; k < M; ++k) p_var += phi[k] * y[k] / sum;
				p_ref = phi[M] * y[M] / sum;
				free(y); free(pdg);
				b->qual = p_var > p_var_thr? -4.343 * log(p_ref) : 4.343 * log(p_var);
				if ((flag&FLAG_VARONLY) && p_var < p_var_thr) continue;
			}
			vcf_write1(out, h1, b);
		}
		hts_itr_destroy(iter);
	}
	free(phi);
	if (idx) hts_idx_destroy(idx);
	hts_close(out);

	bcf_destroy1(b);
	if (imap) {
		for (i = 0; i < n_samples; ++i) free(samples[i]);
		free(samples);
		bcf_hdr_destroy(hsub);
		free(imap);
	}
	bcf_hdr_destroy(h);
	hts_close(in);
	return 0;
}

