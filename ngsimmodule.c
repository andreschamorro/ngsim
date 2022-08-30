/* Extension code.
**
** Andres Chamorro 
** This code is in the public domain.
*/
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>
#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define PACKAGE_VERSION "0.0.1-r1"

const uint8_t nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
/* Simple normal random number generator, copied from genran.c */

double ran_normal()
{ 
	static int iset = 0; 
	static double gset; 
	double fac, rsq, v1, v2; 
	if (iset == 0) {
		do { 
			v1 = 2.0 * drand48() - 1.0;
			v2 = 2.0 * drand48() - 1.0; 
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq); 
		gset = v1 * fac; 
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

/* wgsim */

enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
	int l, m; /* length and maximum buffer size */
	mut_t *s; /* sequence */
} mutseq_t;

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTEND = 0.3;
static double MAX_N_RATIO = 0.05;

void wgsim_mut_diref(const kseq_t *ks, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
{
	int i, deleting = 0;
	mutseq_t *ret[2];

	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
	ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
	ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
	ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
	for (i = 0; i != ks->seq.l; ++i) {
		int c;
		c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (deleting) {
            if (drand48() < INDEL_EXTEND) {
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                continue;
            } else deleting = 0;
        }
		if (c < 4 && drand48() < MUT_RATE) { // mutation
			if (drand48() >= INDEL_FRAC) { // substitution
				double r = drand48();
				c = (c + (int)(r * 3.0 + 1)) & 3;
				if (is_hap || drand48() < 0.333333) { // hom
					ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
				} else { // het
					ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
				}
			} else { // indel
				if (drand48() < 0.5) { // deletion
					if (is_hap || drand48() < 0.333333) { // hom-del
						ret[0]->s[i] = ret[1]->s[i] = DELETE;
                        deleting = 3;
					} else { // het-del
                        deleting = drand48()<0.5?1:2;
						ret[deleting-1]->s[i] = DELETE;
					}
				} else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        num_ins++;
                        ins = (ins << 2) | (int)(drand48() * 4.0);
                    } while (num_ins < 4 && drand48() < INDEL_EXTEND);

					if (is_hap || drand48() < 0.333333) { // hom-ins
						ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
					} else { // het-ins
						ret[drand48()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c;
					}
				}
			}
		}
	}
}

typedef struct {
    PyObject_HEAD
    char *fn;
		kseq_t *ks;
		size_t total_len, num_ref;
} FileState;

/* ReadState - ngs generator instance.
 *
 * sequence: ref to the sequence that's being iterated
 * seq_index: index of the next element in the sequence to yield
 * enum_index: next enumeration index to yield
 *
 * In pseudo-notation, the yielded tuple at each step is:
 *  enum_index, sequence[seq_index]
*/
typedef struct {
    PyObject_HEAD
		FileState fstate;
		Py_ssize_t size, len_l, len_r, std_dev, is_hap, dist;
    size_t buff_size, buff_index, enum_index;
    char* buff_l;
    char* buff_r;

		uint64_t cn_pairs;
		int cr_len;
		kstring_t c_refname;
    mutseq_t rseq[2];
} ReadState;

static void
fill_buff(ReadState *rstate){
	rstate->buff_index = 0;
  mut_t *target;
	int n_sub[2], n_indel[2], n_err[2], ext_coor[2], s[2], i, j, k, c, is_flip;
	size_t d, pos;
	uint8_t *tmp_seq[2];
	double ran;
	size_t max_size = rstate->len_l > rstate->len_r? rstate->len_l: rstate->len_r;
	tmp_seq[0] = (uint8_t*)calloc(rstate->len_l+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(rstate->len_r+2, 1);
	while(rstate->buff_index < rstate->buff_size){
		if (rstate->cn_pairs == 0){
			if((rstate->cr_len = kseq_read(rstate->fstate.ks)) < 0) break;
			// Skip sequence as it is shorter
			if (rstate->cr_len < rstate->dist + 3 * rstate->std_dev) break;

			rstate->cn_pairs = (uint64_t)((long double)rstate->cr_len / rstate->fstate.total_len * rstate->size + 0.5);
			rstate->c_refname = rstate->fstate.ks->name;
			// generate mutations
			wgsim_mut_diref(rstate->fstate.ks, rstate->is_hap, rstate->rseq, rstate->rseq+1);
		}
		do { // avoid boundary failure
			ran = ran_normal() * rstate->std_dev + rstate->dist;
			d = (size_t)(ran + 0.5);
			d = d > max_size? d : max_size;
			pos = (int)((rstate->cr_len - d + 1) * drand48());
		} while (pos < 0 || pos >= rstate->fstate.ks->seq.l || pos + d - 1 >= rstate->fstate.ks->seq.l);

		if (drand48() < 0.5){
			is_flip = 0;
			s[0] = rstate->len_l;
			s[1] = rstate->len_r;
		} else {
			is_flip = 1;
			s[1] = rstate->len_l;
			s[0] = rstate->len_r;
		}

		target = rstate->rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
		n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = 0;

#define __gen_read(x, start, iter) do {	\
	for (i = (start), k = 0, ext_coor[x] = -10; i >= 0 && i < rstate->fstate.ks->seq.l && k < s[x]; iter) {	\
		int c = target[i], mut_type = c & mutmsk;			\
		if (ext_coor[x] < 0) {								\
			if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
			ext_coor[x] = i;								\
		}													\
		if (mut_type == DELETE) ++n_indel[x];				\
		else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
			tmp_seq[x][k++] = c & 0xf;						\
			if (mut_type == SUBSTITUTE) ++n_sub[x];			\
		} else {											\
			int n, ins;										\
			++n_indel[x];									\
			tmp_seq[x][k++] = c & 0xf;						\
			for (n = mut_type>>12, ins = c>>4; n > 0 && k < s[x]; --n, ins >>= 2) \
				tmp_seq[x][k++] = ins & 0x3;				\
		}													\
	}														\
	if (k != s[x]) ext_coor[x] = -10;						\
} while (0)
		__gen_read(0, pos, ++i);
		__gen_read(1, pos + d -1, --i);
		for (k = 0; k < s[1]; ++k) tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
		if (ext_coor[0] < 0 || ext_coor[1] < 0) { // fail to generate the read(s)
			continue;
		}
		j = is_flip? 1: 0;
		for (i = 0; i < rstate->len_l; ++i){
			c = tmp_seq[j][i];
			if (c >= 4)
				c = 4;
			else if (drand48() < ERR_RATE)
				c = (c + 1) & 3; // recurrent sequencing errors
			rstate->buff_l[rstate->buff_index  * rstate->len_l + i] = "ACGTN"[c]; 
		}
		j = is_flip? 0: 1;
		for (i = 0; i < rstate->len_r; ++i){
			c = tmp_seq[j][i];
			if (c >= 4)
				c = 4;
			else if (drand48() < ERR_RATE)
				c = (c + 1) & 3; // recurrent sequencing errors
			rstate->buff_r[rstate->buff_index  * rstate->len_r + i] = "ACGTN"[c]; 
		}
		rstate->cn_pairs--;
		rstate->buff_index++;
	}
	rstate->buff_size = rstate->buff_index + 1;
	rstate->buff_index = 0;
}

static PyObject *
readgen_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    char *fn;
		int size, len_l, len_r, std_dev, is_hap, dist, buff_size;
		len_l = len_r = 70;
		std_dev = 50;
		dist = 500;
		buff_size = 10000;
		static char *kwlist[] = {"fafile", "size", "len_l", "len_r", "std_dev", "dist", "is_hap", "buff_size", NULL};

		if (!PyArg_ParseTupleAndKeywords(args, kwargs, "si|iiiipi:new", kwlist, 
					&fn, &size, &len_l, &len_r, &std_dev, &dist, &is_hap, &buff_size)){
        return NULL;
		}
		gzFile fp_fa = gzopen(fn, "r");
		kseq_t *ks = kseq_init(fp_fa);
		size_t tot_len = 0, n_ref = 0;
		int l = 0;
		while ((l = kseq_read(ks)) >= 0) {
			tot_len += l;
			++n_ref;
		}
		kseq_destroy(ks);
		gzclose(fp_fa);

    /* Create a new ReadState and initialize its state - pointing to the last
     * index in the sequence.
    */
    ReadState *rstate = (ReadState *)type->tp_alloc(type, 0);
    if (!rstate)
        return NULL;
		rstate->fstate.fn = fn;
		rstate->fstate.total_len = tot_len;
		rstate->fstate.num_ref = n_ref;
		rstate->fstate.ks = kseq_init(gzopen(fn, "r"));
		rstate->size = size;
		rstate->len_l = len_l;
		rstate->len_r = len_r;
		rstate->std_dev = std_dev;
		rstate->dist = dist;
		rstate->is_hap = is_hap;
		rstate->buff_size = buff_size;
		rstate->enum_index = 0;
		rstate->buff_index = 0;
		rstate->buff_l = (char *)calloc(len_l*(buff_size+2), sizeof(char)); 
		rstate->buff_r = (char *)calloc(len_r*(buff_size+2), sizeof(char)); 

		fill_buff(rstate);
    return (PyObject *)rstate;
}


static void
readgen_dealloc(ReadState *rstate)
{
    /* We need XDECREF here because when the generator is exhausted,
     * rgstate->sequence is cleared with Py_CLEAR which sets it to NULL.
    */
		free(rstate->buff_l); free(rstate->buff_r);
    Py_TYPE(rstate)->tp_free(rstate);
}


static PyObject *
readgen_next(ReadState *rstate)
{
    /* seq_index < 0 means that the generator is exhausted.
     * Returning NULL in this case is enough. The next() builtin will raise the
     * StopIteration error for us.
    */
  //if ((rstate->size >= 0) & ((l = kseq_read(ks)) >= 0)) {
	if ((rstate->enum_index < rstate->size)){
		if ((rstate->buff_index >= rstate->buff_size)){
			fill_buff(rstate);
		}
		PyObject *result = Py_BuildValue("s#s#", rstate->c_refname.s, rstate->c_refname.l, rstate->buff_l + (rstate->buff_index * rstate->len_l), rstate->len_l);
		rstate->buff_index++;
    rstate->enum_index++;
    return result;
  }

  /* The reference to the sequence is cleared in the first generator call
   * after its exhaustion (after the call that returned the last element).
   * Py_CLEAR will be harmless for subsequent calls since it's idempotent
   * on NULL.
  */
  rstate->size = -1;
  Py_CLEAR(rstate->fstate.fn);
  return NULL;
}


PyTypeObject PyReadgen_Type = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "readgen",                       /* tp_name */
    sizeof(ReadState),            /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)readgen_dealloc,     /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    PyObject_SelfIter,              /* tp_iter */
    (iternextfunc)readgen_next,      		/* tp_iternext */
    0,                              /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    readgen_new, 		                    /* tp_new */
};


static struct PyModuleDef ngsimmodule = {
   PyModuleDef_HEAD_INIT,
   "ngsim",                  /* m_name */
   "",                      /* m_doc */
   -1,                      /* m_size */
};


PyMODINIT_FUNC
PyInit_ngsim(void)
{
    PyObject *module = PyModule_Create(&ngsimmodule);
    if (!module)
        return NULL;

    if (PyType_Ready(&PyReadgen_Type) < 0)
        return NULL;
    Py_INCREF((PyObject *)&PyReadgen_Type);
    PyModule_AddObject(module, "readgen", (PyObject *)&PyReadgen_Type);

    return module;
}
