#ifndef _MAPRED_GP_H
#define _MAPRED_GP_H
#include "nr3/nr3.h"
#include "map_reduce.h"

#undef DEBUG
struct mvas_gp;

struct local_value_t{
	MatDoub *ls_kuu;
	VecDoub *ls_zu;
} ;

struct global_value_t{
	MatDoub *gs_kuu;
	VecDoub *gs_zu;
} ;

struct map_data_t{
	mvas_gp * mvas_pointer;
	MatDoub * D;
	int dsv;
	MatDoub * ast;
	int as;
	MatDoub* tsetk;
	Int tsk;
	Cholesky * chol_kuu;
	Cholesky * chol_sddk;
	Cholesky * chol_suu;
	MatDoub* kuu;
	local_value_t * lvt;
	global_value_t * gvt;

	VecDoub* pmuk;
	VecDoub* pvark;
};

struct mrgp_data_t{
	map_data_t* mdt;
	int size;
	int pos;
};

struct map_regr_data {
	mvas_gp * mvas_pointer;
	MatDoub * u;
	Int ss;
	MatDoub * tset_blk;
	int ts_blk;
	VecDoub * pmu_blk;
	VecDoub * pvar_blk;
	VecDoub * alpha;
	Cholesky * chol_suu;
	Cholesky * chol_kuu;
};

struct regr_data_t {
	map_regr_data * mrd;
	int size;
	int pos;
};

int keycomp( void *key1, void *key2 );

int mvasgp_splitter ( void *data_in, int req_units, map_args_t *out );

void mvasgp_map ( map_args_t *args );

void mvasgp_reduce ( void *key_in, iterator_t *itr );

void mvasgp_start( mrgp_data_t* data, int num_procs, int num_threads);


int regr_splitter ( void *data_in, int req_units, map_args_t *out );

void regr_map ( map_args_t *args );

void regr_reduce ( void *key_in, iterator_t *itr );

void regr_start( regr_data_t * data, int num_procs, int num_threads);


int pic_splitter( void* data_in, int req_units, map_args_t* out );

void pic_map( map_args_t* args );

void pic_reduce( void *key_in, iterator_t *itr );

void picmr_start( mrgp_data_t* data, int num_procs, int num_threads );

#endif
