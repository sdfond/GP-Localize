#include <assert.h>
#include <stdint.h>
#include "mvas_gp.h"
#include "nr3/nr3.h"
#include "map_reduce.h"
#include "mapred_gp.h"
int keycomp( const void *key1, const void *key2 )
{
  return 0;
}

int mvasgp_splitter( void *data_in, int reg_units, map_args_t *out )
{
  assert( data_in );
  mrgp_data_t* data = (mrgp_data_t *) data_in;

  if ( data->pos >= data->size ) {
    return 0;
  }

  //printf( "Split Function Running \n" );
  out-> data = (void *) &data->mdt[data->pos];
  out-> length = sizeof(map_data_t);

  data-> pos ++;
  return 1;
}

void mvasgp_map ( map_args_t *args )
{
  //printf( "Map function Running \n" );
  assert( args );

  map_data_t* key = (map_data_t*) args->data;

  assert(key);
  assert(key->mvas_pointer);
  assert(key->D);
  assert(key->ast);
  assert(key->chol_kuu);

  //printf( "dsv : %d\n", key->dsv );
  key->chol_sddk = key->mvas_pointer->chol_pcov( *(key->D), key->dsv, *(key->ast), key->as, key->chol_kuu );

  //printf( "Error 2" );
  key->mvas_pointer->pitc_prep( *(key->D), key->dsv, key->chol_sddk, *(key->ast), key->as, key->chol_kuu, *(key->lvt->ls_zu), *(key->lvt->ls_kuu));

  //printf( "Error 3" );
  global_value_t* gvt = key->gvt;
  emit_intermediate ( (void *) gvt, (void *) key->lvt, sizeof(gvt) );

}

void mvasgp_reduce( void *key_in, iterator_t *itr )
{
  global_value_t* gvt = (global_value_t*) key_in;
  void * val;

  assert( key_in );
  assert( itr );

  int nrows = gvt->gs_kuu->nrows();
  int ncols = gvt->gs_kuu->ncols();
  while( iter_next( itr, &val ) ) {
    local_value_t * lvt = (local_value_t *) val;

    for( Int r = 0; r < nrows; r++ ) {
      (*gvt->gs_zu)[r] += (*lvt->ls_zu)[r];
      for( Int c = 0; c < ncols; c++ ) {
        (*gvt->gs_kuu)[r][c] += (*lvt->ls_kuu)[r][c];
      }
    }
  }
  emit( (void *) gvt, (void *) 1 );
}

void mvasgp_start( mrgp_data_t* data, int num_procs, int num_threads)
{
  final_data_t result;

  map_reduce_init();

  map_reduce_args_t mapred_args;
  memset( &mapred_args, 0 , sizeof(mapred_args) );
  mapred_args.task_data = data;
  mapred_args.map = mvasgp_map;
  mapred_args.reduce = mvasgp_reduce;
  mapred_args.combiner = NULL;
  mapred_args.splitter = mvasgp_splitter;
  mapred_args.key_cmp = keycomp;
  mapred_args.locator = NULL;
  mapred_args.partition = NULL;
  mapred_args.result = &result;
  mapred_args.data_size = sizeof(mrgp_data_t);
  mapred_args.unit_size = 1;

  mapred_args.L1_cache_size = 1024 * 1024 * 2;
  mapred_args.num_map_threads = num_threads;
  mapred_args.num_reduce_threads = 1;
  mapred_args.num_merge_threads = 1;
  mapred_args.num_procs = num_procs;
  mapred_args.key_match_factor = 1;

  map_reduce( &mapred_args );

  map_reduce_finalize();

  //return (global_value_t*) result.keyval_t->key;
}

int regr_splitter ( void *data_in, int req_units, map_args_t *out )
{
  assert( data_in );
  regr_data_t* data = (regr_data_t *) data_in;

  if ( data->pos >= data->size ) {
    return 0;
  }
#ifdef DEBUG
  printf( "PITC Regr Split Function Running \n" );
#endif

  out-> data = (void *) &data->mrd[data->pos];
  out-> length = sizeof(map_regr_data);

  data-> pos ++;
  return 1;
}

void regr_map ( map_args_t *args )
{
  assert( args );

  map_regr_data* key = (map_regr_data*) args->data;

  assert(key);
  assert(key->mvas_pointer);
  assert(key->u);
  assert(key->alpha);
  assert(key->chol_kuu);
  assert(key->chol_suu);
  key->mvas_pointer->pitc_regr_low2_core( *(key->u), key->ss, *(key->tset_blk), key->ts_blk,
                                          *(key->pmu_blk), *(key->pvar_blk), *(key->alpha), *(key->chol_suu), *(key->chol_kuu));
}

void regr_reduce ( void *key_in, iterator_t *itr )
{
}

void regr_start( regr_data_t * data, int num_procs, int num_threads)
{
  final_data_t result;

  map_reduce_init();

  map_reduce_args_t mapred_args;
  memset( &mapred_args, 0 , sizeof(mapred_args) );
  mapred_args.task_data = data;
  mapred_args.map = regr_map;
  mapred_args.reduce = regr_reduce;
  mapred_args.combiner = NULL;
  mapred_args.splitter = regr_splitter;
  mapred_args.key_cmp = keycomp;
  mapred_args.locator = NULL;
  mapred_args.partition = NULL;
  mapred_args.result = &result;
  mapred_args.data_size = sizeof(regr_data_t);
  mapred_args.unit_size = 1;

  mapred_args.L1_cache_size = 1024 * 1024 * 2;
  mapred_args.num_map_threads = num_threads;
  mapred_args.num_reduce_threads = 1;
  mapred_args.num_merge_threads = 1;
  mapred_args.num_procs = num_procs;
  mapred_args.key_match_factor = 1;

  map_reduce( &mapred_args );

  map_reduce_finalize();
}

int pic_splitter( void* data_in, int req_units, map_args_t* out )
{
  assert( data_in );
  mrgp_data_t* data = (mrgp_data_t *) data_in;

  if ( data->pos >= data->size ) {
    return 0;
  }
#ifdef DEBUG
  printf( "PIC MR Split Function Running \n" );
#endif

  out-> data = (void *) &data->mdt[data->pos];
  out-> length = sizeof(map_data_t);

  data-> pos ++;
  return 1;
}

void pic_map( map_args_t* args )
{
#ifdef DEBUG
  printf( "Map function Running \n" );
#endif
  assert( args );

  map_data_t* key = (map_data_t*) args->data;

  assert(key);
  assert(key->mvas_pointer);
  assert(key->D);
  assert(key->ast);
  assert(key->chol_kuu);
  assert(key->chol_sddk);

  key->mvas_pointer->pic_core( *(key->D), key->dsv,
                               *(key->ast), key->as,
                               *(key->tsetk), key->tsk,
                               *(key->pmuk), *(key->pvark),
                               *(key->kuu), key->chol_sddk,
                               key->chol_kuu, key->chol_suu,
                               *(key->lvt), *(key->gvt) );
#ifdef DEBUG
  printf( "Map function finished \n" );
#endif
}

void pic_reduce( void *key_in, iterator_t *itr )
{
}


void picmr_start( mrgp_data_t* data, int num_procs, int num_threads )
{
  final_data_t result;

  map_reduce_init();

  map_reduce_args_t mapred_args;
  memset( &mapred_args, 0 , sizeof(mapred_args) );
  mapred_args.task_data = data;
  mapred_args.map = pic_map;
  mapred_args.reduce = pic_reduce;
  mapred_args.combiner = NULL;
  mapred_args.splitter = pic_splitter;
  mapred_args.key_cmp = keycomp;
  mapred_args.locator = NULL;
  mapred_args.partition = NULL;
  mapred_args.result = &result;
  mapred_args.data_size = sizeof(mrgp_data_t);
  mapred_args.unit_size = 1;

  mapred_args.L1_cache_size = 1024 * 1024 * 2;
  mapred_args.num_map_threads = num_threads;
  mapred_args.num_reduce_threads = 1;
  mapred_args.num_merge_threads = 1;
  mapred_args.num_procs = num_procs;
  mapred_args.key_match_factor = 1;

  map_reduce( &mapred_args );

  map_reduce_finalize();
}


