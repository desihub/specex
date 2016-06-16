
#ifndef CSPECEX_H
#define CSPECEX_H

#include <specex_desi.h>
#include <specex_merge.h>


#ifdef __cplusplus
extern "C" {
#endif


int cspecex_desi_psf_fit(int argc, char * argv[]);

int cspecex_psf_merge(int argc, char * argv[]);

int cspecex_spot_merge(int argc, char * argv[]);


#ifdef __cplusplus
}
#endif

#endif

