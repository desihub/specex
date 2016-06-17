
#include <cspecex.h>


int cspecex_desi_psf_fit(int argc, char * argv[]) {
    return specex_desi_psf_fit_main(argc, argv);
}

int cspecex_psf_merge(int argc, char * argv[]) {
    return specex_merge_psf_main(argc, argv);
}

int cspecex_spot_merge(int argc, char * argv[]) {
    return specex_merge_spot_main(argc, argv);
}

