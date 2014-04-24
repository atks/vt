/* Use librandom from the julia source tree */

extern double dsfmt_gv_genrand_close_open();
extern double randmtzig_gv_randn (void);
extern double randmtzig_exprnd (void);

double unif_rand(void) {
    return dsfmt_gv_genrand_close_open();
}

double exp_rand(void) {
    return randmtzig_exprnd();
}

double norm_rand(void) {
    return randmtzig_gv_randn();
}
