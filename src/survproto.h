/*
 ** $Id $
 **
 ** Prototypes of all the survival functions
 **  Including this in each routine helps prevent mismatched argument errors
 */

void agexact(int *maxiter,  int *nusedx,   int *nvarx,   double *start, 
             double *stop,   int *event,    double *covar2,double *offset, 
             int   *strata, double *means,  double *beta,  double *u, 
             double *imat2,  double loglik[2], int *flag,  double *work, 
             int   *work2,  double *eps,    double *tol_chol, double *sctest);
             
             void agfit2( int   *maxiter,  int   *nusedx,  int   *nvarx, 
                          double *start,    double *stop,    int   *event, 
                          double *covar2,   double *offset,  double *weights,
                          int   *strata,   double *means,   double *beta, 
                          double *u,        double *imat2,   double loglik[2], 
                                                                          int   *flag,     double *work,    int   *end,
                                                                          double *eps,      double *tol_chol,double *sctest);
             
             void agfit4_a(int *nusedx, int *nvarx, double *yy, 
                           double *covar2, double *offset2,
                           double *weights2, int *strata2,
                           double *means, double *beta, double *u, 
                           double *loglik, 
                           int *methodx, int *ptype2, int *pdiag2,
                           int *nfrail,  int *frail2);
             
             void agfit4_b(int *maxiter, int *nusedx, int *nvarx, 
                           double *beta, double *u,
                           double *imat2,  double *jmat2, double *loglik, 
                           int *flag,  double *eps, double *tolerch, int *methodx, 
                           int *nfrail, double *fbeta, double *fdiag);
             
             void agfit_null(int   *n,      int   *method,   double *start, double *stop, 
                             int   *event,  double * offset,  double *weights,
                             int   *strata, double loglik[2]);
             
             void aghaz2(int   *n,     double *start,   double *stop,   int   *event, 
                         double *score, int   * strata, double *hazard, double * cumhaz);
             
             void agmart(int   *n,     int   *method,  double *start,   double *stop, 
                         int   *event, double *score,   double *wt,      int   *strata, 
                         double *resid);
             
             void agres12(int   *nx,     int   *nvarx,   double *y,    double *covar2, 
                          int   *strata, double *score,   int *method, double *resid2, 
                          double *a);
             
             void agscore(int   *nx,       int   *nvarx,      double *y,
                          double *covar2,   int   *strata,     double *score,
                          double *weights,  int   *method,     double *resid2, double *a);
             
             void agsurv1(int   *sn,     int   *snvar,  double *y,      double *score, 
                          int   *strata, 
                          double *surv,   double *varh,   int   *snsurv,
                          double *xmat,   double *d,      double *varcov, double *yy,
                          int   *shisn,  double *hisy,   double *hisxmat,double *hisrisk, 
                          int   *hisstrat);
             
             void agsurv2(int   *sn,      int   *snvar,    double *y, 
                          double *score,   int   *strata,   double *wt,    double *surv, 
                          double *varh,    double *xmat,     double *varcov, 
                          int   *snsurv,  double *d,        int   *sncurve,
                          double *newx,    double *newrisk);
             
             void chinv2  (double **matrix, int n);
             int cholesky2(double **matrix, int n, double toler);
             void chsolve2(double **matrix, int n, double *y);
             void chinv3(double **matrix , int n, int m, double *fdiag);
             int cholesky3(double **matrix, int n, int m, double *diag, double toler);
             void chsolve3(double **matrix, int n, int m, double *diag, double *y);
             
             void coxdetail(int   *nusedx,   int   *nvarx,    int   *ndeadx, 
                            double *y,        double *covar2,   int   *strata,  
                            double *score,    double *weights,  double *means2, 
                            double *u2,       double *var,      int   *rmat,
                            double *nrisk2,   double *work);
             
             void coxfit2(int   *maxiter,   int   *nusedx,    int   *nvarx, 
                          double *time,      int   *status,    double *covar2, 
                          double *offset,	double *weights,   int   *strata,
                          double *means,     double *beta,      double *u, 
                          double *imat2,     double loglik[2],  int   *flag, 
                          double *work,	double *eps,       double *tol_chol,
                          double *sctest,double *sctest2, double *sctest3);
             
             void coxfit4_a(int *nusedx, int *nvarx, double *yy, 
                            double *covar2, double *offset2,
                            double *weights2, int *strata2,
                            double *means, double *beta, double *u, 
                            double *loglik, 
                            int *methodx, int *ptype2, int *pdiag2,
                            int *nfrail,  int *frail2);
             
             void coxfit4_b(int *maxiter, int *nusedx, int *nvarx, 
                            double *beta, double *u,
                            double *imat2,  double *jmat2, double *loglik, 
                            int *flag,  double *eps, double *tolerch, int *methodx, 
                            int *nfrail, double *fbeta, double *fdiag);
             
             void coxfit4_c (int *nusedx, int *nvar, int *methodx, double *expect);
             
             void coxfit_null(int   *nusedx,    int   *method,   double *time, 
                              int   *status,    double *score,    double *weights, 
                              int   *strata,    double *loglik, double *resid);
             
             void coxhaz2(int   *n,      double *score,   int   *mark, 
                          int   *strata, double *hazard,  double *cumhaz);
             
             void coxmart(int   *sn,     int   *method,    double *time, 
                          int   *status, int   * strata,   double *score, 
                          double *wt,     double *expect);
             
             void coxph_wtest(int *nvar2, int *ntest, double *var, double *b,
                              double *scratch, double *tolerch);
             
             void coxscho(int   *nusedx,    int   *nvarx,    double *y, 
                          double *covar2,    double *score,    int   *strata,  
                          int   *method2,   double *work);
             
             void coxscore(int   *nx,      int   *nvarx,    double *y, 
                           double *covar2,  int   *strata,   double *score, 
                           double *weights, int   *method,   double *resid2,
                           double *scratch);
             
             double coxsafe(double x);
             double **dmatrix(double *array, int ncol, int nrow);
             
             void init_doloop(int min, int max);
             int doloop      (int nloops, int *index);
             
             void pyears1(int   *sn,      int   *sny,      int   *sdoevent, 
                          double *sy,      double *wt,       
                          int   *sedim,   int   *efac, 
                          int   *edims,   double *secut,    double *expect, 
                          double *sedata,  int   *sodim,    int   *ofac, 
                          int   *odims,   double *socut,    int   *smethod, 
                          double *sodata,  double *pyears,   double *pn, 
                          double *pcount,  double *pexpect,  double *offtable);
             
             void pyears2(int   *sn,      int   *sny,   int   *sdoevent, 
                          double *sy,      double *wt,    
                          int   *sodim,   int   *ofac, 
                          int   *odims,   double *socut, double *sodata,
                          double *pyears,  double *pn,    double *pcount, 
                          double *offtable);
             
             double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
                           double *data,  int *fac,    int *dims,     double **cuts, 
                           double step,   int  edge);
             
             void survdiff2(int   *nn,     int   *nngroup,    int   *nstrat, 
                            double *rho,    double *time,       int   *status, 
                            int   *group,  int   *strata,	   double *obs, 
                            double *exp,    double *var,        double *risk, 
                            double *kaplan);
             
             void survfit2(int   *sn,      double *y,       double *wt,
                           int   *strata,  int   *method,  int   *error, 
                           double *mark,    double *surv,	double *varh,
                           double *risksum);
             
             void survfit3(int   *sn,        double *y,               double *wt,
                           int   *strata,    int   *method,          int   *error, 
                           int   *nstrat,    double *ntimes_strata,  
                           double *timelist,  double *weighted_event,  double *surv,
                           double *varh,	 double *risksum,         double *enter,
                           double *exit_censored);
             
             
             SEXP survreg6(SEXP maxiter2,   SEXP nvarx,  SEXP y,
                           SEXP ny2,        SEXP covar2, SEXP wtx,
                           SEXP offset2,    SEXP beta2,  SEXP nstratx,
                           SEXP stratax,    SEXP epsx,   SEXP tolx,       
                           SEXP dist,       SEXP expr,   SEXP rho);
             
             double survregc1(int n,          int nvar,     int nstrat,      int whichcase,
                              double *beta,   int dist,     int *strat,     double *offset,
                              double *time1,  double *time2, double *status, double *wt,
                              double **covar, double **imat, double **JJ,    double *u, 
                              SEXP expr,      SEXP rho,      double *dummy,  int nf,
                              int *frail,    double *fdiag, double *jdiag );
             
             double survregc2(int n,          int nvar,     int nstrat,      int whichcase,
                              double *beta,   int dist,     int *strat,     double *offset,
                              double *time1,  double *time2, double *status, double *wt,
                              double **covar, double **imat, double **JJ,    double *u, 
                              SEXP expr,      SEXP rho,      double *dummy,  int nf,
                              int *frail,    double *fdiag, double *jdiag );
             
             void survpenal(int whichcase, int nfrail,    int  nvar2,    double **hmat, 
                            double **JJ,   double *hdiag, double *jdiag,
                            double *u,     double *beta,  double *loglik,
                            int ptype,     int pdiag,     SEXP pexpr1,   double *cptr1, 
                            SEXP pexpr2,   double *cptr2, SEXP rho);
             
             void coxpenal (int whichcase, int nfrail,   int  nvar,     double **hmat, 
                            double *hdiag, double *u,    double *fbeta, double *beta,  
                            double *penalty,
                            int ptype,     int pdiag,     SEXP pexpr1,   double *cptr1, 
                            SEXP pexpr2,   double *cptr2, SEXP rho);
