/* Use gsl library to fit bspline model to data */

#define BS_T double

void set_vector_from_data (double *data, size_t n, gsl_vector *v) {
  v->size = n;
  v->stride = 1;
  v->block = NULL;
  v->owner = 0;
  v->data = data;
}

void bspline_linear_fit( BS_T *xin, BS_T *yin, int n, int ncoeffs, int k,
                         BS_T *xout, BS_T *yout, BS_T *yout_deriv,
                         BS_T *yout_err, BS_T *yout_deriv_err, int nout) {

  size_t nbreak = ncoeffs + 2 - k;
  int nderiv=1; // number of derivatives to compute
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_bspline_deriv_workspace *dw;
  gsl_vector *B;
  gsl_matrix *dB;
  gsl_vector *c;
  //  gsl_vector *w;
  gsl_vector *x, *y;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;
  double chisq, Rsq, dof, tss;
  double xmin, xmax;
  x = (gsl_vector *) malloc(sizeof(gsl_vector));
  y = (gsl_vector *) malloc(sizeof(gsl_vector));

  bw = gsl_bspline_alloc(k, nbreak);
  if (yout_deriv != NULL) dw = gsl_bspline_deriv_alloc(k);
  B = gsl_vector_alloc(ncoeffs);
  if (yout_deriv != NULL)   dB = gsl_matrix_alloc(ncoeffs,nderiv+1);

  set_vector_from_data(xin, n, x);
  set_vector_from_data(yin, n, y);
  xmin = xin[0];
  xmax = xin[n-1];
  
  X = gsl_matrix_alloc(n, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  //  w = gsl_vector_alloc(n);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);

  mw = gsl_multifit_linear_alloc(n, ncoeffs);
  gsl_bspline_knots_uniform(xmin, xmax, bw);

  /* construct the fit matrix X */
  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(x, i);

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);

      /* fill in row i of X */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
        }
    }
  //  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);
  gsl_multifit_linear(X, y, c, cov, &chisq, mw);

  dof = n - ncoeffs;
  //  tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
  tss = gsl_stats_tss(y->data, 1, y->size);
  Rsq = 1.0 - chisq / tss;

  //  fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
  //                   chisq / dof, Rsq);

  {
    double yi, yerr;

    if ( yout_deriv != NULL ) {
      int nderiv = 1;
      gsl_vector_view v;
      if ( yout != NULL ) {
        for(i=0; i < nout; i++) {
          gsl_bspline_deriv_eval(xout[i], nderiv, dB, bw, dw);
          
          v = gsl_matrix_column(dB,0);
          gsl_multifit_linear_est(&v.vector, c, cov, &yi, &yerr);
          yout[i] = yi;
          if(yout_err != NULL) yout_err[i] = yerr;
          
          v = gsl_matrix_column(dB,1);
          gsl_multifit_linear_est(&v.vector, c, cov, &yi, &yerr);
          yout_deriv[i] = yi;
          if(yout_deriv_err != NULL)  yout_deriv_err[i] = yerr;
        }
      }
      else {
        for(i=0; i < nout; i++) {
          gsl_bspline_deriv_eval(xout[i], nderiv, dB, bw, dw);
          v = gsl_matrix_column(dB,1);
          gsl_multifit_linear_est(&v.vector, c, cov, &yi, &yerr);
          yout_deriv[i] = yi;
          if(yout_deriv_err != NULL)  yout_deriv_err[i] = yerr;
        }
      }
    }
    else {
      if ( yout != NULL ) {
        for(i=0; i < nout; i++) {
          gsl_bspline_eval(xout[i], B, bw);
          gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
          yout[i] = yi;
          if(yout_err != NULL) yout_err[i] = yerr;
        }
      }
    }
  }

  gsl_bspline_free(bw);
  if (yout_deriv != NULL) gsl_bspline_deriv_free(dw);
  if (yout_deriv != NULL) gsl_matrix_free(dB);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  //  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);

}
