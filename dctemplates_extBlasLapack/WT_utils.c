/* functions needed by some part of templates */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"

/* WT: if we set 
#define STATIC 
only we can make many places more thread safe! */
#include "STATICdef.c"


/* Common Block Declarations */
struct {
    doublereal a[40000], m[200];
} system_;
#define system_1 system_

struct {
    integer n, lda;
} matdim_;
#define matdim_1 matdim_

struct {
    char curpform[5];
} forms_;
#define forms_1 forms_

/* Table of constant values */
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;


/* function getbreak_ from Tester.c */

/*     ================================================================ */
doublereal getbreak_()
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    extern doublereal dlamch_();  /* dlamch_ is in lapack */
    STATIC doublereal eps;


/*     Get breakdown parameter tolerance; for the test routine, */
/*     set to machine precision. */


    eps = dlamch_("EPS", 3L);
/* Computing 2nd power */
    d__1 = eps;
    ret_val = d__1 * d__1;

    return ret_val;

} /* getbreak_ */



/* *********************************************************************** */
/* some functions from Utils.c with all Fortran craziness and IO replaced
   by C standard stuff  */
/* *********************************************************************** */
/* I think stuff like:
     lsame_(flag_, "SPLIT", 3L, 5L)
   can be replaced by:
     strncmp(flag_, "SPLIT", 1)==0
*/

/*  This file contains routines used by Jacobi, SOR, and Chebyshev: */

/*  Jacobi/SOR: */
/*     MATSPLIT  calls specific matrix splitting routine */
/*     JACSPLIT */
/*     SORSPLIT */
/*     BACKSOLVE */

/*  Chebyshev: */
/*     GETEIG    computes eigenvalue of iteration matrix. */

/*     =========================================================== */
/* Subroutine */ int jacsplit_(n, a, lda, work, ldw, flag_, flag_len)
integer *n;
doublereal *a;
integer *lda;
doublereal *work;
integer *ldw;
char *flag_;
ftnlen flag_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2;

    ///* Builtin functions */
    //integer s_wsle(), do_lio(), e_wsle();
    ///* Subroutine */ int s_stop();

    /* Local variables */
    STATIC integer i, j;
    extern logical lsame_();

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };

    /* Parameter adjustments */
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    //if (lsame_(flag_, "SPLIT", 1L, 5L))
    if(strncmp(flag_, "SPLIT", 1)==0) {
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    work[i + work_dim1] = 1. / a[i + i * a_dim1];
	    a[i + i * a_dim1] = 0.;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		a[i + j * a_dim1] = -a[i + j * a_dim1];
/* L10: */
	    }
/* L20: */
	}
    }
    //else if (lsame_(flag_, "RECONSTRUCT", 1L, 11L)) {
    else if(strncmp(flag_, "RECONSTRUCT", 1)==0) {
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		a[i + j * a_dim1] = -a[i + j * a_dim1];
/* L30: */
	    }
	    a[i + i * a_dim1] = 1. / work[i + work_dim1];
/* L40: */
	}
    } else {
	//s_wsle(&io___4);
	//do_lio(&c__9, &c__1, "UNKNOWN SPLITTING OPTION. QUITTING...", 37L);
	//e_wsle();
	//s_stop("", 0L);
	printf("jacsplit_: UNKNOWN SPLITTING OPTION. QUITTING...\n");
	exit(911);
    }

    return 0;

} /* jacsplit_ */




/*     =========================================================== */
/* Subroutine */ int sorsplit_(omega, n, a, lda, b, work, ldw, flag_, 
	flag_len)
doublereal *omega;
integer *n;
doublereal *a;
integer *lda;
doublereal *b, *work;
integer *ldw;
char *flag_;
ftnlen flag_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2;

    ///* Builtin functions */
    //integer s_wsle(), do_lio(), e_wsle();
    ///* Subroutine */ int s_stop();

    /* Local variables */
    STATIC integer i, j;
    extern logical lsame_();

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, 0, 0 };

    /* Parameter adjustments */
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    //if (lsame_(flag_, "SPLIT", 3L, 5L))
    if(strncmp(flag_, "SPLIT", 1)==0) {

/*        Set M. */

	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    work[i + i * work_dim1] = a[i + i * a_dim1];
	    i__2 = i - 1;
	    for (j = 1; j <= i__2; ++j) {
		work[i + j * work_dim1] = *omega * a[i + j * a_dim1];
/* L10: */
	    }
/* L20: */
	}

/*        Set NN and B. */

/*        Temporarily store the matrix A in order to reconstruct */
/*        the original matrix. Because the lower triangular portion */
/*        of A must be zeroed, this is the easiest way to deal with it
. */
/*        This causes the requirement that WORK be N x (2N+3). */

	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		work[i + (j + *n + 3) * work_dim1] = a[i + j * a_dim1];
/* L30: */
	    }
/* L40: */
	}

	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    b[i] = *omega * b[i];
	    a[i + i * a_dim1] = (1. - *omega) * a[i + i * a_dim1];
	    i__2 = *n;
	    for (j = i + 1; j <= i__2; ++j) {
		a[i + j * a_dim1] = -(*omega) * a[i + j * a_dim1];
/* L50: */
	    }
/* L60: */
	}

	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
	    i__2 = i - 1;
	    for (j = 1; j <= i__2; ++j) {
		a[i + j * a_dim1] = 0.;
/* L70: */
	    }
/* L80: */
	}

    }
    //else if (lsame_(flag_, "RECONSTRUCT", 3L, 11L))
    else if(strncmp(flag_, "RECONSTRUCT", 1)==0) {
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    b[i] /= *omega;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		a[i + j * a_dim1] = work[i + (j + *n + 3) * work_dim1];
/* L90: */
	    }
/* L100: */
	}

    } else {
	//s_wsle(&io___7);
	//do_lio(&c__9, &c__1, "UNKNOWN SPLITTING OPTION. QUITTING...", 37L);
	//e_wsle();
	//s_stop("", 0L);
        printf("sorsplit_: UNKNOWN SPLITTING OPTION. QUITTING...\n");
        exit(911);
    }
    return 0;

} /* sorsplit_ */


/*     ======================================================== */
/* Subroutine */ int matsplit_(omega, b, work, ldw, method, flag_, method_len,
	 flag_len)
doublereal *omega, *b, *work;
integer *ldw;
char *method, *flag_;
ftnlen method_len;
ftnlen flag_len;
{
    /* System generated locals */
    integer work_dim1, work_offset;

    ///* Builtin functions */
    //integer s_wsle(), do_lio(), e_wsle();
    ///* Subroutine */ int s_stop();

    /* Local variables */
    extern /* Subroutine */ int jacsplit_(), sorsplit_();
    extern logical lsame_();

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };

/*     .. */
/*     .. Parameters .. */

/*     MAXDIM2 = MAXDIM*MAXDIM. */


/*     .. Common Blocks .. */


    /* Parameter adjustments */
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --b;

    /* Function Body */
    //if (lsame_(method, "JACOBI", 1L, 6L))
    if(strncmp(method, "JACOBI", 1)==0) {
	jacsplit_(&matdim_1.n, system_1.a, &matdim_1.lda, &work[work_offset], 
		ldw, flag_, 1L);
    }
    //else if (lsame_(method, "SOR", 1L, 3L))
    else if (strncmp(method, "SOR", 1)==0) {
	sorsplit_(omega, &matdim_1.n, system_1.a, &matdim_1.lda, &b[1], &work[
		work_offset], ldw, flag_, 1L);
    } else {
	//s_wsle(&io___1);
	//do_lio(&c__9, &c__1, "ERROR: UNKNOW METHOD. QUITTING...", 33L);
	//e_wsle();
	//s_stop("", 0L);
	printf("matsplit_: ERROR: UNKNOW METHOD. QUITTING...\n");
	exit(911);
    }

    return 0;

} /* matsplit_ */



/* renamed backsolve_ from Utils.c */
/*     ======================================================== */
/* Subroutine */ int templates_backsolve(n, a, lda, x)
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    extern /* Subroutine */ int dtrsv_();


/*     .. Argument Declarations .. */

/*     Mask to BLAS routine. X overwritten with inv(A)*X */

    /* Parameter adjustments */
    --x;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    dtrsv_("LOWER", "NOTRANS", "NONUNIT", n, &a[a_offset], lda, &x[1], &c__1, 
	    5L, 7L, 7L);

    return 0;

} /* backsolve_ */
