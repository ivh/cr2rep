/*
 * This file is part of the CR2RES Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/
#include <math.h>
#include <cpl.h>
#include "cr2res_slitdec.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int slit_func_vert(
         int         ncols,
         int         nrows,
         int         osample,
         double  *   im,
         int     *   mask,
         double  *   ycen,
         double  *   sL,
         double  *   sP,
         double  *   model,
         double      lambda_sP,
         double      lambda_sL,
         double      sP_stop,
         int         maxiter) ;
static int bandsol(double *, double *, int, int) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_slitdec     Slit Decomposition
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief   Main slit decomposition function
  @param    img_in	    full detector image
  @param    ycen        current order mid-line y-coord
  @param    height      number of pix above and below mid-line
  @param    swath       width per swath
  @param    oversample  factor for oversampling
  @param    smooth_slit
  @param    slit_func   the returned slit function
  @param    spec        the returned spectrum
  @param    model       the returned model
  @return   int 0 if ok, -1 otherwise

  This func takes a single image (contining many orders), and a *single*
  order definition in the form of central y-corrds., plus the height.
  Swath widht and oversampling are passed through.

  The task of this function then is to

    -cut out the relevant pixels of the order
    -shift im in y integers, so that nrows becomes minimal,
        adapt ycen accordingly
    -loop over swaths, in half-steps
    -derive a good starting guess for the spectrum, by median-filter
        average along slit, beware of cosmics
    -run slit_func_vert()
    -merge overlapping swath results by linear weights from swath-width to edge.
    -return re-assembled model image, slit-fu, spectrum, new mask.
    -calculate the errors and return them. This is done by comparing the
        variance of (im-model) to the poisson-statistics of the spectrum.

 */
/*----------------------------------------------------------------------------*/
int cr2res_slitdec_vert(
        cpl_image   *   img_in,
        cpl_vector  *   ycen,
        int             height,
        int             swath,
        int             oversample,
        double          smooth_slit,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model)
{
    int i, j, nswaths;
    int row, col, x, y, ny_os;
    int sw_start, sw_end;
    double pixval, img_median;
    int * ycen_int;
    int badpix;
    double * ycen_rest;
    double * ycen_sw;
    double * img_sw_data;
    double * spec_sw_data;
    double * slitfu_sw_data;
    double * model_sw;
    int * mask_sw;
    cpl_size lenx;
    cpl_image * img_sw;
    cpl_image * tmp;
    cpl_image * img_out;
    cpl_vector * spec_sw;
    cpl_vector * slitfu_sw;
    cpl_vector * spc;
    cpl_vector * slitfu;

    lenx = cpl_image_get_size_x(img_in);
    ny_os = oversample*(height+1) +1; // number of rows after oversampling
    nswaths = (lenx / swath) ; // TODO: Allow last swath be partial

    mask_sw = cpl_malloc(height*swath*sizeof(int));
    model_sw = cpl_malloc(height*swath*sizeof(double));
    img_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    ycen_int = cpl_malloc(lenx*sizeof(int));
    ycen_rest = cpl_malloc(lenx*sizeof(double));
    ycen_sw = cpl_malloc(swath*sizeof(double));

    for (i=0;i<lenx;i++){
        ycen_int[i] = (int)cpl_vector_get(ycen,i) ;
        ycen_rest[i] = fmod(cpl_vector_get(ycen,i), 1.0) ;
    }

    // Local versions of return data
    slitfu = cpl_vector_new(ny_os);
    spc = cpl_vector_new(lenx);
    img_out = cpl_image_new(lenx, cpl_image_get_size_y(img_in), CPL_TYPE_DOUBLE);
    slitfu_sw = cpl_vector_new(ny_os);
    slitfu_sw_data = cpl_vector_get_data(slitfu_sw);

    for (i=0;i<nswaths;i++){
        sw_start = i*swath;
        sw_end = (i+1)*swath;
        cpl_msg_debug(__func__,"Img: x:%d-%d y:%d%d",
            sw_start+1, sw_end, ycen_int[0]-(height/2), ycen_int[0]+(height/2));
        for(col=0; col<swath; col++){      // col is x-index in cut-out
            x = i*swath + col;          // coords in large image
            for(row=0;row<height;row++){   // row is y-index in cut-out
                y = ycen_int[x] - (height/2) + row;
                pixval = cpl_image_get(img_in, x+1, y+1, &badpix);
                cpl_image_set(img_sw, col+1, row+1, pixval);
                if (badpix ==0) mask_sw[row*swath+col] = 1;
                else mask_sw[row*swath+col] = 0;
            }
        }
        cpl_image_save(img_sw, "tmp.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);

        img_median = cpl_image_get_median(img_sw);
        for (j=0;j<ny_os;j++) cpl_vector_set(slitfu_sw,j,img_median);
        img_sw_data = cpl_image_get_data_double(img_sw);
        tmp = cpl_image_collapse_median_create(img_sw, 0, 0, 0);
        spec_sw = cpl_vector_new_from_image_row(tmp,1);
        cpl_image_delete(tmp);
        spec_sw_data = cpl_vector_get_data(spec_sw);
        for (j=sw_start;j<sw_end;j++) ycen_sw[j-sw_start] = ycen_rest[j];

        /* Finally ready to call the slit-decomp */
        slit_func_vert(swath, height, oversample, img_sw_data, mask_sw,
                        ycen_sw, slitfu_sw_data, spec_sw_data, model_sw,
                        0.0, smooth_slit, 1.0e-5, 20);

        for(col=0; col<swath; col++){      // col is x-index in cut-out
            for(row=0;row<height;row++){   // row is y-index in cut-out
                x = i*swath + col;          // coords in large image
                y = ycen_int[x] - (height/2) + row;
                cpl_image_set(img_out,x+1,y+1, model_sw[row*swath+col]);
            }
        }

        if (i==0) cpl_vector_copy(slitfu,slitfu_sw);
        else cpl_vector_add(slitfu,slitfu_sw);

        for (j=sw_start;j<sw_end;j++) {
            cpl_vector_set(spc, j, cpl_vector_get(spec_sw,j-sw_start));
        }
        cpl_vector_delete(spec_sw);
    } // End loop over swaths
    cpl_vector_delete(slitfu_sw);

    // divide by nswaths to make the slitfu into the average over all swaths.
    cpl_vector_divide_scalar(slitfu,nswaths);

    // TODO: Update BPM in img_out
    // TODO: Calculate error and return it.

    // TODO: Deallocate return arrays in case of error, return -1
    cpl_image_delete(img_sw);
    cpl_free(mask_sw) ;
    cpl_free(model_sw) ;
    cpl_free(ycen_int) ;
    cpl_free(ycen_rest);
    cpl_free(ycen_sw);

    *slit_func = slitfu;
    *spec = spc;
    *model = hdrl_image_create(img_out, NULL);

    // If I don't read hrdl wrong, hdrl_image_create() does not copy memory but wrap,
    // so no delete needed here.
    //cpl_image_delete(img_out) ;

    return 0;
}

/** @} */

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param    ncols       Swath width in pixels
  @param    nrows       Extraction slit height in pixels
  @param    osample     Subpixel ovsersampling factor
  @param    im          Image to be decomposed
  @param    mask        int mask of same dimension as image
  @param    ycen        Order centre line offset from pixel row boundary
  @param    sL          Slit function resulting from decomposition, start
                        guess is input, gets overwriteten with result
  @param    sP          Spectrum resulting from decomposition
  @param    model       the model reconstruction of im
  @param    lambda_sP   Smoothing parameter for the spectrum, could be zero
  @param    lambda_sL   Smoothing parameter for the slit function, usually >0
  @param    sP_stop     Fraction of spectyrum change, stop condition
  @param    maxiter     Max number of iterations
  @return
 */
/*----------------------------------------------------------------------------*/
static int slit_func_vert(
        int         ncols,
        int         nrows,
        int         osample,
        double  *   im,
        int     *   mask,
        double  *   ycen,
        double  *   sL,
        double  *   sP,
        double  *   model,
        double      lambda_sP,
        double      lambda_sL,
        double      sP_stop,
        int         maxiter)
{
    int x, y, iy, jy, iy1, iy2, ny, nd, i, j;
	double step, d1, d2, sum, norm, dev, lambda, diag_tot, sP_change, sP_max;
	int info, iter, isum;
    /* Initialise */
    nd=2*osample+1;
	ny=osample*(nrows+1)+1; /* The size of the sf array */
    step=1.e0/osample;
    double * E = cpl_malloc(ncols*sizeof(double)); // double E[ncols];
    double * sP_old = cpl_malloc(ncols*sizeof(double)); // double sP_old[ncols];
    double * Aij = cpl_malloc(ny*ny*sizeof(double)); // double Aij[ny*ny];
    double * bj = cpl_malloc(ny*sizeof(double)); // double bj[ny];
    double * Adiag = cpl_malloc(ncols*3*sizeof(double)); // double Adiag[ncols*3];
    double * omega = cpl_malloc(ny*nrows*ncols*sizeof(double)); // double omega[ny][nrows][ncols];
    // index as: [iy+(y*ny)+(x*ny*nrows)]

    /*
      Construct the omega tensor. Normally it has the dimensionality of
      ny*nrows*ncols.
      The tensor is mostly empty and can be easily compressed to ny*nx, but
      this will complicate matrix operations at later stages. I will keep
      it as it is for now.
      Note, that omega is used in in the equations for sL, sP and for the model
      but it does not involve the data, only the geometry. Thus it can be
      pre-computed once.
      */
    for(x=0; x<ncols; x++) {
		iy2=(1.e0-ycen[x])*osample;
        /*
           The initial offset should be reconsidered.
           It looks fine but needs theory.
         */
		iy1=iy2-osample;
		if(iy2==0) d1=step;
		else if(iy1==0) d1=0.e0;
		else d1=fmod(ycen[x], step);
		d2=step-d1;
		for(y=0; y<nrows; y++)
		{
			iy1+=osample;
			iy2+=osample;
			for(iy=0; iy<ny; iy++)
			{
				if(iy<iy1) omega[iy+(y*ny)+(x*ny*nrows)]=0.;
				else if(iy==iy1) omega[iy+(y*ny)+(x*ny*nrows)]=d1;
				else if(iy>iy1 && iy<iy2) omega[iy+(y*ny)+(x*ny*nrows)]=step;
				else if(iy==iy2) omega[iy+(y*ny)+(x*ny*nrows)]=d2;
				else omega[iy+(y*ny)+(x*ny*nrows)]=0.;
			}
		}
	}

    /* Loop through sL , sP reconstruction until convergence is reached */
	iter=0;
    do
    {
        /* Compute slit function sL */

        /* Fill in SLE arrays */
    	diag_tot=0.e0;
        for(iy=0; iy<ny; iy++)
        {
            bj[iy]=0.e0;
            for(jy=max(iy-osample,0); jy<=min(iy+osample,ny-1); jy++)
            {
            /* printf("iy=%d jy=%d %d\n", iy, jy, iy+ny*(jy-iy+osample)); */
                Aij[iy+ny*(jy-iy+osample)]=0.e0;
                for(x=0; x<ncols; x++)
                {
                    sum=0.e0;
                   for(y=0; y<nrows; y++)
                       sum+=omega[iy+(y*ny)+(x*ny*nrows)]*omega[jy+(y*ny)+(x*ny*nrows)]*mask[y*ncols+x];
                   Aij[iy+ny*(jy-iy+osample)]+=sum*sP[x]*sP[x];
                }
            }
            for(x=0; x<ncols; x++)
           	{
           		sum=0.e0;
                for(y=0; y<nrows; y++) sum+=omega[iy+(y*ny)+(x*ny*nrows)]*mask[y*ncols+x]*im[y*ncols+x];
                bj[iy]+=sum*sP[x];
            }
            diag_tot+=Aij[iy+ny*osample];
        }
        // printf("SUM : %e\n", sum);

        /* Scale regularization parameters */
	    lambda=lambda_sL*diag_tot/ny;

        /* Add regularization parts for the slit function */
        Aij[ny*osample]    +=lambda;           /* Main diagonal  */
        Aij[ny*(osample+1)]-=lambda;           /* Upper diagonal */
        for(iy=1; iy<ny-1; iy++)
        {
            Aij[iy+ny*(osample-1)]-=lambda;      /* Lower diagonal */
            Aij[iy+ny*osample    ]+=lambda*2.e0; /* Main diagonal  */
            Aij[iy+ny*(osample+1)]-=lambda;      /* Upper diagonal */
        }
        Aij[ny-1+ny*(osample-1)]-=lambda;      /* Lower diagonal */
        Aij[ny-1+ny*osample]    +=lambda;      /* Main diagonal  */

        /* Solve the system of equations */
        info=bandsol(Aij, bj, ny, nd);
        if(info) printf("info(sL)=%d\n", info);

        /* Normalize the slit function */
        norm=0.e0;
        for(iy=0; iy<ny; iy++)
        {
            sL[iy]=bj[iy];
            norm+=sL[iy];
        }
        norm/=osample;
        for(iy=0; iy<ny; iy++) sL[iy]/=norm;

        /* Compute spectrum sP */
        for(x=0; x<ncols; x++)
        {
            Adiag[x+ncols]=0.e0;
        	E[x]=0.e0;
        	for(y=0; y<nrows; y++)
            {
            	sum=0.e0;
        	    for(iy=0; iy<ny; iy++)
        	    {
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*sL[iy];
        	    }

                Adiag[x+ncols]+=sum*sum*mask[y*ncols+x];
                E[x]+=sum*im[y*ncols+x]*mask[y*ncols+x];
            }
        }
        if(lambda_sP>0.e0)
        {
        	norm=0.e0;
        	for(x=0; x<ncols; x++)
        	{
        		sP_old[x]=sP[x];
        		norm+=sP[x];
        	}
        	norm/=ncols;
        	lambda=lambda_sP*norm;
            Adiag[0        ] = 0.e0;
            Adiag[0+ncols  ]+= lambda;
            Adiag[0+ncols*2] =-lambda;
        	for(x=1; x<ncols-1; x++)
        	{
        		Adiag[x]=-lambda;
                Adiag[x+ncols  ]+= 2.e0*lambda;
                Adiag[x+ncols*2] =-lambda;
        	}
            Adiag[ncols-1        ] =-lambda;
            Adiag[ncols*2-1+ncols]+= lambda;
            Adiag[ncols*3-1+ncols] = 0.e0;

            info=bandsol(Adiag, E, ncols, 3);
            for(x=0; x<ncols; x++) sP[x]=E[x];
        }
        else
        {
        	for(x=0; x<ncols; x++)
        	{
        	    sP_old[x]=sP[x];
                sP[x]=E[x]/Adiag[x+ncols];
        	}
        }

        /* Compute the model */
  	    for(y=0; y<nrows; y++)
  	    {
            for(x=0; x<ncols; x++)
            {
        	    sum=0.e0;
        	    for(iy=0; iy<ny; iy++) sum+=omega[iy+(y*ny)+(x*ny*nrows)]*sL[iy];
        	    model[y*ncols+x]=sum*sP[x];
            }
        }

        /* Compare model and data */
        sum=0.e0;
        isum=0;
        for(y=0; y<nrows; y++)
        {
        	for(x=0;x<ncols; x++)
        	{
                sum+=mask[y*ncols+x]*(model[y*ncols+x]-im[y*ncols+x]) *
                                (model[y*ncols+x]-im[y*ncols+x]);
                isum+=mask[y*ncols+x];
        	}
        }
        dev=sqrt(sum/isum);

        /* Adjust the mask marking outlyers */
        for(y=0; y<nrows; y++)
        {
        	for(x=0;x<ncols; x++)
        	{
                if(fabs(model[y*ncols+x]-im[y*ncols+x])>6.*dev)
                    mask[y*ncols+x]=0;
                else mask[y*ncols+x]=1;
        	}
        }

        /* Compute the change in the spectrum */
        sP_change=0.e0;
        sP_max=1.e0;
        for(x=0; x<ncols; x++)
        {
            if(sP[x]>sP_max) sP_max=sP[x];
            if(fabs(sP[x]-sP_old[x])>sP_change) sP_change=fabs(sP[x]-sP_old[x]);
        }
        /* Check the convergence */
    } while(iter++ < maxiter && sP_change > sP_stop*sP_max);

    cpl_free(E);
    cpl_free(sP_old);
    cpl_free(omega);
    cpl_free(Aij);
    cpl_free(bj);
    cpl_free(Adiag);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   solve a sparse system of linear equations
  @param    a
  @param    r
  @param    n
  @param    nd
  @return
   bandsol solve a sparse system of linear equations with band-diagonal matrix.
   Band is assumed to be symmetrix relative to the main diaginal.
   Parameters are:
         a is 2D array [n,nd] where n - is the number of equations and nd
           is the width of the band (3 for tri-diagonal system).
           nd must be an odd number. The main diagonal should be in a(*,nd/2)
           The first lower subdiagonal should be in a(1:n-1,nd/2-1), the first
           upper subdiagonal is in a(0:n-2,nd/2+1) etc. For example:
                  / 0 0 X X X \
                  | 0 X X X X |
                  | X X X X X |
                  | X X X X X |
              A = | X X X X X |
                  | X X X X X |
                  | X X X X X |
                  | X X X X 0 |
                  \ X X X 0 0 /
         r is the array of RHS of size n.
   bandsol returns 0 on success, -1 on incorrect size of "a" and -4 on
   degenerate matrix.
 */
/*----------------------------------------------------------------------------*/
static int bandsol(
        double  *   a,
        double  *   r,
        int         n,
        int         nd)
{
    double aa;
    int i, j, k;

    //  if(mod(nd,2)==0) return -1;

    /* Forward sweep */
    for(i=0; i<n-1; i++)
    {
        aa=a[i+n*(nd/2)];
        //    if(aa==0.e0) return -3;
        r[i]/=aa;
        for(j=0; j<nd; j++) a[i+j*n]/=aa;
        for(j=1; j<min(nd/2+1,n-i); j++)
        {
            aa=a[i+j+n*(nd/2-j)];
            //      if(aa==0.e0) return -j;
            r[i+j]-=r[i]*aa;
            for(k=0; k<n*(nd-j); k+=n) a[i+j+k]-=a[i+k+n*j]*aa;
        }
    }

    /* Backward sweep */
    r[n-1]/=a[n-1+n*(nd/2)];
    for(i=n-1; i>0; i--)
    {
        for(j=1; j<=min(nd/2,i); j++) r[i-j]-=r[i]*a[i-j+n*(nd/2+j)];
        //    if(a[i-1+n*(nd/2)]==0.e0) return -5;
        r[i-1]/=a[i-1+n*(nd/2)];
    }

    //  if(a[n*(nd/2)]==0.e0) return -6;
    r[0]/=a[n*(nd/2)];
    return 0;
}
