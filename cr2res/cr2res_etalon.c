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

#include <cpl.h>
#include "cr2res_etalon.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_mask * cr2res_etalon_binary_image(const cpl_image *) ;
static cpl_vector * cr2res_etalon_get_maxpos(const cpl_vector *) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_etalon       Etalon related
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the etalon Image, and detects the blobs
  @param    in      The Etalon image
  @return   The labels int image or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_etalon_computation(const cpl_image * in)
{
    cpl_mask        *   mask ;
    cpl_size            nlabels ;
    cpl_image       *   labels ;
    cpl_apertures   *   aperts ;

    /* Compute Binary image */
    if ((mask = cr2res_etalon_binary_image(in)) == NULL) {
        cpl_msg_error(__func__, "Cannot Binarise the image") ;
        return NULL ;
    }
        
    /* Debugging */
    cpl_mask_save(mask, "mask.fits", NULL, CPL_IO_CREATE) ;

    /* Labelise the different detected apertures */
    if ((labels = cpl_image_labelise_mask_create(mask, &nlabels))==NULL) {
        cpl_msg_error(cpl_func, "Cannot Labelise") ;
        cpl_mask_delete(mask) ;
        return NULL ;
    }
    cpl_mask_delete(mask) ;

    cpl_msg_debug(__func__, "Number of Apertures: %"CPL_SIZE_FORMAT, nlabels) ;

    /* Create the detected apertures list */
    if ((aperts = cpl_apertures_new_from_image(in, labels)) == NULL) {
        cpl_msg_error(cpl_func, "Cannot Compute the apertures") ;
        cpl_image_delete(labels) ;
        return NULL ;
    }

    /* Display Apertures */
    cpl_apertures_dump(aperts, stdout) ;

    /* Free and return */
    cpl_apertures_delete(aperts) ;
    return labels ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a mask identifying the blogs (all separated)
  @param    in      The Etalon image
  @return   The mask
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_etalon_binary_image(const cpl_image * in)
{
    cpl_mask        *   mask ;
    cpl_vector      *   collapsed_vec ;
    cpl_image       *   collapsed_im ;
    cpl_vector      *   local_maxima ;
    double          *   pmax ;
    cpl_image       *   blob_image ;
    cpl_mask        *   blob_mask ;
    double              threshold ;
    int                 i, j, start_x, end_x ;

    /* Check Entries */
    if (in == NULL) return NULL; 

    /*************************/
    /* Get Local Maxima list */
    /*************************/
    /* Collapse in y */
    collapsed_im = cpl_image_collapse_create(in, 0) ;
    
    /* Smooth */
    /* TODO */

    /* Convert as vector */
    collapsed_vec = cpl_vector_new_from_image_row(collapsed_im, 1) ;
    cpl_image_delete(collapsed_im) ;

    /* Debugging */
    cpl_plot_vector(NULL, "w lines", NULL, collapsed_vec) ;

    /* Store local maxima positions of the plot */
    local_maxima = cr2res_etalon_get_maxpos(collapsed_vec) ;
    pmax = cpl_vector_get_data(local_maxima) ;
    cpl_vector_delete(collapsed_vec) ;

    /* Create the output mask  */
    mask = cpl_mask_new(cpl_image_get_size_x(in), cpl_image_get_size_y(in)) ;
   
    /* Loop on the Maxima positions and isolate the blob */
    for (i=0 ; i<cpl_vector_get_size(local_maxima) ; i++) {

        /* Get start_x end_x */
        start_x = 1;
        end_x = cpl_image_get_size_x(in) ;
        if (i>0) 
            start_x = (int)(pmax[i-1] + (pmax[i]-pmax[i-1])/2.) ;
        if (i<cpl_vector_get_size(local_maxima)-1)
            end_x = (int)(pmax[i] + (pmax[i+1]-pmax[i])/2.) ;

        cpl_msg_debug(__func__, "Extract %d -> %d", start_x, end_x) ;

        /* Extact the current blob */
        blob_image = cpl_image_extract(in, start_x, 1, end_x,
                cpl_image_get_size_y(in)) ;

        /* Compute the threshold */
        /* TODO */
        threshold = cpl_image_get_mean(blob_image) ;
        
        cpl_msg_debug(__func__, "Threshold: %g", threshold) ;

        /* Binarise the image */
        blob_mask = cpl_mask_threshold_image_create(blob_image, threshold, 
                DBL_MAX) ;
        cpl_image_delete(blob_image); 

        /* Set the left column to 0 to separate from the neighbor */
        for (j=0 ; j<cpl_mask_get_size_y(blob_mask) ; j++) 
            cpl_mask_set(blob_mask, 1, j+1, CPL_BINARY_0) ;

        /* Fill the Binary with the current blob */
        cpl_mask_copy(mask, blob_mask, start_x, 1); 

        cpl_mask_delete(blob_mask) ;
    }
    /* Free and return */
    cpl_vector_delete(local_maxima) ;

    return mask ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect maxima from a 1d periodic signal and store their positions
  @param    in      The 1d signal as a vector
  @return   The vector with the maxima positions
 */
/*----------------------------------------------------------------------------*/
static cpl_vector * cr2res_etalon_get_maxpos(const cpl_vector * in)
{
    const double    *   pin ;
    cpl_vector      *   maxima_pos ;
    double          *   pmax ;
    int                 i, nb_max ;

    /* Check Entries */
    if (in==NULL) return NULL ;

    /* Initialise */
    pin = cpl_vector_get_data_const(in) ;

    /* Count the number of Max positions */
    nb_max = 0 ;
    for (i=1 ; i<cpl_vector_get_size(in)-1 ; i++) {
        if (pin[i] > pin[i+1] && pin[i] > pin[i-1]) nb_max++ ;
    }

    /* Create the output vector */
    maxima_pos = cpl_vector_new(nb_max) ;
    pmax = cpl_vector_get_data(maxima_pos) ;
    nb_max = 0 ;
    for (i=1 ; i<cpl_vector_get_size(in)-1 ; i++) {
        if (pin[i] > pin[i+1] && pin[i] > pin[i-1]) {
            /* Store the Max X position */
            pmax[nb_max] = i+1 ;
            nb_max++ ;
        }
    }
    return maxima_pos ;
}
