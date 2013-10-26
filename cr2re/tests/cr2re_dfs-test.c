/* $Id: cr2re_dfs-test.c,v 1.3 2013-09-13 08:16:32 yjung Exp $
 *
 * This file is part of the CR2RE Pipeline
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

/*
 * $Author: yjung $
 * $Date: 2013-09-13 08:16:32 $
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include <cpl.h>

#include <cr2re_dfs.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2re_dfs_test  Unit test of cr2re_dfs
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Textual representation of CPL frame group
  @param    group     to convert
  @return   textual representation
 */
/*----------------------------------------------------------------------------*/
static const char *
frame_group_to_string(cpl_frame_group group)
{
    const char * self = "???"; /* Only used if cpl_frame_group is extended */

    switch(group) {
    case CPL_FRAME_GROUP_RAW:
        self = CPL_FRAME_GROUP_RAW_ID;
        break;
    case CPL_FRAME_GROUP_NONE:
        self = "NONE";
        break;
    case CPL_FRAME_GROUP_CALIB:
        self = CPL_FRAME_GROUP_CALIB_ID;
        break;
    case CPL_FRAME_GROUP_PRODUCT:
        self = CPL_FRAME_GROUP_PRODUCT_ID;
        break;
    default:
        break;
    }

    return self;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit test of cr2re_dfs_set_groups
 */
/*----------------------------------------------------------------------------*/
static void test_set_groups(void)
{
    const char *const fctid = "test_set_groups";
    const char *const test_subject = "cr2re_dfs_set_groups";
    cpl_errorstate prestate = cpl_errorstate_get();

    /* Test with invalid input */
    if (cr2re_dfs_set_groups(NULL) == 0) {
        cpl_msg_error(fctid, "Function %s did not fail on NULL input",
                      test_subject);
        cpl_end();
        exit(EXIT_FAILURE);       
    }

    cpl_errorstate_set(prestate);

    /* Test with valid input */
    {
        /* Simulate data */
        const int N = 3;
        const char *const filename[] = {"raw1.fits", 
                                        "raw2.fits", 
                                        "calib.fits"};
        const char *const tag[] = {CR2RE_BIAS_RAW,
                                   CR2RE_BIAS_RAW, 
                                   CR2RE_BIAS_PROCATG};
        cpl_frame_group const expected_group[] = {CPL_FRAME_GROUP_RAW,
                                                  CPL_FRAME_GROUP_RAW,
                                                  CPL_FRAME_GROUP_CALIB};
        cpl_frameset *frames = cpl_frameset_new();
        int i;
    
        for (i = 0; i < N; i++) {
            cpl_frame *frame = cpl_frame_new();

            cpl_frame_set_filename(frame, filename[i]);
            cpl_frame_set_tag(frame, tag[i]);
            cpl_frameset_insert(frames, frame);
        }

        /* Call the function */
        if (cr2re_dfs_set_groups(frames) != 0) {
            cpl_msg_error(fctid, "Function %s failed",
                          test_subject);
            cpl_errorstate_dump(prestate, CPL_FALSE, NULL);
            cpl_frameset_delete(frames);
            cpl_end();
            exit(EXIT_FAILURE);       
        }
        
        /* Verify results */
        for (i = 0; i < N; i++) {
            cpl_frame *frame = cpl_frameset_get_frame(frames, i);

            if (frame == NULL) {
                cpl_msg_error(fctid, "Missing frame number %d",
                              i);
                cpl_errorstate_dump(prestate, CPL_FALSE, NULL);
                cpl_frameset_delete(frames);
                cpl_end();
                exit(EXIT_FAILURE);       
            }

            if (cpl_frame_get_group(frame) != expected_group[i]) {
                cpl_msg_error(fctid, "Frame number %d has group %s, "
                              "%s expected",
                              i,
                              frame_group_to_string(cpl_frame_get_group(frame)),
                              frame_group_to_string(expected_group[i]));
                cpl_errorstate_dump(prestate, CPL_FALSE, NULL);
                cpl_frameset_delete(frames);
                cpl_end();
                exit(EXIT_FAILURE);
            }
        }
        
        cpl_frameset_delete(frames);
    }
    
    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Unit tests of cr2re_dfs module
 */
/*----------------------------------------------------------------------------*/

int main(void)
{
#if defined CPL_VERSION_CODE && CPL_VERSION_CODE >= CPL_VERSION(4, 0, 0)
    cpl_init(CPL_INIT_DEFAULT);
#else
    cpl_init();
#endif

    test_set_groups();

    cpl_end();
    exit(EXIT_SUCCESS);
}

/**@}*/
