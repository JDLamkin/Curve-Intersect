/*=============================================================================
  This file is part of SLV/

  Copyright (C) 2015, Elias Tsigaridas (Elias.Tsigaridas@inria.fr)

  SLV is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  SLV is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  If interested by a copy of the GNU Lesser General Public License, write to
  the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
  MA 02110-1301, USA. 
=============================================================================*/

/*

  Author: Elias Tsigaridas
  
  Comment: Insprired by "Zed's Awesome Debug Macros"

 */ 

#ifndef __DBG_H__
#define __DBG_H__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef NDEBUG

#define SLV_ERROR 

#define debugf(M, ...)
#define debug(M, ...)


#define clean_errno() 

#define log_err(M, ...)

#define log_warn(M, ...)

#define log_info(M, ...)

#define check(A, M, ...)

#define sentinel(M, ...)

#define check_mem(A) 

#define check_debug(A, M, ...)


#else

#define SLV_ERROR error:


#define debugf(M, ...) fprintf(stderr, "DBG %s:%s:%d: " M "\n", __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__)
#define debug(M, ...) fprintf(stderr, "DBG[ %s ]::%d:: " M "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__)


#define clean_errno() (errno == 0 ? "None" : strerror(errno))

#define log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

#define log_warn(M, ...) fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

#define log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define check(A, M, ...) if(!(A)) { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define sentinel(M, ...)  { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define check_mem(A) check((A), "Out of memory.")

#define check_debug(A, M, ...) if(!(A)) { debug(M, ##__VA_ARGS__); errno=0; goto error; }

#endif


#endif /* __DBG_H__ */
