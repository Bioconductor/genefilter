/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001   The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "R.h"
#include "genefilter.h"
#include "R_ext/Rdynload.h"

static const R_CMethodDef CEntries[] = {
/*
  Called with NAOK and DUP which are passed down to give 8 args.
  Fix naokfind()
    {"mm_distance", (DL_FUNC) &mm_distance, 6},  
*/
    {NULL, NULL, 0}
};

void R_init_genefilter(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}
