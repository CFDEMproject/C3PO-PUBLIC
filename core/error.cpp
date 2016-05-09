/*-----------------------------------------------------------------------------*\
                  ___   _____   _____   _____     ___   
                 / ___\/\  __`\/\  __`\/\  __`\  / __`\ 
                /\ \__/\ \ \_\ \ \ \_\ \ \ \_\ \/\ \_\ \
                \ \____\\ \  __/\ \  __/\ \  __/\ \____/
                 \/____/ \ \ \/  \ \ \/  \ \ \/  \/___/ 
                          \ \_\   \ \_\   \ \_\         
                           \/_/    \/_/    \/_/         

         A Compilation for Fluid-Particle Data Post PrOcessing

Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
               2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
---------------------------------------------------------------------------------
License
    CPPPO is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with CPPPO. If not, see <http://www.gnu.org/licenses/lgpl.html>.

	This code is designed for on-the-fly post processing of fluid-particle
	data (e.g., of velocity, pressure, concentration, temperature field).

	Parts of the code were developed in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------*/


#include "error.h"
#include "input.h"
#include "output.h"
#include "comm.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

Error::Error(c3po *ptr) : c3poBase(ptr)
{

}

/* ----------------------------------------------------------------------
   called by one proc
------------------------------------------------------------------------- */

void Error::throw_error_one(const char *file, int line,const char *msg1,
                            const char *msg2,const char *msg3,const char *msg4) const
{

    char errstr[500];
    sprintf(errstr,"\n\nERROR on process %d (in file %s, line %d): ",comm().me(),file,line);


    // just one proc is calling, so call all() write
    output().write_screen_all(errstr);
    output().write_screen_all(msg1);
    if(msg2) output().write_screen_all(msg2);
    if(msg3) output().write_screen_all(msg3);
    if(msg4) output().write_screen_all(msg4);
    output().write_log_all(errstr);
    output().write_log_all(msg1);
    if(msg2) output().write_log_all(msg2);
    if(msg3) output().write_log_all(msg3);
    if(msg4) output().write_log_all(msg4);

    
    comm().abort_one();

}

/* ----------------------------------------------------------------------
   called by all proc
------------------------------------------------------------------------- */

void Error::throw_error_all(const char *file, int line,const char *msg1,
                            const char *msg2,const char *msg3,const char *msg4) const
{

    char errstr[500];
    sprintf(errstr,"\n\nERROR on process %d (in file %s, line %d): ",comm().me(),file,line);


    // all procs are calling, just 0 should write
    output().write_screen_one(errstr);
    output().write_screen_one(msg1);
    if(msg2) output().write_screen_all(msg2);
    if(msg3) output().write_screen_all(msg3);
    if(msg4) output().write_screen_all(msg4);
    output().write_log_one(errstr);
    output().write_log_one(msg1);
    if(msg2) output().write_log_all(msg2);
    if(msg3) output().write_log_all(msg3);
    if(msg4) output().write_log_all(msg4);

    comm().finalize_all();

}
