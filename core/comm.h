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

#ifndef C3PO_COMM_H
#define C3PO_COMM_H

#include "stdio.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include "mpi.h"

#include <QtCore/QJsonArray>
#include <QtCore/QJsonObject>
#include <QtCore/QJsonDocument>
#include <QtCore/QString>
#include <QtCore/QFile>

namespace C3PO_NS
{

class Comm : public c3poBase, public c3poBaseInterface
{
    public:

      Comm(c3po *ptr,MPI_Comm &communicator);
      ~Comm();

      // inherited from c3poBaseInterface

      void init();
      void parse_command(int narg,char const* const* arg);

      // inline access

      inline int me()           const { return me_;}
      inline int nprocs()       const { return nprocs_;}
      inline bool is_proc_0()   const { return (0 == me_)    ? true: false; }
      inline bool is_parallel() const { return (nprocs_ > 1) ? true: false; }
      inline MPI_Comm world()   const { return world_; }            

      //TODO route this to CommModel
      inline void exchange () {}
      
      inline void bcast_int(int  *value) const 
      {
         int size = sizeof(value);
         MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);        
         MPI_Bcast(value,size,MPI_INT,0,MPI_COMM_WORLD);
      };

      inline void bcast_char(char *value) const //Use also for QByteArray!
      {
         int size = sizeof(value);
         MPI_Bcast(&size,1,MPI_INT,0,world_);        
         MPI_Bcast(value,size,MPI_CHAR,0,world_);
      };
      
      inline void allreduce_int(int send,int receive) const {MPI_Allreduce(&send, &receive, 1, MPI_INT, MPI_SUM, world_);};
      
      inline void barrier() const {MPI_Barrier(world_);};
      
   //   MPI::Win* goRMA(double*) const;
   //   MPI::Win* goRMA(int*) const;

      void abort_one() const;
      void finalize_all() const;

    private:

      //TODO
      //class CommModel *commModel_;
      // can be CommModelCartesian or CommModelMany2Many

      int me_, nprocs_;        // number of MPI procs and my rank
      MPI_Comm &world_;        // MPI communicator
      
};

} //end c3po_NS

#endif
