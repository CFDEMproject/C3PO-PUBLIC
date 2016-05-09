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

#ifndef C3PO_TIMER_H
#define C3PO_TIMER_H

#include "c3po_base.h"
#include "c3po_base_interface.h"

enum{TIME_OVERALL,TIME_SELECTOR,TIME_SELECTOR_MPI, TIME_FILTER_MPI,TIME_FILTER, TIME_FILTER_TOTAL,
      TIME_OUTPUT,TIME_SYNC,TIME_SAMPLING,TIME_BINNING, TIME_INPUT,TIME_N};

namespace C3PO_NS {

class Timer  : public c3poBase, public c3poBaseInterface
{
 public:
  mutable double *array;

  Timer(c3po *ptr);
  ~Timer();
  void init();
  void stamp() const;
  void stamp(int) const;
  void synchronized_start(int);
  void synchronized_stop(int);
  double elapsed(int);

  void dumpStats();
  void globalReport();

 private:
  mutable double previous_time;

  std::string     timeDir_;
  bool            timeDirIsCreated_;
};

}

#endif

