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

#include "c3po.h"
#include "control.h"
#include "input.h"
#include "output.h"
#include "comm.h"
#include "error.h"
#include "mesh.h"
#include "operation_container.h"
#include "selector_container.h"
#include "data_storage.h"
#include "c3po_base_interface_vector.h"
#include "qjson_includes.h"
#include <cmath>
#include <ctime>
#include "timer.h"
#include "multiphase_FlowBasic.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   Central c3po Class
    - allocates memory, error, input, etc.
------------------------------------------------------------------------- */

c3po::c3po(int narg, char **arg, MPI_Comm communicator) :
  c3poBaseInterfaceVector_( new c3poBaseInterfaceVector),
  timer_             ( c3poBaseInterfaceVector_->add<Timer>(new Timer(this),"timer")),
  control_           ( c3poBaseInterfaceVector_->add<Control>(new Control(this),"control")),    //the std::string is that in the input script"
  comm_              ( c3poBaseInterfaceVector_->add<Comm>(new Comm(this,communicator),"comm")),    //the std::string is that in the input script"
  input_             ( c3poBaseInterfaceVector_->add<Input>(new Input(this),"input")),  //the std::string is that in the input script"
  output_            ( c3poBaseInterfaceVector_->add<Output>(new Output(this),"output")), //the std::string is that in the input script"
  error_             ( c3poBaseInterfaceVector_->add<Error>(new Error(this),"error")),  //the std::string is that in the input script"
  operationContainer_( c3poBaseInterfaceVector_->add<OperationContainer>(new OperationContainer(this),"operation" ) ), //the std::string is that in the input script"
  selectorContainer_ ( c3poBaseInterfaceVector_->add<SelectorContainer>(new SelectorContainer(this),"selector" ) ), //the std::string is that in the input script"
  dataStorageInternal_(c3poBaseInterfaceVector_->add<DataStorage>(new DataStorage(this),"storage")),  //the std::string is that in the input script"
  mesh_(c3poBaseInterfaceVector_->add<c3poMesh>(new c3poMesh(this),"mesh")),  //the std::string is that in the input script"
  basicMultiphaseQty_(c3poBaseInterfaceVector_->add<multiphaseFlowBasic>(new multiphaseFlowBasic(this),"basicMultiphase"))
{
    input_->process_input_script();
    dataStorageInternal_->init();
}

/* ----------------------------------------------------------------------
   shutdown c3po
   delete top-level classes
------------------------------------------------------------------------- */

c3po::~c3po()
{
    timer_->synchronized_stop(TIME_OVERALL);
    timer_->dumpStats();
    timer_->globalReport();

    delete c3poBaseInterfaceVector_;

    delete operationContainer_;
    delete selectorContainer_;
}

/* **********************************************************************
       MEMBER FUNCTIONS
********************************************************************** */
/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void c3po::init()
{
    printf("\n*******************************************\n");
    printf("\nCreating a c3po object \n");
    printf("Will initialize individual components of c3po now... \n");

    c3poBaseInterfaceVector_->init();


    printf("\n*******************************************\n\n");

}

// *************************************************************

int c3po::selectorCount() const
{
    return selectorContainer_->selectorCount();
}
// *************************************************************
bool c3po::verbose() const
{
    return input_->mainSettings()["verbose"].toBool();
}

// *************************************************************
double c3po::meshCheckTolerance() const
{
    return mesh_->meshCheckTolerance();
}

// *************************************************************
double c3po::meshFilterWidthTolerance() const
{
    return mesh_->meshFilterWidthTolerance();
}

// *************************************************************
bool c3po::meshVerbose() const
{
    return mesh_->meshVerbose();
}

// *************************************************************
void c3po::setCells(int * number_of_cells, int * global_number_of_cells, double * cellsize) const
{
    //Set global cell counts (Structured mesh only)
    mesh_->setCellsize(cellsize);
    mesh_->setcell(number_of_cells);
    mesh_->setcell_global(global_number_of_cells);

}
// *************************************************************
void c3po::preRunOperations() const
{

 dataStorageInternal_->readParticles();

 selectorContainer_->bubble_run();
 dataStorageInternal_->writeRegions();


}

// *************************************************************
void c3po::runFilters(int id) const
{
    timer_->synchronized_start(TIME_FILTER_TOTAL);
   // MPI_Barrier(MPI_COMM_WORLD);

    if(id==0) dataStorageInternal_->writeParticles(); //Write particle data just one time

    if(!input_->mainSettings()["doFiltering"].toBool())
        return;


    int MaxNofCells=*(mesh_->MaxNofCells());

    selectorContainer_->setFilter(id);
    int NofFilteringOp_=operationContainer_->filterCount();
    std::string filtName_(selectorContainer_->filter(id)->name());

    bool doParticle_= false;
    bool doCell_=false;

      for(int s=0;s<NofFilteringOp_;s++)
      {
       if(operationContainer_->filter(s)->particleBased())
        doParticle_= true;

       if(!(operationContainer_->filter(s)->particleBased()))
        doCell_= true;
      }
   //Begin of step
   for(int s=0;s<NofFilteringOp_;s++)
    operationContainer_->filter_begin_of_step(s);

 if(!dataStorageInternal_->useProbes()) doParticle_=false;

 if(doParticle_)
 {
  for(int probe=0;probe<dataStorageInternal_->NofProbes();probe++)
  {
    dataStorageInternal_->setCurrentProbes(probe);
    dataStorageInternal_->clearProbeFields();

    if(!dataStorageInternal_->runProbes(filtName_)) continue;

    int numPar_=dataStorageInternal_-> MaxNumOfParticles();
    selectorContainer_->isParticle();
   //Run particleBased Filtering operations

   for(int n=0;n<numPar_;n++)
   {

     selectorContainer_->setCurrentParticle(n);
      selectorContainer_->begin_of_step();



      for(int s=0;s<NofFilteringOp_;s++)
       if(operationContainer_->filter(s)->particleBased())
         operationContainer_->filter_middle_of_step(s);

   }

  }

 }

 if(doCell_)
 {
  #ifdef C3PO_DEBUG_TRUE
  int perc=int(MaxNofCells/100);
  #endif
  selectorContainer_->resetFSize();
  selectorContainer_->resetFilterVolume();
  selectorContainer_->isCell();
   //Run Eulerian Filtering operations
    for (int n=0;n<MaxNofCells;n++)
     {

      #ifdef C3PO_DEBUG_TRUE
      if(comm_->is_proc_0())
       if(n % perc ==0)
        std::cout <<  "\n cell filter " << filtName_ <<" running =====>  " << n/MaxNofCells*100 << "%";
      #endif


       selectorContainer_->setCurrentCell(n);

       selectorContainer_->begin_of_step();
       timer_->stamp();
       MPI_Barrier(MPI_COMM_WORLD);
       timer_->stamp(TIME_SYNC);

       for(int s=0;s<NofFilteringOp_;s++)
        if(!operationContainer_->filter(s)->particleBased())
         operationContainer_->filter_middle_of_step(s);

     }
  }

  //End of step
   for(int s=0;s<NofFilteringOp_;s++)
    operationContainer_->filter_end_of_step(s);

   selectorContainer_->resetFSize();
   selectorContainer_->resetFilterVolume();

   dataStorageInternal_->writeFields(filtName_);
   dataStorageInternal_->writeParticleFields(filtName_);


 if(input_->verbose())
    std::cout << "\nFiltering operations ended successfully for proc: " << comm_->me() << "\n";

  // MPI_Barrier(MPI_COMM_WORLD);

   timer_->synchronized_stop(TIME_FILTER_TOTAL);

};


// *************************************************************
void c3po::runSampling() const
{
    timer_->synchronized_start(TIME_SAMPLING);
    if(!input_->mainSettings()["doSampling"].toBool())
        return;
  //This function allows to use previously filtered fields for sampling
  dataStorageInternal_->refreshRMAfields();



  for(int s=0;s<operationContainer_->sampleCount();s++)
  {
    operationContainer_->sample_begin_of_step(s);
    operationContainer_->sample_middle_of_step(s);
    operationContainer_->sample_end_of_step(s);
  }
    timer_->synchronized_stop(TIME_SAMPLING);

   if(input_->verbose())
    std::cout << "\nSampling operations ended successfully for proc: " << comm_->me() << "\n";
}

// *************************************************************
void c3po::runBinning() const
{
    timer_->synchronized_start(TIME_BINNING);
    if(!input_->mainSettings()["doBinning"].toBool())
        return;

  for(int s=0;s<operationContainer_->binCount();s++)
  {
    //operationContainer_->bin_begin_of_step(s);
    //operationContainer_->bin_middle_of_step(s);
    operationContainer_->bin_end_of_step(s);
  }


   if(input_->verbose())
    std::cout << "\nBinning operations ended successfully for proc: " << comm_->me() << "\n";

   // MPI_Barrier(MPI_COMM_WORLD);
     timer_->synchronized_stop(TIME_BINNING);
}

// *************************************************************
void c3po::postRunOperations() const
{

 dataStorageInternal_->deleteParticles();

}

// *************************************************************
void c3po::flush() const
{
   for(int iOp=0; iOp<operationContainer_->operationCount(); iOp++)
        operationContainer_->operation(iOp)->flush();
}
/******************************************************************/

int c3po::getVFnamesNumber()
{

    return input_->getnumber_ofVF();

}

int c3po::getSFnamesNumber()
{

    return input_->getnumber_ofSF();

}


std::string c3po::getVFnames(int i)
{

    return input_->getVFname(i);

}

std::string c3po::getSFnames(int i)
{

    return input_->getSFname(i);

}

void c3po::registerVF(std::string name, double* x,double* y,double* z,int sp)
{

  dataStorageInternal_->addfVF(name, x, y, z, sp);


}

void c3po::GlobalVF(std::string name, double* x,double* y,double* z, int sp)
{

  dataStorageInternal_->addRMAvF(name, x, y, z, sp);

}


void c3po::registerSF(std::string name, double* x)
{

  dataStorageInternal_->addfSF(name, x);

}

void c3po::GlobalSF(std::string name, double* x)
{

  dataStorageInternal_->addRMAsF(name, x);

}

void c3po::addParticle(std::string        groupName_,
                       double                      m,
                       double*                   pos,
                       double*                   vel,
                       std::vector< double* >* force,
                       std::vector<double>*  scalars,
                       double*                torque
                      )
{
 dataStorageInternal_->addParticle(groupName_,m,pos,vel,force,scalars, torque);
}

void c3po::deleteParticles()
{
 dataStorageInternal_->deleteParticles();
}

void c3po::resetFields()
{

 dataStorageInternal_->deleteFields();

}

void c3po::resetGlobalFields()
{

 dataStorageInternal_->deleteRMA();

}

int c3po::getVFnumber()
{

 return dataStorageInternal_->numberOfVF();

}

int c3po::getSFnumber()
{

 return dataStorageInternal_->numberOfSF();

}


double c3po::getMaxDomainfromJson(int i) const
{
 return input_->getMaxDomainfromJson(i);
}

double c3po::getMinDomainfromJson(int i) const
{
 return input_->getMinDomainfromJson(i);
}

double c3po::cellSizefromJson(int i) const
{
 return input_->cellSizefromJson(i);
}

const char* c3po::getOpFilterName(int id)
{
 return operationContainer_->filter(id)->name();
}

void c3po::setTime(std::string t)
{
 dataStorageInternal_->setTimeName(t);
}

int c3po::getOpFiltNum()
{
 return operationContainer_->filterCount();
}

int c3po::getOpSampNum()
{
 return operationContainer_->sampleCount();
}

int c3po::getOpBinNum()
{
 return operationContainer_->binCount();
}

void c3po::addFieldToConvert(std::string name_)
{
 dataStorageInternal_->addFieldToConvert(name_);
}

void c3po::refreshFields()
{
 dataStorageInternal_->refreshRMAfields();
}

std::string c3po::getFilterName(int id)
{
 return selectorContainer_->filter(id)->name();
}

std::vector<std::string> c3po::vectorFTF(int id)
{
 return operationContainer_->filter(id)->returnVectorFTFNames();
}

std::vector<std::string> c3po::scalarFTF(int id)
{
 return operationContainer_->filter(id)->returnScalarFTFNames();
}

std::vector<std::string> c3po::vectorFTFVariance(int id)
{
 return operationContainer_->filter(id)->returnVectorFTFVarianceNames();
}

std::vector<std::string> c3po::scalarFTFVariance(int id)
{
 return operationContainer_->filter(id)->returnScalarFTFVarianceNames();
}

int c3po::numberOfFilters()
{
 return selectorContainer_->filterCount();
}

void c3po::registerDomainInfo(double maxDomain[3],double minDomain[3],double maxDomainGlobal[3],double minDomainGlobal[3])
{
 mesh_->registerDomainInfo( maxDomain, minDomain, maxDomainGlobal, minDomainGlobal);
}

void c3po::registerCells(double* Vvector, double* posVector, int NumberOfCells, int meshXYZspacing)
{
 mesh_->registerCells(Vvector,posVector,NumberOfCells,meshXYZspacing);
}

void c3po::setNofCells(const int x) const
{
 mesh_->setNofCells(x);
}

int c3po::getDomainDecomposition() const
{
 return input_->getDomainDecomposition();
}

bool c3po::writeFields() const
{
 return input_->writeInterface();
}

void c3po::csv_columnsToRead( int *columns) const
{
 input_->csv_columnsToRead(columns);
}
void c3po::checkIJK(bool struct_) const
{
 selectorContainer_->checkIJK(struct_);
}
void c3po::clearMesh() const
{
 mesh_->clearCells();
}
bool c3po::registerCFDEMparticles() const
{
 return input_->registerCFDEMparticles();
}
std::vector<std::string> c3po::getGradientScalarList() const
{
 return input_->getGradientScalarList();
}
std::vector<std::string> c3po::getGradientVectorList() const
{
 return input_->getGradientVectorList();
}

std::vector<std::string> c3po::getShearRateList() const
{
 return input_->getShearRateList();
}

/* ----------------------------------------------------------------------
   help message for command line options and styles present in executable
------------------------------------------------------------------------- */

void c3po::help()
{

  printf( //(screen,
          "\nCommand line options:\n\n"
          "-in filename                : read input from file, not stdin NOTIMPLEMENTED (-i)\n"
          "-help                       :  help message NOTIMPLEMENTED (-h)\n"
          "-log none/filename          : where to send log output NOTIMPLEMENTED (-l)\n"
          "-screen none/filename       : where to send screen output NOTIMPLEMENTED (-sc)\n"
          "-other command line Options          : NOTIMPLEMENTED\n\n");

}

/* ----------------------------------------------------------------------
Reporting Function - all the reporting functions should be called by this function
------------------------------------------------------------------------- */

void c3po::report()
{

 selectorContainer_->SCreport();

}
