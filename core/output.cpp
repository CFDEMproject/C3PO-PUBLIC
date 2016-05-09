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


#include "output.h"
#include "comm.h"
#include "version.h"
#include <cstring>
#include <stdlib.h>
#include "error.h"
#include "timer.h"
#include <sys/stat.h>


using namespace C3PO_NS;

struct stat dirLogic;
/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

Output::Output(c3po *ptr) : c3poBase(ptr),
    screen_(0),
    logfile_(0)
{

    version_ = new char[strlen(C3PO_VERSION)+100];
    sprintf(version_,"Version: %s",C3PO_VERSION);

    screen_ = stdout;
    if (0 == screen_)
        printf("Cannot open screen");
    else
        fprintf(screen_,"C3PO (%s)\n",version_);//", screen %d, this %d\n", version_,screen_,this);

    logfile_ = fopen("log.c3po","w");
    if (0 == logfile_)
        write_screen_one("Cannot open logfile");
    else
        fprintf(logfile_,"c3po (%s)\n",version_);
}

Output::~Output()
{

    if (screen_) fclose(screen_);
    screen_ = 0;
    if (logfile_) fclose(logfile_);
    logfile_ = 0;
}

/* ----------------------------------------------------------------------
   write out to screen
------------------------------------------------------------------------- */
void Output::write_screen_one( const char *message) const
{
    if(comm().is_proc_0())
        printf("%s\n",message);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_screen_one_int(const int number, const char *message) const
{
    if(comm().is_proc_0())
        printf("%s %i\n",message, number);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_screen_one_int2(const int number, const char *message,const char *message2) const
{
    if(comm().is_proc_0())
        printf("%s %i %s\n",message, number,message2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_screen_all( const char *message) const
{
    fprintf(screen_,"[Process %d/%d] %s\n",comm().me(),comm().nprocs(),message);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_screen_all_int(const int number, const char *message) const
{
    fprintf(screen_,"[Process %d] %s %i\n",comm().me(),message, number);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_log_one(const char *message) const
{
    if(comm().is_proc_0())
        fprintf(logfile_,"%s\n",message);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_log_all(const char *message) const
{
    fprintf(logfile_,"[Process %d] %s\n",comm().me(),message);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::write_time_one(double *message)
{
    fprintf(logfile_,"%g \n",*message);
}

/* -------------------------------------------------------------------------- */
//              Write to Json
/* --------------------------------------------------------------------------- */ 
void Output::createQJsonArrays(std::string fileName_,std::string mainObject_, std::vector<std::string> _ObjectName,
                               std::vector<double**> _data,
                               int _datanum, int _datarows, int spacingJ,
                               bool overwrite) 
{

 QFile saveFile( fileName_.c_str() );
 if (overwrite) 
 { 
   if (!saveFile.open(QIODevice::WriteOnly))
     printf("Can not open savefile: %s",fileName_.c_str() );
 }
 else
 {
   if (!saveFile.open(QIODevice::ReadWrite))
     printf("Can not open savefile: %s",fileName_.c_str() );
     QByteArray    saveData = saveFile.readAll();
     QJsonDocument loadDoc  = QJsonDocument::fromJson(saveData);
 }
   QJsonObject * Object_ = new QJsonObject;
  QJsonObject * MainObject_= new QJsonObject;
  
  for(unsigned int i = 0; i < _ObjectName.size(); i++)
  {
         double ** currDat = _data[i];
         QJsonArray dataArray;
         QString qstr(_ObjectName[i].c_str());
         for (int j = 0; j < _datanum; j++)
         {
           QJsonArray dataArraySub;
           for(int rowI=0; rowI < _datarows;rowI++)
               dataArraySub << QJsonValue(currDat[rowI][j*spacingJ]);
        
           dataArray << dataArraySub;

         } 
         (*Object_)[qstr] = dataArray;     
  }

  (*MainObject_)[mainObject_.c_str()]=(*Object_);
  QJsonDocument * saveDoc = new QJsonDocument(*MainObject_);

  saveFile.write( saveDoc->toJson() );
  saveFile.close();

  delete Object_;
  delete MainObject_;
  delete saveDoc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void Output::createQJsonArrays(std::string fileName_,std::string mainObject_, std::vector<std::string> ObjectName_,std::vector<double*> data_,int datanum, bool overwrite,std::vector<int> * datanumvec_) 
{

 QFile saveFile( fileName_.c_str() );
 if (overwrite) 
 { 
   if (!saveFile.open(QIODevice::WriteOnly))
     printf("Can not open savefile: %s",fileName_.c_str() );
 }
 else
 {
   if (!saveFile.open(QIODevice::ReadWrite))
     printf("Can not open savefile: %s",fileName_.c_str() );
     QByteArray    saveData = saveFile.readAll();
     QJsonDocument loadDoc  = QJsonDocument::fromJson(saveData);
 }
  QJsonObject * Object_ = new QJsonObject;
  QJsonObject * MainObject_= new QJsonObject;
  
  for(unsigned int i = 0; i < ObjectName_.size(); i++)
      {
         QJsonArray dataArray;
         QString qstr(ObjectName_[i].c_str());

         int datanum_;
         if(datanum==-1 && datanumvec_!=NULL) datanum_=(*datanumvec_)[i];
         else if (datanum > -1) datanum_=datanum;
         else datanum_=0;
         
         for (int j = 0;j < datanum_; j++)
           dataArray << QJsonValue(data_[i][j]);

         (*Object_)[qstr] = dataArray;     
      }

 (*MainObject_)[mainObject_.c_str()]=(*Object_);
 QJsonDocument * saveDoc = new QJsonDocument(*MainObject_);

 saveFile.write( saveDoc->toJson() );
 saveFile.close();
 
 delete Object_;
  delete MainObject_;
  delete saveDoc;
    
}

void Output::createQJsonArrays(std::string fileName_,std::string mainObject_, std::vector<std::string> ObjectName_,std::vector<int*> data_,int datanum, bool overwrite,std::vector<int> * datanumvec_) 
{

 QFile saveFile( fileName_.c_str() );
 if (overwrite) 
 { 
   if (!saveFile.open(QIODevice::WriteOnly))
     printf("Can not open savefile: %s",fileName_.c_str() );
 }
 else
 {
   if (!saveFile.open(QIODevice::ReadWrite))
     printf("Can not open savefile: %s",fileName_.c_str() );
     QByteArray    saveData = saveFile.readAll();
     QJsonDocument loadDoc  = QJsonDocument::fromJson(saveData);
 }
  QJsonObject * Object_ = new QJsonObject;
  QJsonObject * MainObject_= new QJsonObject;
  
  for(unsigned int i = 0; i < ObjectName_.size(); i++)
      {
         QJsonArray dataArray;
         QString qstr(ObjectName_[i].c_str());

         int datanum_;
         if(datanum==-1 && datanumvec_!=NULL) datanum_=(*datanumvec_)[i];
         else if (datanum > -1) datanum_=datanum;
         else datanum_=0;
         
         for (int j = 0;j < datanum_; j++)
           dataArray << QJsonValue(data_[i][j]);

         (*Object_)[qstr] = dataArray;     
      }

 (*MainObject_)[mainObject_.c_str()]=(*Object_);
 QJsonDocument * saveDoc = new QJsonDocument(*MainObject_);

 saveFile.write( saveDoc->toJson() );
 saveFile.close();
 
 delete Object_;
  delete MainObject_;
  delete saveDoc;
    
}


/* ---------------------------------------------------------------------------*/
void Output::generateDir(std::string dir, bool &check) const
{
    if( (!check) && comm().is_proc_0() )
    {

        if(stat(dir.c_str(), &dirLogic) == 0 && S_ISDIR(dirLogic.st_mode))
            std::cout << ""; 
        else
        {
            std::string command("mkdir "+dir);
            system(command.c_str());
        }
        check = true;
    }

    MPI_Barrier(MPI_COMM_WORLD); //must sync here!

}

