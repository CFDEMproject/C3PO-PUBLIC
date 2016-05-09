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


#include "input.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include "stdio.h"
#include "input_properties.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <stdlib.h> 
#include "timer.h"
using std::cout;
using std::endl;

#include <cstring>

const int MAX_CHARS_PER_LINE = 1024;
const int MAX_TOKENS_PER_LINE = 50;


using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

Input::Input(c3po *ptr) : c3poBase(ptr),
    //infile_(0),
    myName("test"),
    writePrecision_(5),
    scope_("none"),
    dumpFormat_("json"),
    storageWriteFields_(false),
    storageWriteParticles_(false)
{
    // parse input switches - NOTIMPLEMENTED
 
    //printf("Input constructor called, I am %d, Comm is %d\n",this,comm());
    //printf("Input constructor called comm refernce to me %d\n",this,comm()->input());
}

Input::~Input()
{
    //if (infile_) infile_->close();
    //infile_ = 0;
}

/* ----------------------------------------------------------------------
   process an input script
------------------------------------------------------------------------- */

void Input::process_input_script()
{
     
    timer().stamp();     
    char filename [400];
    
    sprintf(filename,"./c3po_control/c3po.input");
    
    std::ifstream ifnfile (filename,std::ifstream::in);
    
     if (!ifnfile.is_open())
          output().write_screen_one("Cannot open input script");

    output().write_screen_one("\nc3po is opening the main JSON document...\n");
        
    loadDoc_ = openJsonFile("c3po_control",   "c3po", "mainSettings", mainObj_);
               openJsonFile("c3po_control", "mesh", "mainSettings", meshMainSettings_);
   
    readFieldNames();
    QString Qdf=mainSettings()["dumpFormat"].toString();
    std::string df=Qdf.toUtf8().constData();
    dumpFormat_.assign(df);

    storageWriteFields_    = mainSettings()["storageWriteFields"].toBool();
    storageWriteParticles_ = mainSettings()["storageWriteParticles"].toBool();
    
   #ifndef H5_LIB
    if(dumpFormat_.compare("hdf5")==0) error().throw_error_one(FLERR,"hdf5 libraries are not compiled, please change \"dumpFormat\" to \"json\". Check! \n");
   #endif
   
    if(dumpFormat_.compare("hdf5")==1 &&  dumpFormat_.compare("json")==1) error().throw_error_one(FLERR,"Not a valid \"dumpFormat\"! Valid dump formats are:\njson \nhdf5 \n");
   
  
   readBCfromJson();
   
    std::istream& infile_(ifnfile);    
       
    while (!infile_.eof())
    {
    
            char buf[MAX_CHARS_PER_LINE];
                 
	        infile_.getline(buf, MAX_CHARS_PER_LINE);
	
            int n = 0; // a for-loop index

            const char* token[MAX_TOKENS_PER_LINE] = {}; 
            
            token[0] = strtok(buf, " \t"); // first token
            if (token[0]) // zero if line is blank
            {
              for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
              {
                token[n] = strtok(0, " \t"); // subsequent tokens
                if (!token[n]) break; // no more tokens
              }

              //TODO error if n < 2

              // process the tokens

              std::string str(token[0], strlen(token[0]));
              c3poInterface().parse_command(str,n-1,&(token[1]));
            }

    }

    output().write_screen_one("\n ***   C3PO processed your inputs.   ***\n\n");
    timer().stamp(TIME_INPUT);
    
}


// *************************************************************
QJsonDocument Input::openJsonFile(const char* dirName, const char* fileName, const char* objectName, QJsonObject &json ) const
{
    char jsonfile[200];
    sprintf(jsonfile,"%s/%s.json",dirName,fileName);
	
    QFile    loadFile( jsonfile );
    if(!loadFile.open(QIODevice::ReadOnly)) 
        error().throw_error_one(FLERR,"can not open loadfile ",jsonfile);

    QByteArray    saveData = loadFile.readAll();
    QJsonDocument loadDoc  = QJsonDocument::fromJson(saveData);

    if(loadDoc.isNull())
	{ 
	    printf("Problematic file: %s \n", jsonfile);
    	error().throw_error_one(FLERR,"QJsonDocument is invalid. Check! \n");
	}

	json = loadDoc.object()[objectName].toObject();

    return loadDoc;
}
//*****************************************************************
void Input::readFieldNames()
{
 
  QJsonObject qOb=mainSettings()["FieldsToRegister"].toObject();
 
  QString qsV=qOb["Vectorfields"].toString();
   std::string big=qsV.toUtf8().constData();
   
  for (unsigned int it=0; it<big.size();it++)
  {
        
        if (big[it] != ' ' )
        {  
          
          std::string name;
         
          while (big[it] != ' ')
            {
               name.push_back(big[it]);
               it++;  
               if (it==big.size()) break;           
              
            }  
           
          VF_names.push_back(name); 
           
        }
  }
  

 
  QString qsS=qOb["Scalarfields"].toString();
   std::string bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
           
          SF_names.push_back(name); 
           
        }
  }
  
 if(!qOb["computeScalarGradient"].isNull()) 
 {
  qsS=qOb["computeScalarGradient"].toString();
  bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
           
          SF_grad.push_back(name); 
        }
  }
 }
 
 if(!qOb["computeVectorGradient"].isNull()) 
 {
  qsS=qOb["computeVectorGradient"].toString();
  bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
           
          VF_grad.push_back(name); 
           
        }
  }
 }
 
 if(!qOb["computeShearRate"].isNull()) 
 {
  qsS=qOb["computeShearRate"].toString();
  bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
           
          shearRates_.push_back(name); 
           
        }
  }
 }
 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void Input::readBCfromJson() const
{
   if(mainSettings()["filterBCs"].isNull()) error().throw_error_one(FLERR,"Boundary conditions for filtering operations are not specified in c3po.json! \n You should edit the 'filterBCs' field" );
   QJsonObject qOb=mainSettings()["filterBCs"].toObject();
   
   
   QString qsS=qOb["x"].toString();
   std::string bigS=qsS.toUtf8().constData();
   
   if(bigS.compare("periodic")==0) bound_[0]=-2;
   else if (bigS.compare("wall")==0) bound_[0]=-1;
   else  error().throw_error_one(FLERR,"Boundary conditions for filtering operations are not specified in c3po.json! \n Your 'x' boundary condition is wrong. Only 'wall' and 'periodic' boundary conditions are allowed");
   
   qsS=qOb["y"].toString();
   bigS=qsS.toUtf8().constData();
  
   if(bigS.compare("periodic")==0) bound_[1]=-2;
   else if (bigS.compare("wall")==0) bound_[1]=-1;
   else  error().throw_error_one(FLERR,"Boundary conditions for filtering operations are not specified in c3po.json! \n Your 'y' boundary condition is wrong. Only 'wall' and 'periodic' boundary conditions are allowed");

   qsS=qOb["y"].toString();
   bigS=qsS.toUtf8().constData();
  
   if(bigS.compare("periodic")==0) bound_[2]=-2;
   else if (bigS.compare("wall")==0) bound_[2]=-1;
   else  error().throw_error_one(FLERR,"Boundary conditions for filtering operations are not specified in c3po.json! \n Your 'z' boundary condition is wrong. Only 'wall' and 'periodic' boundary conditions are allowed");
   

}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double Input::getMaxDomainfromJson(int i) const
{
   if(meshMainSettings()["MaxDomain"].isNull()) error().throw_error_one(FLERR, "You have to specify the 'MaxDomain' field in mesh.json");
   QJsonArray array=meshMainSettings()["MaxDomain"].toArray();
    
   double Msize_=(array.at(i)).toDouble();
   
   return Msize_;   

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double Input::getMinDomainfromJson(int i) const
{
   if(meshMainSettings()["MinDomain"].isNull()) error().throw_error_one(FLERR, "You have to specify the 'MinDomain' field in mesh.json");
   QJsonArray array=meshMainSettings()["MinDomain"].toArray();
   
   double msize_=(array.at(i)).toDouble();
   
   return msize_;   

}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int Input::getDomainDecomposition() const
{
  if(meshMainSettings()["domain decomposition"].isNull())
  {
   output().write_screen_one("\n WARNING: because a domain decomposition direction was NOT specified in mesh.json, the k direction will be used\n");
   return 2;
  }

  int dir_ = meshMainSettings()["domain decomposition"].toInt();
  if(dir_>2 || dir_<0)
  {
  output().write_screen_one("\n WARNING: because the domain decomposition direction was ill specified in mesh.json, the k direction will be used\n");
   return 2;
  
  }
  return dir_;
} 

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double Input::cellSizefromJson(int i) const
{
   if(meshMainSettings()["Cell Size"].isNull()) error().throw_error_one(FLERR, "You have to specify the 'Cell Size field in mesh.json");
   QJsonArray array=meshMainSettings()["Cell Size"].toArray();
   
   double Csize_=(array.at(i)).toDouble();
   
   return Csize_;   

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */  
bool Input::writeInterface() const
{

 return mainObj_["interfaceWriteFields"].toBool();
   
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void Input::csv_columnsToRead(int * columns) const
{
 
 if(meshMainSettings()["CsvColumns"].isNull())
 {
  output().write_screen_one("\nNo information regarding the columns to read in mesh.json was given. The first four columns will be read.");
  for(int i=0;i<4;i++)
   columns[i]=i;
 } 
 else
 {
  QJsonObject qOb=meshMainSettings()["CsvColumns"].toObject();
  
  if(qOb["cellx"].isNull()) error().throw_error_one(FLERR, "You have to specify the column corresponding to the cell center x coordinates in mesh.json");
  if(qOb["celly"].isNull()) error().throw_error_one(FLERR, "You have to specify the column corresponding to the cell center y coordinates in mesh.json");
  if(qOb["cellz"].isNull()) error().throw_error_one(FLERR, "You have to specify the column corresponding to the cell center z coordinates in mesh.json");
  if(qOb["cellV"].isNull()) error().throw_error_one(FLERR, "You have to specify the column corresponding to the cell volumes in mesh.json");
  
  columns[0] = qOb["cellx"].toInt();
  columns[1] = qOb["celly"].toInt();
  columns[2] = qOb["cellz"].toInt();
  columns[3] = qOb["cellV"].toInt();
  
 }  

}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

std::string Input::getVFname(int x)
{
 std::string emp("empty");
 int size_= VF_names.size();                                   
 if(x<=size_)
  return VF_names[x];
 else
  return emp;
}

std::string Input::getSFname(int x)
{
 std::string emp("empty");
 int size_= SF_names.size();                                    
 if(x<=size_)
  return SF_names[x];
 else
  return emp;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void Input::readParticles(std::vector<double>* positions_, std::string probesName) const
{
 
 
 positions_->clear();
 
 QJsonDocument loadDoc;
 
 QJsonObject    parObj_;
 
 loadDoc = openJsonFile("c3po_control",probesName.c_str(),"probeData", parObj_ );
 
 if(parObj_["numberOfProbes"].isNull())
 {
  error().throw_error_one(FLERR,"ERROR: can not find \"numberOfprobes\" entry in the probe file. \n If you are not interest in probes or particles please set 'useProbes' to false in c3po.json 'mainSettings' \n"); 
 }
 
 int nPar_ = parObj_["numberOfProbes"].toInt();
 
 if(parObj_["positions"].isNull())
  error().throw_error_one(FLERR,"ERROR: entry \"positions\" not specified in probes file");
 
 
 QJsonObject    posObj_ = parObj_["positions"].toObject();
 
 for(int par=0;par<nPar_;par++) 
 {
  
  std::stringstream ss;
  ss << par;
  std::string str = ss.str();
  
  if(posObj_[str.c_str()].isNull())
   error().throw_error_one(FLERR,"ERROR: found samples are inconsistent with the \"numberOfProbes\" entry or invalid sample id");
   
  QJsonArray dataArray =  posObj_[str.c_str()].toArray();
  
  if(dataArray.count()!=3)
   error().throw_error_one(FLERR,"ERROR: wrong array length in probes file!");
   
  for(int i=0;i<3;i++)
   positions_->push_back((dataArray[i]).toDouble());


 }
 
 int regPar_= positions_->size()/3;
 
 if(regPar_!=nPar_)
   error().throw_error_one(FLERR,"ERROR: found probes are inconsistent with the \"numberOfProbes\" entry!");
 

}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool Input::registerCFDEMparticles() const
{
 if(mainObj_["registerCFDEMparticles"].isNull()) return false;
 else return mainObj_["registerCFDEMparticles"].toBool();
}

