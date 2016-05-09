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
\*-----------------------------------------------------------------------------
    Description
    tool for interpretation and execution of user-defined formulas
\*-----------------------------------------------------------------------------*/

#include "formula.h"
#include "stdlib.h"
#include <cstdio>
#include <cmath>
#include "mpi.h"


using namespace C3PO_NS;

Formula::Formula(const char* formula)
:
raw_formula_(formula)
{
 limitNumeratorActive_       = false;
 limitDenominatorActive_     = false;    absoluteNumerator_   = false;
 limitNumerator_[0]          = -1e99;    limitNumerator_[1]   = 1e99;
 limitNumeratorSymmetric_    = true;     absoluteDenominator_ = false;
 limitDenominator_[0]        = -1e99;    limitDenominator_[1] = 1e99;
 limitDenominatorSymmetric_  = true;
 normalize_                  = false;

 ParseFormula();
 CheckFormula();
}
/*-----------------------------------------------------------------------------*/
Formula::Formula(const char* formula, const char* normalization)
:
raw_formula_(formula),
normalization_(normalization)
{
 limitNumeratorActive_       = false;
 limitDenominatorActive_     = false;    absoluteNumerator_   = false;
 limitNumerator_[0]          = -1e99;    limitNumerator_[1]   = 1e99;
 limitNumeratorSymmetric_    = true;     absoluteDenominator_ = false;
 limitDenominator_[0]        = -1e99;    limitDenominator_[1] = 1e99;
 limitDenominatorSymmetric_  = true;
 normalize_                  = true;

 ParseFormula();
 CheckFormula();
}
/*-----------------------------------------------------------------------------*/
Formula::~Formula()
{
}
/*-----------------------------------------------------------------------------*/
void Formula::throw_error(std::string function ,std::string msg_) const
{
 std::cout << "\nFormula::"<< function << " -> ERROR: " << msg_;
 MPI_Barrier(MPI_COMM_WORLD);
 //MPI_Finalize();
 exit(1);
}
/*-----------------------------------------------------------------------------*/
void Formula::ParseFormula() const
{
  bool foundNumerator_   = false;
  bool foundDenominator_ = false;
  for (unsigned int it=0; it<raw_formula_.size();it++)
  {
        if (raw_formula_[it] == '(' && !foundDenominator_ )
            foundNumerator_ = true;

        if (raw_formula_[it] == '%')
        {
            if(foundDenominator_ == true) throw_error("ParseFormula","Only one numerator/denominator separator is allowed!");
            foundNumerator_   = false;
            foundDenominator_ = true;
        }

        if ( raw_formula_[it] != ' ' )
        {  
          formula_.push_back(raw_formula_[it]);

          if(     foundNumerator_ 
               && !(raw_formula_[it]=='(') 
               && !(raw_formula_[it]==')')  
            )
            numerator_.push_back(raw_formula_[it]);

          if(      foundDenominator_ 
                && !(raw_formula_[it]=='(') 
                && !(raw_formula_[it]==')')
                && !(raw_formula_[it]=='%')  
            )
            denominator_.push_back(raw_formula_[it]);

        }
  }

  //If user did not specify brackets, fill numerator
  if(numerator_.size()==0)
      numerator_.assign(formula_.begin(),formula_.end());

   std::cout << "formula: ";
   std::cout << (getFormula()).c_str() << " \n";

}
/*-----------------------------------------------------------------------------*/
void Formula::interpretFormula(int numVecs, int numScalars) const
{
    //Replace words starting with 'vec' and 'scalar' with ids
    for(unsigned int i=(numerator_.size()-1);i>0;i--) //start at the back
    {
        if(   ( numerator_[i]=='v' && numerator_[i+2]!='c')  
            ||( numerator_[i]=='e' && numerator_[i+1]!='c') 
            ||( numerator_[i]=='s' && numerator_[i+5]!='r') 
            ||( numerator_[i]=='a' && numerator_[i+3]!='r' && numerator_[i+1]!='r')  
            ||( numerator_[i]=='l' && numerator_[i+2]!='r') 
            ||( numerator_[i]=='c' && numerator_[i+4]!='r' && numerator_[i-1]!='e') 
          )
        throw_error("interpretFormula","fields are not correctly defined in your formula numerator! Use 'vec' or 'scalar' to identify them.");
				
		if( numerator_[i]=='c' )
        {
          if( numerator_[i-1]=='e' && numerator_[i-2]=='v')
          {
           if(i==numerator_.size()-1) throw_error("interpretFormula","std::vector field number missing!");          
           numerator_[i+1]=numerator_[i+1]-1; //for std::vectors just decrease the index           
           if(numerator_[i+1]-'0'>9) throw_error("interpretFormula","std::vector field number missing!");                     
           numerator_.erase(numerator_.begin()+i-2,numerator_.begin()+i+1);
		   i = i - 2;
		   if (i == 0) i = 1;
         }
	    else
         {
           throw_error("interpretFormula","std::vector field not correctly defined in your formula numerator!");      
         }
        } 
        
       
        if( numerator_[i]=='r' )
        {
          if( numerator_[i-1]=='a' && numerator_[i-2]=='l'  && numerator_[i-3]=='a' && numerator_[i-4]=='c' && numerator_[i-5]=='s' )
          {
            if(i==numerator_.size()-1) throw_error("interpretFormula","scalar field number missing!");
            numerator_[i+1]=numerator_[i+1]-1+numVecs;
            if(numerator_[i+1]-'0'>9) throw_error("interpretFormula","scalar field number missing!"); 
            numerator_.erase(numerator_.begin()+i-5,numerator_.begin()+i+1);
            i = i - 5;
			if (i == 0) i = 1;
          }
         else
          {
           throw_error("interpretFormula","scalar field not correctly defined in your formula numerator!");      
          }
        }   
    }
    std::cout << "interpreted numerator: " << (getNumerator()).c_str() << " \n";
    

    //Replace words starting with 'vec' and 'scalar' with ids
    if(denominator_.size()>0)
    for(unsigned int i=(denominator_.size()-1);i>0;i--) //start at the back
    {
       if(   ( denominator_[i]=='v' && denominator_[i+2]!='c')  
            ||( denominator_[i]=='e' && denominator_[i+1]!='c') 
            ||( denominator_[i]=='s' && denominator_[i+5]!='r') 
            ||( denominator_[i]=='a' && denominator_[i+3]!='r' && denominator_[i+1]!='r')  
            ||( denominator_[i]=='l' && denominator_[i+2]!='r') 
            ||( denominator_[i]=='c' && denominator_[i+4]!='r' && denominator_[i-1]!='e') 
         )
        throw_error("interpretFormula","fields are not correctly defined in your formula denominator! Use 'vec' or 'scalar' to identify them.");

       if( denominator_[i]=='c' )
       {
        if( denominator_[i-1]=='e' && denominator_[i-2]=='v' )
         {
            if(i==denominator_.size()-1) throw_error("interpretFormula","std::vector field number missing!");           
            denominator_[i+1]=denominator_[i+1]-1; //for std::vectors just decrease the index       
            if(denominator_[i+1]-'0'>9) throw_error("interpretFormula","std::vector field number missing!");           
            denominator_.erase(denominator_.begin()+i-2,denominator_.begin()+i+1);
			i = i - 2;
			if (i == 0) i = 1;
         }
	    else
         {
           throw_error("interpretFormula","std::vector field not correctly defined in your formula denominator!");      
         }
       }
       
       
       if( denominator_[i]=='r' )
       {
         if( denominator_[i-1]=='a' && denominator_[i-2]=='l'  && denominator_[i-3]=='a' && denominator_[i-4]=='c' && denominator_[i-5]=='s' )
         {
            if(i==denominator_.size()-1) throw_error("interpretFormula","scalar field number missing!");
            denominator_[i+1]=denominator_[i+1]-1+numVecs; //for 
            if(denominator_[i+1]-'0'>9) throw_error("interpretFormula","scalar field number missing!");
            denominator_.erase(denominator_.begin()+i-5,denominator_.begin()+i+1);
			i = i - 5;
			if (i == 0) i = 1;
         }
         else
         {
           throw_error("interpretFormula","scalar field not correctly defined in your formula denominator!");      
         }
       }
    }

    std::cout << "interpreted denominator: " << (getDenominator()).c_str() << " \n";
    
    if(normalize_)
     CheckNormalization( numVecs,  numScalars);
}

/*-----------------------------------------------------------------------------*/
void Formula::CheckNormalization(int numVecs, int numScalars) const
{
 
 bool typeFound_=false;


 for(unsigned int id=0; id<normalization_.size();id++)
 {
  
  std::cout << "\n norm:" << normalization_[id];
  
  if(normalization_[id]==' ') continue;
  if(normalization_[id]-'0' < 9) continue;
  
  if(!typeFound_)
  {
   //Add one for every Normalization
   if(normalization_[id]=='D') 
   {
    typeFound_=true;
    normalizationType_= normalization_[id];
   }
   else throw_error("Formula::CheckNormalization()" ,"Incorrect normalization (0)!");
  }
  else 
  {
   if(normalization_[id]=='m')
   {
    
    if(normalization_[id+1]-'0' < 9)
    {
      normalizationId_.push_back(normalization_[id+1]-'0'-1);
      normalizationClass_.push_back('m');
      
      //continue;
    }
    else  throw_error("Formula::CheckNormalization()" ,"Incorrect normalization (1)!");
      
   }
   else if (normalization_[id]=='v')
   {
   
    if(normalization_[id+1]-'0' < 9  )
    {
      normalizationClass_.push_back('v');
      normalizationId_.push_back(normalization_[id+1]-'0'-1);
    }
    else  throw_error("Formula::CheckNormalization()" ,"Incorrect normalization (2)!");
      
   }
   else if (normalization_[id]=='s')
   {
    
    if(normalization_[id+1]-'0' < 9 )
    {
      normalizationId_.push_back(normalization_[id+1]-'0'+numVecs-2);
      normalizationClass_.push_back('s');
    }
    else  throw_error("Formula::CheckNormalization()" ,"Incorrect normalization (3)!");
      
   }
   else throw_error("Formula::CheckNormalization()" ,"Incorrect normalization! (4)");
  
  }
 
 }
 
 if(!typeFound_)
   throw_error("Formula::CheckNormalization()" ,"Incorrect normalization! (6)");

//  std::cout << "\nprinting normalizations...";
//  for(unsigned int i=0;i<normalizationClass_.size();i++)
//   std::cout << "\nnormClass: "<<normalizationClass_[i]<<" normId: "<<normalizationId_[i];

 

}
/*-----------------------------------------------------------------------------*/
void Formula::CheckFormula() const
{
 bool operator_before=true;

 for (unsigned int it=0; it<formula_.size();it++)
  {
   if(formula_[it]=='(' || formula_[it]==')' 
      || formula_[it]=='v'|| formula_[it]=='e'|| formula_[it]=='c'
      || formula_[it]=='s'|| formula_[it]=='a'|| formula_[it]=='l'|| formula_[it]=='r'
     )
       continue;

   if(formula_[it]=='%')
   {
       operator_before=true;
       continue;
   }

   if(operator_before )    
   {
    operator_before=false;
   }
   else
   {
    if(formula_[it]!='*' && formula_[it]!='/' && formula_[it]!='+' && formula_[it]!='-'  )
    {
    
     char buf[512];
     sprintf(buf,"\nFormula::CheckFormula -> ERROR: '%c' is not a valid entry for CPPPO formulas OR the formula is not correct (i.e., does not contain 'scalar' or 'vec' quantities!\n Valid entries for position %i are:\n '+'\n '-'\n '*'\n '/'\n",formula_[it],it);
     std::string msg_(buf);

     throw_error("CheckFormula",msg_); 
    
    }
   
    operator_before=true;
    
   }   
  }
 

 
}

/*-----------------------------------------------------------------------------*/
std::string Formula::getFormula() const
{
 std::string out_;
 for(unsigned int i=0;i<formula_.size();i++)
  out_.push_back(formula_[i]);
  
 return out_;
}


/*-----------------------------------------------------------------------------*/
std::string Formula::getNumerator() const
{
 std::string out_;
 for(unsigned int i=0;i<numerator_.size();i++)
  out_.push_back(numerator_[i]);
  
 return out_;
}


/*-----------------------------------------------------------------------------*/
std::string Formula::getDenominator() const
{
 std::string out_;
 for(unsigned int i=0;i<denominator_.size();i++)
  out_.push_back(denominator_[i]);
  
 return out_;
}

/*-----------------------------------------------------------------------------*/
//normalizes the sample
void Formula::normalize( std::vector< std::vector<double>  > * _sample, 
                         std::vector< std::vector<double>* >*  _markers,
                         std::vector< double > * outputField_ 
                       )
{
 

 //Determine size of sample
 int size_=outputField_->size();
 
 
  //Loop sample vectors
  for(int i=0;i<size_;i++)
  {
   bool haveValidOutput=true; 
   double valueNumerator_=1;
   //This is the drag force calculation
   if(normalizationType_=='D') //'interphaseEC' is the identifier for interphase Exchange Coefficient. TODO; implement more identifiers if necessary
   { 
     //will use marker values to access vol fraction and slip velocity
     if(normalizationClass_.size()!=3) throw_error("Formula::normalize()" ,"Wrong number of arguments!");    
     double valueParticleVolFraction = getValueForNormalization(0,tmpIndex_[i],_sample,_markers);
     double valueGasVelocity        = getValueForNormalization(1,tmpIndex_[i],_sample,_markers);
	 double valueSolidsVelocity        = getValueForNormalization(2,tmpIndex_[i],_sample,_markers);
     
    
   //  std::cout << "\n volFrac : " << valueParticleVolFraction << " slip : "<< valueSlipVelocity ;
     
    valueNumerator_ *= basicMultiphaseQty_->interphaseEC(valueParticleVolFraction, valueGasVelocity, valueSolidsVelocity);
   } 
    
 
   //only insert in case we have a valid output
   if(haveValidOutput && std::abs(valueNumerator_)>1e-15)
   {
       //TODO: normalize ONLY desired values of sample array
       (*outputField_)[i] /= valueNumerator_;
   }
  }
  
  tmpIndex_.clear();
}   


/*-----------------------------------------------------------------------------*/
void Formula::evaluate( std::vector< std::vector<double>  > * inputFields_, 
                        std::vector< std::vector<double>* >*  _markers,
                        std::vector< double > * outputField_,
                        std::vector< std::vector<double> >*    newMarkers_
                        
                      )
{
 int size_=(*inputFields_)[0].size();
 
 //start calculation
 for(int i=0;i<size_;i++)
 {
  bool haveValidOutput=true;

  //Compute numerator
  double valueNumerator_=(*inputFields_)[(numerator_[0]-'0')][i];
  for(int n=1;n<int(numerator_.size());n++)
  {
   if(numerator_[n]=='+') valueNumerator_+= (*inputFields_)[(numerator_[n+1]-'0')][i];
   if(numerator_[n]=='-') valueNumerator_-= (*inputFields_)[(numerator_[n+1]-'0')][i];
   if(numerator_[n]=='*') valueNumerator_*= (*inputFields_)[(numerator_[n+1]-'0')][i];
   if(numerator_[n]=='/') valueNumerator_/= (*inputFields_)[(numerator_[n+1]-'0')][i] + 1e-15;
   n++; //Advance another value because of operator
  }
  
  if(absoluteNumerator_)
   valueNumerator_=std::fabs(valueNumerator_);
  
  if(limitNumeratorActive_)
    haveValidOutput =  haveValidOutput 
                     && testLimits(valueNumerator_, limitNumerator_, limitNumeratorSymmetric_);
  
    

  //Compute denominator
  double valueDenominator_ = 1.0;

  if(denominator_.size()>0)
  {
    valueDenominator_=(*inputFields_)[(denominator_[0]-'0')][i];
    for(int n=1;n<int(denominator_.size());n++)
    {
     if(denominator_[n]=='+') valueDenominator_+= (*inputFields_)[(denominator_[n+1]-'0')][i];
     if(denominator_[n]=='-') valueDenominator_-= (*inputFields_)[(denominator_[n+1]-'0')][i];
     if(denominator_[n]=='*') valueDenominator_*= (*inputFields_)[(denominator_[n+1]-'0')][i];
     if(denominator_[n]=='/') valueDenominator_/= (*inputFields_)[(denominator_[n+1]-'0')][i] + 1e-15;
     n++; //Advance another value because of operator
    }
  }
 
 if(absoluteDenominator_)
   valueDenominator_=std::fabs(valueDenominator_);
      
 if(limitDenominatorActive_)
    haveValidOutput =  haveValidOutput
                    && testLimits(valueDenominator_, limitDenominator_, limitDenominatorSymmetric_);

  //only insert in case we have a valid output
  if(haveValidOutput)
  {
   outputField_->push_back(valueNumerator_ / valueDenominator_); 
   for(unsigned int mark=0;mark< newMarkers_->size();mark++)
       (*newMarkers_)[mark].push_back((*_markers)[0][mark][i]);
       
   if(normalize_) tmpIndex_.push_back(i);
  }    
 
 }
 
 if(normalize_)
    normalize( inputFields_, 
               _markers,
               outputField_ 
             );
}

/*-----------------------------------------------------------------------------*/
bool Formula::testLimits(double& value, double* limits, bool& isSymmetric) const
{
    //check the limits: limits[0] is lower limit, limits[1] is upper limit
    if(isSymmetric) //this is typically the more useful one!
    {
        if( std::fabs(value)>limits[0] && std::fabs(value)<limits[1] )
            return true;
        else
            return false;
    }
    else
    {
        if( value>limits[0] && value<limits[1] )
            return true;
        else
            return false;
    }
}
/*-----------------------------------------------------------------------------*/
double Formula::getValueForNormalization(int id,
                                         int index,
                                         std::vector< std::vector<double> > *  _sample,
                                         std::vector< std::vector<double>* >*  _markers 
                                        ) const
{
 double value_=1;
 
  if(normalizationClass_[id] == 'm')
  {
     value_ = (*_markers)[0][normalizationId_[id]][index];
  }
  else
  {
     value_ = (*_sample)[normalizationId_[id]][index];
  }      
  return value_;
}
/*-----------------------------------------------------------------------------*/
void Formula::setMaxNumberOfFields(int maxMark_, int maxScal_, int maxVec_) const
{
    NofVecSample_  =    maxVec_;
    NofScalSample_  =  maxScal_;
    NofMarkSample_  =  maxMark_;
    
    if(normalize_)
    {
   
     
     for(unsigned int i=0; i<normalizationClass_.size(); i++)
     {
     
    
      if(normalizationClass_[i] == 'm')
      {
       if(normalizationId_[i] > NofMarkSample_ -1  )
        throw_error("Formula::setMaxNumberOfFields()" ,"Markers defined in normalization do not match with markers defined in the sampling operation!!"); 
        
      }
      else if(normalizationClass_[i] == 'v')
      {
       std::cout << "\nnormalizationId_[i]: " << normalizationId_[i] << " NofVecSample_ -1: " << NofVecSample_ -1;
       if(normalizationId_[i] > NofVecSample_ -1  )
        throw_error("Formula::setMaxNumberOfFields()" ,"Vectors defined in normalization do not match with vectors defined in the sampling operation!!"); 
        
      }
      else if(normalizationClass_[i] == 's')
      {
       if(normalizationId_[i] > NofScalSample_ + NofVecSample_ -2  )
        throw_error("Formula::setMaxNumberOfFields()" ,"Scalars defined in normalization do not match with scalars defined in the sampling operation!!"); 
        
      }
 
     
     }
     
     std::cout << "\nFields correctly registered in Formula class.";
         
    }
}
