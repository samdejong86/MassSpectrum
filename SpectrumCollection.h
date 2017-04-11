#include <TString.h>  //ROOT Class: a string implementation
#include <ostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <TMatrixD.h> //ROOT Class: Matrix implementation



/* A collection of mass specta. Includes methods for fitting the collection to data. 
 *
 *  Fitting is done using least squares:
 *
 *   soln = (Xt*X)^-1*Xt*Y
 *   where Y is the input data
 *
 * This uses the Spectrum class defined in Spectrum.h
 */
class SpectrumCollection{
 public:

  //default constructor
  SpectrumCollection(){
    XtEval=false;
    XtXEval=false;
    XtXInvEval=false;
  }

  //add a spectrum to the list using the filename, create the X matrix
  void AddSpectrum(TString fname){
    Spectrum a(fname);

    cout<<"Adding "<<a.getName()<<" to collection\n";
    gases.push_back(a);

    //X is an n x 50 matrix, where n is the number of spectra currently in the collection. 
    int currentNcol = X.GetNcols();
    X.ResizeTo(50, currentNcol+1); //add a new column to the matrix for the spectum being added
    for(int i=0; i<50; i++){
      X[i][currentNcol] = a.getY()[i];
    }
    
    //will need to reevaluate Xt, XtX, and XtX^-1
    XtEval=false;
    XtXEval=false;
    XtXInvEval=false;

  }

  //transpose the X matrix to get Xt
  void Transpose(){
    Xt.ResizeTo(X.GetNcols(), 50);

    for(int i=0; i<50; i++){
      for(int j=0; j<X.GetNcols(); j++){
	Xt[j][i] = X[i][j];
      }
    }
    XtEval=true;  //note that the matrix has been transposed, to avoid re-transposing.
  }

  //multiply Xt by X
  void SelfMultiply(){
    if(!XtEval) this->Transpose();  //transpose if it hasn't been done yet
    
    TMatrixD prod= this->MatrixMultiply(Xt, X);
    
    XtX.ResizeTo(prod);
    XtX=prod;

    XtXEval=true; //note that XtX exists, to avoid recalculating.
  }
  
  void InvertXtX(){
    if(!XtXEval) this->SelfMultiply();  //calculate XtX if it hasn't been calculated yet
    XtXInv.ResizeTo(XtX);
    XtXInv = XtX.Invert(); 
    XtXInvEval=true; //note that XtX^-1 exists, to avoid recalculating
  }

  //Get the solution to (Xt*X)^-1*Xt*Y
  TMatrixD Evaluate(TMatrixD input, bool withError=false){
    //make sure Xt, XtX, and XtX^-1 have been evaluated
    if(!XtEval) this->Transpose();
    if(!XtXEval) this->SelfMultiply();
    if(!XtXInvEval) this->InvertXtX();

    TMatrixD XtY = this->MatrixMultiply(Xt, input);

    TMatrixD soln = this->MatrixMultiply(XtXInv, XtY);   //calculate the least squares solution

    //estimate the uncertainty on the fit parameters
    if(withError){

      TMatrixD appliedSoln(input.GetNrows(), input.GetNcols());
      
      for(int i=0; i<input.GetNrows(); i++){
	for(int j=0; j<(int)gases.size(); j++){
	  appliedSoln[i][0] +=  soln[j][0]*gases[j].getY()[i];
	}
      }
      
      
      TMatrixD difference = input-appliedSoln;
      
      double sum;
      double sumSq;
      int n=difference.GetNrows();
      for(int i=0; i<n; i++){
	sum +=difference[i][0];
	sumSq+=difference[i][0]*difference[i][0];      
      }
      
      double sigmaSq = (sumSq/n)-(sum/n)*(sum/n);
      TMatrixD Var = XtXInv*sigmaSq;
      
      soln.ResizeTo(soln.GetNrows(),2);
      for(int i=0; i<soln.GetNrows(); i++){
	soln[i][1] = sqrt(Var[i][i]);
      }
    }
    
    return soln;

  }
  
  
  //get methods
  vector<Spectrum> getGases(){ return gases;}
  vector<TString> getGasNames(){
    vector<TString> names;
    for(int i=0; i<(int)gases.size(); i++){
      names.push_back(gases[i].getName());
    }
    return names;
  }

  //get the various matrices used
  TMatrixD getX() { return X;}
  TMatrixD getXt() {
    if(!XtEval) this->Transpose();
    return Xt;}
  TMatrixD getXtX() {
    if(!XtXEval) this->SelfMultiply();
    return XtX;
  }

  void PrintXtX(){
    if(!XtXEval) this->SelfMultiply();
    for(int i=0; i<XtX.GetNrows(); i++){
      for(int j=0; j<XtX.GetNcols(); j++){
	cout<<XtX[i][j]<<"  ";
      }
      cout<<endl;
    }
  }
  

 private:

  //the TMatrix multiplier doesn't seem to work...
  TMatrixD MatrixMultiply(TMatrixD a, TMatrixD b){
    
    int r1=a.GetNrows();
    int r2=b.GetNrows();
    int c1=a.GetNcols();
    int c2=b.GetNcols();

    TMatrixD mult(r1,c2);
    
    int i,j,k;
    
    /* Initializing elements of matrix mult to 0.*/
    for(i=0; i<r1; ++i){
      for(j=0; j<c2; ++j){
       mult[i][j]=0;
      }
    }
    /* Multiplying matrix a and b and storing in array mult. */
    for(i=0; i<r1; ++i){
      for(j=0; j<c2; ++j){
	for(k=0; k<c1; ++k){
	  mult[i][j]+=a[i][k]*b[k][j];
	}
      }
    }
    
    return mult;
  }
  



  const int nRows=50;
  TMatrixD X;
  TMatrixD Xt;
  TMatrixD XtX;
  TMatrixD XtXInv;
  vector<Spectrum> gases;
  bool XtEval;
  bool XtXEval;
  bool XtXInvEval;
};

