#include <TString.h>  //ROOT Class: a string implementation
#include <ostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <TH1F.h>      //ROOT Class: 1-d histogram implementation

/*
 *  This class is an implementation of a mass spectrum. This class reads a mass spectrum file from the NIST webbook (http://webbook.nist.gov/chemistry/form-ser.html). 
 *
 *  The class calculates the atomic number of the molecule (the number of protons) and the number of H, C, O, N, Ar, or D atoms in the molecule. It also provides methods
 *  to visualize the mass spectrum, and to get the raw mass spectrum data.
 *
 *  It is meant to be used with the SpectrumCollection class, in SpectrumCollection.h.
 */




//element names and the number of protons they have
TString elementNames[6] = {"H", "C", "O", "N", "Ar", "D"};  
int    elementNumber[6] = { 1,   6,   8,   7,   18,   1}; 

//class for mass spectra
class Spectrum {
 public:
  
  //default constructor - creates an empty mass spectrum
  Spectrum(){
    filename="";
    for(int i=0; i<nMassSpecEntries; i++){  
      MZ.push_back(i+1);
      relInt.push_back(0);
    }
    for(int i=0; i<6; i++){
      contents[i] = 0;
    }
    Z=0;
    name="";
  }

  //constructor that loads a file
  Spectrum(TString file){
    filename=file;
    for(int i=0; i<nMassSpecEntries; i++){
      MZ.push_back(i+1);
      relInt.push_back(0);
    }
    for(int i=0; i<6; i++){
      contents[i] = 0;
    }

    readFile();
    ParseName();

  }
  
  //destructor
  ~Spectrum(){
  }

  //method that reads a jdx file
  void readFile(){
  
    ifstream InFile; 
    InFile.open(filename);

    TString EleName="";

    while(true){

      //this parses the name of the molecule. The name is of the form '##MOLFORM=C O2'
      TString word1;
      InFile>>word1;
      if(word1.Contains("##MOLFORM=")){
	word1.ReplaceAll("##MOLFORM=", "");  //remove the '##MOLFORM=' portion of the string
	EleName = EleName + word1;
	while(true){       
	  TString word2;
	  InFile>>word2;
	  if(word2.Contains("#")) break;  //end the loop when the next data tag is reached
	  EleName = EleName + "_" + word2;
	}	  
      }

      //this parses the mass spectrum of the molecule. The mass spectum data is of the form :
      /*
	##PEAK TABLE=(XY..XY)
	1,210 2,9999
	##END=
      */
      // The first number is m/z, the second number is the relative intensity.

      if(word1=="TABLE=(XY..XY)"){
	while(true){ 
	  TString word3;
	  InFile>>word3;
	  if(word3=="##END=") break; //end the loop when the next data tag is found
	  
	  int from = 0;
	  TString tok;
	  int index=-1;
	  
	  bool first=true;
	  
	  while(word3.Tokenize(tok, from, "[,]")){ //tokenize the string based on ','
	    if(first){
	      index = tok.Atoi()-1; //the first token is the index (m/z)
	      first=false;
	    }
	    else if(!false) {
	      if(index!=-1) relInt[index] = (double)tok.Atoi()/9999; //the second token is the relative intensity. Divide by 9999 to normalize to 1.
	    }
	  }
	}
      }
      if(InFile.eof()) break;
    }

    name = EleName;  //set the name of the molecule to the name read from the jdx file.
    return;

  }

  //reads the name, gets the number of each element and calculates Z of the atom
  void ParseName(){
    int from = 0;
    TString tok;
    while(name.Tokenize(tok, from, "[_]")){  //tokenize the string based on '_'. Splits C_O2 into C and O2
      int n=1; //the multiplicity of this atom in the molecule
      //if the token contains a number, the multiplicity is that number.
      if(tok.Contains("1")||tok.Contains("2")||tok.Contains("3")||tok.Contains("4")||tok.Contains("5")||tok.Contains("6")||tok.Contains("7")||tok.Contains("8")||tok.Contains("9")){
	TString sub = tok(tok.Length()-1,1);
	n = sub.Atoi();
      }
      //fill the 'contents[]' array with the multiplicity of that atom
      for(int i=0; i<6; i++){
	if(tok.Contains(elementNames[i])){
	  contents[i] = n;
	  break;
	}
      }
    }
    //calculate Z of the atom.
    Z=0;
    for(int i=0; i<6; i++){
      Z=Z+contents[i]*elementNumber[i];
    }

  }
  
  //print the non-zero entries in the spectrum
  void printNonZero(){
    for(int i=0; i<50; i++){
      if(y[i]!=0) std::cout<<i+1<<"\t"<<relInt[i]<<std::endl;
    }

  }

  //get a ROOT histogram of the mass spectrum
  TH1F* getHistogram(){
    TH1F *hist = new TH1F(name, "", 150, 0, 50);

    for(int i=0; i<nMassSpecEntries; i++){
      hist->AddBinContent(3*(i+1)+1, relInt[i]);
    }
        
    hist->SetFillColor(1002);
    hist->SetLineColor(1002);
    hist->GetXaxis()->SetTitle("m/z");
    hist->GetYaxis()->SetTitle("Relative Abundance (AU)");
    
    return hist;
  }


  //get methods
  TString getName(){ return name;}
  TString getFile(){ return filename;}
  vector<double> getRelativeIntensity(){ return relInt;}
  vector<double> getMZ(){ return MZ;}
  double getRelativeIntensity(int zm){ 
    if(zm<1||zm>50) return -1; 
    return relInt[zm-1];}
  double* getContents(){ return contents;}
  int getZ(){return Z;}   //the number of protons in this gas
  
  void setFile(TString file) {
    filename = file;
    readFile();
    ParseName();
  }

 private:
  static const int nMassSpecEntries=51;
  TString filename;
  double contents[6];
  vector<double> MZ;
  vector<double> relInt;
  TString name;
  int Z;
};

