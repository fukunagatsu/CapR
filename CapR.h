/*
 * CapR.h
 *
 *  Created on: 2016/6/29
 *      Author: Tsukasa Fukunaga
 */

#ifndef CAPR_H
#define CAPR_H

#include <string>
#include <stdlib.h>
#include <vector>
#include "energy_par.h"
#include "intloops.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

class CapR{
 public:
  CapR(string input_file, string output_file, int maximal_span){
    _input_file_name = input_file;
    _output_file_name = output_file;
    _maximal_span = maximal_span;
     _seq_length = 0;
    set_energy_parameters();
    ofstream of(_output_file_name.c_str(), ios::out);
    if (!of){
      cerr  << "Error: Cannot make " << _output_file_name << "." << endl;
      exit(1);
    }
    of.close();
  }
  void Run();
 private:
  string _input_file_name;
  string _output_file_name;
  int _maximal_span;

  double hairpin[31];
  double mismatchH[7][5][5];
  double mismatchI[7][5][5];
  double stack[7][7];
  double bulge[31];
  double TermAU;
  double int11[8][8][5][5];
  double int21[8][8][5][5][5];
  double int22[8][8][5][5][5][5];
  double internal[31];
  double MLclosing;
  double MLintern;
  double MLbase;
  double dangle5[8][5];
  double dangle3[8][5];
  double ninio[MAXLOOP+1];
  
  vector<int> _int_sequence;
  int _seq_length;

  vector<double> _Alpha_outer;
  vector<vector<double> > _Alpha_stem;
  vector<vector<double> > _Alpha_stemend;
  vector<vector<double> > _Alpha_multi;
  vector<vector<double> > _Alpha_multibif;
  vector<vector<double> > _Alpha_multi1;
  vector<vector<double> > _Alpha_multi2;
  
  vector<double> _Beta_outer;
  vector<vector<double> > _Beta_stem;
  vector<vector<double> > _Beta_stemend;
  vector<vector<double> > _Beta_multi;
  vector<vector<double> > _Beta_multibif;
  vector<vector<double> > _Beta_multi1;
  vector<vector<double> > _Beta_multi2;

  void set_energy_parameters(){
    MLclosing = -ML_closing37*10/kT;
    MLintern = -ML_intern37*10./kT;
    MLbase = -ML_BASE37*10./kT;
    TermAU= -TerminalAU*10/kT;
    
    for (int i=0; i<=30; i++) {
      hairpin[i] = -hairpin37[i]*10./kT;
      bulge[i] = - bulge37[i]*10./kT;
      internal[i] = -internal_loop37[i]*10./kT;
    }
    
    for (int i=0; i< 7; i++){
      for (int j=0; j<5; j++){
	for (int k=0; k<5; k++) {
	  mismatchI[i][j][k] = -mismatchI37[i][j][k]*10.0/kT;
	  mismatchH[i][j][k] = -mismatchH37[i][j][k]*10.0/kT;
	}
      }

      for (int j=0; j<7; j++) {
	stack[i][j] = -stack37[i][j]*10./kT;
      }

      for (int j=0; j<=4; j++) {
	dangle5[i][j] = -dangle5_37[i][j]*10./kT;
	dangle3[i][j] = -dangle3_37[i][j]*10./kT;
	if (i>2){
	  dangle3[i][j] += TermAU;
	}
      }
    }
    
    for (int i=0; i<=7; i++){
      for (int j=0; j<=7; j++){
	for (int k=0; k<5; k++){
	  for (int l=0; l<5; l++){
	    int11[i][j][k][l] = -int11_37[i][j][k][l]*10./kT;
	    for (int m=0; m<5; m++){
	      int21[i][j][k][l][m] = -int21_37[i][j][k][l][m]*10./kT;
	      for (int n=0; n<5; n++){
		int22[i][j][k][l][m][n] = -int22_37[i][j][k][l][m][n]*10./kT;
	      }
	    }
	  }
	}
      }
    }
    
    for (int i=0; i<=MAXLOOP; i++){
      ninio[i]= -min(MAX_NINIO, i*F_ninio37)*10/kT ;
    } 
  }
  void CalcMain(string &sequence, string &name);
  void Initiallize(string &sequence);
  void CalcInsideVariable();
  void CalcOutsideVariable();
  void CalcStructuralProfile(string &name);
  double CalcExteriorProbability(int x);
  void CalcHairpinProbability(vector<double> &hairpin_probability);
  double CalcMultiProbability(int x);
  void CalcBulgeAndInternalProbability(vector<double> &bulge_probability, vector<double> &internal_probability);
  void CalcLogSumBulgeAndInternalProbability(vector<double> &bulge_probability, vector<double> &internal_probability);
  void Clear();

  double CalcDangleEnergy(int type,int a, int b);
  double logsumexp(double x,double y);
  double LoopEnergy(int type, int type2,int i,int j,int p,int q);
  double HairpinEnergy(int type, int i, int j);
};

#endif
