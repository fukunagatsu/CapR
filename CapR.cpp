/*
 * CapR.cpp
 *
 *  Created on: 2016/6/29
 *      Author: Tsukasa Fukunaga
 */

#include "CapR.h"
#include "fastafile_reader.h"

void CapR::Run(){
  vector<string> sequences; sequences.reserve(10);
  vector<string> names; names.reserve(10);

  FastafileReader fastafile_reader;
  fastafile_reader.ReadFastafile(_input_file_name, sequences, names);
  
  for (int i = 0; i < sequences.size() ; i++){
    CalcMain(sequences[i],names[i]);
  }
}

void CapR::CalcMain(string &sequence, string &name){
  Initiallize(sequence);
  CalcInsideVariable();
  CalcOutsideVariable();
  CalcStructuralProfile(name);
  Clear();
}

void CapR::Initiallize(string &sequence){
  _seq_length = sequence.length();
  _int_sequence.resize(_seq_length+1);
  for(int i = 0;i < _seq_length;i++){
    if(sequence[i] == 'A' || sequence[i] == 'a'){
      _int_sequence[i+1] = 1;
    }else if(sequence[i] == 'C' || sequence[i] == 'c'){
      _int_sequence[i+1] = 2;
    }else if(sequence[i] == 'G' || sequence[i] == 'g'){
      _int_sequence[i+1] = 3;
    }else if(sequence[i] == 'T' || sequence[i] == 't' || sequence[i] == 'U' || sequence[i] == 'u'){
      _int_sequence[i+1] = 4;
    }else{
      _int_sequence[i+1] = 0;
    }
  }
  
  _Alpha_outer.resize(_seq_length+1, 0);
  _Alpha_stem.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_stemend.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multi.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multibif.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multi1.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multi2.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  

  _Beta_outer.resize(_seq_length+1, 0);
  _Beta_stem.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_stemend.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multi.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multibif.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multi1.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multi2.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
}

void CapR::CalcInsideVariable(){
  for (int j =TURN+1; j <= _seq_length; j++){
    for (int i=j-TURN; i >= max(0,j-_maximal_span-1); i--){
      //Alpha_stem
      int type = BP_pair[_int_sequence[i+1]][_int_sequence[j]];
      int type2 = BP_pair[_int_sequence[i+2]][_int_sequence[j-1]];
      
      double temp = 0; bool flag = 0;
      if (type != 0) {
	type2 = rtype[type2];
	if(_Alpha_stem[i+1][j-i-2] != -INF){
	  //Stem¨Stem
	  if(type2 != 0){
	    temp = _Alpha_stem[i+1][j-i-2]+ LoopEnergy(type, type2,i+1,j,i+2,j-1);
	  }
	  flag = 1;
	}
    
	if(_Alpha_stemend[i+1][j-i-2] != -INF){
	  //Stem¨StemEnd
	  temp = flag == 1 ? logsumexp(temp,_Alpha_stemend[i+1][j-i-2]) : _Alpha_stemend[i+1][j-i-2];
	  flag = 1;
	}

	_Alpha_stem[i][j-i] = flag == 0 ? -INF : temp;
      }else{
	_Alpha_stem[i][j-i] = -INF;
      }
      
      //Alpha_multiBif
      temp = 0; flag = 0;
      for (int k=i+1; k<=j-1; k++){
	if(_Alpha_multi1[i][k-i] != -INF && _Alpha_multi2[k][j-k] != -INF){
	  temp = flag == 0 ? _Alpha_multi1[i][k-i]+_Alpha_multi2[k][j-k] : logsumexp(temp,_Alpha_multi1[i][k-i]+_Alpha_multi2[k][j-k]);
	  flag = 1;
	}
      }
      _Alpha_multibif[i][j-i] = flag == 0 ? -INF : temp;
      
      //Alpha_multi2
      temp = 0; flag = 0; 
      if (type != 0) {
	if(_Alpha_stem[i][j-i] != -INF){
	  temp = _Alpha_stem[i][j-i]+MLintern+CalcDangleEnergy(type,i,j);
	  flag = 1;
	}
      }
      if(_Alpha_multi2[i][j-i-1] != -INF){
	_Alpha_multi2[i][j-i] = _Alpha_multi2[i][j-i-1]+MLbase;
	if(flag == 1){
	  _Alpha_multi2[i][j-i] = logsumexp(temp,_Alpha_multi2[i][j-i]);
	}
      }else{
	_Alpha_multi2[i][j-i] = flag == 0 ? -INF : temp;
      }
      
      //Alpha_multi1
      if(_Alpha_multi2[i][j-i] != -INF && _Alpha_multibif[i][j-i] != -INF){
	_Alpha_multi1[i][j-i] = logsumexp(_Alpha_multi2[i][j-i],_Alpha_multibif[i][j-i]);
      }else if(_Alpha_multi2[i][j-i] == -INF){
	_Alpha_multi1[i][j-i] = _Alpha_multibif[i][j-i];
      }else if(_Alpha_multibif[i][j-i] == -INF){
	_Alpha_multi1[i][j-i] = _Alpha_multi2[i][j-i];
      }else{
	_Alpha_multi1[i][j-i] = -INF;
      }
      
      //Alpha_multi
      flag = 0;
      if(_Alpha_multi[i+1][j-i-1] != -INF){
	_Alpha_multi[i][j-i] = _Alpha_multi[i+1][j-i-1]+MLbase;
	flag = 1;
      }
      
      if(flag == 1){
	if(_Alpha_multibif[i][j-i] != -INF){
	  _Alpha_multi[i][j-i] = logsumexp(_Alpha_multi[i][j-i],_Alpha_multibif[i][j-i]);
	}
      }else{
	_Alpha_multi[i][j-i] = _Alpha_multibif[i][j-i];
      }
      
      //Alpha_stemend
      if(j != _seq_length){
	temp = 0;
	type = BP_pair[_int_sequence[i]][_int_sequence[j+1]];
	if (type!=0) {
	  //StemEnd¨sn
	  temp = HairpinEnergy(type, i,j+1);
	  
	  //StemEnd¨sm_Stem_sn
	  for (int p =i; p <= min(i+MAXLOOP,j-TURN-2); p++) {
	    int u1 = p-i;
	    for (int q=max(p+TURN+2,j-MAXLOOP+u1); q<=j; q++) {
	      type2 = BP_pair[_int_sequence[p+1]][_int_sequence[q]];
	      if(_Alpha_stem[p][q-p] != -INF){
		if (type2 != 0 && !(p == i && q == j)) {
		  type2 = rtype[type2];
		  temp = logsumexp(temp,_Alpha_stem[p][q-p]+LoopEnergy(type, type2,i,j+1,p+1,q)); 
		}
	      }
	    }
	  }
	  
	  //StemEnd¨Multi
	  int tt = rtype[type];
	  temp = logsumexp(temp,_Alpha_multi[i][j-i]+MLclosing+MLintern+dangle3[tt][_int_sequence[i+1]]+dangle5[tt][_int_sequence[j]]);
	  _Alpha_stemend[i][j-i] = temp;
	}else{
	  _Alpha_stemend[i][j-i] = -INF;
	}
      }
    }
  }
  
  //Alpha_Outer
  for(int i = 1;i <= _seq_length;i++){
    double temp = _Alpha_outer[i-1];
    for(int p = max(0,i-_maximal_span-1); p <i;p++){
      if(_Alpha_stem[p][i-p] != -INF){
	int type = BP_pair[_int_sequence[p+1]][_int_sequence[i]];
	double ao = _Alpha_stem[p][i-p]+CalcDangleEnergy(type,p,i);
	temp = logsumexp(temp,ao+_Alpha_outer[p]);
      }
    }
    _Alpha_outer[i] = temp;
  }
}

void CapR::CalcStructuralProfile(string &name){
  vector<double> bulge_probability; bulge_probability.resize(_seq_length, 0.0);
  vector<double> internal_probability; internal_probability.resize(_seq_length, 0.0);
  vector<double> hairpin_probability; hairpin_probability.resize(_seq_length, 0.0);
  vector<double> multi_probability; multi_probability.resize(_seq_length, 0.0);
  vector<double> exterior_probability; exterior_probability.resize(_seq_length, 0.0);
  vector<double> stem_probability; stem_probability.resize(_seq_length, 0.0);

  double pf = _Alpha_outer[_seq_length];
  if(pf >= -690 && pf <= 690){
    CalcBulgeAndInternalProbability(bulge_probability, internal_probability);
  }else{
    CalcLogSumBulgeAndInternalProbability(bulge_probability, internal_probability);}
  
  CalcHairpinProbability(hairpin_probability);
  
  for(int i = 1; i <=_seq_length;i++){
    exterior_probability[i-1] = CalcExteriorProbability(i);
    multi_probability[i-1] = CalcMultiProbability(i);
    stem_probability[i-1] = 1.0 - bulge_probability[i-1] - exterior_probability[i-1] - hairpin_probability[i-1] - internal_probability[i-1] - multi_probability[i-1];
  }
  ofstream of(_output_file_name.c_str(), ios::out | ios::app);
  of << ">" + name << endl;
  of << "Bulge ";
  for(int i = 0; i <_seq_length;i++){
    of << bulge_probability[i]<<" ";
  }
  of << endl;
  of << "Exterior ";
  for(int i = 0; i <_seq_length;i++){
    of << exterior_probability[i]<<" ";
  }
  of << endl;
  of << "Hairpin ";
  for(int i = 0; i <_seq_length;i++){
    of << hairpin_probability[i]<<" ";
  }
  of << endl;
  of << "Internal ";
  for(int i = 0; i <_seq_length;i++){
    of << internal_probability[i]<<" ";
  }
  of << endl;
  of << "Multibranch ";
  for(int i = 0; i <_seq_length;i++){
    of << multi_probability[i]<<" ";
  }
  of << endl;
  of << "Stem ";
  for(int i = 0; i <_seq_length;i++){
    of << stem_probability[i]<<" ";
  }
  of << endl;
  of << endl;
  of.close();
}

double CapR::CalcExteriorProbability(int x){
  double probability = exp(_Alpha_outer[x-1]+_Beta_outer[x]-_Alpha_outer[_seq_length]);
  return(probability);
}

void CapR::CalcHairpinProbability(vector<double> &hairpin_probability){
  for(int x = 1; x <=_seq_length;x++){
    double temp = 0.0;
    int type = 0;
    bool flag = 0;
    double h_energy = 0.0;
    
    for(int i = max(1,x-_maximal_span);i<x ;i++){
      for(int j = x+1; j<=min(i+_maximal_span,_seq_length);j++){
	type = BP_pair[_int_sequence[i]][_int_sequence[j]];
	if(_Beta_stemend[i][j-i-1] != -INF){
	  h_energy = _Beta_stemend[i][j-i-1] + HairpinEnergy(type, i,j);
	  temp = flag == 1 ? logsumexp(temp, h_energy) : h_energy;
	  flag = 1;
	}
      }
    }

    if(flag == 1){
      hairpin_probability[x-1] = exp(temp-_Alpha_outer[_seq_length]);
    }else{
      hairpin_probability[x-1] = 0.0;
    }
  }
}

double CapR::CalcMultiProbability(int x){
  double probability = 0.0;
  double temp = 0.0;
  bool flag = 0;
  
  for(int i = x; i<=min(x+_maximal_span,_seq_length);i++){
    if(_Beta_multi[x-1][i-x+1] != -INF && _Alpha_multi[x][i-x] != -INF){
      temp = flag == 0 ? _Beta_multi[x-1][i-x+1] + _Alpha_multi[x][i-x] : logsumexp(temp,_Beta_multi[x-1][i-x+1] + _Alpha_multi[x][i-x]);
      flag = 1;
    }
  }
  
  for(int i = max(0,x-_maximal_span); i<x;i++){
    if(_Beta_multi2[i][x-i] != -INF && _Alpha_multi2[i][x-i-1] != -INF){
      temp = flag == 0 ? _Beta_multi2[i][x-i] + _Alpha_multi2[i][x-i-1] : logsumexp(temp,_Beta_multi2[i][x-i] + _Alpha_multi2[i][x-i-1]);
      flag = 1;
    }
  }
  if(flag == 1){ probability = exp(temp-_Alpha_outer[_seq_length]); }
  return(probability);
}

void CapR::CalcBulgeAndInternalProbability(vector<double> &bulge_probability, vector<double> &internal_probability){
  double temp = 0;
  int type = 0;
  int type2 = 0;
 
  for(int i = 1; i<_seq_length-TURN-2;i++){
    for(int j = i+TURN+3; j<=min(i+_maximal_span,_seq_length);j++){
      type = BP_pair[_int_sequence[i]][_int_sequence[j]];
      if (type!=0) {
	for (int p =i+1; p <= min(i+MAXLOOP+1,j-TURN-2); p++) {
	  int u1 = p-i-1;
	  for (int q=max(p+TURN+1,j-MAXLOOP+u1-1); q<j; q++) {
	    type2 = BP_pair[_int_sequence[p]][_int_sequence[q]];
	    if (type2 != 0 && !(p == i+1 && q == j-1)) {
	      type2 = rtype[type2];
	      if(_Beta_stemend[i][j-i-1] != -INF && _Alpha_stem[p-1][q-p+1] != -INF){
		temp = exp(_Beta_stemend[i][j-i-1] + LoopEnergy(type, type2,i,j,p,q)+_Alpha_stem[p-1][q-p+1]);
		
		for(int k = i+1; k <= p-1;k++){
		  if(j == q+1){
		    bulge_probability[k-1] += temp;		   
		  }else{
		    internal_probability[k-1] += temp;		  
		  }
		}
		
		for(int k = q+1; k <= j-1;k++){
		  if(i == p-1){
		    bulge_probability[k-1] += temp;		   
		  }else{
		    internal_probability[k-1] += temp;		  
		  }
		}
	      } 
	    }
	  }
	}
      }
    }
  }
  
  for(int i=0;i<_seq_length;i++){
    if(bulge_probability[i] != 0){
      bulge_probability[i] /= exp(_Alpha_outer[_seq_length]);
    }
    if(internal_probability[i] != 0){
      internal_probability[i] /= exp(_Alpha_outer[_seq_length]);
    }
  }
}

void CapR::CalcLogSumBulgeAndInternalProbability(vector<double> &bulge_probability, vector<double> &internal_probability){
  double probability = 0;
  double temp = 0;
  int type = 0;
  int type2 = 0;
  vector<bool> b_flag_array; b_flag_array.resize(_seq_length,0);
  vector<bool> i_flag_array; i_flag_array.resize(_seq_length,0);
  
  for(int i = 1; i<_seq_length-TURN-2;i++){
    for(int j = i+TURN+3; j<=min(i+_maximal_span,_seq_length);j++){
      type = BP_pair[_int_sequence[i]][_int_sequence[j]];
      if (type!=0) {
	for (int p =i+1; p <= min(i+MAXLOOP+1,j-TURN-2); p++) {
	  int u1 = p-i-1;
	  for (int q=max(p+TURN+1,j-MAXLOOP+u1-1); q<j; q++) {
	    type2 = BP_pair[_int_sequence[p]][_int_sequence[q]];
	    if (type2 != 0 && !(p == i+1 && q == j-1)) {
	      type2 = rtype[type2];
	      if(_Beta_stemend[i][j-i-1] != -INF && _Alpha_stem[p-1][q-p+1] != -INF){
		temp = _Beta_stemend[i][j-i-1] + LoopEnergy(type, type2,i,j,p,q)+_Alpha_stem[p-1][q-p+1];
		
		for(int k = i+1; k <= p-1;k++){
		  if(j == q+1){
		    bulge_probability[k-1] = b_flag_array[k-1] == 1 ? logsumexp(bulge_probability[k-1], temp) : temp;
		    b_flag_array[k-1] = 1;		   
		  }else{
		    internal_probability[k-1] = i_flag_array[k-1] == 1 ? logsumexp(internal_probability[k-1], temp) : temp;
		    i_flag_array[k-1] = 1;		  
		  }
		}
		
		for(int k = q+1; k <= j-1;k++){
		  if(i == p-1){
		    bulge_probability[k-1] = b_flag_array[k-1] == 1 ? logsumexp(bulge_probability[k-1], temp) : temp;
		    b_flag_array[k-1] = 1;
		  }else{
		    internal_probability[k-1] = i_flag_array[k-1] == 1 ? logsumexp(internal_probability[k-1], temp) : temp;
		    i_flag_array[k-1] = 1;
		  }
		}
	      } 
	    }
	  }
	}
      }
    }
  }
  
  for(int i=0;i<_seq_length;i++){
    if(b_flag_array[i]==1){
      bulge_probability[i] = exp(bulge_probability[i]-_Alpha_outer[_seq_length]);
    }
    if(i_flag_array[i]==1){
      internal_probability[i] = exp(internal_probability[i]-_Alpha_outer[_seq_length]);
    }
  }
}

double CapR::CalcDangleEnergy(int type, int a, int b){
  double x = 0;
  if (type != 0) {
    if (a>0) x += dangle5[type][_int_sequence[a]];
    if (b<_seq_length) x += dangle3[type][_int_sequence[b+1]];
    if( b == _seq_length && type>2){
      x += TermAU;
    }
  }
  return(x);
}

void CapR::CalcOutsideVariable(){
  //Beta_outer
  for(int i = _seq_length-1;i >= 0;i--){
    double temp = _Beta_outer[i+1];
    for(int p = i+1; p <= min(i+_maximal_span+1,_seq_length);p++){
      if(_Alpha_stem[i][p-i] != -INF){
	int type = BP_pair[_int_sequence[i+1]][_int_sequence[p]];
	double bo = _Alpha_stem[i][p-i] + CalcDangleEnergy(type,i,p);
	temp = logsumexp(temp,bo+_Beta_outer[p]);
      }
    }
    _Beta_outer[i] = temp;	
  }
  
  for (int q=_seq_length; q>=TURN+1; q--) {
    for (int p=max(0,q-_maximal_span-1); p<= q-TURN; p++) {
      int type = 0;
      int type2 = 0;

      double temp = 0; bool flag = 0;
      if(p != 0 && q != _seq_length){
	//Beta_stemend
	_Beta_stemend[p][q-p] = q-p >= _maximal_span ? -INF : _Beta_stem[p-1][q-p+2];
	
	//Beta_Multi
	flag = 0;
	if(q-p+1 <= _maximal_span+1){
	  if(_Beta_multi[p-1][q-p+1] != -INF){
	    temp = _Beta_multi[p-1][q-p+1] + MLbase;
	    flag = 1;
	  }
	}
	
	type = BP_pair[_int_sequence[p]][_int_sequence[q+1]];
	int tt = rtype[type];
	if(flag == 1){
	  if(_Beta_stemend[p][q-p] != -INF){
	    temp = logsumexp(temp,_Beta_stemend[p][q-p]+MLclosing+MLintern+ dangle3[tt][_int_sequence[p+1]]+dangle5[tt][_int_sequence[q]]);
	  }
	}else{
	  if(_Beta_stemend[p][q-p] != -INF){
	    temp = _Beta_stemend[p][q-p]+MLclosing+MLintern+dangle3[tt][_int_sequence[p+1]]+dangle5[tt][_int_sequence[q]];
	  }else{
	    temp = -INF;
	  }
	}
	_Beta_multi[p][q-p] = temp;
	
	//Beta_Multi1
	temp = 0; flag = 0;
	for(int k = q+1 ; k<= min(_seq_length,p+_maximal_span);k++){
	  if(_Beta_multibif[p][k-p] != -INF && _Alpha_multi2[q][k-q] != -INF){
	    temp = flag == 0 ? _Beta_multibif[p][k-p]+_Alpha_multi2[q][k-q] : logsumexp(temp,_Beta_multibif[p][k-p]+_Alpha_multi2[q][k-q]) ;
	    flag = 1;
	  }
	}
	_Beta_multi1[p][q-p] = flag == 1 ? temp: -INF;
	
	//Beta_Multi2
	temp = 0; flag = 0;
	if(_Beta_multi1[p][q-p] != -INF){
	  temp = _Beta_multi1[p][q-p];
	  flag = 1;
	}
	if(q-p <= _maximal_span){
	  if(_Beta_multi2[p][q-p+1] != -INF){
	    temp = flag == 1 ? logsumexp(temp,_Beta_multi2[p][q-p+1]+MLbase) : _Beta_multi2[p][q-p+1]+MLbase;
	    flag = 1;
	  }
	}
	
	for(int k = max(0,q-_maximal_span); k < p ;k++){
	  if(_Beta_multibif[k][q-k] != -INF && _Alpha_multi1[k][p-k] != -INF){
	    temp = flag == 0 ? _Beta_multibif[k][q-k]+_Alpha_multi1[k][p-k] : logsumexp(temp,_Beta_multibif[k][q-k]+_Alpha_multi1[k][p-k]);
	    flag = 1;
	  }
	}
	_Beta_multi2[p][q-p] = flag == 0 ? -INF : temp;
	
	//Beta_multibif
	if(_Beta_multi1[p][q-p] != -INF && _Beta_multi[p][q-p] != -INF){
	  _Beta_multibif[p][q-p] = logsumexp(_Beta_multi1[p][q-p],_Beta_multi[p][q-p]);
	}else if(_Beta_multi[p][q-p] == -INF){
	  _Beta_multibif[p][q-p] = _Beta_multi1[p][q-p];
	}else if(_Beta_multi1[p][q-p] == -INF){
	  _Beta_multibif[p][q-p] = _Beta_multi[p][q-p];
	}else{
	  _Beta_multibif[p][q-p] = -INF;
	}
	
      }
      
      //Beta_stem
      type2 = BP_pair[_int_sequence[p+1]][_int_sequence[q]];
      if(type2 != 0){
	temp = _Alpha_outer[p]+_Beta_outer[q]+CalcDangleEnergy(type2,p,q);
	
	type2 = rtype[type2];
	for (int i=max(1,p-MAXLOOP); i<=p; i++){
	  for (int j=q; j<=min(q+ MAXLOOP -p+i,_seq_length-1); j++) {
	    type = BP_pair[_int_sequence[i]][_int_sequence[j+1]];
	    if (type != 0 && !(i == p && j == q)) {
	      if(j-i <= _maximal_span+1 && _Beta_stemend[i][j-i] != -INF){
		temp = logsumexp(temp,_Beta_stemend[i][j-i]+LoopEnergy(type,type2,i,j+1,p+1,q));
	      }
	    }
	  }
	}
	
	if(p != 0 && q != _seq_length){
	  type = BP_pair[_int_sequence[p]][_int_sequence[q+1]];
	  if(type != 0){
	    if(q-p+2 <= _maximal_span+1 && _Beta_stem[p-1][q-p+2] != -INF){
	      temp = logsumexp(temp,_Beta_stem[p-1][q-p+2]+LoopEnergy(type,type2,p,q+1,p+1,q));
	    }
	  }
	}
	_Beta_stem[p][q-p] = temp;
	
	if(_Beta_multi2[p][q-p] != -INF){
	  type2 = rtype[type2];
	  temp = _Beta_multi2[p][q-p] + MLintern + CalcDangleEnergy(type2,p,q);
	  _Beta_stem[p][q-p] = logsumexp(temp,_Beta_stem[p][q-p]);
	}
      }else{
	_Beta_stem[p][q-p] = -INF;
      }
    }
  }
}

double CapR::logsumexp(double x,double y){
  double temp = x > y ? x + log(exp(y-x) + 1.0) : y + log(exp(x-y) + 1.0) ;
  return(temp);
}

double CapR::LoopEnergy(int type, int type2,int i,int j,int p,int q){
  double z=0;
  int u1 = p-i-1;
  int u2 = j-q-1;
  
  if ((u1==0) && (u2==0)){
    z = stack[type][type2];
  }else{
    if ((u1==0)||(u2==0)) {
      int u;
      u = u1 == 0 ? u2 : u1;
      z = u <=30 ? bulge[u] : bulge[30] - lxc37*log( u/30.)*10./kT;
      
      if (u == 1){
	z += stack[type][type2];
      }else {
	if (type>2){ z += TermAU;}
	if (type2>2){ z += TermAU;}
      }
    }else{     
      if (u1+u2==2) {
	z = int11[type][type2][_int_sequence[i+1]][_int_sequence[j-1]];
      }else if ((u1==1) && (u2==2)){
	z = int21[type][type2][_int_sequence[i+1]][_int_sequence[q+1]][_int_sequence[j-1]];
      }else if ((u1==2) && (u2==1)){
	z = int21[type2][type][_int_sequence[q+1]][_int_sequence[i+1]][_int_sequence[p-1]];
      }else if ((u1==2) && (u2==2)){
	z = int22[type][type2][_int_sequence[i+1]][_int_sequence[p-1]][_int_sequence[q+1]][_int_sequence[j-1]];
      }else{
	z = internal[u1+u2]+mismatchI[type][_int_sequence[i+1]][_int_sequence[j-1]]+mismatchI[type2][_int_sequence[q+1]][_int_sequence[p-1]];
	z += ninio[abs(u1-u2)];
      }
    }
  }
  return z;
}

double CapR::HairpinEnergy(int type, int i, int j) {
  int d = j-i-1;
  double q = 0;

  q = d <= 30 ? hairpin[d] : hairpin[30] - lxc37*log( d/30.) *10./kT;  
  if(d!= 3){
    q += mismatchH[type][_int_sequence[i+1]][_int_sequence[j-1]];
  }else{
    if(type > 2){q += TermAU;}
  }
  return q;
}

void CapR::Clear(){
  for(int i = 0; i <= _seq_length;i++){
    _Alpha_stem[i].clear();
    _Alpha_stemend[i].clear();
    _Alpha_multi[i].clear();
    _Alpha_multibif[i].clear();
    _Alpha_multi1[i].clear();
    _Alpha_multi2[i].clear();
    _Beta_stem[i].clear();
    _Beta_stemend[i].clear();
    _Beta_multi[i].clear();
    _Beta_multibif[i].clear();
    _Beta_multi1[i].clear();
    _Beta_multi2[i].clear();
  }

  _int_sequence.clear();
  _seq_length = 0;
  _Alpha_outer.clear();
  _Alpha_stem.clear();
  _Alpha_stemend.clear();
  _Alpha_multi.clear();
  _Alpha_multibif.clear();
  _Alpha_multi1.clear();
  _Alpha_multi2.clear();
  
  _Beta_outer.clear();
  _Beta_stem.clear();
  _Beta_stemend.clear();
  _Beta_multi.clear();
  _Beta_multibif.clear();
  _Beta_multi1.clear();
  _Beta_multi2.clear();
}
