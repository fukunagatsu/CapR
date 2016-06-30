/*
 * fastafile_reader.h
 *
 *  Created on: 2016/6/29
 *      Author: Tsukasa Fukunaga
 */

#ifndef FASTAFILE_READER_H
#define FASTAFILE_READER_H

#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;

class FastafileReader {
 public:
  FastafileReader() {}
  void ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names);
};

#endif
