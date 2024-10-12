//
// Created by caronkey on 28/1/2021.
//

#ifndef JACCARD_BARCODE_H
#define JACCARD_BARCODE_H


#include <iterator>
#include <map>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fstream>
#include <set>
#include <algorithm>
#include "htslib/sam.h"
#include "htslib/faidx.h"

typedef std::unordered_map<std::string , std::vector<std::string>> BarcodeMap;
typedef std::pair<std::string, std::string> ContigPair;
typedef std::map<ContigPair, int> ContigMap;
typedef std::unordered_map<std::string, std::unordered_map<std::string, int>> Interactions;

void readBAM(htsFile *in, std::string &out_file, int readLen, std::vector<std::string>& fais);
#endif //JACCARD_BARCODE_H