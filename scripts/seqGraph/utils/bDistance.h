//
// Created by caronkey on 27/11/2020.
//

#ifndef JACCARD_BAM_H
#define JACCARD_BAM_H

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
#include <iostream>
#include <cstring>
#include "htslib/sam.h"
#include "htslib/faidx.h"

/* A contig end: (FASTA ID, head?) */
typedef std::pair<std::string, bool> CI;
typedef std::pair <std::string, std::string> PairString;

/* ScafMap: <pair(scaffold id, bool), count>, cout =  # times index maps to scaffold (c), bool =
 * true-head, false-tail*/
typedef std::map<CI, int> ScafMap;
/* IndexMap: key = index sequence, value = ScafMap */
typedef std::map <std::string, std::set<std::string>*> IndexMap;

/** maps contig FASTA ID to contig length (bp) */
typedef std::pair<std::string, bool> ScaffoldEnd;

typedef std::map<std::string, int> Jaccard;
typedef std::pair<std::string, int> SI;
typedef std::map <std::string, std::vector<SI>> NearScaf;


/**
 *
 * @param in
 * @param imap
 * @param min_size contig Min size
 * @param end_length
 */
void
readBAM(
        htsFile *in,
        IndexMap *imap,
        int min_size, int end_length);
void
nearScaf(IndexMap *imap, const char* out_file);


#endif //JACCARD_BAM_H
