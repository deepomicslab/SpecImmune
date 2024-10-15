//
// Created by caronkey on 27/11/2020.
//

#include "../include/cxxopts.hpp"
#include "e_result_reads.h"
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
// arg1 输入bam，arg2输出jaccard
float DEPTH = -1;
int CUTOFF = 300;
bool SELFLOOP = false;
bool ISTGS = false;
int THRESHOLD = 0;
int INSERT_SIZE = 500;
int INSERT_STD = 300;
int MIN_COUNT = 2;
int MAX_GAP = 500;

std::string RESULT;

void parse_fai(const std::string& contig_fai, std::unordered_map<std::string, int>& contigs) {
    std::ifstream file(contig_fai);
    std::string line;

//    std::vector<std::string> tokens;

    while(std::getline(file, line)) {     // '\n' is the default delimiter
//        std::vector<std::string> contigs;
        std::istringstream iss(line);
        std::string token;
        int pos = 0;
        std::string prev_token = "";
        while ((pos = line.find('\t')) != std::string::npos) {
            token = line.substr(0, pos);
//            std::cout << token << std::endl;
            if (prev_token != "") {
                contigs[prev_token + token] = 0;
            }
            prev_token = token;
            line.erase(0, pos + 1);
        }
//        std::getline(iss, token, '\t');
//        while(std::getline(iss, token, '\t'))   // but we can specify a different one
//        contigs.push_back(token);
    }
}

int main(int argc, char **argv) {

    //    parse options
    cxxopts::Options options("Reads overlap", "Extract reads and barcodes overlap from bam");

    options.add_options()
            ("fai", "contigs result", cxxopts::value<std::string>())
            ("b,bam", "Bam file", cxxopts::value<std::string>())
            ("o,out","Out bam", cxxopts::value<std::string>())
            ("insert_size","average insertion size", cxxopts::value<int>())
            ("insert_std","insertion size standard deviation", cxxopts::value<int>())
            ("min_count","contig pairs less min count will be break", cxxopts::value<int>())
//            ("max-gap","max gaps", cxxopts::value<int>())
            ("r,result","result N", cxxopts::value<std::string>())
            ("h,help", "Print usage");
    auto result = options.parse(argc,argv);
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    std::string bamF = result["bam"].as<std::string>();
    std::string resultF = result["out"].as<std::string>();
    std::string fai_file;
    if(result.count("fai")) {
        fai_file = result["fai"].as<std::string>();
    }
    if(result.count("result")) {
        RESULT = result["result"].as<std::string>();
    }

    if(result.count("insert_size")) {
        INSERT_SIZE = result["insert_size"].as<int>();
    }
    if(result.count("insert_std")) {
        INSERT_STD = result["insert_std"].as<int>();
    }

    if(result.count("min_count")) {
        MIN_COUNT = result["min_count"].as<int>();
    }

    htsFile *in;

    if ((in = hts_open(bamF.c_str(), "r")) == nullptr) {
        fprintf(stderr, "Error opening '%s'\n", bamF.c_str());
        return -1;
    }
    std::unordered_map<std::string, int> contigs;
    if (!fai_file.empty()) {
        parse_fai(fai_file,contigs);
    }
    readBAM(in, resultF, 100,contigs);

}

void initIMap (sam_hdr_t *hdr,Interactions& iMap, std::map<std::string, float>& refCopys) {
    int nRef = sam_hdr_nref(hdr);
    std::string refName, revRef;
    for(int i = 0; i< nRef; i++) {
        refName = sam_hdr_tid2name(hdr, i);
        std::stringstream ssf(refName);
        std::string item;
        while (std::getline(ssf,item,'_'));
        int idx = item.find(':');
        float cov = 0;
        if (idx != -1) {
            cov = std::stod(item.substr(0,idx));
        } else {
            cov = std::stod(item);
        }
        refCopys.emplace(refName, cov);
        revRef = refName+" -";
        std::unordered_map<std::string, std::vector<bam1_t*>> tmp;
        std::unordered_map<std::string, std::vector<bam1_t*>> rvtmp;
        iMap[refName] = tmp;
        iMap[revRef]= rvtmp;
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
        // elems.push_back(std::move(item)); // if C++11 (based on comment from @mchiasson)
    }
    return elems;
}

int splitstrcigar(char* cigar, bool first) {
    std::string cstr = std::string(cigar);
    if (first) {
        int idx = cstr.find('S');
        if (idx > 0)
            return std::stoi(cstr.substr(0,idx));
    } else {
        int idx = cstr.find_last_of('M');
        if (idx > 0) {
            auto t = cstr.substr(idx+1, cstr.size() - 2 - idx);
            if ( t != "")
                return std::stoi(t);
        }
    }
    return 0;
}

int splitcigar(unsigned int * cigar,int n_cigar, bool first) {
    bam_cigar_op(cigar[0]);
    if (first) {
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
            return bam_cigar_oplen(cigar[0]);
        }
    } else {
        if (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CSOFT_CLIP) {
            return bam_cigar_oplen(cigar[n_cigar - 1]);
        }
    }
    return 0;
}


int nSize(std::vector<bam1_t*>& recs, bam_hdr_t* hdr) {
//    int negative_size = int (((float )INSERT_STD/(float)INSERT_SIZE) * (float )INSERT_STD);
    int negative_size = int (((float )INSERT_STD/(float)INSERT_SIZE) * (float )INSERT_STD) + INSERT_STD;
    int refLen = sam_hdr_tid2len(hdr, recs[0]->core.tid);
    int mrefLen = sam_hdr_tid2len(hdr, recs[0]->core.mtid);
    int total_gap = 0;
    int count = 0;
    for (auto item : recs) {
        int refpadding = 0;
        int mrefPadding = 0;
        int refunMatchs = 0;
        int mrefunMatchs = 0;
        auto mcigar = bam_aux2Z(bam_aux_get(item,"MC"));

//        auto op_len = bam_cigar_oplen(mcigar[1]);
        if (item->core.pos > refLen/2) {
            int pos_end = item->core.pos + item->core.l_qseq;
            pos_end = std::min(refLen, pos_end);
            refpadding = refLen - pos_end;
            refunMatchs = splitcigar(bam_get_cigar(item), item->core.n_cigar, false);
        } else {
            int pos_end = item->core.pos;
            refpadding = pos_end;
            refunMatchs = splitcigar(bam_get_cigar(item), item->core.n_cigar, true);
        }

        if (item->core.mpos > mrefLen/2) {
            int pos_end = item->core.mpos + item->core.l_qseq;
            pos_end = std::min(mrefLen, pos_end);
            mrefPadding = mrefLen - pos_end;
            mrefunMatchs = splitstrcigar(mcigar, false);
        } else {
            int pos_end = item->core.mpos;
            mrefPadding = pos_end;
            mrefunMatchs = splitstrcigar(mcigar, true);
        }
//        if (refpadding > INSERT_SIZE/2 || mrefPadding > INSERT_SIZE/2) continue;

//        int gapsC = INSERT_SIZE + 50 - (refpadding + mrefPadding) - (item->core.l_qseq - refunMatchs) - (item->core.l_qseq - mrefunMatchs);
        int gapsC = (INSERT_SIZE - 2*item->core.l_qseq) + (refunMatchs + mrefunMatchs) - (refpadding + mrefPadding);
        if (gapsC > INSERT_SIZE + INSERT_STD || gapsC < -negative_size) {
            continue;
        }
//        if (gapsC > 0)
        total_gap += gapsC;
        count++;
    }
    if (count < MIN_COUNT) {
        return -1000000;
    }
    if (count != 0) {
        int gap_avg = total_gap/count;
        return gap_avg;
    }
//        return ;
}

void readBAM(htsFile *in, std::string& out_file, int readsLen, std::unordered_map<std::string, int>& fais) {
    Interactions iMap;
    sam_hdr_t *hdr;
    bam1_t *b;
    int ret;
    std::map<std::string, float> refCopys;
    if ((hdr = sam_hdr_read(in)) == nullptr) {
        fprintf(stderr, "[E::%s] couldn't read file header \n", __func__);
        return;
    }
    if ((b = bam_init1()) == nullptr) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    auto fpout = sam_open(out_file.c_str(), "w");
    auto state = sam_hdr_write(fpout,hdr);
    std::string readName, refName, mRefName;
    initIMap(hdr, iMap, refCopys);
    int pos, mpos, refLen, mRefLen;
//    int times = 1;
    while ((ret = sam_read1(in, hdr, b)) >= 0) {
//        if (times > 1000000) {
//            break;
//        } times++;
        if (ret < -1) {
            fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
            if (b) bam_destroy1(b);
            if (hdr) sam_hdr_destroy(hdr);
        }
        pos = b->core.pos;
        mpos = b->core.mpos;
        auto flags = b->core.flag;
        if (b->core.tid == -1 || b->core.mtid == -1) continue;
//        if (flags & 0x8) continue;
//        if (flags & 0x4) continue;
        if (ISTGS && !(flags & 0x800)) continue;
        auto rev = flags & 0x10;
        auto mrev = flags & 0x20;
        if (ISTGS) {
            auto p = bam_aux_get(b, "SA");
            std::string pstring = std::string(reinterpret_cast<const char *>(p));
            auto pvector = split(pstring,',');
//            refName = pvector[0].substr(1);
//            mRefName = std::string(sam_hdr_tid2name(hdr, b->core.tid));
//            rev = pvector[2] == "+" ? 0:1 ;
            auto cigars = bam_get_cigar(b);
            auto op=bam_cigar_op(cigars[0]);
            if (op == BAM_CHARD_CLIP) {
                continue;
            }
            auto op_len = bam_cigar_oplen(cigars[0]);
            int indexS = pvector[3].find('S');
            int indexH = pvector[3].find('H');
//            if (indexH < indexS) continue;
            int minidx = std::min(indexS, indexH);
            int mop_len = std::stoi(pvector[3].substr(0,minidx));
            mpos =  std::stoi(pvector[1]);
            if ((pos > 2000 && mpos > 2000) || (pos < 2000 && mpos < 2000) ) continue;
            if (op_len <= mop_len) {
                refName = std::string(sam_hdr_tid2name(hdr, b->core.tid));
                mRefName = pvector[0].substr(1);
                rev = 0;
                mrev  = 1 ;
            } else {
                mRefName = std::string(sam_hdr_tid2name(hdr, b->core.tid));
                refName = pvector[0].substr(1);
                rev  = 0;
                mrev = 1;
            }
//            mrev =
        } else {
            if ((b->core.tid == -1 || b->core.mtid == -1))
                continue;
            refLen = sam_hdr_tid2len(hdr, b->core.tid);
            mRefLen = sam_hdr_tid2len(hdr, b->core.mtid);
            refName = std::string(sam_hdr_tid2name(hdr, b->core.tid));
            mRefName = std::string(sam_hdr_tid2name(hdr, b->core.mtid));
        }
//        refLen = sam_hdr_tid2len(hdr, b->core.tid);
//        mRefLen = sam_hdr_tid2len(hdr, b->core.mtid);
//        if (refName == mRefName)continue;
//        if (fais.find())
//        if (!fais.empty() && (std::find(fais.begin(),fais.end(),refName) == fais.end() || std::find(fais.begin(),fais.end(),mRefName) == fais.end())) continue;
        if (rev)
            refName = refName.append("-");
        else
            refName = refName.append("+");
        if (!mrev)
            mRefName = mRefName.append("-");
        else
            mRefName = mRefName.append("+");
        if (fais.find(refName+mRefName) == fais.end())
            continue;
//        if (refLen <= readsLen) {
//            if (pos == 0)
//                refName = refName.append(" -");
//        } else if (pos < refLen/2 - readsLen/2) {
//            refName = refName.append(" -");
//        }
//        if (mRefLen <= readsLen) {
//            if (pos != 0)
//                mRefName = mRefName.append(" -");
//        } else if (pos >= mRefLen/2 - readsLen/2) {
//            mRefName = mRefName.append(" -");
//        }
        auto& refMap = iMap[refName];
        if (refMap.find(mRefName) == refMap.end()) {
            std::vector<bam1_t*> breads;
            bam1_t* tmp = bam_init1();
            bam_copy1(tmp,b);
            breads.push_back(tmp);
            refMap[mRefName] = breads;
        } else {
            bam1_t* tmp = bam_init1();
            bam_copy1(tmp,b);
            refMap[mRefName].push_back(tmp);
//            int ncount =  nSize(refMap[mRefName], hdr);
        }
    }
    std::ofstream result_out(RESULT);
    std::string prev;
    std::string next;
//    for (const auto& item : refCopys ){
//        if (!fais.empty() && std::find(fais.begin(),fais.end(),item.first) == fais.end()) continue;
//        int copy;
//        if (DEPTH != -1) {
//            auto hapd = item.second/DEPTH;
//            auto hapdF = std::floor(hapd);
//            if (hapd - hapdF > 0.7) {
//                copy = hapdF + 1;
//            } else {
//                copy = hapdF;
//            }
//            if (copy == 0) copy = 1;
////            int copy = int(cov/DEPTH) == 0 ? 1 : int(cov/DEPTH);
//            refCopys.emplace(refName, copy);
//        } else {
//            copy = 1;
//        }
//        if (ISTGS) {
//            copy = 1;
//        }
//        fout<<"SEG "<<item.first<<" "<<item.second<<" "<<copy<<"\n";
//    }
    for (auto& it: iMap) {
        for(auto& it2 : it.second) {
//            if (fais.find(it.first+it2.first) == fais.end()) {
//                continue;
//            }
            int ncount =  nSize(it2.second, hdr);
            result_out<<it.first<<"\t"<<it2.first<<"\t"<<ncount<<std::endl;
            for (auto item: it2.second) {
                state = sam_write1(fpout,hdr,item);
                int temp = 4;
            }
        }
//        if (it.second > 50) fout<<it.first.first<<"\t"<<it.first.second<<"\t"<<it.second<<"\t"<<"\n";
    }
    sam_close(fpout);
}



