//
// Created by caronkey on 27/11/2020.
//

#include "reads_overlap.h"
#include "../include/cxxopts.hpp"
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
// arg1 输入bam，arg2输出jaccard
float DEPTH = -1;
int CUTOFF = 300;
bool SELFLOOP = false;
bool ISTGS = false;
int THRESHOLD = 0;

void parse_fai(const std::string& contig_fai, std::vector<std::string>& contigs) {
    std::ifstream file(contig_fai);
    std::string line;

//    std::vector<std::string> tokens;

    while(std::getline(file, line)) {     // '\n' is the default delimiter

        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, '\t');
//        while(std::getline(iss, token, '\t'))   // but we can specify a different one
        contigs.push_back(token);
    }
}

int main(int argc, char **argv) {

    //    parse options
    cxxopts::Options options("Reads overlap", "Extract reads and barcodes overlap from bam");

    options.add_options()
            ("fai", "contigs fai file", cxxopts::value<std::string>())
            ("b,bam", "Bam file", cxxopts::value<std::string>())
            ("o,out","Out file", cxxopts::value<std::string>())
            ("c,cut_off", "Cute off size", cxxopts::value<int>()->default_value("300"))
            ("d,depth", "Bam depth", cxxopts::value<float>())
            ("t,threshold", "Threshold for graph weigh filter", cxxopts::value<int>()->default_value("0"))
            ("s,self_l", "If self_loop take into consideration", cxxopts::value<bool>()->default_value("false"))
            ("tgs", "If self_loop take into consideration", cxxopts::value<bool>()->default_value("false"))
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

    if(result.count("cut_off")) {
        CUTOFF = result["cut_off"].as<int>();
    }
    if(result.count("depth")) {
        DEPTH = result["depth"].as<float>();
    }
    if(result.count("self_l")) {
        SELFLOOP = result["self_l"].as<bool>();
    }
    if(result.count("self_l")) {
        SELFLOOP = result["self_l"].as<bool>();
    }
    if(result.count("tgs")) {
        ISTGS = result["tgs"].as<bool>();
    }
    if(result.count("threshold")) {
        THRESHOLD = result["threshold"].as<int>();
    }



    htsFile *in;

    if ((in = hts_open(bamF.c_str(), "r")) == nullptr) {
        fprintf(stderr, "Error opening '%s'\n", bamF.c_str());
        return -1;
    }
    std::vector<std::string> fais;
    if (!fai_file.empty()) {
        parse_fai(fai_file,fais);
    }
    readBAM(in, resultF, 100,fais);

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
        std::unordered_map<std::string, int> tmp;
        std::unordered_map<std::string, int> rvtmp;
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


void readBAM(htsFile *in, std::string& out_file, int readsLen, std::vector<std::string>& fais) {
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
    std::string readName, refName, mRefName;
    initIMap(hdr, iMap, refCopys);
    int pos, mpos, refLen, mRefLen;
    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        if (ret < -1) {
            fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
            if (b) bam_destroy1(b);
            if (hdr) sam_hdr_destroy(hdr);
        }
        pos = b->core.pos;
        mpos = b->core.mpos;
        auto flags = b->core.flag;
        if (flags & 0x80) continue;
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
            if ((pos > 1000 && mpos > 1000) || (pos < 1000 && mpos < 1000) ) continue;
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
            if (!ISTGS && (b->core.tid == -1 || b->core.mtid == -1))
                continue;
            refLen = sam_hdr_tid2len(hdr, b->core.tid);
            mRefLen = sam_hdr_tid2len(hdr, b->core.mtid);
            refName = std::string(sam_hdr_tid2name(hdr, b->core.tid));
            if (!ISTGS)
                mRefName = std::string(sam_hdr_tid2name(hdr, b->core.mtid));
            if (!ISTGS && CUTOFF != -1) {
                if (pos < refLen - CUTOFF && mpos > CUTOFF) continue;
                if (mpos < mRefLen - CUTOFF && pos > CUTOFF) continue;
            }
        }
//        refLen = sam_hdr_tid2len(hdr, b->core.tid);
//        mRefLen = sam_hdr_tid2len(hdr, b->core.mtid);
//        if (refName == mRefName) {
//            if (!SELFLOOP) continue;
//            else if(rev == mrev) continue;
//        }
        if (refName == mRefName)continue;
        if (!fais.empty() && (std::find(fais.begin(),fais.end(),refName) == fais.end() || std::find(fais.begin(),fais.end(),mRefName) == fais.end())) continue;
        if (rev)
            refName = refName.append(" -");
        if (!mrev)
            mRefName = mRefName.append(" -");
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
            refMap[mRefName] = 1;
        } else {
            refMap[mRefName]+=1;
        }
    }
    std::ofstream fout(out_file);
    std::string prev;
    std::string next;
    for (const auto& item : refCopys ){
        if (!fais.empty() && std::find(fais.begin(),fais.end(),item.first) == fais.end()) continue;
        int copy;
        if (DEPTH != -1) {
            auto hapd = item.second/DEPTH;
            auto hapdF = std::floor(hapd);
            if (hapd - hapdF > 0.7) {
                copy = hapdF + 1;
            } else {
                copy = hapdF;
            }
            if (copy == 0) copy = 1;
//            int copy = int(cov/DEPTH) == 0 ? 1 : int(cov/DEPTH);
            refCopys.emplace(refName, copy);
        } else {
            copy = 1;
        }
        if (ISTGS) {
            copy = 1;
        }
        fout<<"SEG "<<item.first<<" "<<item.second<<" "<<copy<<"\n";
    }
    for (auto& it: iMap) {
        for(auto& it2 : it.second) {
            if (it2.second < THRESHOLD) continue;
            if (it.first.back() != '-') prev = it.first + " +";
            else prev = it.first;
            if (it2.first.back() != '-') next = it2.first + " +";
            else next = it2.first;
            fout<<"JUNC "<<prev<<" "<<next<<" "<<it2.second<<"\n";
        }
//        if (it.second > 50) fout<<it.first.first<<"\t"<<it.first.second<<"\t"<<it.second<<"\t"<<"\n";
    }
    fout.close();
}

