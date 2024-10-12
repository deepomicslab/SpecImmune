//
// Created by caronkey on 27/11/2020.
//

#include "bDistance.h"
void print_help() {
    std::cout<<"Extract the barcode distance information from bam file"<<std::endl;
    std::cout<<"usage:\n\tbDistance input.bam graph.txt 500"<<std::endl;
    std::cout<<"parameters:\n";
    std::cout<<"\tinput.bam: bam file\n";
    std::cout<<"\tgraph.txt: output graph file\n";
    std::cout<<"\t500: length of the two ends for barcode counting\n";
}
// arg1 输入bam，arg2输出jaccard, arg3 end_length
int main(int argc, char **argv) {
    auto md = argv[1];
    if (strcmp(md, "help") == 0) {
        print_help();
        exit(0);
    }
    auto imap = new IndexMap();
    htsFile *in;

    if ((in = hts_open(argv[1], "r")) == nullptr) {
        fprintf(stderr, "Error opening '%s'\n", argv[1]);
        return -1;
    }
    int end_length = std::stoi(argv[3]);
    readBAM(in, imap, 0, end_length);
    std::cout<<"read done"<<std::endl;
    nearScaf(imap, argv[2]);
    return 0;
}
void
readBAM(
        htsFile *in,
        IndexMap *imap,
        int min_size, int end_length) {
    sam_hdr_t *hdr;
    bam1_t *b;
    int ret;

    std::string prevRN, readyToAddIndex, prevRef, readyToAddRefName;
    int prevPos = -1, readyToAddPos = -1;
    int ct = 1;
    int readRefSize = 0;
    std::string index;
    std::string readName, scafName, seq, tp_scaf;
    std::size_t found;
//    std::vector<std::string> *tp;

    // Number of unpaired reads.
    size_t countUnpaired = 0;
    /* Read each line of the BAM file */
    if ((hdr = sam_hdr_read(in)) == nullptr) {
        fprintf(stderr, "[E::%s] couldn't read file header \n", __func__);
        return;
    }
    if ((b = bam_init1()) == nullptr) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    while ((ret = sam_read1(in, hdr, b)) >= 0) {
//        std::cout<<b->id<<std::endl;
        if (ret < -1) {
            fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
            if (b) bam_destroy1(b);
            if (hdr) sam_hdr_destroy(hdr);
        }
        int pos;

        readName = bam_get_qname(b);
        if (b->core.tid == -1) {
            scafName = "*";
        } else {
            scafName = sam_hdr_tid2name(hdr, b->core.tid);
        }
        pos = b->core.pos;

        /* Parse the index from the readName */
        found = readName.rfind('#');
        if (found != std::string::npos) {
            index = readName.substr(found + 1);
        }
        if (index == "0_0_0") continue;

        if (ct == 2 && readName != prevRN) {
            if (countUnpaired == 0)
                std::cerr << "Warning: Skipping an unpaired read. Read pairs should be "
                             "consecutive in the SAM/BAM file.\n"
                             "  Prev read: "
                          << prevRN
                          << "\n"
                             "  Curr read: "
                          << readName << std::endl;
            ++countUnpaired;
            if (countUnpaired % 1000000 == 0)
                std::cerr << "Warning: Skipped " << countUnpaired << " unpaired reads."
                          << std::endl;
            ct = 1;
        }

        if (ct >= 3)
            ct = 1;
        if (ct == 1) {
            if (readName != prevRN) {
                prevRN = readName;
                prevRef = scafName;
                prevPos = pos;

                /*
                 * Read names are different so we can add the previous index and scafName as
                 * long as there were only two mappings (one for each read)
                 */
                if (!readyToAddIndex.empty() && !readyToAddRefName.empty() &&
                    readyToAddRefName != "*" && readyToAddPos != -1 && readRefSize != 0) {

                    int size = readRefSize;
                    if (size >= min_size) {

                        /*
                         * If length of sequence is less than 2 x end_length, split
                         * the sequence in half to determing head/tail
                         */
                        int cutOff = end_length;
                        if (cutOff == 0 || size <= cutOff * 2)
                            cutOff = size / 2;

                        /*
                         * pair <X, true> indicates read pair aligns to head,
                         * pair <X, false> indicates read pair aligns to tail
                         */
//                        ScaffoldEnd key(readyToAddRefName, true);
//                        ScaffoldEnd keyR(readyToAddRefName, false);

                        /* Aligns to head */
                        if (readyToAddPos <= cutOff) {
                            tp_scaf = readyToAddRefName.append("\t+");
                            if (imap->find(readyToAddIndex) == imap->end()) {
                                auto tp = new std::set<std::string>();
//                                imap->insert(std::make_pair("xxx", tp));
                                imap->insert(std::make_pair(readyToAddIndex, tp));
//                                (*imap[readyToAddIndex] = tp;
                            }
                            (*imap)[readyToAddIndex]->insert(tp_scaf);
//                            (*imap)[readyToAddIndex][key]++;

                            /* Aligns to tail */
                        } else  {
                            tp_scaf = readyToAddRefName.append("\t-");
                            if (imap->find(readyToAddIndex) == imap->end()) {
                                auto tp = new std::set<std::string>();
                                imap->insert(std::make_pair(readyToAddIndex, tp));
                            }
                            (*imap)[readyToAddIndex]->insert(tp_scaf);
//                            (*imap)[readyToAddIndex]->push_back(tp_scaf);
                        }
                    }
                    readyToAddIndex = "";
                    readyToAddRefName = "";
                    readyToAddPos = -1;
                    readRefSize = 0;
                }
            } else {
                ct = 0;
                readyToAddIndex = "";
                readyToAddRefName = "";
                readyToAddPos = -1;
                readRefSize = 0;
            }
        } else if (ct == 2) {
            if (prevRef == scafName && scafName != "*" &&
                !scafName.empty() && !index.empty()) {

                readyToAddIndex = index;
                readyToAddRefName = scafName;
                readyToAddPos = (prevPos + pos) / 2;
                readRefSize = sam_hdr_tid2len(hdr, b->core.tid);
            }
        }
        ct++;
//        bam_destroy1(b);
    }

    /* Close BAM file */


    if (countUnpaired > 0)
        std::cerr << "Warning: Skipped " << countUnpaired
                  << " unpaired reads. Read pairs should be consecutive in the SAM/BAM file.\n";
}

void nearScaf(IndexMap *imap, const char *out_file) {
    Jaccard jaccard; // 0: hh, 1:ht, 2: th, 3:tt
    NearScaf near_scaf;
    std::string index;
    std::string scafA, scafB;
    std::string tmp;
    bool scafAflag, scafBflag;
    for (auto &it : *imap) {
        if (it.second->size() == 1) continue;
        for (auto &o : *it.second) {
            if (o.back() == '+'){
                scafA = o;
                scafA[o.size()-1] = '-';
            } else {
                scafA = o;
                scafA[o.size()-1] = '+';
            }
            auto orig = scafA;
            for (auto &p : *it.second) {
                if(p<o) continue;
                scafA = orig;
                tmp = scafA.append("\t" + p);
                if (jaccard.find(tmp) == jaccard.end()) {
                    jaccard[tmp] =1;
                } else{
                    jaccard[tmp] +=1;
                }
//                jaccard[]
//                std::tie(scafA, scafAflag) = o->first;
//                std::tie(scafB, scafBflag) = p.first;

                /* Only insert into pmap if scafA < scafB to avoid duplicates */
//                if (scafA < scafB) {
//                    if (scafAflag && scafBflag) {
//                        jaccard[scafA + "\t" + scafB][0]++;
//                    } else if (scafAflag) {
//                        jaccard[scafA + "\t" + scafB][1]++;
//                    } else if (scafBflag) {
//                        jaccard[scafA + "\t" + scafB][2]++;
//                    } else {
//                        jaccard[scafA + "\t" + scafB][3]++;
//                    }
//                }
            }
        }
    }
    std::ofstream fout(out_file);
//    统计barcode数量
//    std::string k1;
//    std::string k2;
//    std::istringstream iss(line);
//    std::string k1, k2;
//    while(std::getline(iss, token, '\t'))   // but we can specify a different one
//        tokens.push_back(token);
    for (auto &it: jaccard) {
//         if (it.second[0] <= 10) continue;
//        std::istringstream iss(it.first);
//        iss >> k1 >> k2;
////         k1 = it.first.substr(0, 16);
////         k2 = it.first.substr(17, 34);
        fout << it.first << "\t" << it.second << "\n";
//        fout << k1 << "\t" << "+\t" << k2 << "\t"<< "-\t" << it.second[1] << "\n";
//        fout << k1 << "\t" << "-\t" << k2 << "\t"<< "+\t" << it.second[2] << "\n";
//        fout << k1 << "\t" << "-\t" << k2 << "\t"<< "-\t" << it.second[3] << "\n";

    }
//char * token;
//   for (auto &it : jaccard) {
//       char *dup = strdup(it.first.c_str());
//       token = std::strtok(dup, "\t");
//       std::string k1 = token;
//       std::string k2 = strtok(nullptr, "\t");
//       SI s1;
//       SI s2;
//       bool flag = near_scaf.find(k1) != near_scaf.end();
//       if (it.second[0] <= 5 &&
//           it.second[1] <= 5 &&
//           it.second[2] <= 5 &&
//           it.second[3] <= 5)
//           continue;
//       if (flag) {
//           if (it.second[0] > near_scaf[k1][0].second && it.second[0] >= it.second[1])
//               s1 = SI(k2 + "\t1", it.second[0]);
//           else if (it.second[1] >= near_scaf[k1][0].second && it.second[1] >= it.second[0])
//               s1 = SI(k2 + "\t0", it.second[1]);
//
//           if (it.second[2] > near_scaf[k1][1].second && it.second[2] >= it.second[3])
//               s2 = SI(k2 + "\t1", it.second[2]);
//           else if (it.second[3] > near_scaf[k1][1].second && it.second[3] >= it.second[2])
//               s2 = SI(k2 + "\t0", it.second[3]);
//
//           if (s1.second != 0) near_scaf[k1][0] = s1;
//           if (s2.second != 0) near_scaf[k1][1] = s2;
//       } else {
//           std::vector <SI> ss;
//           if (it.second[0] > it.second[1])
//               s1 = SI(k2 + "\t1", it.second[0]);
//           else if (it.second[1] >= it.second[0])
//               s1 = SI(k2 + "\t0", it.second[1]);
//
//           if (it.second[2] > it.second[3])
//               s2 = SI(k2 + "\t1", it.second[2]);
//           else if (it.second[3] >= it.second[2])
//               s2 = SI(k2 + "\t0", it.second[3]);
//
//           ss.push_back(s1);
//           ss.push_back(s2);
//           near_scaf[k1] = ss;
//       }
//       free(dup);
//   }
//    for (auto &it : near_scaf) {
//        if(it.second[0].second >= 5 && it.second[0].second >= it.second[1].second) {
//            fout << it.first << "\t1\t" << it.second[0].first << "\t" << it.second[0].second << "\n";
//        } else if (it.second[1].second >= 5 && it.second[0].second <= it.second[1].second) {
//            fout << it.first << "\t0\t" << it.second[1].first << "\t" << it.second[1].second << "\n";
//        }
//    }
    fout.close();
}
