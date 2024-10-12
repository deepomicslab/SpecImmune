//
// Created by caronkey on 2022/3/18.
//

#ifndef SEQGRAPH_UTIL_H
#define SEQGRAPH_UTIL_H
namespace seqGraph {
   inline int conjugateIdx(int idx) {
        int cI;
        if(idx == 0) return idx;
        cI = idx%2 == 0? idx-1:idx+1;
        return cI;
    }
    inline void conjugatePath(std::vector<int>* path) {
        for (auto &it : *path) {
            it = conjugateIdx(it);
        }
        std::reverse(path->begin(), path->end());
   }
}
#endif //SEQGRAPH_UTIL_H

