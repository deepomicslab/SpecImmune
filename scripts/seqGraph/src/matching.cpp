//
// Created by caronkey on 28/12/2021.
//

#include "../include/matching.h"
#include <cstring>
#include <algorithm>
#include <queue>
#include <iostream>
#include "limits"
#include <fstream>

const float ZERO = 0.000001;
//
bool cmp(std::pair<int, std::vector<int>* >& a,
         std::pair<int, std::vector<int>* >& b)
{
    return a.second->size() < b.second->size();
}
//TODO, 需要在迭代过程中加入这个过程

int Partition(matching* m,float uCov, std::vector<int> &v, int start, int end){

    int pivot = end;
    int j = start;
    for(int i=start;i<end;++i){
        if( uCov - m->idx2VertexInCurrentGraph(v[i])->getWeight()->getCoverage() < uCov - m->idx2VertexInCurrentGraph(v[pivot])->getWeight()->getCoverage()){
            std::swap(v[i],v[j]);
            ++j;
        }
    }
    std::swap(v[j],v[pivot]);
    return j;

}

void Quicksort(matching* m ,float uCov, std::vector<int> &v, int start, int end ){

    if(start<end){
        int p = Partition(m,uCov, v,start,end);
        Quicksort(m,uCov, v,start,p-1);
        Quicksort(m,uCov, v,p+1,end);
    }

}


//bool matching::cmpVertex (int i,int j) {
//    auto v1 = this->idx2VertexInCurrentGraph(i)->getWeight()->getCoverage();
//    auto v2 = this->idx2VertexInCurrentGraph(j)->getWeight()->getCoverage();
//    return v1 < v2;
//}
// Function to sort the map according
// to value in a (key-value) pairs
void sort(std::map<int, std::vector<int>*>& M, std::vector<std::pair<int, std::vector<int>*>>& A)
{
    // Copy key-value pair from Map
    // to vector of pairs
    A.reserve(M.size());
    for (auto& it : M) {
        A.emplace_back(it);
    }

    // Sort using comparator function
    std::sort(A.begin(), A.end(), cmp);
}

void matching::printM(int i){
    std::cout<<i<<"|";
    for (int j = 0; j <  N; j++){
        std::cout<<this->getIJ(i,j)<<"\t";
    }
}
const float INF = std::numeric_limits<float>::max();
matching::matching(seqGraph::Graph* graph1) {
    this->graph = graph1;
    this->graph->initRowMax();
//    this->graph->removeByGeneAndScore();
    N = 2 * graph1->getVCount();
    this->matched = new int[N + 1];
    std::fill_n(this->matched,N+1, -1);
//    for (int i = 0 ; i < N + 1; i++) this->matched[i] = -1;
//    currentMatrix = this->graph->getConjugateMatrix();
//    this->originalGraph = graph1;
//    this->originalMatrix = new float*[N+1];
//    for(int i = 0; i < N+1; ++i)
//        originalMatrix[i] = new float [N+1];
//    std::memcpy(originalMatrix, currentMatrix, sizeof(int)*(N+1)*(N+1));
//    this->originalVertices = new std::vector<seqGraph::Vertex*>();
    this->originalGraph = graph1;
//    for (auto item : *this->graph->getVertices()) {
//        this->originalVertices->push_back(new seqGraph::Vertex(*item));
//    }
//    this->originalJunctions = new std::vector<seqGraph::Junction*>();
//    for (auto item : *this->graph->getJunctions()) {
//        this->originalJunctions->push_back(new seqGraph::Junction(*item));
//    }
//    this->originalVertices = this->graph->getVertices();
    if (VERBOSE == 2) {
        this->graph->getConjugateMatrix().debugPrint();
    }
}

matching::~matching() {
//    delete this->matched;
    if (graph != nullptr)
        free(graph);
//    for (auto item: *originalVertices) {
//        delete item;
//    }
//    originalVertices->clear();
//    originalVertices->shrink_to_fit();
//    for (int i = 0; i < N+1; ++i) {
//        free(originalMatrix[i]);
//    }
//    free(originalMatrix);
}



void matching::resetGraph(seqGraph::Graph* g) {
    if (this->graph->isReconstructed)
        free(this->graph);
    std::cout<<"free done"<<std::endl;
    this->graph =  g;
    this->graph->initRowMax();
    N = 2 *  g->getVCount();
    this->matched = new int[N + 1];
    std::fill_n(this->matched,N+1, -1);
//    for (int i = 0 ; i < N + 1; i++) this->matched[i] = -1;
//    currentMatrix = this->graph->getConjugateMatrix();
    this->getMatrix();
    this->cyclePaths.clear();
    this->cyclePaths.shrink_to_fit();
    this->mergedPaths.clear();
}

bool matching::kmDfs(int u, bool visity[],bool visitx[], std::set<int>* pre, float ex[], float ey[], float slack[]) {
//    auto matrix = getMatrix();
    visity[seqGraph::conjugateIdx(u)] = true;
    visitx[u] = true;
    std::vector<int> candidates;
    for (int i = 1; i < N+1; i++) {
        if(i == u) continue;
//        auto t1 = ex[u];
//        auto t2 = ey[i];
//        auto t3 = matrix[u][i];
//        auto t4 = visity[i];
//        auto t5 = matched[i];
//        auto t6 = slack[i];
//        auto p = std::find(pre->begin(), pre->end(), i) != pre->end();
        if (visity[i]) continue;
//        if ( u == 19 && i == 1615) {
//            std::cout<<matched[i];
//            printM(u);
//            std::cout<<std::endl;
//             int km = 0;
//            std::cout<<std::endl;
//        }
        auto mIJ = this->getIJ(u,i);
        if(std::abs(ex[u] + ey[i] - mIJ) <= ZERO) {
            visity[i] = true;
            candidates.push_back(i);
            if((matched[i] == -1) || (kmDfs(matched[i], visity, visitx, pre, ex, ey, slack))) {
                pre->insert(u);
                pre->insert(seqGraph::conjugateIdx(i));
                matched[i] = u;
                matched[seqGraph::conjugateIdx(u)] = seqGraph::conjugateIdx(i);
                return true;
            }
        } else if(slack[i] > ex[u] + ey[i] - mIJ) {
//            auto tt1 = ex[u];
//            auto tt2 = ey[i];
//            auto tt3 = matrix[u][i];
//            if (ex[u] + ey[i]  - matrix[u][i] >= ZERO)
                slack[i] = ex[u] + ey[i]  - mIJ;
//            auto tts = slack[i];
//            auto ttttt = 99;
        }
    }
//    sort(candidates.begin(),candidates.end(), this->cmpVertex());
//    Quicksort(this, this->idx2VertexInCurrentGraph(u)->getOutCoverage(), candidates, 0, candidates.size() - 1);
//    for (auto i : candidates) {
//        if((matched[i] == -1) || (kmDfs(matched[i], visity, visitx, pre, ex, ey, slack))) {
//            pre->insert(u);
//            pre->insert(seqGraph::conjugateIdx(i));
//            matched[i] = u;
//            matched[seqGraph::conjugateIdx(u)] = seqGraph::conjugateIdx(i);
//            return true;
//        }
//    }

    return false;
}


bool matching::dfs(int u, bool visity[], std::vector<int>* pre) {
//    auto matrix = getMatrix();
    visity[seqGraph::conjugateIdx(u)] = true;
    for (int i = 1; i < N+1; i++) {
//        auto p = std::find(pre->begin(), pre->end(), i) != pre->end();
        if(this->getIJ(u,i) != 0 && !visity[i]) {
            visity[i] = true;
            if((matched[i] == 0) || (dfs(matched[i], visity, pre))) {
                pre->push_back(u);
                pre->push_back(seqGraph::conjugateIdx(i));
                matched[i] = u;
                matched[seqGraph::conjugateIdx(u)] = seqGraph::conjugateIdx(i);
                return true;
            }
        }
    }
    return false;
}

void matching::bfs(int u, float ex[], float ey[], bool visity[], int pre[], std::set<int>& skipped,float slack[]) {
//    auto matrix = getMatrix();
    int x,cY,y=0, yy=0;
    float d = 0;
    std::fill_n(slack, N+1, INF);
    matched[y] = u;

    while (true){
        x = matched[y], d = INF, visity[y] = true;
        visity[seqGraph::conjugateIdx(x)] = true;
//        visity[x] = true;
        for(int i = 1; i < N + 1; i++){
//            if i is visited or i == x or (A1 -- B1 is in match already ,A2 -- B2, and model == 2)
            if(visity[i] or i == x or this->vertexLookup(x,i)) continue;
//            if(visity[i] or i == x) continue;
            auto mIJ = this->getIJ(x,i);
            if(slack[i] >= ex[x] + ey[i] - mIJ){
                slack[i] = ex[x] + ey[i] - mIJ;
                pre[i] = y; //表示 y对应的 点 y 需要减小的权值
            }
            if(slack[i] <= d){
                d = slack[i],yy = i; //找出减少最小的那条边
                //            如果没有matched，并且是zero就break
                if ((matched[i] == -1 || getIJ(matched[i],i) == 0) && d < ZERO)
                    break;
            }
        }
        if (d < ZERO) {
            d = 0;
        }
        if (d != 0) {
            for(int i = 0; i < N+1; i++){
                if(visity[i]) ex[matched[i]] -= d,ey[i] += d;
                else slack[i] -= d;
            }
        }

        y = yy;

        if(matched[y] == -1 || getIJ(matched[y],y) == 0) {
            break;
        }
    }
    while(y){
        matched[y] = matched[pre[y]];
        matched[seqGraph::conjugateIdx(matched[pre[y]])] = seqGraph::conjugateIdx(y);
        skipped.emplace(matched[pre[y]]);
        if (getIJ(matched[y], y) != 0)
            skipped.emplace(seqGraph::conjugateIdx(y));
        y = pre[y]; // 更新每个点对应的点
    }
}


//void matching::bfs(int u, float * ex, float* ey, bool* visity, int* pre, float* slack) {
//    auto matrix = this->graph->getConjugateMatrix();
//    int x,cY,y=0, yy=0;
//    float delta;
//
//    for (int i = 1; i < N + 1; i++) slack[i] = INF;
//
////    当前的match以及共轭的match
//    matched[y] = u;
////    matched[seqGraph::conjugateIdx(u)] = y;
//    while (matched[y]) {
//        x = matched[y], delta = INF;
////        cY = matched[seqGraph::conjugateIdx(u)];
////      当前点y和x共轭点同时被visit
//        visity[y] = true;
////        visity[seqGraph::conjugateIdx(x)] = true;
//        for (int i = 1; i < N + 1; i++) {
//            if(!visity[i]) { // visty and its conjugate will always sync.
//                if(slack[i] > ex[x] + ey[i] - matrix[x][i]){
//                    slack[i] = ex[x] + ey[i] - matrix[x][i];
//                    pre[i] = y;
//                }
//                if(slack[i] < delta) {
//                    delta = slack[i], yy = i;
//                }
//            }
//        }
//        for(int i = 0; i < N + 1; i++){
//            if(visity[i]) ex[matched[i]] -= delta,ey[i] += delta;
//            else slack[i] -= delta;
//        }
//        y = yy;
//    }
//    while(y){
//        matched[y] = matched[pre[y]];
////        auto t1 = seqGraph::conjugateIdx(matched[pre[y]]);
////        auto t2 = seqGraph::conjugateIdx(y);
////        matched[t1] = seqGraph::conjugateIdx(y);
//        y = pre[y];
//    }
//}

void matching::hungarian() {
//    auto matrix = getMatrix();
    bool visity[N + 1];
    bool visitx[N + 1];
    float slack[N + 1], ex[N + 1], ey[N + 1];

    std::fill_n(ex,N+1,0);
    std::fill_n(ey,N+1,0);
//    std::fill_n(ex,N+1,0);

//    memset(ey,0,sizeof ey);
//    memset(ex,0,sizeof ex);
//    memset(visity,0,sizeof visity);

// TODO, get row max
    for( int i = 1 ; i < N +1 ; i++ ){
        ex[i] = this->getIRowMax(i);
//        for( int j = 1 ; j <= N ; j++ ){
//            if( ex[i] < matrix[i][j] ) ex[i] = matrix[i][j];
//        }
    }
    auto* pre = new std::set<int>();
    if(VERBOSE >= 2) {
        std::cout<<"         ";
        for (int i = 1; i < N+1;i++) {
            std::cout << this->idx2StrDir(i) << "\t";
        }
        std::cout<<std::endl;
    }
    for (int i = 1; i < N + 1; i++) {
        std::fill_n(slack, N+1, 10000000);
//        for( int l = 1 ; l <=N ; l++ ) slack[l] = 1000;
        if (pre->find(i) != pre->end()) {
            continue;
        }
        while (true) {
            memset(visity, 0, sizeof visity);
            memset(visitx, 0, sizeof visitx);
            auto r = kmDfs(i, visity,visitx, pre,ex, ey, slack);
            if (r) {
                break;
            } else {
                float delta = INF;
                for( int j = 1 ; j <= N ; j++ ){
                    if( !visity[j] && delta > slack[j] ){
                        delta = slack[j];
                    }
                }

                for( int k = 1 ; k <= N; k++ ){
                    if( visitx[k] ) ex[k] -= delta;
                    if( visity[k] ) ey[k] += delta;
                    else slack[k] -= delta;
                }
            }
        }
        if (VERBOSE >= 2) {
            std::cout << i << " " << this->idx2StrDir(i) << "\t|";
            for(int k = 1; k < N + 1; k++) {
                if(this->matched[k] != -1)
                    std::cout << this->idx2StrDir(matched[k]) << "\t";
                else
                    std::cout<< this->matched[k]<<"\t";

            }
            std::cout<<std::endl;
        }
        if (VERBOSE >= 1) {
            std::cout << i << " " << this->idx2StrDir(i) <<std::endl;
        }
    }
    free(pre);
}
void matching::main_steps() {
//    auto matrix = getMatrix();
    bool visity[N + 1];
//    bool visitx[N + 1];
    float slack[N + 1], ex[N + 1], ey[N + 1];
    int pre[N + 1];
    std::fill_n(ex,N+1,0);
    std::fill_n(ey,N+1,0);
//    std::fill_n(ex,N+1,0);

//    memset(ey,0,sizeof ey);
//    memset(ex,0,sizeof ex);
//    memset(visity,0,sizeof visity);

//TODO get row max
    for( int i = 1 ; i < N +1 ; i++ ){
        ex[i] = this->getIRowMax(i);
//        for( int j = 1 ; j <= N ; j++ ){
//            if( ex[i] < matrix[i][j] ) ex[i] = matrix[i][j];
//        }
    }
    std::set<int> skipped;
    for( int i = 1 ; i < N+1 ; i++ ){
        if (i == 19) {
            auto tmpp = 33;
        }
        if (skipped.find(i) != skipped.end()) continue;
//        if (ex[i] == 0){
//            matched[i] = i;
//            continue;
//        }
        memset( visity , false , sizeof(visity) );
        std::fill_n(pre, N+1, 0);
        std::fill_n(slack, N+1, 10000000);

        bfs(i,ex,ey, visity,pre, skipped, slack);
        if (VERBOSE >= 1) {
            std::cout << i << " " << this->idx2StrDir(i) <<std::endl;
        }
    }

    for(int k = 1; k < N + 1; k++) {
        if(this->matched[k] != -1)
            std::cout << this->idx2StrDir(matched[k]) << "\t";
        else
            std::cout<< this->matched[k]<<"\t";

    }
    std::cout<<std::endl;
}

int  matching::inDegree(int idx) {
    int now = idx;
    int vIdx = (now + 1) / 2;
    char dir = now % 2 == 0 ? '-':'+';
    if ((*this->originalGraph->getVertices())[vIdx - 1]->getId() == "NZ_CP028684.1:62109-92725_0") {
        auto tmp = 44;
    }
//    int edge_count = (*this->originalGraph->getVertices())[vIdx - 1]->getPrevJuncCountIgnoreCopy();

    if (dir == '+')
        return (*this->originalGraph->getVertices())[vIdx - 1]->getPrevJuncCountIgnoreCopyPositive();
    else return  (*this->originalGraph->getVertices())[vIdx - 1]->getPrevJuncCountIgnoreCopyNegative();
//    for (auto item : *(*this->originalGraph->getVertices())[vIdx - 1]->getEp3()->getInEdges()) {
//        if (item.)
//    }
//    if (edge_count >1 ) {
//        auto tmp = 22;
//    }
//    return (*this->originalGraph->getVertices())[vIdx - 1]->getEp3()->getInEdges()->size();
}

int  matching::outDegree(int idx) {
    int now = idx;
    int vIdx = (now + 1) / 2;
    char dir = now % 2 == 0 ? '-':'+';
//    return (*this->originalGraph->getVertices())[vIdx - 1]->getEp5()->getOutEdges()->size();
    if ((*this->originalGraph->getVertices())[vIdx - 1]->getId() == "NZ_CP028684.1:62109-92725_0") {
        auto tmp = 44;
    }
//    return (*this->originalGraph->getVertices())[vIdx - 1]->getNextJuncCountIgnoreCopy();
    if (dir == '+')
        return (*this->originalGraph->getVertices())[vIdx - 1]->getNextJuncCountIgnoreCopyPositive();
    else return  (*this->originalGraph->getVertices())[vIdx - 1]->getNextJuncCountIgnoreCopyNegative();
}

std::string matching::idx2StrDir(int idx, const std::string& token) {
    int now = idx;
    int vIdx = (now + 1) / 2;
    char dir = now % 2 == 0 ? '-':'+';
    auto idStr = (*this->originalGraph->getVertices())[vIdx - 1]->getId();
    auto len = idStr.size();
    auto pos = idStr.find_last_of('_');
//    if (idStr[len - 2] == '_') {
//        idStr.pop_back();
//        idStr.pop_back();
//    }
    return idStr.substr(0,pos)+token+dir;
}

std::string matching::idx2Str(int idx) {
    int now = idx;
    int vIdx = (now + 1) / 2;
    auto idStr = (*this->originalGraph->getVertices())[vIdx - 1]->getId();
    return idStr;
}

seqGraph::Vertex *matching::idx2VertexInCurrentGraph(int idx) {
    int vIdx = (idx + 1) / 2;
    auto idStr = (*this->graph->getVertices())[vIdx - 1]->getId();
    auto v = this->graph->getVertexByIdQ(idStr);
    return v;
}
seqGraph::Vertex *matching::idx2VertexInOriginalGraph(int idx) {
//    int vIdx = (idx + 1) / 2;
    auto idStr = idx2Str(idx);
    auto v = this->graph->getVertexByIdQ(idStr);
    return v;
}
std::map<int, std::vector<int>*>* matching::resolvePath(std::map<int, std::vector<int>*>* prevPaths) {
//    auto matrix = graph->getConjugateMatrix();
    bool visited[N+1];
    memset(visited, 0, sizeof visited);
    auto* resolvedPath = new std::map<int, std::vector<int> *>();
    auto* resPath = new std::map<int, std::vector<int> *>();
//    auto checkedC = checkConjugateMatch();
//    if(checkedC!= 0)
//        std::cout<<"Conjugate not checked\n";
//    else
//        std::cout<<"check conjugate done\n";
    std::vector<int> cyclePaths;
    for(int i = 1; i < N + 1; i++) {
        if(i == 93) {
            int m = 88;
        }
        if (visited[i]) continue;
        auto* currentPath = new std::vector<int>();

//        if (matrix[matched[i]][i] == 0){
//            if (matrix[matched[seqGraph::conjugateIdx(i)]][seqGraph::conjugateIdx(i)] == 0) {
//                currentPath->push_back(i); resolvedPath->emplace(i,currentPath);
//                visited[i] = true;
//                visited[seqGraph::conjugateIdx(i)] = true;
//            } else {
//                visited[i] = true;
//            }
//            continue;
////            currentPath->push_back(i); resolvedPath->emplace(i,currentPath); continue;
//        }
        int now = i;
        int vIdx = (now + 1) / 2;
        char dir = now % 2 == 0 ? '-':'+';
//        currentPath->push_back((*this->graph->getVertices())[vIdx - 1]->getId() + dir);
        currentPath->push_back(now);
//        if ((*this->graph->getVertices())[vIdx - 1]->getId() == "52_0")
//            int tm = 99;

//        std::cout<< (*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
        visited[now] = true;
//        a+ b+ ..., a- ..., can merge two
        if (this->getIJ(matched[seqGraph::conjugateIdx(now)], seqGraph::conjugateIdx(now)) == 0) {
            visited[seqGraph::conjugateIdx(now)] = true;
        }
//        visited[seqGraph::conjugateIdx(now)] = true;
        bool currentInsert = true;
        bool isCycle = false;
        std::deque<int> zeroBreakPoint;
        int timesi = 0;
        while (true) {
            timesi++;
            if (timesi >= 10* this->getN()) {
                break;
            }
//            auto tmp = this->idx2StrDir(now);
//            auto tmp2 = this->idx2StrDir(matched[now]);
            auto mIJ = this->getIJ(matched[now],now);
            if (visited[matched[now]]) {
                if (mIJ == 0) {
                    zeroBreakPoint.push_back(currentPath->front());
                } else {
                if (matched[now] != currentPath->front()) {
                    zeroBreakPoint.push_back(currentPath->front());
                    currentPath->push_back(matched[now]);
                }
//                    zeroBreakPoint.push_back(matched[now]);
                }
                break;
            }
//            如果连上之前断裂的路径
//            if (resolvedPath->find(matched[now]) != resolvedPath->end()) {
//                auto oldPath = (*resolvedPath)[matched[now]];
////                +"_"+ std::to_string(now)
//                for (auto it = currentPath->rbegin(); it != currentPath->rend(); it++){
//                    oldPath->insert(oldPath->begin(), *it);
//                }
//                resolvedPath->erase(matched[now]);
////                auto extendPath = this->addPrevPath(prevPaths, oldPath);
//                resolvedPath->emplace(i, oldPath);
//                currentInsert = false;
//                break;
//            }
//            环状路径
//            if(visited[matched[now]]){
//                isCycle = true;
////                if (matched[now] == i) currentPath->push_back(-1);
//                currentPath = breakResolvedPaths(currentPath);
//                break;
//            }
//            visited[matched[now]] = true;
//            visited[seqGraph::conjugateIdx(matched[now])] = true;
//            vIdx = (matched[now] + 1) / 2;
//            dir = matched[now] % 2 == 0 ? '-':'+';
//            std::cout<<(*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
            currentPath->push_back(matched[now]);
            if (mIJ == 0) {
                zeroBreakPoint.push_back(matched[now]);
            }
            now = matched[now];
            visited[now] = true;
            visited[seqGraph::conjugateIdx(now)] = true;
        }
        if(currentInsert) {
            //  seems bug here, a cycle path want to connect with a cycle path 126 127 , 125 128, fix with add it conjugate into visited
            if (zeroBreakPoint.empty() || zeroBreakPoint.back() != currentPath->front()) {
                visited[seqGraph::conjugateIdx(currentPath->front())] = true;
            }
//            break zero edge and connect 1 2 3 4 and 4 5 6 7
            breakResolvedPaths(currentPath, zeroBreakPoint, resolvedPath);
            zeroBreakPoint.clear();
//            resolvedPath->emplace((*currentPath)[0], currentPath);
//            if (isCycle) {
//                this->cyclePaths.push_back((*currentPath)[0]);
//            }
        }
    }
    if (BREAK_C) {
        breakAndMergeCycle(resolvedPath);
    }
//    auto tmp = this->idx2VertexInCurrentGraph(this->cyclePaths.front());
    for (auto path : *resolvedPath) {
        auto extendPath = this->addPrevPath(prevPaths, path.second);
        resPath->emplace(extendPath->front(),extendPath);
//        path.second->clear();
//        path.second->shrink_to_fit();
    }
    free(resolvedPath);
//    if(prevPaths != nullptr) {
//        for(auto item: *prevPaths) {
//            item.second->clear();
//            item.second->shrink_to_fit();
//        }
//    }
    free(prevPaths);

    free(matched);
    return resPath;
}

bool matching::vertexLookup(int i, int j)  {
    auto isMatched = false;
    auto iVertex = this->idx2VertexInCurrentGraph(i);
    auto jVertex = this->idx2VertexInCurrentGraph(j);

    int jStart = j - (jVertex->getCopyIdx()) * 2;
    int jEnd = j + (jVertex->getWeight()->getCopyNum() - jVertex->getCopyIdx()) * 2;

    int iStart = i - (iVertex->getCopyIdx())  * 2;
    int iEnd = i + (iVertex->getWeight()->getCopyNum() - 1 - iVertex->getCopyIdx()) * 2;
    //    check if i have matched all path already
    auto maxMatched = iVertex->getPrevJuncCountIgnoreCopyNegative() + iVertex->getPrevJuncCountIgnoreCopyPositive();
    int matchedCount = 0;
    bool isMaxMatched = false;
    for (int p = iStart; p <= iEnd; p+=2) {
        if (matched[seqGraph::conjugateIdx(p)] != -1 && this->getIJ(matched[seqGraph::conjugateIdx(p)] ,seqGraph::conjugateIdx(p)) != 0){
            matchedCount++;
            if (matchedCount >= maxMatched) {
                isMaxMatched = true;
                break;
            }
        }
    }
//    if max matched, can use this
    if (isMaxMatched) {
        return false;
    }
//  else  check if have matched[j1] = i1
    for (int p = jStart; p <= jEnd; p+=2) {
        for (int q = iStart; q <= iEnd; q+=2) {
            if (this->matched[p] == q) {
                isMatched = true;
                break;
            }
        }
    }

    if (isMatched) {
        return true;
    }
    return false;
}

float* matching::mergePath(std::vector<int>* p1, std::vector<int>* p2, float* result) {
//    auto result = new float[4];
//    for(int i = 0 ; i< 4 ; i ++) result[i] = 0;
//            TODO, only add end node
    bool isCycle_1 = this->isCycle(p1->front());
    bool isCycle_2 = this->isCycle(p2->front());
    auto front_1 = p1->front();
    auto front_2 = p2->front();
    auto back_1 = p1->back();
    auto back_2 = p2->back();
    result[0] = this->getIJFromOrigG(front_2,back_1);
    result[1] = this->getIJFromOrigG(seqGraph::conjugateIdx(back_2),back_1);
    result[2] = this->getIJFromOrigG(front_2, seqGraph::conjugateIdx(front_1));
    result[3] = this->getIJFromOrigG(seqGraph::conjugateIdx(back_2),seqGraph::conjugateIdx(front_1));
    float max_r = result[0];
    if (max_r < result[1]) {
        max_r = result[1];
    }
    if (max_r < result[2]) {
        max_r = result[2];
    }
    if (max_r < result[3]) {
        max_r = result[3];
    }

    if (isCycle_1) {
        front_1 = seqGraph::conjugateIdx(front_1);
        back_1 = seqGraph::conjugateIdx(back_1);
//        front_2 = seqGraph::conjugateIdx(front_2);
        auto tmp0 = this->getIJFromOrigG(front_2,back_1);
        auto tmp1 = this->getIJFromOrigG(seqGraph::conjugateIdx(back_2),back_1);
        auto tmp2 = this->getIJFromOrigG(front_2, seqGraph::conjugateIdx(front_1));
        auto tmp3 = this->getIJFromOrigG(seqGraph::conjugateIdx(back_2),seqGraph::conjugateIdx(front_1));

        if (tmp0 > max_r || tmp1 > max_r || tmp2 > max_r || tmp3 > max_r) {
            result[0] = tmp0;
            result[1] = tmp1;
            result[2] = tmp2;
            result[3] = tmp3;
        }
    }

    if (isCycle_2) {
        front_2 = seqGraph::conjugateIdx(front_2);
        back_2 = seqGraph::conjugateIdx(back_2);
//        front_2 = seqGraph::conjugateIdx(front_2);
        auto tmp0 = this->getIJFromOrigG(front_2,back_1);
        auto tmp1 = this->getIJFromOrigG(seqGraph::conjugateIdx(back_2),back_1);
        auto tmp2 = this->getIJFromOrigG(front_2, seqGraph::conjugateIdx(front_1));
        auto tmp3 = this->getIJFromOrigG(seqGraph::conjugateIdx(back_2),seqGraph::conjugateIdx(front_1));

        if (tmp0 > max_r || tmp1 > max_r || tmp2 > max_r || tmp3 > max_r) {
            result[0] = tmp0;
            result[1] = tmp1;
            result[2] = tmp2;
            result[3] = tmp3;
        }
    }
//    int max_p1 = 5 > p1->size() ? p1->size():5;
//    int max_p2 = 5 > p2->size() ? p2->size():5;
//
//    for (int i = 0 ; i < max_p1; i ++) {
//        int ip1 = (*p1)[i];
//        if (ip1 == -1) continue;
//
//        for (int j = 0 ; j < max_p2; j ++) {
//            int ip2 = (*p2)[j];
//            if (ip2 == -1) continue;
////            result[0] += matrix[ip2][ip1];
////            result[1] += matrix[seqGraph::conjugateIdx(ip2)][ip1];
////            result[2] += matrix[ip2][seqGraph::conjugateIdx(ip1)];
////            result[3] += matrix[seqGraph::conjugateIdx(ip2)][seqGraph::conjugateIdx(ip1)];
//            if (ip1 % 2 == 1) {
//                if (ip2 % 2 == 1) {//            ++
//                    result[0] += this->getIJ(ip1,ip2);
//                } else {
//                    result[1] += this->getIJ(seqGraph::conjugateIdx(ip2), ip1);
//                }
//            } else {
//                if (ip2 % 2 == 1) {//            -+
//                    result[2] += this->getIJ(ip2, seqGraph::conjugateIdx(ip1));
//                } else {
//                    result[3] += this->getIJ(seqGraph::conjugateIdx(ip2), seqGraph::conjugateIdx(ip1));
//                }
//            }
//        }
//    }
    return result;
}

void matching::reconstructMatrix(std::map<int, std::vector<int>*>* paths, seqGraph::Graph* originGraph) {
    auto* resultG = new seqGraph::Graph();
    resultG->isReconstructed = true;
//    auto tm = resultG->getConjugateMatrix() == nullptr;
    auto* values = new float[4];
    for (auto iPath: *paths) {
        if (this->isCycle(iPath.first) && !BREAK_C) continue;
//        if (this->isCycle((*iPath.second)[0])) continue;
        resultG->addVertex(std::to_string(iPath.second->front()),"xx",1,2,1,1,2,0);
    }
    seqGraph::Vertex* v1;
    seqGraph::Vertex* v2;
//    for(auto junc : *originalJunctions) {
//        int i = junc->getSource()->getIdx();
//        int j = junc->getTarget()->getIdx();
//        int sDir = junc->getSourceDir();
//        int tDir = junc->getTargetDir();
//        float weightValue = junc->getWeight()->getCopyNum();
//        if (sDir == '+') {
//            if (tDir == '+') {
//                this->ConjugateMatrix[2*j + 1][2*i + 1] = weightValue;
//                this->ConjugateMatrix[2*(i+1)][2*(j+1)] = weightValue;
//            } else {
//                this->ConjugateMatrix[2*(j + 1)][2*i+1] = weightValue;
//                this->ConjugateMatrix[2*(i+1)][2*j+1] = weightValue;
//            }
//        } else {
//            if (tDir == '+') {
//                this->ConjugateMatrix[2*j+1][2*(i + 1)] = weightValue;
//                this->ConjugateMatrix[2*i+1][2*(j+1)] = weightValue;
//            } else {
//                this->ConjugateMatrix[2*i+1][2*j+1] = weightValue;
//                this->ConjugateMatrix[2*(j+1)][2*(i+1)] = weightValue;
//            }
//        }
//    }
//    auto matrix = originGraph->getConjugateMatrix();
    for (auto iPath: *paths) {
        if (!BREAK_C){
            if (this->isCycle(iPath.first)) continue;
        }
        v1 = resultG->getVertexByIdQ(std::to_string(iPath.second->front()));
        if (idx2Str(iPath.second->front()) == "EDGE_5532997_length_1440_cov_2.471753_0")
            int iii = 9;
        for (auto jPath: *paths) {
            if (!BREAK_C){
                if (this->isCycle(jPath.first)) continue;
            }
            if (idx2Str(jPath.second->front()) == "EDGE_62388_length_682_cov_5.543802_0")
                int idi = 9;
            if (iPath.first == jPath.first) continue;
            v2 = resultG->getVertexByIdQ(std::to_string(jPath.second->front()));
            mergePath(iPath.second, jPath.second, values);
//            if (iPath.second->size() == 1 && jPath.second->size() == 1) {
//                int i = iPath.second->front();
//                int j = jPath.second->front();
//                values[0] = this->getIJFromOrigG(j,i);
//                values[1] = this->getIJFromOrigG(seqGraph::conjugateIdx(j), i);
//                values[2] = this->getIJFromOrigG(j, seqGraph::conjugateIdx(i));
//                values[3] = this->getIJFromOrigG(seqGraph::conjugateIdx(j), seqGraph::conjugateIdx(i));
//            } else {
//                mergePath(iPath.second, jPath.second, values);
//            }
//            if (values[0] == 0 && values[1]==0 && values[2]==0 && values[3]==0) continue;
//            std::string v1Str, v2Str;
//            for (auto item : *iPath.second) {
//                v1Str+=idx2StrDir(item);
//            }
//            for (auto item : *jPath.second) {
//                v2Str+= idx2StrDir(item);
//            }
            if (values[0] != 0) {
                resultG->addJunction(v1, v2, '+', '+', values[0], 1 , 1);
            }

            if (values[1] != 0) {
                resultG->addJunction(v1, v2, '+', '-', values[1], 1 , 1);
            }
            if (values[2] != 0) {
                resultG->addJunction(v1, v2, '-', '+', values[2], 1 , 1);
            }
            if (values[3] != 0) {
                resultG->addJunction(v1, v2, '-', '-', values[3], 1 , 1);
            }
        }
    }
    delete[] values;
    std::cout<<"start reset"<<std::endl;
    resetGraph(resultG);
}

int matching::checkConjugateMatch() {
    int r = 0;
//    auto matrix = this->getMatrix();
    for (int i = 1; i < N+1; i++) {
        auto left = matched[i];
        if (this->getIJ(left, i) == 0) continue;
        auto conjugateI = seqGraph::conjugateIdx(i);
        auto conjugateLeft = seqGraph::conjugateIdx(left);
        if (matched[conjugateLeft] != conjugateI) {
            r++;
            std::cout<<i<<'\t'<<matched[i]<<'\t'<<conjugateLeft<<'\t'<<conjugateI<<'\t'<<matched[conjugateLeft]<<'\n';
        }
    }
    std::cout<<std::endl;
    return r;
}

void matching::checkConjugateMatrix() {
    for(int i = 1; i < N+1; i++) {
        for(int j = 1; j < N+1; j++) {
            if(i == 5 && j == 1260)
                auto oo = 33;
            auto t1 = this->getIJ(i,j);
            auto t2 = this->getIJ(seqGraph::conjugateIdx(j), seqGraph::conjugateIdx(i));
            if (std::abs(t1 - t2) > ZERO) {
                std::cout<<"Error "<<i<<" "<<j<<std::endl;
                exit(0);
            }
//            std::cout<<matrix[i][j]<<" ";
        }
    }
}

float matching::getIJ(int i, int j) {
    if (this->graph->isSparse())
        return this->getMatrix().getIJ(i,j);
    else {
        auto iIdx  = (i + 1) / 2;
        auto jIdx = (j + 1) / 2;
        auto iDir = i % 2 == 0 ? '-' : '+';
        auto jDir = j % 2 == 0 ? '-' : '+';
//    get junction j to i
        return this->graph->getIJ(jIdx - 1, iIdx - 1, jDir, iDir);
    }
}

float matching::getIJFromOrigG(int i, int j) {
    if (this->originalGraph->isSparse())
        return this->originalGraph->getConjugateMatrix().getIJ(i,j);
    else {
        auto iIdx  = (i + 1) / 2;
        auto jIdx = (j + 1) / 2;
        auto iDir = i % 2 == 0 ? '-' : '+';
        auto jDir = j % 2 == 0 ? '-' : '+';
//    get junction j to i
        return this->originalGraph->getIJ(jIdx - 1, iIdx - 1, jDir, iDir);
    }
}


std::vector<int>* matching::breakCycle(std::vector<int> * cyclePath) {
    auto res = new std::vector<int>();
//    auto matrix = this->getMatrix();
    int minI = cyclePath->size() - 1;
    float minW = INF;
    for (int i = 0; i < cyclePath->size(); i++) {
        int idx = (*cyclePath)[i];
        int nextIdx = (*cyclePath)[i+1];
        if (i == cyclePath->size() - 1) {
            if(this->getIJ(cyclePath->front(),idx) < minW) {
                minI = cyclePath->size() - 1;
            }
            continue;
        }
        auto mIJ = this->getIJ(nextIdx, idx);
        if (mIJ < minW) {
            minW = mIJ;
            minI = i;
        }
    }
    if (minI == cyclePath->size() - 1) return cyclePath;
    for (int i = minI + 1; i < cyclePath->size(); i++) {
        res->push_back((*cyclePath)[i]);
    }
    for (int i = 0; i <= minI; ++i) {
        res->push_back((*cyclePath)[i]);
    }
    return res;
}

void matching::breakAndMergeCycle(std::map<int,std::vector<int>*> *result) {
    if (this->cyclePaths.empty()) return;
    std::vector<std::pair<int, std::vector<int>*>> A;
    sort(*result,A);
    for (auto item: this->cyclePaths) {
        if (item == 125) {
            auto tmp = 33;
        }
//        Here we only consider the break point of the cycle path, have copy (a bug here may be, we don't consider the copied node at other place)
        auto v1 = idx2VertexInCurrentGraph(item);
        auto cPath = (*result)[item];
        for (auto& p : *result) {
//           avoid A merged into B and B merged into A again
            if (p.first == item or (this->mergedPaths.find(p.first) != this->mergedPaths.end())) continue;
            bool hit = false;
            for (int i = 0 ; i < p.second->size(); i++) {
                auto idx = (*p.second)[i];
                auto v2 = idx2VertexInCurrentGraph(idx);
                if(v1->sameVertex(*v2)) {
//                    merge a cycle into this path.
                    this->mergedPaths.insert({item, p.first});
                    auto pos = p.second->begin() + i;
//                    keep the original first to be first
                    p.second->insert(pos, cPath->begin(), cPath->end());
                    hit = true;
//                    process A1xxxxx insert to A0xxxxx, this will make the front diff from the key.
                    if (i == 0) {
                        auto cPathLen = cPath->size();
                        int swapV = p.second->front();
                        (*p.second)[0] = (*p.second)[cPathLen];
                        (*p.second)[cPathLen] = swapV;
                    }
                    break;
                }
            }
            if (hit) break;
        }
    }
}

void matching::breakResolvedPaths(std::vector<int> *cur, std::deque<int> & zereBK, std::map<int,std::vector<int>*> *result) {
//    环状路,找到有拷贝数的地方断开
    if (zereBK.empty()) {
        bool isNoCopyCycle = true;
        auto cSize = cur->size();
        for (int i = 0; i < cSize;i++){
            if (this->idx2VertexInOriginalGraph((*cur)[0])->getWeight()->getCopyNum() >= 2) {
                isNoCopyCycle = false;
//                this->cyclePaths.push_back((*cur)[0]);
//                result->emplace((*cur)[0], cur);
                break;
            }
            cur->push_back(cur->front());
            cur->erase(cur->begin());
        }
//        如果环中没有任何一个点拷贝数大于1那就是一个单纯环，无法merge到任何上面，直接再权重最小处断开
        if(BREAK_C && isNoCopyCycle) {
            auto tmp = breakCycle(cur);
            this->cyclePaths.push_back((*tmp)[0]);
            result->emplace((*tmp)[0],tmp);
//            cur->clear();
//            cur->shrink_to_fit();
        } else {
            this->cyclePaths.push_back((*cur)[0]);
            result->emplace((*cur)[0], cur);
        }
        return;
    }

//    如果最后断开
    auto lastCfirst = false;
    if (zereBK.back() == cur->front()) zereBK.pop_back();
    else{ //如果最后一个链接第一个
        lastCfirst = true;
    }
    auto * tmp = new std::vector<int>();
//    check the first break path, situation: 1- 2+ ... can reversed merged with 1+ 3+...
    auto isTheFirstPath = true;
    for(auto &item : *cur) {
        if (!zereBK.empty() && item == zereBK.front()) {
//          1- 2+ ... can reversed merged with 1+ 3+...
            if (isTheFirstPath && !tmp->empty()) {
                if(!lastCfirst) {
                    // connect 1 2 3 4 and 4 5 6 7
                    for (auto oldPath: *result) {
//                        the cycle path will merge at next method
                        if(this->isCycle(oldPath.first)) {
                            continue;
                        }
                        if (oldPath.first == seqGraph::conjugateIdx(tmp->front())) {
                            seqGraph::conjugatePath(tmp);
                            tmp->pop_back();
                            for (auto i: *oldPath.second) {
                                tmp->push_back(i);
                            }
                            result->erase(oldPath.first);
                            break;
                        }
                    }
                }
                isTheFirstPath = false;
            }
            //            if front v same back v
            if (tmp->size() > 1 && this->idx2VertexInCurrentGraph(tmp->front())->sameVertex(*this->idx2VertexInCurrentGraph(
                    tmp->back()))) {
                this->cyclePaths.push_back((*tmp)[0]);
                tmp->pop_back();
            }
            result->emplace((*tmp)[0], tmp);
            zereBK.pop_front();
            tmp = new std::vector<int>();
            tmp->push_back(item);
        } else {
            tmp->push_back(item);
        }
    }
    if (!tmp->empty()) {
        if(!lastCfirst) {
//  if the path only break at the last,  1- 2+ ... view not reversed merged with 1+ 3+. we need check if this path might be a "middle path",
//  which means can connect two other path. xxx curr_path xxx
// connect 1 2 3 4 and 4 5 6 7
            for (auto item : *result) {
                if(this->isCycle(item.first)) {
                    continue;
                }
                if (item.first == tmp->back()) {
                    tmp->pop_back();
                    for(auto i : *item.second){
                        tmp->push_back(i);
                    }
                    result->erase(item.first);
                    break;
                }
            }
//  1- 2+ ... can reversed merged with 1+ 3+...
            if (isTheFirstPath) {
                if(!lastCfirst) {
                    // connect 1 2 3 4 and 4 5 6 7
                    for (auto oldPath: *result) {
                        if(this->isCycle(oldPath.first)) {
                            continue;
                        }
                        if (oldPath.first == seqGraph::conjugateIdx(tmp->front())) {
                            seqGraph::conjugatePath(tmp);
                            tmp->pop_back();
                            for (auto i: *oldPath.second) {
                                tmp->push_back(i);
                            }
                            result->erase(oldPath.first);
                            break;
                        }
                    }
                }
                isTheFirstPath = false;
            }
            //            if front v same back v
            if (tmp->size() > 1 && this->idx2VertexInCurrentGraph(tmp->front())->sameVertex(*this->idx2VertexInCurrentGraph(
                    tmp->back()))) {
                this->cyclePaths.push_back((*tmp)[0]);
                tmp->pop_back();
            }

            result->emplace((*tmp)[0], tmp);
        }
        else{
            auto oldPath = (*result)[cur->front()];
            for (auto item : *oldPath) {
                tmp->push_back(item);
            }
            result->erase(cur->front());
            //            if front v same back v
            if (tmp->size() > 1 && this->idx2VertexInCurrentGraph(tmp->front())->sameVertex(*this->idx2VertexInCurrentGraph(
                    tmp->back()))) {
                this->cyclePaths.push_back((*tmp)[0]);
                tmp->pop_back();
            }
            result->emplace((*tmp)[0], tmp);
        }
    }
}

// This function will expand the curPath with prevPath, at the same time, it will modify this.cyclePaths and this.mergedPath corresponding
std::vector<int>* matching::addPrevPath(std::map<int, std::vector<int>*>* prevPaths, std::vector<int>* curPath) {
    if (prevPaths == nullptr) return curPath;
    auto res = new std::vector<int>();
    for(int & item : *curPath) {
        int vIdx = (item + 1) / 2;
        int dir = item % 2;
        auto prevIdx = std::stoi((*this->graph->getVertices())[vIdx - 1]->getId());
//        if (prevIdx == 2980)
//            int mkh = 9;
        if(prevPaths->find(prevIdx) != prevPaths->end()) {
            auto pPath = (*prevPaths)[prevIdx];
//            if +
            if (dir == 1) {
                for (auto pItem: *pPath) {
//                    if (pItem == 0)
//                        int kh = 0;
                    res->push_back(pItem);
                }
            } else { // if -, reverse add
                for (auto it = pPath->rbegin(); it != pPath->rend(); ++it) {
//                    if (*it == -1)
//                        int kh = 0;
                    res->push_back(seqGraph::conjugateIdx(*it));
                }
            }
        } else {
            std::cout<<"error, prev path not found, when parse cycle "<<prevIdx<<std::endl;
        }
    }

//    modify the cyclePaths and mergedPath
    for(int & item : cyclePaths) {
        int vIdx = (item + 1) / 2;
        int dir = item % 2;
        auto prevIdx = std::stoi((*this->graph->getVertices())[vIdx - 1]->getId());
        if (prevPaths->find(prevIdx) != prevPaths->end()) {
            item = prevIdx;
        } else {
            std::cout << "error, prev path not found, when modify cycle " << prevIdx << std::endl;
        }
    }

    std::map<int, int> tmpMergedPaths;
    for(auto & item : mergedPaths) {
//        first
        auto first = item.first;
        int vIdx = (first + 1) / 2;
        auto firstPrevIdx = std::stoi((*this->graph->getVertices())[vIdx - 1]->getId());
//        second
        auto second = item.second;
        vIdx = (second + 1) / 2;
        auto secondPrevIdx = std::stoi((*this->graph->getVertices())[vIdx - 1]->getId());
        if (prevPaths->find(firstPrevIdx) != prevPaths->end() || prevPaths->find(second) != prevPaths->end()) {
            tmpMergedPaths[firstPrevIdx] = secondPrevIdx;
        } else {
            std::cout << "error, prev path not found, when modify cycle " << firstPrevIdx << std::endl;
        }
    }
    this->mergedPaths = tmpMergedPaths;
    return res;
}

void matching::writeMatchResult(std::ofstream& outS) {
//    std::ofstream outF(outFile);
    for (int i = 1; i < this->N - 1;i++) {
        outS << "JUNC " << this->idx2StrDir(i, " ") << " " << this->idx2StrDir(this->matched[i], " ") << " " << this->getIJ(this->matched[i], i) << "\n";
    }
}

bool matching::needMatch() {
    if (this->graph->getJuncSize()==0) {
        return false;
    }
//    if (this->graph->getVertices()->size() == 2 && this->graph->getJuncSize() ==1) {
//        auto sV = graph->getJunctions()->front()->getSource()->getOriginId();
//        auto tV = graph->getJunctions()->front()->getTarget()->getOriginId();
//        auto sD = graph->getJunctions()->front()->getSourceDir();
//        auto tD = graph->getJunctions()->front()->getTargetDir();
//        if (sD == '+') {
//            if (tD == '+') {
//                this->matched[1] = 3;
//                this->matched[2] = 0;
//            } else {
//                this->matched[1] = 2;
//                this->matched[3] = 0;
//            }
//        } else {
//            if (tD == '+') {
//                this->matched[0] = 3;
//                this->matched[2] = 1;
//            } else {
//                this->matched[0] = 2;
//                this->matched[3] = 1;
//            }
//        }
//        return false;
//    }
    return true;
}

seqGraph::Graph *matching::getGraph() const {
    return graph;
}
