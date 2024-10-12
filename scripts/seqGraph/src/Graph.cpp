//
// Created by caronkey on 10/5/2021.
//

#include "../include/Graph.h"
#include "../include/Exceptions.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <fstream>
#include "../include/util.h"
#include "exception"

using namespace seqGraph;

SparseMatrix::SparseMatrix() {
    this->IA.push_back(0);
}

SparseMatrix::~SparseMatrix() {

}

float SparseMatrix::getIJ(int i, int j) {
    int cacheSize = 1000;
    std::string tmp = std::to_string(i).append("-").append(std::to_string(j));
    std::string ctmp = std::to_string(seqGraph::conjugateIdx(j)).append("-").append(std::to_string(conjugateIdx(i)));

    if (IA[i + 1] - IA[i] > cacheSize && cache.find(tmp) != cache.end())
        return cache[tmp];
    if (IA[seqGraph::conjugateIdx(j) + 1] - IA[seqGraph::conjugateIdx(j)] > cacheSize && cache.find(ctmp) != cache.end()) {
        return cache[ctmp];
    }
    auto moveC = IA[i];
    auto countI = IA[i + 1];
    for (int k = moveC; k < moveC + (countI - moveC); k++) {
        if (JA[k] == j){
            if (countI - moveC > cacheSize) {
                std::string tmp = std::to_string(i).append("-").append(std::to_string(j));
                cache.emplace(tmp, values[k]);
            }
            return values[k];
        }
    }

//  conjugate
    moveC = IA[seqGraph::conjugateIdx(j)];
    countI = IA[seqGraph::conjugateIdx(j) + 1];
    for (int k = moveC; k < moveC + (countI - moveC); k++) {
        if (JA[k] == conjugateIdx(i)) {
            if (countI - moveC > cacheSize) {
                std::string tmp = std::to_string(i).append("-").append(std::to_string(j));
                cache.emplace(tmp, values[k]);
            }
            return values[k];
        }
    }
    return 0;
}

void SparseMatrix::debugPrint() {
    std::cout<<"values\n";
    for (auto i : values) {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"IA\n";
    for (auto i : IA) {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"JA\n";
    for (auto i : JA) {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
}

float SparseMatrix::getIRowMax(int i) {
//    auto moveC = IA[i];
//    auto countI = IA[i + 1];
//    float maxV = 0;
//    for (int k = moveC; k < moveC + (countI - moveC); k++) {
//        if (values[k] > maxV) maxV = values[k];
//    }
//    return maxV;
    if(rowMaxV.count(i) == 0) {
        return 0;
    }
    return rowMaxV[i];
}


Graph::Graph() {
    this->mAvgCoverage = 0;
    this->vertices = new std::vector<Vertex *>();
    this->junctions = new std::vector<Junction *>();
    verticesIdx = new std::map<std::string, int>();
    junctionIdx = new std::unordered_map<std::string, int>();
    this->isolatedVs = new std::vector<Vertex *>();
    ConjugateMatrix = nullptr;
    source = nullptr;
    sink = nullptr;
}

Graph::~Graph() {
    if (ConjugateMatrix != nullptr) {
        for (int i = 0; i < 2*vertices->size(); ++i) {
            delete[] ConjugateMatrix[i];
        }
    }
    delete[] ConjugateMatrix;
    for (auto item : *vertices) {
        delete item;
    }
    for (auto item : *junctions) {
        delete item;
    }
    vertices->clear();
    vertices->shrink_to_fit();
    junctions->clear();
    junctions->shrink_to_fit();
    verticesIdx->clear();
    junctionIdx->clear();
}

//void Graph::originalGraphFree(){
//    for (auto item : *junctions) {
//        delete item;
//    }
//    junctions->clear();
//    junctions->shrink_to_fit();
//    verticesIdx->clear();
//    junctionIdx->clear();
//}

Graph* Graph::getSubgraph(int i) {
    auto subG = new Graph();
    auto idxs = connectedJunctionsIdx[i];
    for(auto j : *idxs) {
        auto junc = (*this->junctions)[j];
        subG->addJunction(junc);
    }
    subG->isReconstructed = false;
    return subG;
}
Vertex *Graph::addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,float aCoverage, float mCredibility, int aCopyNum, int copy_idx) {
//    create vertex add push
    auto v1 = this->getVertexById(mId);
    if(v1 != nullptr) return v1;
    auto *vertex = new Vertex(mId, aChrom, aStart, aEnd, aCoverage, mCredibility, aCopyNum, this->vertices->size(),copy_idx);
//    vertex->setIdx(this->vertices->size());
    this->vertices->push_back(vertex);
    this->verticesIdx->emplace(mId, vertex->getIdx());
    return vertex;
}

Junction *Graph::addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, float copyNum,
                             float coverage, bool aIsBounded) {
//    auto *junction = new Junction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, aIsBounded);
if (targetVertex->getId() == "EDGE_2331069_length_139_cov_10.154762_0") {
    int kk = 33;
}
auto jun = this->doesJunctionExist(*sourceVertex, *targetVertex, sourceDir, targetDir);
float plasCopy = copyNum;
    if (MODEL == 0)
        plasCopy = (sourceVertex->getScore() + targetVertex->getScore()) / 2 * copyNum;
    if (jun == nullptr) {
//        auto  k = sourceVertex->getId()+ targetVertex->getId()+sourceDir+targetDir;
        auto k = std::to_string(sourceVertex->getIdx()) + sourceDir + std::to_string(targetVertex->getIdx()) +targetDir;
//        throw DuplicateJunctionException(junction);

//        junc with plasmid score
        auto *junction = new Junction(sourceVertex, targetVertex, sourceDir, targetDir, plasCopy, coverage, aIsBounded);
//        junction->junctionToEdge();
        junction->setIdx(this->junctions->size());

        this->junctions->push_back(junction);
        this->junctionIdx->emplace(k, junctions->size());
    //    set vertex not orphan
        sourceVertex->setOrphan(false);
        targetVertex->setOrphan(false);
        sourceVertex->setNextJunc(junction);
        targetVertex->setPrevJunc(junction);
        return junction;
    } else {
        if (plasCopy > jun->getWeight()->getCopyNum()) {
            jun->getWeight()->setCopyNum(plasCopy);
        }
    }
    return jun;
}

Junction *
Graph::addJunction(std::string& sourceId, std::string& targetId, char sourceDir, char targetDir, float copyNum,
                   float coverage ,bool aIsBounded) {
    Vertex *sourceVertex = this->getVertexById(sourceId);
    Vertex *targetVertex = this->getVertexById(targetId);
    return this->addJunction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, aIsBounded);
}

Junction * Graph::addJunction(EndPoint *ep3, EndPoint *ep5, float copyNum, float coverage, bool isBounded) {
    Vertex *sourceVertex = ep3->getVertex();
    Vertex *targetVertex = ep5->getVertex();

    char sourceDir = ep3->getType() == _RIGHT_TOP_ ? '+' : '-';
    char targetDir = ep3->getType() == _LEFT_TOP_ ? '+' : '-';
    return this->addJunction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, isBounded);
}

Vertex *Graph::getVertexById(std::string Id) {
    if (verticesIdx->find(Id) != verticesIdx->end()){
        auto idx = (*verticesIdx)[Id];
        return (*this->vertices)[idx];
    }
//    for (auto *vertex : *this->vertices) {
//        if (vertex->getId() == Id) return vertex;
//    }
//    throw VertexDoesNotExistException(Id);
    return nullptr;
}

Vertex *Graph::getVertexByIdQ(std::string Id) {
//    if (verticesIdx->find(Id) != verticesIdx->end()){
//        auto idx = (*verticesIdx)[Id];
//    }
    return (*this->vertices)[(*verticesIdx)[Id]];

//    for (auto *vertex : *this->vertices) {
//        if (vertex->getId() == Id) return vertex;
//    }
//    throw VertexDoesNotExistException(Id);
//    return nullptr;
}


int Graph::BFS_EndPoint(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint) {
    std::queue<EndPoint *> EPQueue;
    EPQueue.push(sourceEndpoint);

    while (!EPQueue.empty()) {
        EndPoint *currentEP = EPQueue.front();
        EPQueue.pop();

        for (Edge *e : *(currentEP->getOutEdges())) {
            // if (e->hasCopy()) {
            EndPoint *nextEP = e->getTarget();
            if (!nextEP->isVisited()) {
                nextEP->setVisited(true);
                nextEP->setShortestPrevEdge(e);
                EPQueue.push(nextEP);
            }
        }
    }
    this->resetVertexVisitFlag();

    Edge *prevEdge;
    EndPoint *prevEP;
    EndPoint *currentEP = sinkEndpoint;
    int found;
    while (true) {
        prevEdge = currentEP->getShortestPrevEdge();
        if (prevEdge == nullptr) {
            found = -1;
            break;
            // return -1;  // no path from aStartVertex to aTargetVertex
        }

        prevEP = prevEdge->getSource();
        if (prevEP == sourceEndpoint) {
            found = 0;
            break;
            // return 0;
        }
        currentEP = prevEP;
    }
    this->resetShortestPrevEdge();
    return found;
}

void Graph::BFS_Vertices(Vertex* vertex, std::set<int>* connectedIdx){
    std::queue<Vertex*> VQueuen;
    VQueuen.push(vertex);

    while (!VQueuen.empty()) {
        auto currentV = VQueuen.front();
        VQueuen.pop();
        for (auto e : currentV->getNextJuncs()) {
            if (e->isRemoved()) continue;
            // if (e->hasCopy()) {
            auto *nextV = e->getTarget();
            connectedIdx->insert(e->getIdx());
            if (!nextV->isVisited()) {
                nextV->setIsVisited(true);
                VQueuen.push(nextV);
            }
        }
        for (auto e : currentV->getPrevJuncs()) {
            if (e->isRemoved()) continue;
            // if (e->hasCopy()) {
            auto *nextV = e->getSource();
            if (!nextV->isVisited()) {
                nextV->setIsVisited(true);
                VQueuen.push(nextV);
            }
        }
    }
//    this->resetVertexVisitFlag();
}



Junction* Graph::doesJunctionExist(Junction *aJunction) {
    for(auto item: *this->junctions){
        if(!item->isRemoved() && *item == *aJunction) {
            return item;
        }
    }
    return nullptr;
}


Junction* Graph::doesJunctionExist(Vertex& v1, Vertex& v2, char v1d, char v2d) {
    auto  k = std::to_string(v1.getIdx()) +v1d + std::to_string(v2.getIdx())+v2d;
    if (v1d == '-') v1d = '+';
    else v1d = '-';
    if (v2d == '-') v2d = '+';
    else v2d = '-';
//    auto rk = v2.getId()+v1.getId()+v2d+v1d;
    auto  rk = std::to_string(v2.getIdx()) +v2d + std::to_string(v1.getIdx())+v1d;
    if (junctionIdx->find(k) != junctionIdx->end()) {
        auto idx = (*junctionIdx)[k];
        return (*junctions)[idx - 1];
    } else if (junctionIdx->find(rk) != junctionIdx->end()) {
        auto idx = (*junctionIdx)[rk];
        return (*junctions)[idx - 1];
    }
//    for(auto item: *this->junctions){
//        auto& s1 = *item->getSource();
//        auto& s2 = *item->getTarget();
//        auto s1d = item->getSourceDir();
//        auto s2d = item->getTargetDir();
//        if ((s1 == v1 && s2 == v2 && v1d == s1d && v2d == s2d) || (s1 == v2 && s2 == v1 && v1d != s1d && v2d != s2d) )
//            return item;
//    }
    return nullptr;
}
bool Graph::doesPathExists(EndPoint *sourceEndPoint, EndPoint *sinkEndpoint) {
    bool isReach = false;
    EndPoint *startEP = sourceEndPoint;
    EndPointPath EPStack;
    EPStack.push_back(startEP);
    Edge *nextEdge = startEP->getOneNextEdge();
    std::vector<Junction*> js;
    while (true) {
        if (nextEdge == NULL) {
            EPStack.pop_back();
            if (EPStack.empty()) {
                break;
            } else {
                startEP = EPStack.back();
                nextEdge = startEP->getOneNextEdge();
            }
        } else {
            EndPoint *nextEP = nextEdge->getTarget();
            if (nextEP == sinkEndpoint || nextEP == sinkEndpoint->getMateEp()) {
                isReach = true;
                break;
            }
            js.push_back(nextEdge->getJunction());
            nextEdge->setVisited(true);
            EPStack.push_back(nextEP);
            startEP = nextEP;
            nextEdge = startEP->getOneNextEdge();
        }
    }
    this->resetJunctionVisitFlag();
    if(isReach) {
        for (auto j : js) {
            j->getWeight()->setCopyNum(100000);
        }
    }
    return isReach;
}

void Graph::removeByGeneAndScore(std::ofstream& cycleFile) {
//    for vertex, gene score.
    bool  flag = false;
    for (auto v : *this->getVertices()) {
        if (!v->isGeneAndScoreOk()) continue;
        for (auto junc : v->getNextJuncs()) {
            auto sDir = junc->getSourceDir();
            auto tDir = junc->getTargetDir();
            if (sDir == '+') {
                if (tDir == '+') {
                    flag =doesPathExists(junc->getTarget()->getEp3(), v->getEp5());
                } else {
                    flag =doesPathExists(junc->getTarget()->getRep3(), v->getEp5());
                }
            } else {
                if (tDir == '+') {
                    flag =doesPathExists(junc->getTarget()->getEp3(), v->getRep5());
                } else {
                    flag =doesPathExists(junc->getTarget()->getRep3(), v->getRep5());
                }
            }
            if (! flag) {
                this->removeJunc(junc);
            }
        }

        for (auto junc : v->getPrevJuncs()) {
            if (junc->isRemoved()) continue;
            auto sDir = junc->getSourceDir();
            auto tDir = junc->getTargetDir();
            if (sDir == '+') {
                if (tDir == '+') {
                    flag =doesPathExists(junc->getTarget()->getEp3(), v->getEp5());
                } else {
                    flag =doesPathExists(junc->getTarget()->getRep3(), v->getEp5());
                }
            } else {
                if (tDir == '+') {
                    flag =doesPathExists(junc->getTarget()->getEp3(), v->getRep5());
                } else {
                    flag =doesPathExists(junc->getTarget()->getRep3(), v->getRep5());
                }
            }
            if (! flag) {
                this->removeJunc(junc);
            }
        }

    }

}

float Graph::getIJ(int i, int j,char sDir,char tDir) {
    auto k = std::to_string(i) + sDir + std::to_string(j)  + tDir;
//    try {
//        auto idx = (*junctionIdx)[k];
//        return (*this->junctions)[idx - 1]->getWeight()->getCopyNum();
//    } catch (std::exception & e) {
//        if (sDir == '-') sDir = '+';
//        else sDir = '-';
//        if (tDir == '-') tDir = '+';
//        else tDir = '-';
//        auto  rk = std::to_string(j) + tDir + std::to_string(i)  + sDir;
//        try {
//            auto idx = (*junctionIdx)[rk];
//            return (*this->junctions)[idx - 1]->getWeight()->getCopyNum();
//        } catch (std::exception& e1) {
//            return 0;
//        }
//    }
//    if (junctionIdx->find(k) != junctionIdx->end()) {
    std::string test = "sewwe";
    if (junctionIdx->find(k) != junctionIdx->end()) {
        auto idx = (*junctionIdx)[k];
        return (*this->junctions)[idx - 1]->getWeight()->getCopyNum();
    }
//    }
    if (sDir == '-') sDir = '+';
    else sDir = '-';
    if (tDir == '-') tDir = '+';
    else tDir = '-';
//    auto rk = v2.getId()+v1.getId()+v2d+v1d;
    auto  rk = std::to_string(j) + tDir + std::to_string(i)  + sDir;
//    if (junctionIdx->find(rk) != junctionIdx->end()) {
    if (junctionIdx->find(rk) != junctionIdx->end()) {
        auto idx = (*junctionIdx)[rk];
        return (*this->junctions)[idx - 1]->getWeight()->getCopyNum();
    }
//    }
    return 0;
}

void Graph::resetJunctionVisitFlag() {
    for (Junction *junction : *this->junctions) {
        if (junction->isRemoved()) continue;
        junction->getCEdge()->setVisited(false);
        junction->getOEdge()->setVisited(false);
    }
}

void Graph::resetVertexVisitFlag() {
    for (auto *vertex : *this->vertices) {
        vertex->getEp3()->setVisited(false);
        vertex->getRep3()->setVisited(false);
    }
}

void Graph::resetShortestPrevEdge() {
    for (auto *vertex : *this->vertices) {
        vertex->getEp3()->setShortestPrevEdge(nullptr);
        vertex->getRep3()->setShortestPrevEdge(nullptr);
    }
}

float Graph::getMAvgCoverage() const {
    return mAvgCoverage;
}

void Graph::setMAvgCoverage(float mAvgCoverage) {
    Graph::mAvgCoverage = mAvgCoverage;
}

Vertex *Graph::getSource() const {
    return source;
}

void Graph::setSource(std::string& Id) {
    this->source = this->getVertexById(Id);
}

Vertex *Graph::getSink() const {
    return sink;
}

void Graph::setSink(std::string& Id) {
    this->sink = this->getVertexById(Id);
}

std::vector<Vertex *> *Graph::getVertices() const {
    return vertices;
}

std::vector<Junction *> *Graph::getJunctions() const {
    return junctions;
}

SparseMatrix&  Graph::getConjugateMatrix(){
    return this->sparseMatrix;
//    if(!this->sparseMatrix.isEmpty())
//        return this->sparseMatrix;

//    int n = this->getVCount();
//    this->ConjugateMatrix = new float *[2*n + 1];
//    for(int i = 0; i < 2*n + 1; i++) {
//
//        auto t = new float[2*n + 1];
//        std::fill_n(t,(2*n + 1), 0);
//        this->sparseMatrix.addValue(i] = t;
//    }
//    float initV = 0;
//    for(auto junc : *junctions) {
//        int i = junc->getSource()->getIdx();
//        int j = junc->getTarget()->getIdx();
//        int sDir = junc->getSourceDir();
//        int tDir = junc->getTargetDir();
//        float weightValue = junc->getWeight()->getCopyNum();
//        if (sDir == '+') {
//            if (tDir == '+') {
//                this->sparseMatrix.addValue(2*j + 1, 2*i + 1, weightValue);
//                this->sparseMatrix.addValue(2*j + 1, 2*i + 1, weightValue);
//                this->sparseMatrix.addValue(2*(i+1), 2*(j+1), weightValue);
//            } else {
//                this->sparseMatrix.addValue(2*(j + 1), 2*i+1, weightValue);
//                this->sparseMatrix.addValue(2*(i+1), 2*j+1, weightValue);
//            }
//        } else {
//            if (tDir == '+') {
//                this->sparseMatrix.addValue(2*j+1, 2*(i + 1), weightValue);
//                this->sparseMatrix.addValue(2*i+1, 2*(j+1), weightValue);
//            } else {
//                this->sparseMatrix.addValue(2*i+1, 2*j+1, weightValue);
//                this->sparseMatrix.addValue(2*(j+1), 2*(i+1), weightValue);
//            }
//        }
//    }

//   vertices are sorted by idx already
//    this->sparseMatrix.IA.push_back(0);
//    for (auto v : *this->vertices) {
////        v + row
//        auto juncs = v->getPrevJuncs();
//        sort(juncs.begin(), juncs.end(), [](const Junction *lhs, const Junction *rhs) {
//            return lhs->getSource()->getIdx() < rhs->getSource()->getIdx();
//        });
//        std::vector<float> toPositiveVs;
//        std::vector<int> toPositiveJAs;
//        std::vector<float> toNegativeVs;
//        std::vector<int> toNegativeJAs;
//        for (auto junc: juncs) {
//            int j = 2 * (junc->getSource()->getIdx() + 1);
//            if (junc->getSourceDir() == '+')
//                j = 2 * junc->getSource()->getIdx() + 1;
//            if ( junc->getTargetDir() == '+') {
//                toPositiveVs.push_back(junc->getWeight()->getCopyNum());
//                toPositiveJAs.push_back(j);
//            } else {
//                toNegativeVs.push_back(junc->getWeight()->getCopyNum());
//                toNegativeJAs.push_back(j);
//            }
//        }
//
////          value
//        this->sparseMatrix.values.insert(this->sparseMatrix.values.end(), toPositiveVs.begin(),toPositiveVs.end());
//        this->sparseMatrix.values.insert(this->sparseMatrix.values.end(), toNegativeVs.begin(),toNegativeVs.end());
////          value col
//        this->sparseMatrix.JA.insert(this->sparseMatrix.JA.end(), toPositiveJAs.begin(), toPositiveJAs.end());
//        this->sparseMatrix.JA.insert(this->sparseMatrix.JA.end(), toNegativeJAs.begin(), toNegativeJAs.end());
//
//
////      IA, for value count
//        this->sparseMatrix.IA.push_back(this->sparseMatrix.IA.back() + toPositiveVs.size());
//        this->sparseMatrix.IA.push_back(this->sparseMatrix.IA.back() + toNegativeVs.size());
//
//
//        juncs = v->getNextJuncs();
//        for (auto junc : juncs) {
//            int j = junc->getSource()->getIdx();
//            if ( junc->getSourceDir() == '-') {
//                toPositiveVs.push_back(junc->getWeight()->getCopyNum());
//            } else {
//                toNegativeVs.push_back(junc->getWeight()->getCopyNum());
//            }
//        }
//
//        float pM = toPositiveVs.empty() ? 0 : *std::max_element(toPositiveVs.begin(),toPositiveVs.end());
//        float nM = toNegativeVs.empty() ? 0 : *std::max_element(toNegativeVs.begin(),toNegativeVs.end());
//        this->rowMaxV.emplace(2 * v->getIdx() + 1, pM);
//        this->rowMaxV.emplace(2 * (v->getIdx() + 1), nM);
//        this->sparseMatrix.rowMaxV.emplace(2 * v->getIdx() + 1, pM);
//        this->sparseMatrix.rowMaxV.emplace(2 * (v->getIdx() + 1), nM);



//   vertices are sorted by idx already
//    this->sparseMatrix.IA.push_back(0);
//    for (auto v : *this->vertices) {
////        v + row
//        auto ep5InEdges = v->getEp5()->getInEdges();
//        sort(ep5InEdges->begin(), ep5InEdges->end(), [](const Edge* lhs, const Edge* rhs) {
//            return lhs->getSource()->getIdx() < rhs->getSource()->getIdx();
//        });
//        for (auto e: *ep5InEdges) {
//            int j = e->getSource()->getIdx();
////          value
//            this->sparseMatrix.values.push_back(e->getWeight()->getCopyNum());
////          value col
//            this->sparseMatrix.JA.push_back(e->getSource()->getIdx());
//        }
////      IA, for value count
//        this->sparseMatrix.IA.push_back(this->sparseMatrix.IA.back() + ep5InEdges->size());
////        v - row
//        auto rEp5InEdges = v->getRep5()->getInEdges();
//
//        sort(rEp5InEdges->begin(), rEp5InEdges->end(), [](const Edge* lhs, const Edge* rhs) {
//            return lhs->getSource()->getIdx() < rhs->getSource()->getIdx();
//        });
//        for (auto e: *rEp5InEdges) {
//            int j = e->getSource()->getIdx();
////          value
//            this->sparseMatrix.values.push_back(e->getWeight()->getCopyNum());
////          value col
//            this->sparseMatrix.JA.push_back(e->getSource()->getIdx());
//        }
////      IA, for value count
//        this->sparseMatrix.IA.push_back(this->sparseMatrix.IA.back() + rEp5InEdges->size());
//
//    }
}

void Graph::parseConnectedComponents() {
    int i = 0;
    for(auto & vertice : *vertices) {
        if (i == vertices->size()) break;
        if (!vertice->isVisited()) {
            auto * connectedIdx = new std::set<int>;
            BFS_Vertices(vertice, connectedIdx);
            i += connectedIdx->size();
            if (!connectedIdx->empty()){
                this->connectedJunctionsIdx.push_back(connectedIdx);
            } else {
                this->isolatedVs->push_back(vertice);
            }
        }
    }
    if (!SUB_ONLY.empty()) {
        for(i= 0;i < this->connectedJunctionsIdx.size(); i++) {
            std::set<std::string> write_segs;
            std::set<std::string> write_juncs;
            auto item = connectedJunctionsIdx[i];
            if (item->size() < 500) continue;
            auto outfile = SUB_ONLY + std::to_string(i) + ".subg";
            std::ofstream out(outfile);
//            auto idxs = connectedJunctionsIdx[i];
            for(auto j : *item) {
                auto junc = (*this->junctions)[j];
                write_segs.emplace("SEG "+junc->getSource()->getOriginId()+" "+ std::to_string(junc->getSource()->getWeight()->getCoverage())+" " + std::to_string(junc->getSource()->getWeight()->getCopyNum()));
                write_segs.emplace("SEG "+junc->getTarget()->getOriginId()+" "+ std::to_string(junc->getTarget()->getWeight()->getCoverage()) + " "+ std::to_string(junc->getTarget()->getWeight()->getCopyNum()));
                write_juncs.emplace("JUNC "+ junc->getSource()->getOriginId()+" "+junc->getSourceDir()+" "+junc->getTarget()->getOriginId()+" "+junc->getTargetDir() + " "+ std::to_string(junc->getWeight()->getCopyNum()));
            }
            for(const auto& j : write_segs) {
                out<<j<<"\n";
            }
            for(const auto& j : write_juncs) {
                out<<j<<"\n";
            }
            out.close();
        }
    }
}

Vertex *Graph::addVertex(Vertex *v) {
    auto v1 = this->getVertexById(v->getId());
    if(v1 != nullptr) return v1;
    v->setIdx(this->vertices->size());
    this->vertices->push_back(v);
    this->verticesIdx->emplace(v->getId(), v->getIdx());
    return v;
}

Junction *Graph::addJunction(Junction *j) {
    this->addVertex(j->getSource());
    this->addVertex(j->getTarget());
//    auto k = j->getSource()->getId()+ j->getTarget()->getId()+j->getSourceDir()+j->getTargetDir();
    auto k = std::to_string(j->getSource()->getIdx()) + j->getSourceDir() + std::to_string(j->getTarget()->getIdx()) +j->getTargetDir();

    this->junctions->push_back(j);
    j->setIdx(this->junctions->size());
    this->junctionIdx->emplace(k, junctions->size());
    return j;
}

void Graph::removeJunc(Junction *junc) {
    junc->setRemoved(true);
}
void Graph::initRowMax() {

    if (this->junctions->size() == 0 || this->vertices->size() / this->junctions->size() > 30) this->sparse = false;
    else this->sparse = true;
    //   vertices are sorted by idx already
    this->sparseMatrix.IA.push_back(0);
    for (auto v : *this->vertices) {
//        v + row
        auto juncs = v->getPrevJuncs();
        sort(juncs.begin(), juncs.end(), [](const Junction *lhs, const Junction *rhs) {
            return lhs->getSource()->getIdx() < rhs->getSource()->getIdx();
        });
        std::vector<float> toPositiveVs;
        std::vector<int> toPositiveJAs;
        std::vector<float> toNegativeVs;
        std::vector<int> toNegativeJAs;
        for (auto junc: juncs) {
            int j = 2 * (junc->getSource()->getIdx() + 1);
            if (junc->getSourceDir() == '+')
                j = 2 * junc->getSource()->getIdx() + 1;
            if ( junc->getTargetDir() == '+') {
                toPositiveVs.push_back(junc->getWeight()->getCopyNum());
                toPositiveJAs.push_back(j);
            } else {
                toNegativeVs.push_back(junc->getWeight()->getCopyNum());
                toNegativeJAs.push_back(j);
            }
        }

//          value
        this->sparseMatrix.values.insert(this->sparseMatrix.values.end(), toPositiveVs.begin(),toPositiveVs.end());
        this->sparseMatrix.values.insert(this->sparseMatrix.values.end(), toNegativeVs.begin(),toNegativeVs.end());
//          value col
        this->sparseMatrix.JA.insert(this->sparseMatrix.JA.end(), toPositiveJAs.begin(), toPositiveJAs.end());
        this->sparseMatrix.JA.insert(this->sparseMatrix.JA.end(), toNegativeJAs.begin(), toNegativeJAs.end());


//      IA, for value count
        this->sparseMatrix.IA.push_back(this->sparseMatrix.IA.back() + toPositiveVs.size());
        this->sparseMatrix.IA.push_back(this->sparseMatrix.IA.back() + toNegativeVs.size());


        juncs = v->getNextJuncs();
        for (auto junc : juncs) {
            int j = junc->getSource()->getIdx();
            if ( junc->getSourceDir() == '-') {
                toPositiveVs.push_back(junc->getWeight()->getCopyNum());
            } else {
                if (junc->getWeight()->getCopyNum() >= 1000) {
                    auto tmppp = 99;
                }
                toNegativeVs.push_back(junc->getWeight()->getCopyNum());
            }
        }

        float pM = toPositiveVs.empty() ? 0 : *std::max_element(toPositiveVs.begin(),toPositiveVs.end());
        float nM = toNegativeVs.empty() ? 0 : *std::max_element(toNegativeVs.begin(),toNegativeVs.end());
        if (2 * (v->getIdx() + 1) == 54) {
            auto tmppp = 99;
        }
        this->rowMaxV.emplace(2 * v->getIdx() + 1, pM);
        this->rowMaxV.emplace(2 * (v->getIdx() + 1), nM);
        this->sparseMatrix.rowMaxV.emplace(2 * v->getIdx() + 1, pM);
        this->sparseMatrix.rowMaxV.emplace(2 * (v->getIdx() + 1), nM);

    }


//    for (auto v : *this->vertices) {
//        auto juncs = v->getPrevJuncs();
//        std::vector<float> toPositiveVs;
//        std::vector<float> toNegativeVs;
//        for (auto junc: juncs) {
//            int j = 2 * (junc->getSource()->getIdx() + 1);
//            if (junc->getSourceDir() == '+')
//                j = 2 * junc->getSource()->getIdx() + 1;
//            if ( junc->getTargetDir() == '+') {
//                toPositiveVs.push_back(junc->getWeight()->getCopyNum());
//            } else {
//                toNegativeVs.push_back(junc->getWeight()->getCopyNum());
//            }
//        }
//
//        this->sparseMatrix.values.insert(this->sparseMatrix.values.end(), toPositiveVs.begin(),toPositiveVs.end());
//        this->sparseMatrix.values.insert(this->sparseMatrix.values.end(), toNegativeVs.begin(),toNegativeVs.end());
//        juncs = v->getNextJuncs();
//        for (auto junc : juncs) {
//            int j = junc->getSource()->getIdx();
//            if ( junc->getSourceDir() == '-') {
//                toPositiveVs.push_back(junc->getWeight()->getCopyNum());
//            } else {
//                toNegativeVs.push_back(junc->getWeight()->getCopyNum());
//            }
//        }
//
//        float pM = toPositiveVs.empty() ? 0 : *std::max_element(toPositiveVs.begin(),toPositiveVs.end());
//        float nM = toNegativeVs.empty() ? 0 : *std::max_element(toNegativeVs.begin(),toNegativeVs.end());
//        this->rowMaxV.emplace(2 * v->getIdx() + 1, pM);
//        this->rowMaxV.emplace(2 * (v->getIdx() + 1), nM);
}