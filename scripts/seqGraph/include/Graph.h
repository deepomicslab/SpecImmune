//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_GRAPH_H
#define SEQMAP_GRAPH_H

#include "Vertex.h"
#include "Junction.h"
#include  <set>
#include "map"
#include "unordered_map"
extern std::string SUB_ONLY;
extern int MODEL;

namespace seqGraph {
    typedef std::vector<EndPoint *> EndPointPath;


    //    CRS for sparse matrix
    class SparseMatrix {
    public:
        int N;

        std::vector<float> values;
        std::vector<int> IA;
        std::vector<int> JA;
        std::map<int, float> rowMaxV;
        std::unordered_map<std::string, float> cache;


        explicit SparseMatrix ();
        ~SparseMatrix();

        float getIJ(int i, int j);
//        void addValue(int i, int j, float v);
        void debugPrint();
        float getIRowMax(int i);
        inline bool isEmpty() const {
            return JA.empty();
        }
    };


    class Graph {
    protected:
        float mAvgCoverage;
        std::vector<Vertex *> *vertices;
        std::vector<Junction *> *junctions;
        Vertex *source;
        Vertex *sink;
        float ** ConjugateMatrix;
        std::map<std::string, int>* verticesIdx;
        std::unordered_map<std::string, int>* junctionIdx;
        std::vector<std::set<int>* > connectedJunctionsIdx;
        std::vector<Vertex *> *isolatedVs;
        std::map<int, float> rowMaxV;
        SparseMatrix sparseMatrix;
        bool sparse;
    public:
        bool isReconstructed;
        Graph();

        ~Graph();
        Graph* getSubgraph(int i);
        void originalGraphFree();
        inline int getVCount() const {
            return this->vertices->size();
        }

        inline int subGraphCount() {
            return this->connectedJunctionsIdx.size();
        }

        inline int getJuncSize() {
            return this->junctionIdx->size();
        }

        inline float getIRowMaxV(int i ) {
            if (rowMaxV.count(i) == 0){
                return 0;
            }
            return rowMaxV[i];
        }

        inline bool isSparse() const {
//            return false;
            return this->sparse;
        }

        inline std::vector<Vertex*> * getIsolatedVs() const {
            return this->isolatedVs;
        }

        void removeJunc(Junction* junc);
        Vertex *addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,float aCoverage, float mCredibility, int aCopyNum, int copy_idx);

        Junction *
        addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, float copyNum,
                    float coverage, bool aIsBounded);

        Junction *
        addJunction(std::string& sourceId, std::string& targetId, char sourceDir, char targetDir, float copyNum,
                    float coverage, bool aIsBounded);

        Junction * addJunction(EndPoint* ep3, EndPoint* ep5, float copyNum, float converage, bool isBounded);

        Junction *
        addJunction(Junction* j);

        Vertex *addVertex(Vertex *v);

        Vertex *getVertexById(std::string Id);
        Vertex *getVertexByIdQ(std::string Id);
        Vertex *getVertexByIdx(int idx);

        void removeByGeneAndScore(std::ofstream& cycleFile);

        float getIJ(int i, int j,char sDir,char tDir);

        bool doesPathExists(EndPoint *sourceEndPoint, EndPoint *sinkEndpoint);

        int BFS_EndPoint(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint);

        void BFS_Vertices(Vertex* vertex,std::set<int>* connectedIdx);

        Junction* doesJunctionExist(Junction *junction);

        Junction* doesJunctionExist(Vertex& v1, Vertex& v2, char v1d, char v2d);

        void resetJunctionVisitFlag();

        void parseConnectedComponents();

        void resetVertexVisitFlag();

        void resetShortestPrevEdge();

        float getMAvgCoverage() const;

        void setMAvgCoverage(float mAvgCoverage);

        Vertex *getSource() const;

        void setSource(std::string& Id);

        Vertex *getSink() const;

        void setSink(std::string& Id);

        std::vector<Vertex *> *getVertices() const;

        std::vector<Junction *> *getJunctions() const;

//        float ** getConjugateMatrix();

        SparseMatrix& getConjugateMatrix();

        void initRowMax();
    };
}

#endif //SEQMAP_GRAPH_H
