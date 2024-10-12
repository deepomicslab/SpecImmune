//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_EDGE_H
#define SEQMAP_EDGE_H
#define _INNER_EDGE_ 0
#define _OUTER_EDGE_ 1

#include "EndPoint.h"
#include <string>
#include "Junction.h"
namespace seqGraph {
    class EndPoint;
    class Junction;

    class Edge {
    protected:
        EndPoint *source;
        EndPoint *target;
        Weight *weight;
        Junction* junction;
        int type;
        bool visited;
    public:
        Edge(EndPoint *source, EndPoint *target, float aCoverage, float copyNum, int type);

        Edge(EndPoint *source, EndPoint *target, Weight *weight, int type);

        ~Edge();

        std::string getInfo();

        int getType() const;

        bool isVisited() const;

        void setVisited(bool visited);

        bool hasCopy();
        void traverse();

        EndPoint *getSource() const;

        EndPoint *getTarget() const;

        Weight *getWeight() const;

        Junction *getJunction() const;

        void setJunction(Junction *junction);
    };
}


#endif //SEQMAP_EDGE_H
