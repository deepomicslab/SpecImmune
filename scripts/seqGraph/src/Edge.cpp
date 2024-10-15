//
// Created by caronkey on 10/5/2021.
//

#include <string>
#include "../include/Edge.h"

using namespace seqGraph;

Edge::Edge(EndPoint *source, EndPoint *target, float aCoverage,float copyNum, int type) {
    this->source = source;
    this->target = target;
    this->weight = new Weight(aCoverage);
    this->weight->setCopyNum(copyNum);
    this->type = type;
    this->visited = false;
}

Edge::Edge(EndPoint *source, EndPoint *target, Weight *weight, int type) {
    this->source = source;
    this->target = target;
    source->addEdge(this);
    target->addEdge(this);
    this->weight = weight;
    this->type = type;
    this->visited = false;
}

Edge::~Edge() {

}

std::string Edge::getInfo() {
    return this->source->getInfo() + this->target->getInfo();
}

bool Edge::isVisited() const {
    return visited;
}

void Edge::setVisited(bool visited) {
    this->visited = visited;
}

int Edge::getType() const {
    return type;
}

bool Edge::hasCopy() {
    return this->weight->getCopyNum() >= 1;
}

void Edge::traverse() {
    this->weight->decreaseCopyNum();
}

EndPoint *Edge::getSource() const {
    return source;
}

EndPoint *Edge::getTarget() const {
    return target;
}

Weight *Edge::getWeight() const {
    return weight;
}

Junction *Edge::getJunction() const {
    return junction;
}

void Edge::setJunction(Junction *junction) {
    Edge::junction = junction;
}
