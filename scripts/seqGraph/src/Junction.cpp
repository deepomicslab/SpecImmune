//
// Created by caronkey on 12/5/2021.
//

#include "../include/Junction.h"

using namespace seqGraph;

Junction::Junction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, float copyNum,
                   float coverage, bool aIsBounded) {
    this->sourceDir = sourceDir;
    this->targetDir = targetDir;
    this->weight = new Weight(coverage);
    weight->setCopyNum(copyNum);
    this->source = sourceVertex;
    this->target = targetVertex;
    this->hasLowerBoundLimit = aIsBounded;
    this->removed = false;
    this->junctionToEdge();
}

Junction::~Junction() {
    delete oEdge;
    delete cEdge;
    delete weight;
}

void Junction::junctionToEdge() {
    if (this->sourceDir == _POSITIVE_DIR_ && this->targetDir == _POSITIVE_DIR_) {
        this->oEdge = new Edge(this->source->getEp3(), this->target->getEp5(), this->weight, _OUTER_EDGE_);
        this->cEdge = new Edge(this->target->getRep3(), this->source->getRep5(), this->weight,
                               _OUTER_EDGE_);
    } else if (this->sourceDir == _NEGATIVE_DIR_ && this->targetDir == _NEGATIVE_DIR_) {
        this->oEdge = new Edge(this->source->getRep3(), this->target->getRep5(), this->weight,
                               _OUTER_EDGE_);
        this->cEdge = new Edge(this->target->getEp3(), this->source->getEp5(), this->weight, _OUTER_EDGE_);
    } else if (this->sourceDir == _POSITIVE_DIR_ && this->targetDir == _NEGATIVE_DIR_) {
        this->oEdge = new Edge(this->source->getEp3(), this->target->getRep5(), this->weight, _OUTER_EDGE_);
        this->cEdge = new Edge(this->target->getEp3(), this->source->getRep5(), this->weight, _OUTER_EDGE_);
    } else if (this->sourceDir == _NEGATIVE_DIR_ && this->targetDir == _POSITIVE_DIR_) {
        this->oEdge = new Edge(this->source->getRep3(), this->target->getEp5(), this->weight, _OUTER_EDGE_);
        this->cEdge = new Edge(this->target->getRep3(), this->source->getEp5(), this->weight, _OUTER_EDGE_);
    }
    this->oEdge->setJunction(this);
    this->cEdge->setJunction(this);
}

bool Junction::operator==(const Junction &otherJunc) const {
    return((*this->source == *otherJunc.target && *this->target == *otherJunc.source)  ||
            (*this->source == *otherJunc.source && *this->target == *otherJunc.target && this->sourceDir == otherJunc.sourceDir &&
            this->targetDir != otherJunc.targetDir));
//    return (this->oEdge->getInfo() == otherJunc.oEdge->getInfo() &&
//            this->cEdge->getInfo() == otherJunc.cEdge->getInfo()) ||
//           (this->oEdge->getInfo() == otherJunc.cEdge->getInfo() &&
//            this->cEdge->getInfo() == otherJunc.oEdge->getInfo());
}

char Junction::getSourceDir() const {
    return sourceDir;
}

char Junction::getTargetDir() const {
    return targetDir;
}

Weight *Junction::getWeight() const {
    return weight;
}

Edge *Junction::getOEdge() const {
    return oEdge;
}

Edge *Junction::getCEdge() const {
    return cEdge;
}

std::string Junction::getInfo() {
    return this->source->getId() + this->sourceDir + "=>" + this->target->getId() + this->targetDir;
}

void Junction::restoreCopy() {
    this->weight->restore();
}

void Junction::backupCopy() {
    this->weight->backup();
}

Vertex *Junction::getSource() const {
    return source;
}

Vertex *Junction::getTarget() const {
    return target;
}

bool Junction::isHasLowerBoundLimit() const {
    return hasLowerBoundLimit;
}

void Junction::setHasLowerBoundLimit() { hasLowerBoundLimit = true; }
void Junction::resetHasLowerBoundLimit() { hasLowerBoundLimit = false; }
void Junction::checkLowerBound() { hasLowerBoundLimit = this->weight->getCopyNum() >= 1; }

bool Junction::hasCopy(){
    return this->oEdge->hasCopy();
}

bool Junction::isInferred() const {
    return inferred;
}

bool Junction::isRemoved() {
    return this->removed;
}
void Junction::setRemoved(bool r) {
    this->removed = r;
}

void Junction::setInferred(bool inferred) {
    Junction::inferred = inferred;
}

int Junction::getIdx() const {
    return idx;
}

void Junction::setIdx(int idx) {
    Junction::idx = idx;
}
