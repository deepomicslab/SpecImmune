//
// Created by caronkey on 10/5/2021.
//

#include "../include/Vertex.h"

using namespace seqGraph;

bool Vertex::sameVertex(const seqGraph::Vertex &v2) {
    auto idx1 = this->Id.find_last_of('_');
    auto idx2 = v2.getId().find_last_of('_');
    return this->Id.substr(0, idx1) == v2.getId().substr(0,idx2);
}

Vertex::Vertex(std::string mId, std::string aChrom, int aStart, int aEnd,float aCoverage, float aCredibility, int aCopyNum, int idx, int copy_idx) {
    this->copy_idx = copy_idx;
    this->Id = mId;
    this->orphan = true;
    this->chrom = aChrom;
    this->start = aStart;
    this->end = aEnd;
    this->length = this->end - this->start;
    this->depth = aCoverage;
    this->credibility = aCredibility;
    this->visited = false;
    this->idx = idx;
//    create weight
    weight = new Weight(aCoverage);
    weight->setCopyNum(aCopyNum);
//    create endpoints
    this->EP5 = new EndPoint(weight, _LEFT_TOP_, mId);
    this->rEP5 = new EndPoint(weight, _RIGHT_BOTTOM_, mId);
    this->EP3 = new EndPoint(weight, _RIGHT_TOP_, mId);
    this->rEP3 = new EndPoint(weight, _LEFT_BOTTOM_, mId);

    this->EP3->setVertex(this);
    this->rEP3->setVertex(this);
    this->EP5->setVertex(this);
    this->rEP5->setVertex(this);
    this->setMateEP();

    //    create inner edges, inner edges has same weight with vertex weigh
//    auto *topEdge = new Edge(this->EP5, this->EP3, aCoverage, _INNER_EDGE_);
//    auto *bottomEdge = new Edge(this->rEP5, this->rEP3, aCoverage, _INNER_EDGE_);
}


Vertex::~Vertex() {
    delete EP3;
    delete EP5;
    delete rEP3;
    delete rEP5;
    prevJuncs.clear();
    prevJuncs.shrink_to_fit();
    nextJuncs.clear();
    nextJuncs.shrink_to_fit();
    delete weight;
}

bool Vertex::operator==(const Vertex &rhs) const {
    return this->Id == rhs.Id;
};

bool Vertex::operator>(const Vertex &rhs) const {
    return this->getId() > rhs.getId();
}

bool Vertex::operator<(const Vertex &rhs) const {
    return this->getId() > rhs.getId();
}

const std::string Vertex::getId() const {
    return Id;
}

const std::string Vertex::getOriginId() const {
    auto idx1 = this->Id.find_last_of('_');
    return this->Id.substr(0, idx1);
}

void Vertex::setId(const int mId) {
    Vertex::Id = mId;
}

void Vertex::setIdx(const int idx) {
    this->idx = idx;
    this->getEp3()->setIdx(2 * idx + 1);
    this->getEp5()->setIdx(2 * idx + 1);

    this->getRep3()->setIdx(2 * (idx + 1));
    this->getRep5()->setIdx(2 * (idx + 1));
}

int Vertex::getIdx() const {
    return this->idx;
}
Weight *Vertex::getWeight() const {
    return weight;
}

void Vertex::setWeight(Weight *weight) {
    Vertex::weight = weight;
}

EndPoint *Vertex::getEp3() const {
    return EP3;
}

void Vertex::setEp3(EndPoint *ep3) {
    EP3 = ep3;
}

EndPoint *Vertex::getEp5() const {
    return EP5;
}

void Vertex::setEp5(EndPoint *ep5) {
    EP5 = ep5;
}

EndPoint *Vertex::getRep3() const {
    return rEP3;
}

void Vertex::setRep3(EndPoint *rEp3) {
    rEP3 = rEp3;
}

EndPoint *Vertex::getRep5() const {
    return rEP5;
}

void Vertex::setRep5(EndPoint *rEp5) {
    rEP5 = rEp5;
}


bool Vertex::isOrphan() const {
    return orphan;
}

void Vertex::setOrphan(bool orphan) {
    Vertex::orphan = orphan;
}

void Vertex::restoreCopy() {
    this->weight->restore();
}

void Vertex::backupCopy() {
    this->weight->backup();
}

void Vertex::setHasLowerBoundLimit() { hasLowerBoundLimit = true; }
void Vertex::resetHasLowerBoundLimit() { hasLowerBoundLimit = false; }
void Vertex::checkLowerBound() { hasLowerBoundLimit = this->weight->getCopyNum() >= 1; }

const std::string Vertex::getChrom() const {
    return chrom;
}

void Vertex::setChrom(const std::string mChrom) {
    Vertex::chrom = mChrom;
}

int Vertex::getStart() const {
    return start;
}

void Vertex::setStart(int mStart) {
    Vertex::start = mStart;
}

int Vertex::getEnd() const {
    return end;
}

void Vertex::setEnd(int mEnd) {
    Vertex::end = mEnd;
}

float Vertex::getCredibility() const {
    return credibility;
}

void Vertex::setCredibility(float mCredibility) {
    Vertex::credibility = mCredibility;
}

bool Vertex::isHasLowerBoundLimit() {
    return this->hasLowerBoundLimit;
}

bool Vertex::hasCopy() {
    return this->EP3->hasCopy();
}

float Vertex::getInCoverage(){
    return this->EP5->getInCoverage();
}
float Vertex::getOutCoverage(){
    return this->EP3->getOutCoverage();
}

void Vertex::setMateEP() {
    this->EP5->setMateEp(this->EP3);
    this->EP3->setMateEp(this->EP5);
    this->rEP5->setMateEp(this->rEP3);
    this->rEP3->setMateEp(this->rEP5);
}

const std::vector<Junction *> &Vertex::getNextJuncs() const {
    return nextJuncs;
}

const std::vector<Junction *> &Vertex::getPrevJuncs() const {
    return prevJuncs;
}

void Vertex::setNextJunc(Junction* v) {
    this->nextJuncs.push_back(v);
    if(v->getTarget()->getCopyIdx() == 0){
        if (v->getTargetDir() == '+') this->nextJuncCountIgnoreCopyPositive++;
        else this->nextJuncCountIgnoreCopyNegative++;
    }
}
void Vertex::setPrevJunc(Junction* v){
    this->prevJuncs.push_back(v);
    if (v->getSource()->getCopyIdx() == 0) {
        if (v->getTargetDir() == '+')
        this->prevJuncCountIgnoreCopyPositive++;
        else this->prevJuncCountIgnoreCopyNegative++;
    }
}

void Vertex::setGeneAndScore(bool gene, float s) {
    this->containsGene = gene;
    this->score = s;
}

bool Vertex::isGeneAndScoreOk() {
    return this->containsGene && score > 0.7;
}
bool Vertex::isVisited() const {
    return visited;
}

void Vertex::setIsVisited(bool visited) {
    Vertex::visited = visited;
}

float Vertex::getScore() const {
    return score;
}

int Vertex::getCopyIdx() const {
    return copy_idx;
}

int Vertex::getNextJuncCountIgnoreCopyPositive() const {
    return nextJuncCountIgnoreCopyPositive;
}

int Vertex::getNextJuncCountIgnoreCopyNegative() const {
    return nextJuncCountIgnoreCopyNegative;
}

int Vertex::getPrevJuncCountIgnoreCopyPositive() const {
    return prevJuncCountIgnoreCopyPositive;
}

int Vertex::getPrevJuncCountIgnoreCopyNegative() const {
    return prevJuncCountIgnoreCopyNegative;
}
