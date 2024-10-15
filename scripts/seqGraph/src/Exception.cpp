//
// Created by caronkey on 17/5/2021.
//

#include "../include/Exceptions.h"

using namespace std;
using namespace seqGraph;

/* duplicate junction */
DuplicateJunctionException::DuplicateJunctionException(Junction *aJunction) {
    mJunction = aJunction;
    whatMsg = "DuplicateJunctionException: " + mJunction->getInfo();
}

const char *DuplicateJunctionException::what() const throw() {
    return whatMsg.c_str();
}

/* cannot find Vertex */
VertexDoesNotExistException::VertexDoesNotExistException(std::string& aSegId) {
    Id = aSegId;
    whatMsg = "VertexDoesNotExistException: Vertex with ID " + Id + " does not exist";
}

const char *VertexDoesNotExistException::what() const throw() {
    return whatMsg.c_str();
}
