//
// Created by caronkey on 17/5/2021.
//

#ifndef SEQMAP_EXCEPTIONS_H
#define SEQMAP_EXCEPTIONS_H


#include "Junction.h"
#include "EndPoint.h"
#include "Edge.h"

namespace seqGraph {
    class DuplicateJunctionException : public std::exception {
    private:
        Junction *mJunction;
        std::string whatMsg;
    public:
        DuplicateJunctionException(Junction *aJunction);

        virtual const char *what() const throw();
    };

    class VertexDoesNotExistException : public std::exception {
    private:
        std::string Id;
        std::string whatMsg;
    public:
        VertexDoesNotExistException(std::string& Id);

        virtual const char *what() const throw();
    };
}

#endif //SEQMAP_EXCEPTIONS_H
