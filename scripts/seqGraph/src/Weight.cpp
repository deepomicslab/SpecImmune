//
// Created by caronkey on 10/5/2021.
//

#include "../include/Weight.h"

using namespace seqGraph;

Weight::Weight(float aCoverage) {
    mCoverage = aCoverage;
    mCopyNum = 0;
    mCopyNumOriginal = 0;
    mCopyNumBackup = 0;

    mIsInferred = false;
}

Weight::~Weight() { ; }

float Weight::getCoverage() { return mCoverage; }

// float Weight::getOriginalCoverage() { return mCoverageOriginal; }
// float Weight::getAdjustedCoverage() { return mCoverageAdjusted; }
float Weight::getCopyNum() { return mCopyNum; }

float Weight::getCopyNumBackup() { return mCopyNumBackup; }

void Weight::setCoverage(float aCoverage) { mCoverage = aCoverage; }

// void Weight::setOriginalCoverage(float aCoverage) { mCoverageOriginal = aCoverage; }
// void Weight::setAdjustedCoverage(float aCoverage) { mCoverageAdjusted = aCoverage; }
void Weight::setCopyNum(float aCopyNum) { mCopyNum = mCopyNumOriginal = mCopyNumBackup = aCopyNum; }

void Weight::backup() { mCopyNumBackup = mCopyNum; }

void Weight::restore() { mCopyNum = mCopyNumBackup; }

void Weight::increaseCopyNum(float aIncrement) { mCopyNum += aIncrement; }

void Weight::decreaseCopyNum(float aDecrement) { mCopyNum -= aDecrement; }

bool Weight::isInferred() { return mIsInferred; }

void Weight::setInferred() { mIsInferred = true; }

void Weight::resetInferred() { mIsInferred = false; }

