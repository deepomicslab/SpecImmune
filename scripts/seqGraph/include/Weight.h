//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_WEIGHT_H
#define SEQMAP_WEIGHT_H


namespace seqGraph {
    class Weight {
    protected:
        float mCoverage;           // coverage as adjusted by LP
        float mCopyNum;             // the copy number, coverage / average coverage of one copy
        float mCopyNumOriginal;
        float mCopyNumBackup;

        bool mIsInferred;

    public:
        Weight(float aCoverage);

        ~Weight();

        float getCoverage();

        float getCopyNum();

        float getCopyNumBackup();

        void setCoverage(float aCoverage);

        void setCopyNum(float aCopyNum);

        void backup();

        void restore();

        void increaseCopyNum(float aIncrement = 1);

        void decreaseCopyNum(float aDecrement = 1);

        bool isInferred();

        void setInferred();

        void resetInferred();

        void print();
    };
}


#endif //SEQMAP_WEIGHT_H
