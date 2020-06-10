//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "CommonEdge.h"

namespace component {

struct Block {
    Block() { inBlk.reserve(100); }
    std::vector<int>    inBlk;
    std::vector<int>    exBlk;
};

class CBF_Block {
    typedef std::vector<std::vector<std::vector<Block>>> Mtr3D;
    typedef std::vector<Block*>                          VecBlk;
public:
    CBF_Block();
    ~CBF_Block();
public:
    void            buildBlock(std::shared_ptr<Mesh> _pMesh, std::shared_ptr<CommonEdge> _pCE, int _maxNum, value_t _lam);
    const Block&    block(int x, int y, int z) const { return blkArray[x][y][z]; }
    const Block&    block(int _id) const { return *(vecBlock[_id]); }
    void            getxyzBlkNum(int& _x, int& _y, int& _z) const;
    int             getValidBlkNum() const { return static_cast<int>(vecBlock.size()); }
    void            reportInfo(Qostream& strm) const;
    bool            writeBlockInfo(const Qstring& dir) const;
    void            clear();
    //
    void            debug(std::shared_ptr<Mesh> _pMesh, std::shared_ptr<CommonEdge> _pCE, Qostream& strm) const;
private:
    bool            inBlock(const VectorR3& _point, const VectorR3& _min, const VectorR3& _max) const;
    bool            exInBlock(const VectorR3& _point, const VectorR3& _min, const VectorR3& _max) const;
private:
    value_t         blkSize, delta;
    int             xBlkNum, yBlkNum, zBlkNum;
    Mtr3D           blkArray;
    VecBlk          vecBlock;
};

} // namespace component