#include "CBF_Block.h"

using namespace component;

CBF_Block::CBF_Block()
{
}


CBF_Block::~CBF_Block()
{
}

void CBF_Block::buildBlock(std::shared_ptr<Mesh> _pMesh, std::shared_ptr<CommonEdge> _pCE, int _maxNum, value_t _lam)
{
    delta = 0.1f * _lam;
    value_t xSize, ySize, zSize;
    _pMesh->getSize(xSize, ySize, zSize);

    VectorR3 min_Blk, max_Blk;
    _pMesh->getBoundary(min_Blk, max_Blk);

    value_t maxSize = (xSize > ySize) ? (xSize > zSize ? xSize : zSize) : (ySize > zSize ? ySize : zSize);
    blkSize = maxSize / _maxNum;

    xBlkNum = static_cast<int>(xSize / blkSize + 0.5);
    yBlkNum = static_cast<int>(ySize / blkSize + 0.5);
    zBlkNum = static_cast<int>(zSize / blkSize + 0.5);

    if (xBlkNum == 0)   ++xBlkNum;
    if (yBlkNum == 0)   ++yBlkNum;
    if (zBlkNum == 0)   ++zBlkNum;

    blkArray.resize(xBlkNum);
    for (int x = 0; x < xBlkNum; ++x)
    {
        blkArray[x].resize(yBlkNum);
        for (int y = 0; y < yBlkNum; ++y)
            blkArray[x][y].resize(zBlkNum);
    }

    auto unknowns = _pCE->getCommonEdgeNum();
    int nv1, nv2;
    VectorR3 v1, v2, mid, tmpMid;
    for (int n = 0; n < unknowns; ++n)
    {
        _pCE->getCommonEdge(n, nv1, nv2);
        v1 = _pMesh->getVertex(nv1);
        v2 = _pMesh->getVertex(nv2);
        mid = (v1 + v2) / 2;
        tmpMid = mid - min_Blk;

        int xID = static_cast<int>(tmpMid.x / blkSize);
        int yID = static_cast<int>(tmpMid.y / blkSize);
        int zID = static_cast<int>(tmpMid.z / blkSize);

        xID = (xID == xBlkNum) ? xID - 1 : xID;
        yID = (yID == yBlkNum) ? yID - 1 : yID;
        zID = (zID == zBlkNum) ? zID - 1 : zID;

        blkArray[xID][yID][zID].inBlk.push_back(n);

        VectorR3 minBlk(xID * blkSize, yID * blkSize, zID * blkSize);
        VectorR3 maxBlk       = minBlk + VectorR3(blkSize, blkSize, blkSize);
        VectorR3 insideMinBlk = minBlk + VectorR3(delta, delta, delta);
        VectorR3 insideMaxBlk = maxBlk - VectorR3(delta, delta, delta);

        if (inBlock(tmpMid, insideMinBlk, insideMaxBlk))
            continue;

        // process extend block
        int xset = (tmpMid.x - delta < xID*blkSize) ? -1 : ((tmpMid.x + delta >(xID + 1)*blkSize) ? 1 : 0);
        int yset = (tmpMid.y - delta < yID*blkSize) ? -1 : ((tmpMid.y + delta >(yID + 1)*blkSize) ? 1 : 0);
        int zset = (tmpMid.z - delta < zID*blkSize) ? -1 : ((tmpMid.z + delta >(zID + 1)*blkSize) ? 1 : 0);

        xset = (xID == 0 && xset == -1) ? 0 : ((xID == xBlkNum - 1 && xset == 1) ? 0 : xset);
        yset = (yID == 0 && yset == -1) ? 0 : ((yID == yBlkNum - 1 && yset == 1) ? 0 : yset);
        zset = (zID == 0 && zset == -1) ? 0 : ((zID == zBlkNum - 1 && zset == 1) ? 0 : zset);

        std::vector<int> ix = { 0 }, iy = { 0 }, iz = { 0 };
        if (xset != 0)  ix.push_back(xset);
        if (yset != 0)  iy.push_back(yset);
        if (zset != 0)  iz.push_back(zset);
    
        for (int i = 0; i < ix.size(); ++i)
        {
            for (int j = 0; j < iy.size(); ++j)
            {
                for (int k = 0; k < iz.size(); ++k)
                {
                    if (i == 0 && j == 0 && k == 0)
                        continue;
                    auto xindex = xID + ix[i];
                    auto yindex = yID + iy[j];
                    auto zindex = zID + iz[k];

                    VectorR3 minBlk(xindex * blkSize, yindex * blkSize, zindex * blkSize);
                    VectorR3 maxBlk       = minBlk + VectorR3(blkSize, blkSize, blkSize);
                    VectorR3 insideMinBlk = minBlk - VectorR3(delta, delta, delta);
                    VectorR3 insideMaxBlk = maxBlk + VectorR3(delta, delta, delta);

                    if (exInBlock(tmpMid, insideMinBlk, insideMaxBlk))
                        blkArray[xindex][yindex][zindex].exBlk.push_back(n);
                }
            }
        }
    }

    vecBlock.reserve(xBlkNum*yBlkNum*zBlkNum);
    for (int x = 0; x < xBlkNum; ++x)
    {
        for (int y = 0; y < yBlkNum; ++y)
        {
            for (int z = 0; z < zBlkNum; ++z)
            {
                if (!blkArray[x][y][z].inBlk.empty())
                    vecBlock.push_back(&blkArray[x][y][z]);
            }
        }
    }
    vecBlock.shrink_to_fit();
}

void CBF_Block::getxyzBlkNum(int & _x, int & _y, int & _z) const
{
    _x = xBlkNum;
    _y = yBlkNum;
    _z = zBlkNum;
}

void CBF_Block::reportInfo(Qostream & strm) const
{
    int validBlocks = 0, validUnknowns = 0;
    for (int x = 0; x < xBlkNum; ++x)
    {
        for (int y = 0; y < yBlkNum; ++y)
        {
            for (int z = 0; z < zBlkNum; ++z)
            {
                if (!blkArray[x][y][z].inBlk.empty())
                {
                    ++validBlocks;
                    validUnknowns += static_cast<int>(blkArray[x][y][z].inBlk.size());
                }
            }
        }
    }
    auto oldState = strm.flags();
    strm << std::left << HEADING "CBFM Block Information:\n" TRAILING;
    strm << LEVEL1 "Block size:" << blkSize << "m\n"
         << LEVEL1 "Extend size:" << delta << "m\n"
         << LEVEL1 "Total blocks:" << xBlkNum * yBlkNum * zBlkNum << '\n'
         << LEVEL1 "Valid blocks:" << validBlocks << '\n'
         << LEVEL1 "Valid unknowns:" << validUnknowns << '\n';
    strm.flags(oldState);
}

bool CBF_Block::writeBlockInfo(const Qstring & dir) const
{
    Qofstream output(dir + "/cbf_block.txt");
    if (output.fail())
        return false;

    for (int x = 0; x < xBlkNum; ++x)
        for (int y = 0; y < yBlkNum; ++y)
            for (int z = 0; z < zBlkNum; ++z)
                if (!blkArray[x][y][z].inBlk.empty())
                    output << "[" << x << "][" << y << "][" << z << "] inside:" << std::setw(5) 
                            << blkArray[x][y][z].inBlk.size() <<"  extend:" << std::setw(5)
                            << blkArray[x][y][z].exBlk.size() << "\n";


    int coupleBlocks = 0;
    int validBlocks = 0;
    for (int fx = 0; fx < xBlkNum; ++fx)
    {
        for (int fy = 0; fy < yBlkNum; ++fy)
        {
            for (int fz = 0; fz < zBlkNum; ++fz)
            {
                if (blkArray[fx][fy][fz].inBlk.empty())
                    continue;
                ++validBlocks;
                for (int sx = 0; sx < xBlkNum; ++sx)
                {
                    for (int sy = 0; sy < yBlkNum; ++sy)
                    {
                        for (int sz = 0; sz < zBlkNum; ++sz)
                        {
                            if (blkArray[sx][sy][sz].inBlk.empty())
                                continue;
                            if ((abs(fx - sx) < 1) && (abs(fy - sy) < 1) && (abs(fz - sz) < 1))
                                continue;

                            output << "[" << fx << "][" << fy << "][" << fz << "]<->["
                                << sx << "][" << sy << "][" << sz << "]\n";
                            ++coupleBlocks;
                        }
                    }
                }
            }
        }
    }
    output  << "-total blocks:" << xBlkNum * yBlkNum * zBlkNum << '\n'
            << "-valid blocks:" << validBlocks << '\n'
            << "-number of coupling blocks: " << coupleBlocks << std::endl;
    return true;
}

void CBF_Block::clear()
{
    blkArray.clear();
}

void CBF_Block::debug(std::shared_ptr<Mesh> _pMesh, std::shared_ptr<CommonEdge> _pCE, Qostream & strm) const
{
    VectorR3 min_Blk, max_Blk;
    _pMesh->getBoundary(min_Blk, max_Blk);

    int unknowns = _pCE->getCommonEdgeNum();
    int nv1, nv2;
    VectorR3 v1, v2, mid, tmpMid;
    std::vector<int> blkUnknowns;
    for (int x = 0; x < xBlkNum; ++x)
    {
        for (int y = 0; y < yBlkNum; ++y)
        {
            for (int z = 0; z < zBlkNum; ++z)
            {
                VectorR3 minBlk(x * blkSize, y * blkSize, z * blkSize);
                VectorR3 maxBlk = minBlk + VectorR3(blkSize, blkSize, blkSize);
                VectorR3 exMinBlk = minBlk - VectorR3(delta, delta, delta);
                VectorR3 exMaxBlk = maxBlk + VectorR3(delta, delta, delta);
                int curUnknowns = 0;
                for (int u = 0; u < unknowns; ++u)
                {
                    _pCE->getCommonEdge(u, nv1, nv2);
                    v1 = _pMesh->getVertex(nv1);
                    v2 = _pMesh->getVertex(nv2);
                    mid = (v1 + v2) / 2;
                    tmpMid = mid - min_Blk;
                    if (exInBlock(tmpMid, exMinBlk, exMaxBlk))
                        ++curUnknowns;
                }
                blkUnknowns.push_back(curUnknowns);
            }
        }
    }

    for (int n = 0; n < blkUnknowns.size(); ++n)
    {
        strm << n << ": " << blkUnknowns[n] << '\n';
    }
    strm << std::endl;
}


bool CBF_Block::inBlock(const VectorR3 & _point, const VectorR3 & _min, const VectorR3 & _max) const
{
    bool xBool = _point.x < _min.x ? false : (_point.x > _max.x ? false : true);
    bool yBool = _point.y < _min.y ? false : (_point.y > _max.y ? false : true);
    bool zBool = _point.z < _min.z ? false : (_point.z > _max.z ? false : true);

    return xBool && yBool && zBool;
}

bool CBF_Block::exInBlock(const VectorR3 & _point, const VectorR3 & _min, const VectorR3 & _max) const
{
    bool xBool = _point.x <= _min.x ? false : (_point.x >= _max.x ? false : true);
    bool yBool = _point.y <= _min.y ? false : (_point.y >= _max.y ? false : true);
    bool zBool = _point.z <= _min.z ? false : (_point.z >= _max.z ? false : true);

    return xBool && yBool && zBool;
}
