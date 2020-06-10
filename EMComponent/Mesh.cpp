#include "Mesh.h"
#include "Mathdef.h"
#include "MyDebugNew.h"
#include "ErrorHandle.h"
#include "tools.h"

#define THROW_ERROR(str) throw component::FileError(str)

using namespace component;
using std::setw;

Mesh::Mesh()
{
}

Mesh::~Mesh()
{
    if (!triangles_.empty())
        clear();
}

void Mesh::loadMeshFile(const Qstring & filename)
{
    std::string suffix = filename.substr(filename.find_last_of('.') + 1);
    if ("rpt" == suffix)
    {
        if (!loadRptFile(filename))
            THROW_ERROR("fail loading rpt mesh file: " + filename);

        if (!readRptTmpFile())
            THROW_ERROR("fail reading rpt temp file");
    }
    else if ("stl" == suffix)
    {
        if (!loadStlFile(filename))
            THROW_ERROR("fail loading stl mesh file: " + filename);
        if (!readStlTmpFile())
            THROW_ERROR("fail reading stl temp file");
    }
    else if ("nas" == suffix)
    {
        if (!loadNasFile(filename))
            THROW_ERROR("fail loading nas mesh file: " + filename);
        if (!readNasTmpFile())
            THROW_ERROR("fail reading nas temp file");
    }
    else
        THROW_ERROR("cannot identify mesh file type: " + suffix + " (1.rpt, 2.stl(text), 3.nas)");
}

void Mesh::loadMeshForEFIE_PMCHW(const Qstring & filename)
{
    std::string map_file = filename.substr(0, filename.find_last_of('.'));
    map_file += ".map";
    if (!prepareForEFIE_PMCHW(map_file))
        THROW_ERROR("fail loading map(EFIE_PMCHW) file: " + map_file);
    if (!loadNasFileForEFIE_PMCHW(filename))
        THROW_ERROR("fail loading nas(EFIE_PMCHW) mesh file: " + filename);
    if (!readNasTmpFileForEFIE_PMCHW())
        THROW_ERROR("fail reading nas(EFIE_PMCHW) temp file");
}

void Mesh::loadDDMMesh(const Qstring & filename)
{
    std::string map_file = filename.substr(0, filename.find_last_of('.')) + ".map";
    if (!loadDDMMapFile(map_file))
        THROW_ERROR("fail loading map(DDM) file: " + map_file);
    if (!loadDDMNasFile(filename))
        THROW_ERROR("fail loading nas(DDM) mesh file: " + filename);
    if (!readDDMNasTmpFile())
        THROW_ERROR("fail reading nas(DDM) temp file");
}



void Mesh::loadVIEMesh(const Qstring & filename)
{
	std::string map_file = filename.substr(0, filename.find_last_of('.')) + ".map";
	if(!loadVIEMapFile(map_file))
		THROW_ERROR("fail loading map(VIE) file: " + map_file);
	if (!loadVIENasFile(filename))
		THROW_ERROR("fail loading nas(VIE) mesh file: " + filename);
	if (!readVIENasTmpFile())
		THROW_ERROR("fail reading nas(VIE) temp file");
}

void Mesh::loadVSIEMesh(const Qstring& filename)
{
	std::string map_file = filename.substr(0, filename.find_last_of('.')) + ".map";
	if (!loadVSIEMapFile(map_file))
		THROW_ERROR("fail loading map(VIE) file: " + map_file);
	if (!loadVSIENasFile(filename))
		THROW_ERROR("fail loading nas(VIE) mesh file: " + filename);
	if (!readVSIENasTmpFile())
		THROW_ERROR("fail reading nas(VIE) temp file");
}

void Mesh::loadIBCMesh(const Qstring& filename)
{
	std::string map_file = filename.substr(0, filename.find_last_of('.')) + ".map";
	if (!loadIBCMapFile(map_file))
		THROW_ERROR("fail loading map(TDS) file: " + map_file);
	if (!loadIBCNasFile(filename))
		THROW_ERROR("fail loading nas(TDS) mesh file: " + filename);
	if (!readIBCNasTmpFile())
		THROW_ERROR("fail reading nas(TDS) temp file");
}


void Mesh::reportInfo(Qostream& strm) const
{
    VectorR3 min    = boundary_.getBoundaryMin();
    VectorR3 max    = boundary_.getBoundaryMax();
    VectorR3 bSize  = boundary_.getBoundarySize();
    size_t vmem = sizeof(VectorR3) * vertices_.size() / (1024 * 1024);
    size_t tmem = sizeof(Triangle) * triangles_.size() / (1024 * 1024);
    auto oldState = strm.flags();
    strm << std::left << std::fixed << std::setprecision(2);
    strm << HEADING "Mesh Information\n" TRAILING;
    strm << LEVEL1 "Vertices: " << vertices_.size() << '\n';
    strm << LEVEL1 "Triangles: " << triangles_.size() << '\n';
    strm << LEVEL1 "Minimum position:"    << '(' << min.x << ',' << min.y
                                                << ',' << min.z << ')' << '\n';
    strm << LEVEL1 "Maximum position:" << '(' << max.x << ',' << max.y
                                                << ',' << max.z << ')' << '\n';
    strm << LEVEL1 "Mesh size: " << bSize.x << 'x' << bSize.y << 'x' << bSize.z << '\n';
    strm << LEVEL1 "Memory detail:" << '\n'
        << LEVEL2 "vertex:" << vmem << " MB" << FORMAT_MEMORY(vmem) << '\n'
        << LEVEL2 "triangle:" << tmem << " MB" << FORMAT_MEMORY(tmem) << '\n';
    strm.flush();
    strm.flags(oldState);
}

void Mesh::clear()
{
    common_face_.clear();
    subdomain_id_.clear();
    triangles_.clear();
    vertices_.clear();
}

void Mesh::Debug()
{
    auto count = triangles_.size();
    std::vector<int> sub_count(subdomain_num_, 0);
    for (size_t i = 0; i < count; ++i)
    {
        ++sub_count[triangles_[i].getID()];
    }
    Qcout << "Number of subdomain : " << subdomain_num_ << '\n';
    Qcout << "Number of triangles in each domain\n";
    for (int i = 0;i<sub_count.size();++i)
    {
        Qcout << "  --subdomain " << i + 1 << " : " << std::setw(10) << sub_count[i] << '\n';
    }
    Qcout.flush();
}

bool Mesh::prepareForEFIE_PMCHW(const Qstring & filename)
{
    common_face_.clear();
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    Qifstream   input_stream(filename);
    if (input_stream.fail())
    {
        return false;
    }
    int faceID;
    while (input_stream >> faceID)
    {
        common_face_.insert(faceID);
        input_stream.ignore(maxSize, '\n');
    }
    return true;
}


bool Mesh::loadRptFile(const Qstring & filename)
{
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    std::stringstream tmpInfo;
    Qifstream input_stream(filename, std::ios::in);
    if (input_stream.fail())
        return false;

    input_stream.ignore(maxSize, '\n');
    input_stream.ignore(maxSize, '\n');
    
    int  triangleNum = 0, nodeNum = 0;
    int  triID, pos1, pos2, pos3;
    while (input_stream >> triID >> pos1 >> pos2 >> pos3)
    {
        tmpInfo << triID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\n';
        ++triangleNum;
    }

    input_stream.clear();
    input_stream.ignore(maxSize, '\n');
    input_stream.ignore(maxSize, '\n');

    int         nodeID;
    value_t     coord1, coord2, coord3;
    while (input_stream >> nodeID >> coord1 >> coord2 >> coord3 >> pos1 >> pos1)
    {
        tmpInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
        ++nodeNum;
    }
    
    if (!input_stream.eof())
        return false;

    Qofstream tmpFile("temp.dat", std::ios::out);
    if (tmpFile.fail())
        return false;
    tmpFile << triangleNum << '\t' << nodeNum << '\n';
    tmpFile << tmpInfo.rdbuf();

    return true;
}

bool Mesh::loadStlFile(const Qstring & filename)
{
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    std::stringstream tmpTriInfo;
    Qifstream input_stream(filename, std::ios::in);
    if (input_stream.fail())
        return false;
    input_stream.ignore(maxSize, '\n');

    int triNum = 0;
    value_t coord1, coord2, coord3;
    
    while (true)
    {
        input_stream.ignore(maxSize, 'l');
        input_stream >> coord1 >> coord2 >> coord3;
        if (input_stream.fail())   break;
        ++triNum;
        tmpTriInfo << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';

        for (int v = 0; v < 3; ++v)
        {
            input_stream.ignore(maxSize, 'x');
            input_stream >> coord1 >> coord2 >> coord3;
            tmpTriInfo << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
        }

        tmpTriInfo << '\n';
        input_stream.ignore(maxSize, 'l');
    }

    Qofstream tmpFile("temp.dat", std::ios::out);
    if (tmpFile.fail())
        return false;
    tmpFile << triNum << '\n';
    tmpFile << tmpTriInfo.rdbuf();

    return true;
}

bool Mesh::loadNasFile(const Qstring & filename)
{
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    std::stringstream tmpNodeInfo, tmpTriInfo;
    Qifstream input_stream(filename, std::ios::in);
    if (input_stream.fail())
        return false;

    for (int t = 0; t < 9; ++t)
        input_stream.ignore(maxSize, '\n');
    int     nodeID, nodeNum = 0;
    value_t coord1, coord2, coord3;
    int triID, triNum = 0;
    int pos1, pos2, pos3, groupID;

    while (input_stream.peek() != 'E')
    {
        while (true)
        {
            if (input_stream.peek() != 'G')    break;
            input_stream.ignore(maxSize, '*');
            input_stream >> nodeID >> coord1 >> coord2;
            input_stream.ignore(maxSize, '*');
            input_stream >> nodeID >> coord3;
            tmpNodeInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
            if (nodeID != ++nodeNum)
                return false;
            input_stream.get();
        }

        while (true)
        {
            if (input_stream.peek() != 'C') break;
            input_stream.ignore(maxSize, ' ');
            input_stream >> triID >> groupID >> pos1 >> pos2 >> pos3;
            tmpTriInfo << triID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\n';
            if (triID != ++triNum)
                return false;
            input_stream.get();
        }
    }

    Qofstream tmpFile("temp.dat", std::ios::out);
    if (tmpFile.fail())
        return false;
    tmpFile << triNum << '\t' << nodeNum << '\n';
    tmpFile << tmpTriInfo.rdbuf() << tmpNodeInfo.rdbuf();
    return true;
}

bool Mesh::loadNasFileForEFIE_PMCHW(const Qstring & filename)
{
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    std::stringstream tmpNodeInfo, tmpTriInfo;
    Qifstream input_stream(filename, std::ios::in);
    if (input_stream.fail())
        return false;

    for (int t = 0; t < 9; ++t)
        input_stream.ignore(maxSize, '\n');
    int     nodeID, nodeNum = 0;
    value_t coord1, coord2, coord3;
    while (true)
    {
        if (input_stream.peek() == 'C')    break;
        input_stream.ignore(maxSize, '*');
        input_stream >> nodeID >> coord1 >> coord2;
        input_stream.ignore(maxSize, '*');
        input_stream >> nodeID >> coord3;
        tmpNodeInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
        if (nodeID != ++nodeNum)
            return false;
        input_stream.get();
    }

    int triID, triNum = 0;
    int pos1, pos2, pos3, groupID;
    input_stream.ignore(maxSize, ' ');
    while (input_stream >> triID >> groupID >> pos1 >> pos2 >> pos3)
    {
        int isCmface = (common_face_.find(groupID) != common_face_.cend()) ? 1 : 0;
        tmpTriInfo << triID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\t' << isCmface << '\n';
        if (triID != ++triNum)
            return false;
        input_stream.ignore(maxSize, ' ');
    }

    if (!input_stream.eof())
        return false;

    Qofstream tmpFile("temp.dat", std::ios::out);
    if (tmpFile.fail())
        return false;
    tmpFile << triNum << '\t' << nodeNum << '\n';
    tmpFile << tmpTriInfo.rdbuf() << tmpNodeInfo.rdbuf();
    return true;
}

bool Mesh::loadDDMMapFile(const Qstring & filename)
{
    Qifstream input(filename, std::ios::in);
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    int id = 0, subNo = 0;
    if (input.fail())
        return false;
    if (input.peek() == '*')
    {
        input.ignore(maxSize, '\n');
    }
    while (input >> subNo)
    {
        subdomain_id_[subNo] = id;
        input.ignore(maxSize, '\n');
        if (input.peek() == '*')
        {
            input.ignore(maxSize, '\n');
            ++id;
        }
    }
    subdomain_num_ = id + 1;
    return true;
}

bool Mesh::loadDDMNasFile(const Qstring & filename)
{
    auto maxSize = std::numeric_limits<std::streamsize>::max();
    std::stringstream tmpNodeInfo, tmpTriInfo;
    Qifstream input_stream(filename, std::ios::in);
    if (input_stream.fail())
        return false;

    for (int t = 0; t < 9; ++t)
        input_stream.ignore(maxSize, '\n');
    int     nodeID, nodeNum = 0;
    value_t coord1, coord2, coord3;

    int triID, triNum = 0;
    int pos1, pos2, pos3, groupID;

    while (input_stream.peek() != 'E')
    {
        while (true)
        {
            if (input_stream.peek() != 'G')    break;
            input_stream.ignore(maxSize, '*');
            input_stream >> nodeID >> coord1 >> coord2;
            input_stream.ignore(maxSize, '*');
            input_stream >> nodeID >> coord3;
            tmpNodeInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
            if (nodeID != ++nodeNum)
                return false;
            input_stream.get();
        }

        while (true)
        {
            if (input_stream.peek() != 'C') break;
            input_stream.ignore(maxSize, ' ');
            input_stream >> triID >> groupID >> pos1 >> pos2 >> pos3;
            tmpTriInfo << triID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\t' << subdomain_id_[groupID] << '\n';
            if (triID != ++triNum)
                return false;
            input_stream.get();
        }
    }

    Qofstream tmpFile("temp.dat", std::ios::out);
    if (tmpFile.fail())
        return false;
    tmpFile << triNum << '\t' << nodeNum << '\n';
    tmpFile << tmpTriInfo.rdbuf() << tmpNodeInfo.rdbuf();
    return true;
}

bool Mesh::loadVIENasFile(const Qstring& filename)
{
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	std::stringstream tmpNodeInfo, tmpTetInfo;
	Qifstream input_stream(filename, std::ios::in);
	if (input_stream.fail())
		return false;

	for (int t = 0; t < 9; ++t)
		input_stream.ignore(maxSize, '\n');
	int     nodeID, nodeNum = 0;
	value_t coord1, coord2, coord3;
	int tetID, tetNum = 0;
	int pos1, pos2, pos3, pos4, groupID;

	while (input_stream.peek() != 'E')
	{
		while (true)
		{
			if (input_stream.peek() != 'G')    break;
			input_stream.ignore(maxSize, '*');
			input_stream >> nodeID >> coord1 >> coord2;
			input_stream.ignore(maxSize, '*');
			input_stream >> nodeID >> coord3;
			tmpNodeInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
			if (nodeID != ++nodeNum)
				return false;
			input_stream.get();
		}

		while (true)
		{
			if (input_stream.peek() != 'C') break;
			input_stream.ignore(maxSize, ' ');
			input_stream >> tetID >> groupID >> pos1 >> pos2 >> pos3 >> pos4;
			tmpTetInfo << tetID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\t' << pos4 <<'\t'<<subdomain_eps_[groupID]<< '\n';
			if (tetID != ++tetNum)
				return false;
			input_stream.get();
		}
	}

	Qofstream tmpFile("temp.dat", std::ios::out);
	if (tmpFile.fail())
		return false;
	tmpFile << tetNum << '\t' << nodeNum << '\n';
	tmpFile << tmpTetInfo.rdbuf() << tmpNodeInfo.rdbuf();
	return true;
}

bool Mesh::loadVIEMapFile(const Qstring& filename)
{
	Qifstream input(filename, std::ios::in);
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	Qstring domainname; int groupID; value_t eps;
	
	/*if (input.peek() == 'B')
	{
		input.ignore(maxSize, '\n');
	}*/
	while (input.peek() != 'E')
	{
		input >> groupID >> domainname >> eps;
		input.get();
		subdomain_eps_[groupID] = eps;
	}
	return true;
}

bool Mesh::loadVSIEMapFile(const Qstring& filename)
{
	Qifstream input(filename, std::ios::in);
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	Qstring domainname; int groupID; value_t eps;

	/*if (input.peek() == 'B')
	{
	input.ignore(maxSize, '\n');
	}*/
	while (input.peek() != 'E')
	{
		input >> groupID >> domainname >> eps;
		input.get();
		subdomain_eps_[groupID] = eps;
	}
	return true;
}

bool Mesh::loadVSIENasFile(const Qstring& filename)
{
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	std::stringstream tmpNodeInfo, tmpTetInfo, tmpTriInfo;
	Qifstream input_stream(filename, std::ios::in);
	if (input_stream.fail())
		return false;

	for (int t = 0; t < 6; ++t)
		input_stream.ignore(maxSize, '\n');

	input_stream.ignore(maxSize, ':');
	input_stream >> tri_number;
	input_stream.get();
	input_stream.ignore(maxSize, '\n');
	input_stream.ignore(maxSize, ':');
	input_stream >> tet_number;
	input_stream.get();

	int     nodeID, nodeNum = 0;
	value_t coord1, coord2, coord3;
	int tetID=0, tetNum = 0;
	int triID=0, triNum = 0;
	int pos1, pos2, pos3, pos4, groupID;

	while (input_stream.peek() != 'E')
	{
		while (true)
		{
			if (input_stream.peek() != 'G')    break;
			input_stream.ignore(maxSize, '*');
			input_stream >> nodeID >> coord1 >> coord2;
			input_stream.ignore(maxSize, '*');
			input_stream >> nodeID >> coord3;
			tmpNodeInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
			if (nodeID != ++nodeNum)
				return false;
			input_stream.get();
		}

		while (triID != tri_number)
		{
			if (input_stream.peek() != 'C') break;
			input_stream.ignore(maxSize, ' ');
			input_stream >> triID >> groupID >> pos1 >> pos2 >> pos3;
			tmpTriInfo << triID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\t' << groupID << '\n';
			//if (triID == tri_number) break;
			if (triID != ++triNum)
				return false;
			input_stream.get();
		}
		//input_stream.get();
		while (tetID!=tet_number)
		{
			if (input_stream.peek() != 'C') break;
			input_stream.ignore(maxSize, ' ');
			input_stream >> tetID >> groupID >> pos1 >> pos2 >> pos3 >> pos4;
			tmpTetInfo << tetID << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\t' << pos4 << '\t' << subdomain_eps_[groupID] << '\n';
			if (tetID != ++tetNum)
				return false;
			input_stream.get();
		}
	}

	Qofstream tmpFile("temp.dat", std::ios::out);
	if (tmpFile.fail())
		return false;
	tmpFile << tri_number <<'\t'<< tet_number << '\t' << nodeNum << '\n';
	if (tri_number == 0)
	{
		tmpFile << tmpTetInfo.rdbuf() << tmpNodeInfo.rdbuf();
	}
	else if (tet_number == 0)
	{
		tmpFile << tmpTriInfo.rdbuf() << tmpNodeInfo.rdbuf();
	}
	else
	{
		tmpFile << tmpTriInfo.rdbuf() << tmpTetInfo.rdbuf() << tmpNodeInfo.rdbuf();
	}
	return true;
}

bool Mesh::loadIBCMapFile(const Qstring& filename)
{
	Qifstream input(filename, std::ios::in);
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	Qstring domainname; int groupID; Qstring fla;

	while (input.peek() != 'E')
	{
		input >> groupID >> domainname >> fla;
		input.get();
		subdomain_fla_[groupID] = fla;
	}

	return true;
}

bool Mesh::loadIBCNasFile(const Qstring& filename)
{
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	std::stringstream tmpNodeInfo, tmpTriInfo;
	Qifstream input_stream(filename, std::ios::in);
	if (input_stream.fail())
		return false;

	for (int t = 0; t < 9; ++t)
		input_stream.ignore(maxSize, '\n');
	int     nodeID, nodeNum = 0;
	value_t coord1, coord2, coord3;
	int triID, triNum = 0;
	int pos1, pos2, pos3, groupID;

	while (input_stream.peek() != 'E')
	{
		while (true)
		{
			if (input_stream.peek() != 'G')    break;
			input_stream.ignore(maxSize, '*');
			input_stream >> nodeID >> coord1 >> coord2;
			input_stream.ignore(maxSize, '*');
			input_stream >> nodeID >> coord3;
			tmpNodeInfo << nodeID << '\t' << coord1 << '\t' << coord2 << '\t' << coord3 << '\n';
			if (nodeID != ++nodeNum)
				return false;
			input_stream.get();
		}

		while (true)
		{
			if (input_stream.peek() != 'C') break;
			input_stream.ignore(maxSize, ' ');
			input_stream >> triID >> groupID >> pos1 >> pos2 >> pos3;
			tmpTriInfo << triID  << '\t' << pos1 << '\t' << pos2 << '\t' << pos3 << '\t' << subdomain_fla_[groupID] << '\n';
			if (triID != ++triNum)
				return false;
			input_stream.get();
		}
	}

	Qofstream tmpFile("temp.dat", std::ios::out);
	if (tmpFile.fail())
		return false;
	tmpFile << triNum << '\t' << nodeNum << '\n';
	tmpFile << tmpTriInfo.rdbuf() << tmpNodeInfo.rdbuf();

	return true;
}

bool Mesh::readRptTmpFile()
{
    Qifstream input("temp.dat", std::ios::in);
    if (input.fail())
        return false;
    int triNum, nodeNum;
    input >> triNum >> nodeNum;

    triangles_.reserve(triNum);
    vertices_.reserve(nodeNum);

    int  triID, pos1, pos2, pos3;
    Triangle tri;
    for (int i = 0; i < triNum; ++i)
    {
        input >> triID >> pos1 >> pos2 >> pos3;
        tri.setVertex(pos1 - 1, pos2 - 1, pos3 - 1);
        triangles_.push_back(tri);
    }
    
    VectorR3 node;
    VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
    VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

    int nodeID;
    for (int i = 0; i < nodeNum; ++i)
    {
        input >> nodeID >> node.x >> node.y >> node.z;
        vertices_.push_back(std::move(node));
        
        minBox.x = Min(minBox.x, node.x);
        minBox.y = Min(minBox.y, node.y);
        minBox.z = Min(minBox.z, node.z);

        maxBox.x = Max(maxBox.x, node.x);
        maxBox.y = Max(maxBox.y, node.y);
        maxBox.z = Max(maxBox.z, node.z);
    }
    
    input >> nodeID;
    if (!input.eof())
        return false;

    boundary_.setBoundary(minBox, maxBox);

    for (int t = 0; t < triNum; ++t)
    {
        triangles_[t].getVertex(pos1, pos2, pos3);
        triangles_[t].setVertex(vertices_[pos1], vertices_[pos2], vertices_[pos3]);
    }

    return true;
}

bool Mesh::readNasTmpFile()
{
    Qifstream input("temp.dat", std::ios::in);
    if (input.fail())
        return false;
    int triNum, nodeNum;
    input >> triNum >> nodeNum;

    triangles_.reserve(triNum);
    vertices_.reserve(nodeNum);

    int  triID, pos1, pos2, pos3;
    Triangle tri;
    for (int i = 0; i < triNum; ++i)
    {
        input >> triID >> pos1 >> pos2 >> pos3;
        tri.setVertex(pos1 - 1, pos2 - 1, pos3 - 1);
        triangles_.push_back(tri);
    }

    VectorR3 node;
    VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
    VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

    int nodeID;
    for (int i = 0; i < nodeNum; ++i)
    {
        input >> nodeID >> node.x >> node.y >> node.z;
        vertices_.push_back(std::move(node));

        minBox.x = Min(minBox.x, node.x);
        minBox.y = Min(minBox.y, node.y);
        minBox.z = Min(minBox.z, node.z);

        maxBox.x = Max(maxBox.x, node.x);
        maxBox.y = Max(maxBox.y, node.y);
        maxBox.z = Max(maxBox.z, node.z);
    }

    input >> nodeID;
    if (!input.eof())
        return false;

    boundary_.setBoundary(minBox, maxBox);

    for (int t = 0; t < triNum; ++t)
    {
        triangles_[t].getVertex(pos1, pos2, pos3);
        triangles_[t].setVertex(vertices_[pos1], vertices_[pos2], vertices_[pos3]);
    }

    return true;
}

bool Mesh::readStlTmpFile()
{
    Qifstream input("temp.dat", std::ios::in);
    if (input.fail())
        return false;
    int triNum;
    input >> triNum;
    triangles_.reserve(triNum);

    VectorR3 unitNormal, nodeA, nodeB, nodeC;
    VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
    VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);
    Triangle tri;
    for (int t = 0; t < triNum; ++t)
    {
        input >> unitNormal.x >> unitNormal.y >> unitNormal.z;
        input >> nodeA.x >> nodeA.y >> nodeA.z;
        minBox.x = Min(minBox.x, nodeA.x);
        minBox.y = Min(minBox.y, nodeA.y);
        minBox.z = Min(minBox.z, nodeA.z);

        maxBox.x = Max(maxBox.x, nodeA.x);
        maxBox.y = Max(maxBox.y, nodeA.y);
        maxBox.z = Max(maxBox.z, nodeA.z);

        input >> nodeB.x >> nodeB.y >> nodeB.z;
        minBox.x = Min(minBox.x, nodeB.x);
        minBox.y = Min(minBox.y, nodeB.y);
        minBox.z = Min(minBox.z, nodeB.z);

        maxBox.x = Max(maxBox.x, nodeB.x);
        maxBox.y = Max(maxBox.y, nodeB.y);
        maxBox.z = Max(maxBox.z, nodeB.z);

        input >> nodeC.x >> nodeC.y >> nodeC.z;
        minBox.x = Min(minBox.x, nodeC.x);
        minBox.y = Min(minBox.y, nodeC.y);
        minBox.z = Min(minBox.z, nodeC.z);

        maxBox.x = Max(maxBox.x, nodeC.x);
        maxBox.y = Max(maxBox.y, nodeC.y);
        maxBox.z = Max(maxBox.z, nodeC.z);

        tri.setVertex(nodeA, nodeB, nodeC, unitNormal);
        triangles_.push_back(tri);
    }

    input >> triNum;
    if (!input.eof())
        return false;

    boundary_.setBoundary(minBox, maxBox);

    return true;
}

bool Mesh::readNasTmpFileForEFIE_PMCHW()
{
    Qifstream input("temp.dat", std::ios::in);
    if (input.fail())
        return false;
    int triNum, nodeNum;
    input >> triNum >> nodeNum;

    triangles_.reserve(triNum);
    vertices_.reserve(nodeNum);

    int  triID, pos1, pos2, pos3, isCmFace;
    Triangle tri;
    for (int i = 0; i < triNum; ++i)
    {
        input >> triID >> pos1 >> pos2 >> pos3 >> isCmFace;
        tri.setVertex(pos1 - 1, pos2 - 1, pos3 - 1);
        tri.setIsCommonFace(isCmFace);
        triangles_.push_back(tri);
    }

    VectorR3 node;
    VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
    VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

    int nodeID;
    for (int i = 0; i < nodeNum; ++i)
    {
        input >> nodeID >> node.x >> node.y >> node.z;
        vertices_.push_back(std::move(node));

        minBox.x = Min(minBox.x, node.x);
        minBox.y = Min(minBox.y, node.y);
        minBox.z = Min(minBox.z, node.z);

        maxBox.x = Max(maxBox.x, node.x);
        maxBox.y = Max(maxBox.y, node.y);
        maxBox.z = Max(maxBox.z, node.z);
    }

    input >> nodeID;
    if (!input.eof())
        return false;

    boundary_.setBoundary(minBox, maxBox);

    for (int t = 0; t < triNum; ++t)
    {
        triangles_[t].getVertex(pos1, pos2, pos3);
        triangles_[t].setVertex(vertices_[pos1], vertices_[pos2], vertices_[pos3]);
    }
    return true;
}

bool Mesh::readDDMNasTmpFile()
{
    Qifstream input("temp.dat", std::ios::in);
    if (input.fail())
        return false;
    int triNum, nodeNum;
    input >> triNum >> nodeNum;

    triangles_.reserve(triNum);
    vertices_.reserve(nodeNum);

    int  triID, pos1, pos2, pos3, subID;
    Triangle tri;
    for (int i = 0; i < triNum; ++i)
    {
        input >> triID >> pos1 >> pos2 >> pos3 >> subID;
        tri.setVertex(pos1 - 1, pos2 - 1, pos3 - 1);
        tri.setID(subID);
        triangles_.push_back(tri);
    }

    VectorR3 node;
    VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
    VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

    int nodeID;
    for (int i = 0; i < nodeNum; ++i)
    {
        input >> nodeID >> node.x >> node.y >> node.z;
        vertices_.push_back(std::move(node));

        minBox.x = Min(minBox.x, node.x);
        minBox.y = Min(minBox.y, node.y);
        minBox.z = Min(minBox.z, node.z);

        maxBox.x = Max(maxBox.x, node.x);
        maxBox.y = Max(maxBox.y, node.y);
        maxBox.z = Max(maxBox.z, node.z);
    }

    input >> nodeID;
    if (!input.eof())
        return false;

    boundary_.setBoundary(minBox, maxBox);

    for (int t = 0; t < triNum; ++t)
    {
        triangles_[t].getVertex(pos1, pos2, pos3);
        triangles_[t].setVertex(vertices_[pos1], vertices_[pos2], vertices_[pos3]);
    }

    return true;
}

bool Mesh::readVIENasTmpFile()
{
	Qifstream input("temp.dat", std::ios::in);
	if (input.fail())
		return false;
	int tetNum, nodeNum;
	input >> tetNum >> nodeNum;
	
	tetrahedron_.reserve(tetNum);
	vertices_.reserve(nodeNum);

	int  tetID, pos1, pos2, pos3, pos4;
	value_t eps;
	Tetrahedron tet;

	for (int i = 0; i < tetNum; ++i)
	{
		input >> tetID >> pos1 >> pos2 >> pos3 >> pos4 >> eps;
		tet.SetVertices(pos1 - 1, pos2 - 1, pos3 - 1, pos4 - 1);
		tet.SetEpsilon(eps);
		tetrahedron_.push_back(tet);
	}

	VectorR3 node;
	VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
	VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

	int nodeID;
	for (int i = 0; i < nodeNum; ++i)
	{
		input >> nodeID >> node.x >> node.y >> node.z;
		vertices_.push_back(std::move(node));

		minBox.x = Min(minBox.x, node.x);
		minBox.y = Min(minBox.y, node.y);
		minBox.z = Min(minBox.z, node.z);

		maxBox.x = Max(maxBox.x, node.x);
		maxBox.y = Max(maxBox.y, node.y);
		maxBox.z = Max(maxBox.z, node.z);
	}

	input >> nodeID;
	if (!input.eof())
		return false;

	boundary_.setBoundary(minBox, maxBox);
	
	for (int t = 0; t < tetNum; ++t)
	{
		tetrahedron_[t].GetVertices(pos1, pos2, pos3, pos4);
		tetrahedron_[t].SetVertices(vertices_[pos1], vertices_[pos2], vertices_[pos3], vertices_[pos4]);
	}

	return true;
}

bool Mesh::readVSIENasTmpFile()
{
	Qifstream input("temp.dat", std::ios::in);
	if (input.fail())
		return false;
	int triNum, tetNum, nodeNum;
	input >> triNum >> tetNum >> nodeNum;
	triangles_.reserve(triNum);
	tetrahedron_.reserve(tetNum);
	vertices_.reserve(nodeNum);

	int  triID, tetID, pos1, pos2, pos3, pos4, groupID;
	value_t eps;
	Triangle tri;
	Tetrahedron tet;

	for (int i = 0; i < triNum; ++i)
	{
		input >> triID >> pos1 >> pos2 >> pos3 >> groupID;
		tri.setVertex(pos1 - 1, pos2 - 1, pos3 - 1);
		tri.setgroupID(groupID);
		triangles_.push_back(tri);
	}

	for (int i = 0; i < tetNum; ++i)
	{
		input >> tetID >> pos1 >> pos2 >> pos3 >> pos4 >> eps;
		tet.SetVertices(pos1 - 1, pos2 - 1, pos3 - 1, pos4 - 1);
		tet.SetEpsilon(eps);
		tetrahedron_.push_back(tet);
	}

	VectorR3 node;
	VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
	VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

	int nodeID;
	for (int i = 0; i < nodeNum; ++i)
	{
		input >> nodeID >> node.x >> node.y >> node.z;
		vertices_.push_back(std::move(node));

		minBox.x = Min(minBox.x, node.x);
		minBox.y = Min(minBox.y, node.y);
		minBox.z = Min(minBox.z, node.z);

		maxBox.x = Max(maxBox.x, node.x);
		maxBox.y = Max(maxBox.y, node.y);
		maxBox.z = Max(maxBox.z, node.z);
	}

	input >> nodeID;
	if (!input.eof())
		return false;

	boundary_.setBoundary(minBox, maxBox);

	for (int t = 0; t < triNum; ++t)
	{
		triangles_[t].getVertex(pos1, pos2, pos3);
		triangles_[t].setVertex(vertices_[pos1], vertices_[pos2], vertices_[pos3]);
	}

	for (int t = 0; t < tetNum; ++t)
	{
		tetrahedron_[t].GetVertices(pos1, pos2, pos3, pos4);
		tetrahedron_[t].SetVertices(vertices_[pos1], vertices_[pos2], vertices_[pos3], vertices_[pos4]);
	}

	return true;
}

bool Mesh::readIBCNasTmpFile()
{
	Qifstream input("temp.dat", std::ios::in);
	if (input.fail())
		return false;
	int triNum, nodeNum;
	input >> triNum >> nodeNum;

	triangles_.reserve(triNum);
	vertices_.reserve(nodeNum);

	int  triID, pos1, pos2, pos3;
	Qstring fla;
	Triangle tri;
	for (int i = 0; i < triNum; ++i)
	{
		input >> triID >> pos1 >> pos2 >> pos3 >> fla;
		tri.setVertex(pos1 - 1, pos2 - 1, pos3 - 1);
		tri.setFla(fla);
		triangles_.push_back(tri);
	}

	VectorR3 node;
	VectorR3 minBox(MAX_VALUE, MAX_VALUE, MAX_VALUE);
	VectorR3 maxBox(MIN_VALUE, MIN_VALUE, MIN_VALUE);

	int nodeID;
	for (int i = 0; i < nodeNum; ++i)
	{
		input >> nodeID >> node.x >> node.y >> node.z;
		vertices_.push_back(std::move(node));

		minBox.x = Min(minBox.x, node.x);
		minBox.y = Min(minBox.y, node.y);
		minBox.z = Min(minBox.z, node.z);

		maxBox.x = Max(maxBox.x, node.x);
		maxBox.y = Max(maxBox.y, node.y);
		maxBox.z = Max(maxBox.z, node.z);
	}

	input >> nodeID;
	if (!input.eof())
		return false;

	boundary_.setBoundary(minBox, maxBox);

	for (int t = 0; t < triNum; ++t)
	{
		triangles_[t].getVertex(pos1, pos2, pos3);
		triangles_[t].setVertex(vertices_[pos1], vertices_[pos2], vertices_[pos3]);
	}

	return true;
}

bool Mesh::preprocessMesh()
{
	return true;
}
