//**********************************************************
// Author: Zhang Qiang & Xiang Wei
// License: MIT
//**********************************************************

#ifndef _MESH_H_
#define _MESH_H_
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Boundary.h"

namespace component {

class Mesh {
public:
    Mesh();
    ~Mesh();
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

public:
    void loadMeshFile(const Qstring& filename);
    void loadMeshForEFIE_PMCHW(const Qstring& filename);
    void loadDDMMesh(const Qstring& filename);
	void loadVIEMesh(const Qstring& filename);
	void loadVSIEMesh(const Qstring& filename);
	void loadIBCMesh(const Qstring& filename);

    const Triangle& getTriangleRef(size_t i) const;
    Triangle& getTriangleRef(size_t i);
	Tetrahedron& getTetrahedronRef(size_t i);
    const VectorR3& getVertex(size_t i)const;
    int getTriangleNum() const;
	int getNodeNum() const;
    int getSubdomainNum() const;
	int getTetrahedronNum() const;
    void getBoundary(VectorR3 &_min_boundary, VectorR3 &_max_boundary) const;
    VectorR3 getSize() const;
	VectorR3 getCenterBox() const;
    void getSize(value_t &x, value_t &y, value_t &z) const;
    void reportInfo(Qostream& strm) const;
    void clear();

    void Debug();
private:
    bool prepareForEFIE_PMCHW(const Qstring& filename);
    bool loadRptFile(const Qstring& filename);
    bool loadStlFile(const Qstring& filename);
    bool loadNasFile(const Qstring& filename);
    bool loadNasFileForEFIE_PMCHW(const Qstring& filename);
    bool loadDDMMapFile(const Qstring& filename);
	bool loadDDMNasFile(const Qstring& filename);

	bool loadVIEMapFile(const Qstring& filename);
	bool loadVIENasFile(const Qstring& filename);

	bool loadVSIEMapFile(const Qstring& filename);
	bool loadVSIENasFile(const Qstring& filename);

	bool loadIBCMapFile(const Qstring& filename);
	bool loadIBCNasFile(const Qstring& filename);

    bool readRptTmpFile();
    bool readNasTmpFile();
    bool readStlTmpFile();
    bool readNasTmpFileForEFIE_PMCHW();
    bool readDDMNasTmpFile();
	bool readVIENasTmpFile();
	bool readVSIENasTmpFile();
	bool readIBCNasTmpFile();
	bool preprocessMesh();

private:
    int                     subdomain_num_;
	int						tri_number;
	int						tet_number;
    Boundary                boundary_;
    std::set<int>           common_face_;
    std::map<int, int>      subdomain_id_;
	std::map<int, value_t>  subdomain_eps_;
	std::map<int, Qstring>  subdomain_fla_;
    std::vector<VectorR3>   vertices_;
    std::vector<Triangle>   triangles_;
	std::vector<Tetrahedron> tetrahedron_;
};

inline const Triangle & Mesh::getTriangleRef(size_t i) const
{
    return triangles_[i];
}

inline Triangle & Mesh::getTriangleRef(size_t i)
{
    return triangles_[i];
}

inline Tetrahedron & Mesh::getTetrahedronRef(size_t i)
{
	return tetrahedron_[i];
}

inline const VectorR3& Mesh::getVertex(size_t i) const
{
    return vertices_[i];
}

inline int Mesh::getTriangleNum() const
{
    return static_cast<int>(triangles_.size());
}

inline int Mesh::getNodeNum() const
{
	return static_cast<int>(vertices_.size());
}
inline int Mesh::getSubdomainNum() const
{
    return subdomain_num_;
}

inline int Mesh::getTetrahedronNum() const
{
	return static_cast<int>(tetrahedron_.size());
}

inline void Mesh::getBoundary(VectorR3 & _min_boundary, VectorR3 & _max_boundary) const
{
    _min_boundary = boundary_.getBoundaryMin();
    _max_boundary = boundary_.getBoundaryMax();
}

inline void Mesh::getSize(value_t &x, value_t &y, value_t &z) const
{
    x = (boundary_.getBoundaryMax() - boundary_.getBoundaryMin()).x;
    y = (boundary_.getBoundaryMax() - boundary_.getBoundaryMin()).y;
    z = (boundary_.getBoundaryMax() - boundary_.getBoundaryMin()).z;
}

inline VectorR3 Mesh::getSize() const
{
    return VectorR3(boundary_.getBoundaryMax() - boundary_.getBoundaryMin());
}

inline VectorR3 Mesh::getCenterBox() const
{
	return VectorR3((boundary_.getBoundaryMax() + boundary_.getBoundaryMin()) / 2.0f);
}

} // namespace component

#endif