#include "BaseModel.h"
#include "Parameters.h"
#include <windows.h>
#include "gl\gl.h"
#include "gl\glu.h"
#include <float.h>
#include <iostream>
#pragma comment(lib, "OPENGL32.LIB")
#pragma comment(lib, "GLU32.LIB")
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
using namespace std;

namespace
{
	constexpr double kEpsilonTolerance = 1e-12;
	constexpr double kMinDeterminantThreshold = 1e-18;

	bool NormalizeSafe(CPoint3D& v)
	{
		const double len = v.Len();
		if (len <= kEpsilonTolerance)
		{
			return false;
		}
		v /= len;
		return true;
	}

	void BuildTangentFrameFromNormal(const CPoint3D& normalUnit, CPoint3D& t1, CPoint3D& t2)
	{
		CPoint3D axis = (std::fabs(normalUnit.x) < 0.9) ? CPoint3D(1.0, 0.0, 0.0) : CPoint3D(0.0, 1.0, 0.0);
		t1 = axis * normalUnit;
		if (!NormalizeSafe(t1))
		{
			axis = CPoint3D(0.0, 0.0, 1.0);
			t1 = axis * normalUnit;
			if (!NormalizeSafe(t1))
			{
				t1 = CPoint3D(1.0, 0.0, 0.0);
			}
		}

		t2 = normalUnit * t1;
		if (!NormalizeSafe(t2))
		{
			t2 = CPoint3D(0.0, 1.0, 0.0);
		}
	}
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CBaseModel::CBaseModel(const string& filename) : m_filename(filename)
{
}

void CBaseModel::BuildVertexToFaceMapping()
{
	m_vertToFaces.clear();
	const int numVerts = GetNumOfVerts();
	
	for (int vertIndex = 0; vertIndex < numVerts; ++vertIndex)
	{
		m_vertToFaces[vertIndex] = vector<int>();
	}

	const int numFaces = GetNumOfFaces();
	for (int faceIndex = 0; faceIndex < numFaces; ++faceIndex)
	{
		const CBaseModel::CFace& face = Face(faceIndex);
		for (int i = 0; i < 3; ++i)
		{
			const int vertIndex = face[i];
			m_vertToFaces[vertIndex].push_back(faceIndex);
		}
	}
}

const std::vector<int>& CBaseModel::GetFacesFromVertex(int vertIndex) const
{
	static const std::vector<int> emptyVector;
	auto it = m_vertToFaces.find(vertIndex);
	if (it != m_vertToFaces.end())
	{
		return it->second;
	}
	
	return emptyVector;
}

const double CBaseModel::GetAreaOfTriangle(int faceIndex) const
{
	return m_TriangleAreas[faceIndex];
}
void CBaseModel::Render() const
{	
	glPushMatrix();
	GLint shadeModel;
	glGetIntegerv(GL_SHADE_MODEL, &shadeModel);
	if (shadeModel == GL_SMOOTH)
	{
		for (int i = 0; i < GetNumOfFaces(); ++i)
		{
			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; ++j)
			{			
				const CPoint3D &pt = Vert(Face(i)[j]);
				if (!m_NormalsToVerts.empty())
				{
					const CPoint3D &normal = Normal(Face(i)[j]);
					glNormal3f((float)normal.x, (float)normal.y, (float)normal.z);
				}
				glVertex3f((float)pt.x, (float)pt.y, (float)pt.z);
			}	
			glEnd();
		}	
	}
	else
	{
		for (int i = 0; i < GetNumOfFaces(); ++i)
		{
			glBegin(GL_POLYGON);
			CPoint3D normal = VectorCross(Vert(Face(i)[0]), Vert(Face(i)[1]), Vert(Face(i)[2]));
			normal.Normalize();
			glNormal3f((float)normal.x, (float)normal.y, (float)normal.z);
			for (int j = 0; j < 3; ++j)
			{			
				const CPoint3D &pt = Vert(Face(i)[j]);
				glVertex3f((float)pt.x, (float)pt.y, (float)pt.z);
			}	
			glEnd();
		}	
	}
	
	glPopMatrix();
}

void CBaseModel::SetVertex(double x, double y, double z) 
{
	m_Verts.emplace_back(CPoint3D(x, y, z));
}

void CBaseModel::ComputeScaleAndNormals()
{
	if (m_Verts.empty())
		return;
	m_NormalsToVerts.clear();
	m_NormalsToFaces.clear();
	m_NormalsToVerts.resize(m_Verts.size(), CPoint3D(0, 0, 0));
	CPoint3D center(0, 0, 0);
	double sumArea(0);
	CPoint3D sumNormal(0, 0, 0);
	double deta(0);
	for (int i = 0; i < (int)m_Faces.size(); ++i)
	{
		CPoint3D normal = VectorCross(Vert(Face(i)[0]),
			Vert(Face(i)[1]),
			Vert(Face(i)[2]));
		double area = normal.Len();
		m_TriangleAreas.push_back(area);
		CPoint3D gravity3 = Vert(Face(i)[0]) + Vert(Face(i)[1]) + Vert(Face(i)[2]);
		center += area * gravity3;
		sumArea += area;
		sumNormal += normal;
		deta += gravity3 ^ normal;
		normal.x /= area;
		normal.y /= area;
		normal.z /= area;
		m_NormalsToFaces.emplace_back(normal);
		for (int j = 0; j < 3; ++j)
		{
			m_NormalsToVerts[Face(i)[j]] += normal * area;
		}
	}
	center /= sumArea * 3;
	deta -= 3 * (center ^ sumNormal);
	if (true)//deta > 0)
	{
		for (int i = 0; i < GetNumOfVerts(); ++i)
		{
			if (fabs(m_NormalsToVerts[i].x)
				+ fabs(m_NormalsToVerts[i].y)
				+ fabs(m_NormalsToVerts[i].z) >= FLT_EPSILON)
			{					
				m_NormalsToVerts[i].Normalize();
			}
		}
	}
	else
	{
		for (int i = 0; i < GetNumOfFaces(); ++i)
		{
			int temp = m_Faces[i][0];
			m_Faces[i][0] = m_Faces[i][1];
			m_Faces[i][1] = temp;
		}
		for (int i = 0; i < GetNumOfVerts(); ++i)
		{
			if (fabs(m_NormalsToVerts[i].x)
				+ fabs(m_NormalsToVerts[i].y)
				+ fabs(m_NormalsToVerts[i].z) >= FLT_EPSILON)
			{					
				double len = m_NormalsToVerts[i].Len();
				m_NormalsToVerts[i].x /= -len;
				m_NormalsToVerts[i].y /= -len;
				m_NormalsToVerts[i].z /= -len;
			}
		}
	}

	CPoint3D ptUp(m_Verts[0]);
	CPoint3D ptDown(m_Verts[0]);
	for (int i = 1; i < GetNumOfVerts(); ++i)
	{
		if (m_Verts[i].x > ptUp.x)
			ptUp.x = m_Verts[i].x;
		else if (m_Verts[i].x < ptDown.x)
			ptDown.x = m_Verts[i].x;
		if (m_Verts[i].y > ptUp.y)
			ptUp.y = m_Verts[i].y;
		else if (m_Verts[i].y < ptDown.y)
			ptDown.y = m_Verts[i].y;
		if (m_Verts[i].z > ptUp.z)
			ptUp.z = m_Verts[i].z;
		else if (m_Verts[i].z < ptDown.z)
			ptDown.z = m_Verts[i].z;
	}	
	m_boundingBox = std::make_pair(ptUp, ptDown);

	double maxEdgeLenOfBoundingBox = -1;
	if (ptUp.x - ptDown.x > maxEdgeLenOfBoundingBox)
		maxEdgeLenOfBoundingBox = ptUp.x - ptDown.x;
	if (ptUp.y - ptDown.y > maxEdgeLenOfBoundingBox)
		maxEdgeLenOfBoundingBox = ptUp.y - ptDown.y;
	if (ptUp.z - ptDown.z > maxEdgeLenOfBoundingBox)
		maxEdgeLenOfBoundingBox = ptUp.z - ptDown.z;
	m_scale = maxEdgeLenOfBoundingBox / 2;
	//m_center = center;
	//m_ptUp = ptUp;
	//m_ptDown = ptDown;

	ComputePrincipalCurvaturesAndDirections();
}

const CPoint3D& CBaseModel::GetFaceNormal(int faceIndex) const
{
	return m_NormalsToFaces[faceIndex];
}

void CBaseModel::ComputePrincipalCurvaturesAndDirections()
{
	const int numVerts = GetNumOfVerts();
	m_PrincipalKMin.assign(numVerts, 0.0);
	m_PrincipalKMax.assign(numVerts, 0.0);
	m_PrincipalDirMin.assign(numVerts, CPoint3D(1.0, 0.0, 0.0));
	m_PrincipalDirMax.assign(numVerts, CPoint3D(0.0, 1.0, 0.0));

	if (numVerts == 0)
	{
		return;
	}

	if (m_vertToFaces.empty())
	{
		BuildVertexToFaceMapping();
	}

	for (int i = 0; i < numVerts; ++i)
	{
		CPoint3D n = Normal(i);
		if (!NormalizeSafe(n))
		{
			n = CPoint3D(0.0, 0.0, 1.0);
		}

		CPoint3D t1, t2;
		BuildTangentFrameFromNormal(n, t1, t2);

		std::set<int> neighbors;
		const std::vector<int>& incidentFaces = GetFacesFromVertex(i);
		for (int faceId : incidentFaces)
		{
			const CFace& f = Face(faceId);
			for (int k = 0; k < 3; ++k)
			{
				if (f[k] != i)
				{
					neighbors.insert(f[k]);
				}
			}
		}

		if (neighbors.size() < 3)
		{
			m_PrincipalDirMin[i] = t1;
			m_PrincipalDirMax[i] = t2;
			continue;
		}

		Eigen::Matrix3d ata = Eigen::Matrix3d::Zero();
		Eigen::Vector3d atb = Eigen::Vector3d::Zero();
		int validCount = 0;

		const CPoint3D& vi = Vert(i);
		for (int nbr : neighbors)
		{
			const CPoint3D d = Vert(nbr) - vi;
			const double u = d ^ t1;
			const double v = d ^ t2;
			const double w = d ^ n;
			if (std::fabs(u) + std::fabs(v) <= kEpsilonTolerance)
			{
				continue;
			}

			Eigen::Vector3d row(0.5 * u * u, u * v, 0.5 * v * v);
			ata += row * row.transpose();
			atb += row * w;
			++validCount;
		}

		if (validCount < 3 || std::fabs(ata.determinant()) <= kMinDeterminantThreshold)
		{
			m_PrincipalDirMin[i] = t1;
			m_PrincipalDirMax[i] = t2;
			continue;
		}

		// ata is symmetric normal-equation matrix; LDLT is stable and efficient for this case.
		const Eigen::Vector3d coeff = ata.ldlt().solve(atb);
		if (!coeff.allFinite())
		{
			m_PrincipalDirMin[i] = t1;
			m_PrincipalDirMax[i] = t2;
			continue;
		}

		Eigen::Matrix2d hessian;
		hessian << coeff[0], coeff[1],
			coeff[1], coeff[2];
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(hessian);
		if (eig.info() != Eigen::Success)
		{
			m_PrincipalDirMin[i] = t1;
			m_PrincipalDirMax[i] = t2;
			continue;
		}

		const Eigen::Vector2d evals = eig.eigenvalues();
		const Eigen::Matrix2d evecs = eig.eigenvectors();
		m_PrincipalKMin[i] = evals[0];
		m_PrincipalKMax[i] = evals[1];

		CPoint3D vMin = t1 * evecs(0, 0) + t2 * evecs(1, 0);
		CPoint3D vMax = t1 * evecs(0, 1) + t2 * evecs(1, 1);

		vMin -= n * (vMin ^ n);
		if (!NormalizeSafe(vMin))
		{
			vMin = t1;
		}

		vMax -= n * (vMax ^ n);
		vMax -= vMin * (vMax ^ vMin);
		if (!NormalizeSafe(vMax))
		{
			vMax = n * vMin;
			if (!NormalizeSafe(vMax))
			{
				vMax = t2;
			}
		}

		m_PrincipalDirMin[i] = vMin;
		m_PrincipalDirMax[i] = vMax;
	}
}

void CBaseModel::LoadModel()
{
	ReadFile(m_filename);
	ComputeScaleAndNormals();
	BuildVertexToFaceMapping();
}

string CBaseModel::GetFileShortName() const
{
	int pos = (int)m_filename.size() - 1;
	while (pos >= 0)
	{
		if (m_filename[pos] == '\\')
			break;
		--pos;
	}	
	++pos;
	string str(m_filename.substr(pos));
	return str;
}

string CBaseModel::GetFileFullName() const
{
	return m_filename;
}

void CBaseModel::ReadObjFile(const string& filename)
{
	ifstream in(filename.c_str());
	if (in.fail())
	{
		throw "fail to read file";
	}
	char buf[256];
	while (in.getline(buf, sizeof buf))
	{
		istringstream line(buf);
		string word;
		line >> word;
		if (word == "v")
		{
			CPoint3D pt;
			line >> pt.x;
			line >> pt.y;
			line >> pt.z;

			m_Verts.push_back(pt);
		}
		else if (word == "f")
		{
			CFace face;
			int tmp;
			vector<int> polygon;
			polygon.reserve(4);
			while (line >> tmp)
			{
				polygon.push_back(tmp);
				char tmpBuf[256];
				line.getline(tmpBuf, sizeof tmpBuf, ' ');
			}
			for (int j = 1; j < (int)polygon.size() - 1; ++j)
			{
				face[0] = polygon[0] - 1;
				face[1] = polygon[j] - 1;
				face[2] = polygon[j + 1] - 1;
				m_Faces.push_back(face);
			}
		}
		else
		{
			continue;
		}
	}
	m_Verts.swap(vector<CPoint3D>(m_Verts));
	m_Faces.swap(vector<CFace>(m_Faces));
	in.close();
}

void CBaseModel::ReadFile(const string& filename)
{
	int nDot = (int)filename.rfind('.');
	if (nDot == -1)
	{
		throw "File name doesn't contain a dot!";
	}
	string extension = filename.substr(nDot + 1);
	
	if (extension == "obj")
	{
		ReadObjFile(filename);
	}
	else if (extension == "off")
	{
		ReadOffFile(filename);
	}
	else if (extension == "m")
	{
		ReadMFile(filename);
	}
	else
	{
		throw "This format can't be handled!";
	}
}

void CBaseModel::ReadOffFile(const string& filename)
{
	ifstream in(filename.c_str());
	if (in.fail())
	{
		throw "fail to read file";
	}
	char buf[256];
	in.getline(buf, sizeof buf);
	int vertNum, faceNum, edgeNum;
	in >> vertNum >> faceNum >> edgeNum;

	for (int i = 0; i < vertNum; ++i)
	{
		CPoint3D pt;
		in >> pt.x;
		in >> pt.y;
		in >> pt.z;
		m_Verts.push_back(pt);
	}
	m_Verts.swap(vector<CPoint3D>(m_Verts));
	
	int degree;
	while (in >> degree)
	{
		int first, second;
		in >> first >> second;

		for (int i = 0; i < degree - 2; ++i)
		{
			CFace f;
			f[0] = first;
			f[1] = second;
			in >> f[2];
			m_Faces.push_back(f);			
			second = f[2];
		}
	}

	in.close();
	m_Faces.swap(vector<CFace>(m_Faces));
}

void CBaseModel::ReadMFile(const string& filename)
{
	ifstream in(filename.c_str());
	if (in.fail())
	{
		throw "fail to read file";
	}
	char buf[256];
	while (in.getline(buf, sizeof buf))
	{
		istringstream line(buf);
		if (buf[0] == '#')
			continue;
		string word;
		line >> word;
		if (word == "Vertex")
		{
			int tmp;
			line >> tmp;
			CPoint3D pt;
			line >> pt.x;
			line >> pt.y;
			line >> pt.z;

			m_Verts.push_back(pt);
		}
		else if (word == "Face")
		{
			CFace face;
			int tmp;
			line >> tmp;
			vector<int> polygon;
			polygon.reserve(4);
			while (line >> tmp)
			{
				polygon.push_back(tmp);
			}
			for (int j = 1; j < (int)polygon.size() - 1; ++j)
			{
				face[0] = polygon[0] - 1;
				face[1] = polygon[j] - 1;
				face[2] = polygon[j + 1] - 1;
				m_Faces.push_back(face);
			}
		}
		else
		{
			continue;
		}
	}
	m_Verts.swap(vector<CPoint3D>(m_Verts));
	m_Faces.swap(vector<CFace>(m_Faces));
	in.close();
}


void CBaseModel::SaveMFile(const string& filename) const
{
	ofstream outFile(filename.c_str());
	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << "Vertex " << i + 1 << " " << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	int cnt(0);
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile <<"Face " << ++cnt << " " << Face(i)[0] + 1 << " " << Face(i)[1] + 1 << " " << Face(i)[2] + 1 << endl;
	}
	outFile.close();
}

void CBaseModel::SaveOffFile(const string& filename) const
{
	ofstream outFile(filename.c_str());
	outFile << "OFF" << endl;
	outFile << GetNumOfVerts() << " " << GetNumOfFaces() << " " << 0 << endl;
	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile << 3 << " " << Face(i)[0]<< " " << Face(i)[1] << " " << Face(i)[2] << endl;
	}
	outFile.close();
}

void CBaseModel::SaveObjFile(const string& filename) const
{
	ofstream outFile(filename.c_str());
	outFile << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << "v " << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile << "f " << Face(i)[0] + 1 << " " << Face(i)[1] + 1<< " " << Face(i)[2] + 1<< endl;
	}
	outFile.close();
}

void CBaseModel::SaveScalarFieldObjFile(const vector<double>& vals, const string& filename) const
{
	ofstream outFile(filename.c_str());
	outFile << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
	outFile << "# maxDis: " << *max_element(vals.begin(), vals.end()) << endl;
	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << "v " << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	for (int i = 0; i < (int)vals.size(); ++i)
	{
		outFile << "vt " << vals[i] << " " << 0 << endl;
	}
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile << "f " << Face(i)[0] + 1 << "/" << Face(i)[0] + 1
			<< " " << Face(i)[1] + 1 << "/" << Face(i)[1] + 1
			<< " " << Face(i)[2] + 1 << "/" << Face(i)[2] + 1 << endl;
	}
	outFile.close();
}

void CBaseModel::SaveScalarFieldObjFile(const vector<double>& vals, const string& comments, const string& filename) const
{
	ofstream outFile(filename.c_str());
	outFile << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
	outFile << comments << endl;
	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << "v " << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	for (int i = 0; i < (int)vals.size(); ++i)
	{
		outFile << "vt " << vals[i] << " " << 0 << endl;
	}
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile << "f " << Face(i)[0] + 1 << "/" << Face(i)[0] + 1
			<< " " << Face(i)[1] + 1 << "/" << Face(i)[1] + 1
			<< " " << Face(i)[2] + 1 << "/" << Face(i)[2] + 1 << endl;
	}
	outFile.close();
}

void CBaseModel::SaveScalarFieldObjFile(const vector<double>& vals, double maxV, const string& filename) const
{
	ofstream outFile(filename.c_str());
	outFile << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
	outFile << "# maxValue = " << *max_element(vals.begin(), vals.end()) / maxV << endl;

	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << "v " << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	for (int i = 0; i < (int)vals.size(); ++i)
	{
		outFile << "vt " << vals[i] / maxV << " " << 0 << endl;
	}
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile << "f " << Face(i)[0] + 1 << "/" << Face(i)[0] + 1
			<< " " << Face(i)[1] + 1 << "/" << Face(i)[1] + 1
			<< " " << Face(i)[2] + 1 << "/" << Face(i)[2] + 1 << endl;
	}
	outFile.close();
}

void CBaseModel::SavePamametrizationObjFile(const vector<pair<double, double>>& uvs, const string& filename) const
{
	ofstream outFile(filename.c_str());
	outFile << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
	for (int i = 0; i < (int)GetNumOfVerts(); ++i)
	{
		outFile << "v " << Vert(i).x << " " << Vert(i).y << " " << Vert(i).z << endl;
	}
	for (int i = 0; i < (int)uvs.size(); ++i)
	{
		outFile << "vt " << uvs[i].first << " " << uvs[i].second << endl;
	}
	for (int i = 0; i < (int)GetNumOfFaces(); ++i)
	{
		if (m_UselessFaces.find(i) != m_UselessFaces.end())
			continue;
		outFile << "f " << Face(i)[0] + 1 << "/" << Face(i)[0] + 1
			<< " " << Face(i)[1] + 1 << "/" << Face(i)[1] + 1
			<< " " << Face(i)[2] + 1 << "/" << Face(i)[2] + 1 << endl;
	}
	outFile.close();
}

void CBaseModel::PrintInfo(ostream& out) const
{
	out << "Model info is as follows.\n";
	out << "Name: " << GetFileShortName() << endl;
	out << "VertNum = " << GetNumOfVerts() << endl;
	out << "FaceNum = " << GetNumOfFaces() << endl;
	out << "Scale = " << m_scale << endl;
}

CPoint3D CBaseModel::GetShiftVertex(int indexOfVert) const
{
	return Vert(indexOfVert) + Normal(indexOfVert) * RateOfNormalShift * GetScale();
}

//CPoint3D CBaseModel::ShiftVertex(int indexOfVert, double epsilon) const
//{
//	return Vert(indexOfVert) +  Normal(indexOfVert) * epsilon;
//}

int CBaseModel::GetVertexID(const CPoint3D& pt) const
{
	double dis = DBL_MAX;
	int id;
	for (int i = 0; i < GetNumOfVerts(); ++i)
	{
		if ((Vert(i) - pt).Len() < dis)
		{
			id = i;
			dis = (Vert(id)-pt).Len();
		}
	}
	return id;
}

string CBaseModel::GetComments(const char* filename)
{
	ifstream in(filename);
	char buf[256];
	string result;
	while (in.getline(buf, sizeof buf))
	{
		if (buf[0] == '#')
		{
			result += buf;
			result += "\n";
		}
	}
	in.close();
	return result;
}

vector<double> CBaseModel::GetScalarField(string filename)
{
	vector<double> scalarField;
	ifstream in(filename);
	char buf[256];
	while (in.getline(buf, sizeof buf))
	{
		istringstream line(buf);
		string word;
		line >> word;
		if (word == "vt")
		{
			double value;
			line >> value;
			scalarField.push_back(value);
		}
	}

	in.close();
	return scalarField;
}

void CBaseModel::SetFaces(const vector<CBaseModel::CFace>& faces)
{
	m_Faces = faces;
}

const vector<CBaseModel::CFace>& CBaseModel::GetFaces() const
{
	return m_Faces;
}