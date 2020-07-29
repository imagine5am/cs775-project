#pragma once

#include "basics.h"

class DuneSediment
{
private:
	float tanThresholdAngleSediment = 0.60f;		// ~33°
	float tanThresholdAngleWindShadowMin = 0.08f;	// ~5°
	float tanThresholdAngleWindShadowMax = 0.26f;	// ~15°
	float tanThresholdAngleBedrock = 2.5f;			// ~68°

	bool vegetationOn = false;
	bool abrasionOn = false;

protected:
	ScalarField2D bedrock;			//!< Bedrock elevation layer, in meter.
	ScalarField2D sediments;		//!< Sediment elevation layer, in meter.
	ScalarField2D vegetation;		//!< Vegetation presence in [0, 1].

	Box2D box;						//!< World space bounding box.
	int nx, ny;						//!< Grid resolution.
	float matterToMove;				//!< Amount of sand transported by the wind, in meter.
	float cellSize;					//!< Size of one cell in meter, squared. Stored to speed up the simulation.
	Vector2 wind;					//!< Base wind direction.

public:
	std::vector<std::vector<Quad*>> roots;
	int maxlevels;
	DuneSediment();
	DuneSediment(const Box2D& bbox, float rMin, float rMax, const Vector2& w);
	DuneSediment(const Box2D& bbox, float rMin, float rMax, const Vector2& w, const int resx, const int resy);
	DuneSediment(const DuneSediment lowDune, const int factor);
	//QT compatible:
	DuneSediment(const DuneSediment lowDune, std::function<bool(const Vector2i&, const DuneSediment& )> cellQuery, const int numLevels);
	
	~DuneSediment();

	// Simulation
	int ToIndex1D(const Vector2i& q) const;
	int ToIndex1D(int i, int j) const;
	void SimulationStepMultiThreadAtomic();
	//QT compatible:
	void SimulationStepQTMultiThreadAtomic(int factor, bool stabilizeBR);

	void EndSimulationStep();
	void SimulationStepWorldSpace();
	//QT compatible:
	void SimulationStepQT();

	void PerformReptationOnCell(int i, int j, int bounce);
	//QT compatible:
	void PerformReptationOnCell(CellData* cd, int bounce);

	void ComputeWindAtCell(int i, int j, Vector2& windDir) const;
	//QT compatible:
	void ComputeWindAtCell(CellData* cd, Vector2& windDir) const;
	float IsInShadow(int i, int j, const Vector2& wind) const;
	//QT compatible:
	float IsInShadow(CellData* cd, const Vector2& wind) const;

	void SnapWorld(Vector2& p) const;
	int CheckSedimentFlowRelative(const Vector2i& p, float tanThresholdAngle, Vector2i* nei, float* nslope) const;
	//QT compatible:
	int CheckSedimentFlowRelative(CellData* cd, float tanThresholdAngle, std::vector<CellData*> &neighbors, std::vector<float> &nslope) const;

	int CheckBedrockFlowRelative(const Vector2i& p, float tanThresholdAngle, Vector2i* nei, float * nslope) const;
	//QT compatible:
	int CheckBedrockFlowRelative(CellData* cd, float tanThresholdAngle, std::vector<CellData*>& neighbors, std::vector<float>& nslope) const;

	void StabilizeSedimentRelative(int i, int j);
	//QT compatible:
	void StabilizeSedimentRelative(CellData* cd);

	bool StabilizeBedrockRelative(int i, int j);
	//QT compatible:
	bool StabilizeBedrockRelative(CellData* cd);

	void StabilizeBedrockAll();
	//QT Compatible
	void StabilizeBedrockAllQT();

	void PerformAbrasionOnCell(int i, int j, const Vector2& windDir);
	//QT Compatible
	void PerformAbrasionOnCell(CellData* cd, const Vector2& windDir);

	// Exports & Meshing
	void ExportObj(const std::string& file) const;
	//QT Compatible
	void ExportQTObj(const std::string& file) const;

	// Inlined functions and query

	float Height(int i, int j) const;
	//QT Compatible
	float Height(float i, float j) const;

	float Height(const Vector2& p) const;
	//QT Compatible
	float Height(Quad* qd) const;

	float Bedrock(int i, int j) const;
	float Sediment(int i, int j) const;
	void SetAbrasionMode(bool c);
	void SetVegetationMode(bool c);
	void SetBedrockData(const ScalarField2D& f);
	void SetSedimentData(const ScalarField2D& f);

	Box2D getBox() const;
	int getNx() const;
	int getNy() const;
	Vector2 getWind() const;
};



inline Box2D DuneSediment::getBox() const
{
	return box;
}

inline int DuneSediment::getNx() const
{
	return nx;
}

inline int DuneSediment::getNy() const
{
	return ny;
}

inline Vector2 DuneSediment::getWind() const
{
	return wind;
}

/*!
\brief Compute the 1D index from a given grid vertex.
\param q grid vertex.
*/
inline int DuneSediment::ToIndex1D(int i, int j) const
{
	return bedrock.ToIndex1D(i, j);
}

/*!
\brief Compute the 1D index from a given grid vertex.
\param q grid vertex.
*/
inline int DuneSediment::ToIndex1D(const Vector2i& q) const
{
	return bedrock.ToIndex1D(q);
}

/*!
\brief
*/
inline float DuneSediment::Height(int i, int j) const
{
	return bedrock.Get(i, j) + sediments.Get(i, j);
}

//QT compatible:
inline float DuneSediment::Height(float i, float j) const
{
	//Quad* roots[16][16];
	CellData* cd =  roots[int(i)][int(j)]->search(Vector2(i,j) );
	return cd->vals['b'] + cd->vals['s'];
	//return bedrock.Get(i, j) + sediments.Get(i, j);
}

//QT compatible:
inline float DuneSediment::Height(Quad* qd) const
{
	//Quad* roots[16][16];
	//CellData* cd = roots[int(i)][int(j)]->search(Vector2(i, j));
	return qd->fetchData()->vals['b'] + qd->fetchData()->vals['s'];
	//return bedrock.Get(i, j) + sediments.Get(i, j);
}

/*!
\brief
*/
inline float DuneSediment::Height(const Vector2& p) const 
{
	return bedrock.GetValueBilinear(p) + sediments.GetValueBilinear(p);
}

/*!
\brief
*/
inline float DuneSediment::Bedrock(int i, int j) const
{
	return bedrock.Get(i, j);
}

/*!
\brief
*/
inline float DuneSediment::Sediment(int i, int j) const
{
	return sediments.Get(i, j);
}

/*!
\brief
*/
inline void DuneSediment::SetAbrasionMode(bool c)
{
	abrasionOn = c;
}

/*!
\brief
*/
inline void DuneSediment::SetVegetationMode(bool c)
{
	vegetationOn = c;
}

/*!
\brief
*/
inline void DuneSediment::SetBedrockData(const ScalarField2D& f)
{
	bedrock = f;
}

/*!
\brief
*/
inline void DuneSediment::SetSedimentData(const ScalarField2D& f)
{
	sediments = f;
}
