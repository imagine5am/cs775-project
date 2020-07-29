#include "desert.h"
#include "noise.h"

#include <iostream>
#include <fstream>
#include <queue>
#include <string>

std::vector<std::vector<Quad*>>* Quad::qroots;
int Quad::nx;
int Quad::ny;

/*!
\brief Todo
*/
DuneSediment::DuneSediment()
{
	nx = ny = 256;
	box = Box2D(Vector2(0), 1);
	wind = Vector2(1, 0);

	bedrock = ScalarField2D(nx, ny, box, 0.0);
	vegetation = ScalarField2D(nx, ny, box, 0.0);
	sediments = ScalarField2D(nx, ny, box, 0.0);

	matterToMove = 0.1f;
	Vector2 celldiagonal = Vector2((box.TopRight()[0] - box.BottomLeft()[0]) / (nx - 1), (box.TopRight()[1] - box.BottomLeft()[1]) / (ny - 1));
	cellSize = Box2D(box.BottomLeft(), box.BottomLeft() + celldiagonal).Size().x; // We only consider squared heightfields
}

/*!
\brief Todo
*/
DuneSediment::DuneSediment(const Box2D& bbox, float rMin, float rMax, const Vector2& w)
{
	box = bbox;
	nx = ny = 256;
	wind = w;

	bedrock = ScalarField2D(nx, ny, box, 0.0);
	vegetation = ScalarField2D(nx, ny, box, 0.0);
	sediments = ScalarField2D(nx, ny, box, 0.0);
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Vector2 p = bedrock.ArrayVertex(i, j);

			// Vegetation
			// Arbitrary clamped 2D noise - but you can use whatever you want.
			float v = PerlinNoise::fBm(Vector3(i * 7.91247f, j * 7.91247f, 0.0f), 1.0f, 0.002f, 3) / 1.75f;
			if (v > 0.58f)
				vegetation.Set(i, j, v);

			// Sand
			sediments.Set(i, j, Random::Uniform(rMin, rMax));
		}
	}

	// Debug code to write a ppm file showing the vegetation.
	// @Todo: could be refactored in the scalarfield2 class.
	//int i, j;
	//FILE* fp;
	//fopen_s(&fp, "first.ppm", "wb"); /* b - binary mode */
	//(void)fprintf(fp, "P6\n%d %d\n255\n", nx, ny);
	//for (j = 0; j < nx; ++j)
	//{
	//	for (i = 0; i < ny; ++i)
	//	{
	//		static unsigned char color[3];
	//		int v = (int)(vegetation.Get(i, j) * 256);
	//		color[0] = v % 256;  /* red */
	//		color[1] = v % 256;  /* green */
	//		color[2] = v % 256;  /* blue */
	//		(void)fwrite(color, 1, 3, fp);
	//	}
	//}
	//(void)fclose(fp);

	Vector2 celldiagonal = Vector2((box.TopRight()[0] - box.BottomLeft()[0]) / (nx - 1), (box.TopRight()[1] - box.BottomLeft()[1]) / (ny - 1));
	cellSize = Box2D(box.BottomLeft(), box.BottomLeft() + celldiagonal).Size().x; // We only consider squared heightfields
	
	matterToMove = 0.1f;
}


DuneSediment::DuneSediment(const Box2D& bbox, float rMin, float rMax, const Vector2& w, const int resx, const int resy)
{
	box = bbox;
	nx = resx;
	ny = resy;
	wind = w;

	bedrock = ScalarField2D(nx, ny, box, 0.0);
	vegetation = ScalarField2D(nx, ny, box, 0.0);
	sediments = ScalarField2D(nx, ny, box, 0.0);
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Vector2 p = bedrock.ArrayVertex(i, j);

			// Vegetation
			// Arbitrary clamped 2D noise - but you can use whatever you want.
			float v = PerlinNoise::fBm(Vector3(i * 7.91247f, j * 7.91247f, 0.0f), 1.0f, 0.002f, 3) / 1.75f;
			if (v > 0.58f)
				vegetation.Set(i, j, v);

			// Sand
			sediments.Set(i, j, Random::Uniform(rMin, rMax));
		}
	}

	// Debug code to write a ppm file showing the vegetation.
	// @Todo: could be refactored in the scalarfield2 class.
	//int i, j;
	//FILE* fp;
	//fopen_s(&fp, "first.ppm", "wb"); /* b - binary mode */
	//(void)fprintf(fp, "P6\n%d %d\n255\n", nx, ny);
	//for (j = 0; j < nx; ++j)
	//{
	//	for (i = 0; i < ny; ++i)
	//	{
	//		static unsigned char color[3];
	//		int v = (int)(vegetation.Get(i, j) * 256);
	//		color[0] = v % 256;  /* red */
	//		color[1] = v % 256;  /* green */
	//		color[2] = v % 256;  /* blue */
	//		(void)fwrite(color, 1, 3, fp);
	//	}
	//}
	//(void)fclose(fp);

	Vector2 celldiagonal = Vector2((box.TopRight()[0] - box.BottomLeft()[0]) / (nx - 1), (box.TopRight()[1] - box.BottomLeft()[1]) / (ny - 1));
	cellSize = Box2D(box.BottomLeft(), box.BottomLeft() + celldiagonal).Size().x; // We only consider squared heightfields

	matterToMove = 0.1f;
}


DuneSediment::DuneSediment( const DuneSediment lowDune, const int factor)
{
	//box = Box2D(Vector2(lowDune.get_box(0)*factor), Vector2(lowDune.get_box(1) * factor)),
	//box = bbox;
	//nx = resx*factor;
	//ny = resy*factor;
	//wind = w;

	box = lowDune.getBox();
	nx = lowDune.getNx() * factor;
	ny = lowDune.getNy() * factor;
	wind = lowDune.getWind();

	bedrock = ScalarField2D(nx, ny, box, 0.0);
	vegetation = ScalarField2D(nx, ny, box, 0.0);
	sediments = ScalarField2D(nx, ny, box, 0.0);
	for (int i = 0; i < nx; i++)
	{
		int li = floor(i / factor);
		for (int j = 0; j < ny; j++)
		{	
			int lj = floor(j / factor);

			Vector2 p = bedrock.ArrayVertex(i, j);


			// Vegetation
			// Arbitrary clamped 2D noise - but you can use whatever you want.
			float v = PerlinNoise::fBm(Vector3(li * 7.91247f, lj * 7.91247f, 0.0f), 1.0f, 0.002f, 3) / 1.75f;
			if (v > 0.58f)
				vegetation.Set(i, j, v);

		
			bedrock.Set(i, j, lowDune.Bedrock(li, lj));
			sediments.Set(i, j, lowDune.Sediment(li,lj));
			//vegetation.Set(i, j, lowDune.Sediment(li, lj));
			
			//sediments.Set(i, j, Random::Uniform(rMin, rMax));
		}
	}



	Vector2 celldiagonal = Vector2((box.TopRight()[0] - box.BottomLeft()[0]) / (nx - 1), (box.TopRight()[1] - box.BottomLeft()[1]) / (ny - 1));
	cellSize = Box2D(box.BottomLeft(), box.BottomLeft() + celldiagonal).Size().x; // We only consider squared heightfields

	matterToMove = 0.1f;
}


DuneSediment::DuneSediment(const DuneSediment lowDune, std::function<bool(const Vector2i&, const DuneSediment&)> cellQuery, const int numLevels) {
	//cellQuery(Vector2i(1, 1));
	
	maxlevels = numLevels;
	box = lowDune.getBox();
	nx = lowDune.getNx();// *factor;
	ny = lowDune.getNy();// *factor;
	wind = lowDune.getWind();



	bedrock = ScalarField2D(nx, ny, box, 0.0);
	vegetation = ScalarField2D(nx, ny, box, 0.0);
	sediments = ScalarField2D(nx, ny, box, 0.0);

	bedrock = lowDune.bedrock;
	sediments = lowDune.sediments;
	vegetation = lowDune.vegetation;

	//std::cerr << "chkpt1" << std::endl;
	Quad* q1 = new Quad(Vector2(0,0), Vector2(0, 0), NULL, false, -1);
	q1->qroots = &roots;
	q1->nx = nx;
	q1->ny = ny;

	for (int i = 0; i < nx; i++) {
		roots.push_back(std::vector<Quad*>());
		for (int j = 0; j < ny; j++)
			roots[i].push_back(NULL);
	}
#pragma omp parallel num_threads(8)
	{
#pragma omp for
		for (int i = 0; i < nx; i++)
		{
			//std::cerr << "chkpt2" << std::endl;
			//roots.push_back(std::vector<Quad*>());
			//int li = floor(i / factor);
			for (int j = 0; j < ny; j++)
			{
				//std::cerr << "chkpt3" << std::endl;
				Vector2 botL, topR;
				botL = Vector2(float(i), float(j));
				topR = Vector2(float(i + 1), float(j + 1));
				int lvl = 0;
				std::map<char, float> vals;
				vals['b'] = lowDune.Bedrock(i, j);
				vals['s'] = lowDune.Sediment(i, j);
				vals['v'] = PerlinNoise::fBm(Vector3(i * 7.91247f, j * 7.91247f, 0.0f), 1.0f, 0.002f, 3) / 1.75f;
				//lowDune.vegetation(i, j);


				//cellQuery(Vector2i(i, j))
				if (cellQuery(Vector2i(i, j), lowDune)) {
					Quad* qd = new Quad(botL, topR, NULL, false, 0);
					//roots[i].push_back(&Quad(botL, topR, NULL, false, 0)  );
					//roots[i].push_back(qd);
					roots[i][j] = qd;
					roots[i][j]->insert(numLevels, vals);


				}
				else {
					//roots[i].push_back(&Quad(botL, topR, NULL, true, 0));
					Quad* qd = new Quad(botL, topR, NULL, true, 0);
					//roots[i].push_back(qd);
					roots[i][j] = qd;
					roots[i][j]->fetchData()->vals = vals;
					roots[i][j]->fetchData()->qd = roots[i][j];
				}
			}//end for ny
		} ////end for nx
	}


	Vector2 celldiagonal = Vector2((box.TopRight()[0] - box.BottomLeft()[0]) / (nx - 1), (box.TopRight()[1] - box.BottomLeft()[1]) / (ny - 1));
	cellSize = Box2D(box.BottomLeft(), box.BottomLeft() + celldiagonal).Size().x; // We only consider squared heightfields

	matterToMove = 0.1f;

}


/*!
\brief Todo
*/
DuneSediment::~DuneSediment()
{

}

/*!
\brief Todo
*/
void DuneSediment::ExportObj(const std::string& url) const
{
	// Clear old data
	
	std::vector<Vector3> vertices;
	std::vector<Vector3> normals;
	std::vector<int> indices;

	// Vertices & UVs & Normals
	normals.resize(nx * ny, Vector3(0));
	vertices.resize(nx * ny, Vector3(0));
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			int id = ToIndex1D(i, j);
			normals[id] = -Normalize(Vector2(bedrock.Gradient(i, j) + sediments.Gradient(i, j)).ToVector3(-2.0f));
			vertices[id] = Vector3(
				box[0][0] + i * (box[1][0] - box[0][0]) / (nx - 1),
				Height(i, j),
				box[0][1] + j * (box[1][1] - box[0][1]) / (ny - 1)
			);
		}
	}

	// Triangles
	int c = 0;
	int vertexArrayLength = ny * nx;
	std::queue <std::string> mtlqueue;

	while (c < vertexArrayLength - nx - 1)
	{
		if (c == 0 || (((c + 1) % nx != 0) && c <= vertexArrayLength - nx))
		{
			indices.push_back(c + nx + 1);
			indices.push_back(c + nx);
			indices.push_back(c);

			float query = sediments.Get(c + nx + 1) + sediments.Get(c + nx) + sediments.Get(c) / 3;
			if (query < 0.45)
				mtlqueue.push("desert_brown");
			else {
				//if (query < 0.97)
					mtlqueue.push("desert_yellow");
				//else
					//mtlqueue.push("highlight_blue");
			}

			indices.push_back(c);
			indices.push_back(c + 1);
			indices.push_back(c + nx + 1);

			query = sediments.Get(c + nx + 1) + sediments.Get(c + nx) + sediments.Get(c) / 3;
			if (query < 0.45)
				mtlqueue.push("desert_brown");
			else {
				//if (query < 0.97)
					mtlqueue.push("desert_yellow");
				//else
					//mtlqueue.push("highlight_blue");
			}
		}
		c++;
	}

	// Export as .obj file
	std::ofstream out;	
	out.open(url);
	if (out.is_open() == false)
		return;
	out << "mtllib materials.mtl" << std::endl;
	out << "g " << "Obj" << std::endl;

	for (int i = 0; i < vertices.size(); i++)
		out << "v " << vertices.at(i).x << " " << vertices.at(i).y << " " << vertices.at(i).z << '\n';
	for (int i = 0; i < normals.size(); i++)
		out << "vn " << normals.at(i).x << " " << normals.at(i).z << " " << normals.at(i).y << '\n';
	
	
	for (int i = 0; i < indices.size(); i += 3)
	{
		

		out << "usemtl " + mtlqueue.front() << std::endl;
		mtlqueue.pop();
		out << "f " << indices.at(i) + 1 << "//" << indices.at(i) + 1
			<< " " << indices.at(i + 1) + 1 << "//" << indices.at(i + 1) + 1
			<< " " << indices.at(i + 2) + 1 << "//" << indices.at(i + 2) + 1
			<< '\n';
	}
	out.close();
}

void DuneSediment::ExportQTObj(const std::string& url) const
{
	// Clear old data
	std::vector<Vector3> vertices;
	std::vector<Vector3> normals;
	std::vector<int> indices;
	std::vector<float> cqVals; //color query values

	//int maxLevels = max;
	std::vector<std::vector<int>> indexMap((nx)*(std::pow(2,maxlevels) ), std::vector<int>((ny) * (std::pow(2, maxlevels)) ));


	int id=-1;
	std::vector<CellData*> cdList;

	//int NX = 4, NY = 4;
	for (int i = 0; i < nx; i++)
	{

		for (int j = 0; j < ny; j++)
		{
			//std::cerr << "chkpt2.2" << std::endl;
			cdList = roots[i][j]->allLeafCellData();

			for (int k = 0; k < cdList.size(); k++) {


				//ScalarField2D sf;
				//std::cerr << "grad check  " << bedrock.Gradient(1, 1) <<std::endl;
				cqVals.push_back(cdList[k]->vals['s']);
				
				normals.push_back(Normalize(Vector2(bedrock.Gradient(cdList[k]->qd, 'b') + sediments.Gradient(cdList[k]->qd, 's')).ToVector3(-2.0f)));
				//normals.push_back(-Normalize(Vector2(bedrock.Gradient(cdList[k]->qd, 'b') + sediments.Gradient(cdList[k]->qd, 's')).ToVector3(-2.0f)));

				vertices.push_back(Vector3(
					//box[0][0] + i * (box[1][0] - box[0][0]) / (nx - 1),
					box[0][0] + (cdList[k]->centre.x) * (box[1][0] - box[0][0]) / (nx - 1),
					Height(cdList[k]->qd),
					//box[0][1] + j * (box[1][1] - box[0][1]) / (ny - 1)
					box[0][0] + (cdList[k]->centre.y) * (box[1][1] - box[0][1]) / (ny - 1)
				));

				id += 1;

				float fracx = cdList[k]->botLeft.x - float(i);
				float fracy = cdList[k]->botLeft.y - float(j);

				long basex = int( (std::pow(2, maxlevels-1) * fracx)  + 0.00000001);
				basex += (std::pow(2, maxlevels-1) * i);
				long basey = int( (std::pow(2, maxlevels-1) * fracy)  + 0.00000001);
				basey += (std::pow(2, maxlevels-1) * j);
				int pseudoCells = std::pow(2, maxlevels -1 - cdList[k]->qd->level);

				for (long p = 0; p < pseudoCells; p++) {
					for (long q = 0; q < pseudoCells; q++) {
						indexMap[basex + p][basey + q] = id;
					}
				}

				//int lvlDiff = maxlevels - cdList[k]->qd->level;



			}

			

		}
	}

	// Triangles
	
	//int c = 0;
	float colQuery;
	int vertexArrayLength = vertices.size();
	std::queue <std::string> mtlqueue;
	int a, b, c, d;
	for (long i = 0; i < nx * std::pow(2, maxlevels -1 ) - 1; i++) {
		for (long j = 0; j < ny * std::pow(2, maxlevels -1 ) -1; j++) {
			a = indexMap[i][j];
			b = indexMap[i][j+1];
			c = indexMap[i+1][j+1];

			d = indexMap[i+1][j];

			
			if ((a - b) * (b - c) * (c - a) != 0) { //ensures triangle is niether point nor line 
				indices.push_back(c); //ccw
				indices.push_back(b);
				indices.push_back(a);
				//mtlqueue.push("highlight_blue");

				colQuery = (cqVals[a]+ cqVals[b]+ cqVals[c]) / 3;
				if (colQuery < 0.45)
					mtlqueue.push("desert_brown");
				else {
						mtlqueue.push("desert_yellow");

				}

			}


			if ((d - a) * (a - c) * (c - d) != 0) {
				indices.push_back(a);
				indices.push_back(d);
				indices.push_back(c);

				colQuery = (cqVals[b] + cqVals[c] + cqVals[d]) / 3;
				if (colQuery < 0.45)
					mtlqueue.push("desert_brown");
				else {
						mtlqueue.push("desert_yellow");
				}
			}
		}
	}

	/*int c = 0;
	while (c < vertexArrayLength - nx - 1)
	{
		if (c == 0 || (((c + 1) % nx != 0) && c <= vertexArrayLength - nx))
		{
			indices.push_back(c + nx + 1);
			indices.push_back(c + nx);
			indices.push_back(c);

			float query = sediments.Get(c + nx + 1) + sediments.Get(c + nx) + sediments.Get(c) / 3;
			if (query < 0.20)
				mtlqueue.push("desert_brown");
			else {
				if (query < 0.97)
					mtlqueue.push("desert_yellow");
				else
					mtlqueue.push("highlight_blue");
			}

			indices.push_back(c);
			indices.push_back(c + 1);
			indices.push_back(c + nx + 1);

			query = sediments.Get(c + nx + 1) + sediments.Get(c + nx) + sediments.Get(c) / 3;
			if (query < 0.20)
				mtlqueue.push("desert_brown");
			else {
				if (query < 0.97)
					mtlqueue.push("desert_yellow");
				else
					mtlqueue.push("highlight_blue");
			}
		}
		c++;
	}*/

	// Export as .obj file
	std::ofstream out;
	out.open(url);
	if (out.is_open() == false)
		return;
	out << "mtllib materials.mtl" << std::endl;
	out << "g " << "Obj" << std::endl;

	for (int i = 0; i < vertices.size(); i++)
		out << "v " << vertices.at(i).x << " " << vertices.at(i).y << " " << vertices.at(i).z << '\n';
	for (int i = 0; i < normals.size(); i++)
		out << "vn " << normals.at(i).x << " " << normals.at(i).z << " " << normals.at(i).y << '\n';


	for (int i = 0; i < indices.size(); i += 3)
	{


		out << "usemtl " + mtlqueue.front() << std::endl;
		mtlqueue.pop();
		out << "f " << indices.at(i) + 1 << "//" << indices.at(i) + 1
			<< " " << indices.at(i + 1) + 1 << "//" << indices.at(i + 1) + 1
			<< " " << indices.at(i + 2) + 1 << "//" << indices.at(i + 2) + 1
			<< '\n';
	}
	out.close();
}