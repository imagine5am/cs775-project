/*
	This is an example implementation of some of the results described in the paper "Desertscapes Simulation".
	No real time visualization in order to reduce dependencies. Running the program will output 3 .obj files which 
	can then be visualized in another application (Blender, MeshLab).

	If you have any questions, you can contact me at:
	axel(dot)paris(at)liris(dot)cnrs(dot)fr
*/

#include <cstdio>

//#include "ctime"
#include "desert.h"

#include "logging.h"

#include <queue>
//#include <string>



void create_mtl_file_deprecated()
{
	std::ofstream out;
	out.open("materials.mtl");
	if (out.is_open() == false)
		return;
	std::cout << "Writing materials.mtl" << std::endl << std::endl;
	out << "newmtl desert_yellow" << std::endl;
	out << "Ka 1 1 0" << std::endl; // AMBIENT COLOR - R G B
	out << "Kd 1 1 0" << std::endl; // DIFFUSE COLOR - R G B
	out << "d 1" << std::endl; // TRANSPARENCY - 0 (transparent) to 1 (opaque)
	out << "illum 1" << std::endl; // ILLUMINATION MODE
	/*
	ILLUMINATION MODES
	0. Color on and Ambient off
	1. Color on and Ambient on
	2. Highlight on
	3. Reflection on and Ray trace on
	4. Transparency: Glass on, Reflection: Ray trace on
	5. Reflection: Fresnel on and Ray trace on
	6. Transparency: Refraction on, Reflection: Fresnel off and Ray trace on
	7. Transparency: Refraction on, Reflection: Fresnel on and Ray trace on
	8. Reflection on and Ray trace off
	9. Transparency: Glass on, Reflection: Ray trace off
	10. Casts shadows onto invisible surfaces
	*/


	std::cout << "Writing materials.mtl" << std::endl << std::endl;
	out << "newmtl highlight_blue" << std::endl;
	out << "Ka 1 1 0" << std::endl; // AMBIENT COLOR - R G B
	out << "Kd 0 0 1" << std::endl; // DIFFUSE COLOR - R G B
	out << "d 1" << std::endl; // TRANSPARENCY - 0 (transparent) to 1 (opaque)
	out << "illum 1" << std::endl; // ILLUMINATION MODE

	out.close();
	//out << "g " << "Obj" << std::endl;
}

void create_mtl_file()
{
	std::ofstream out;
	out.open("materials.mtl");
	if (out.is_open() == false)
		return;

	std::cout << "Writing materials.mtl" << std::endl << std::endl;
	out << "newmtl desert_yellow" << std::endl;
	out << "Ka 0.94 0.784 0.553" << std::endl; // AMBIENT COLOR - R G B
	out << "Kd 0.94 0.784 0.553" << std::endl; // DIFFUSE COLOR - R G B
	out << "d 1" << std::endl; // TRANSPARENCY - 0 (transparent) to 1 (opaque)
	out << "illum 1" << std::endl; // ILLUMINATION MODE

	out << "newmtl desert_brown" << std::endl;
	out << "Ka 0.94 0.784 0.553" << std::endl; // AMBIENT COLOR - R G B
	out << "Kd 0.98 0.62 0.37" << std::endl; // DIFFUSE COLOR - R G B
	out << "d 1" << std::endl; // TRANSPARENCY - 0 (transparent) to 1 (opaque)
	out << "illum 1" << std::endl; // ILLUMINATION MODE

	out << "newmtl highlight_blue" << std::endl;
	out << "Ka 1 1 0" << std::endl; // AMBIENT COLOR - R G B
	out << "Kd 0 0 1" << std::endl; // DIFFUSE COLOR - R G B
	out << "d 0.5" << std::endl; // TRANSPARENCY - 0 (transparent) to 1 (opaque)
	out << "illum 1" << std::endl; // ILLUMINATION MODE

	/*
	ILLUMINATION MODES
	0. Color on and Ambient off
	1. Color on and Ambient on
	2. Highlight on
	3. Reflection on and Ray trace on
	4. Transparency: Glass on, Reflection: Ray trace on
	5. Reflection: Fresnel on and Ray trace on
	6. Transparency: Refraction on, Reflection: Fresnel off and Ray trace on
	7. Transparency: Refraction on, Reflection: Fresnel on and Ray trace on
	8. Reflection on and Ray trace off
	9. Transparency: Glass on, Reflection: Ray trace off
	10. Casts shadows onto invisible surfaces
	*/
	out.close();
	//out << "g " << "Obj" << std::endl;
}

bool cellQueryFunction(Vector2i v, DuneSediment lowdune) {
	
	return lowdune.Sediment(v.x, v.y) >= 0.45;
}

/*!
\brief Running this program will export some
meshes similar to the ones seen in the paper.
*/
int main()
{
	time_t start;
	time(&start);

	int backSteps = 10;
	int qsteps = 1;
	std::queue <DuneSediment> duneQue;
	int iter = 300;
	//assertiopn: backSteps < iter

	create_mtl_file();
	std::string logmsg;
	
	time_t t1;
	time_t t2;
	time_t finish;
	
	Logger("--- Fresh Run begins ----");

	//// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes" << std::endl;
	//DuneSediment dune = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0));
	//for (int i = 0; i < 300; i++)
	//	dune.SimulationStepMultiThreadAtomic();
	//dune.ExportObj("transverse.obj");
	//std::cout << "Done 1/4" << std::endl << std::endl;

	//// Barchan dunes appears under similar wind conditions, but lower sand supply
	//std::cout << "Barchan dunes" << std::endl;
	//dune = DuneSediment(Box2D(Vector2(0), Vector2(256)), 0.5, 0.5, Vector2(3, 0));
	//for (int i = 0; i < 300; i++)
	//	dune.SimulationStepMultiThreadAtomic();
	//dune.ExportObj("barchan.obj");
	//std::cout << "Done 2/4" << std::endl;

	//// Yardangs are created by abrasion, activated with a specific flag in our simulation.
	//std::cout << "Yardangs" << std::endl;
	//dune = DuneSediment(Box2D(Vector2(0), Vector2(256)), 0.5, 0.5, Vector2(6, 0));
	//dune.SetAbrasionMode(true);
	//for (int i = 0; i < 600; i++)
	//	dune.SimulationStepMultiThreadAtomic();
	//dune.ExportObj("yardangs.obj");
	//std::cout << "Done 3/4" << std::endl;

	// Nabkha are created under the influence of vegetation
	//std::cout << "Nabkha" << std::endl;
	//DuneSediment dune = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0));
	//dune.SetVegetationMode(true);
	//for (int i = 0; i < 300; i++)
	//	dune.SimulationStepMultiThreadAtomic();
	//dune.ExportObj("nabkha.obj");
	//std::cout << "Done 4/4" << std::endl;


	//// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes 256 by 256" << std::endl;
	//DuneSediment dune1 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0));
	//for (int i = 0; i < 300; i++)
	//	dune1.SimulationStepMultiThreadAtomic();
	//dune1.ExportObj("transverse256.obj");
	//std::cout << "Done 1/2" << std::endl << std::endl;

	//// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes 512 by 512" << std::endl;
	//DuneSediment dune2 = DuneSediment(Box2D(Vector2(0), Vector2(512)), 1.0, 3.0, Vector2(3, 0));
	//for (int i = 0; i < 300; i++)
	//	dune2.SimulationStepMultiThreadAtomic();
	//dune2.ExportObj("transverse512.obj");
	//std::cout << "Done 2/2" << std::endl << std::endl;

	// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes res: 256, size 256" << std::endl;
	//DuneSediment dune1 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0),256,256);
	//for (int i = 0; i < 300; i++)
	//	dune1.SimulationStepMultiThreadAtomic();
	//dune1.ExportObj("transverseRes256.obj");
	//std::cout << "Done 1/3" << std::endl << std::endl;


	//// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes res: 512, size 256" << std::endl;
	//DuneSediment dune2 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0),512,512);
	//for (int i = 0; i < 300; i++)
	//	dune2.SimulationStepMultiThreadAtomic();
	//dune2.ExportObj("transverseRes512.obj");
	//std::cout << "Done 2/3" << std::endl << std::endl;

	// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes res: 1024, size 256" << std::endl;
	//DuneSediment dune2 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0), 1024, 1024);
	//for (int i = 0; i < 300; i++)
	//	dune2.SimulationStepMultiThreadAtomic();
	//dune2.ExportObj("transverseRes1024.obj");
	//std::cout << "Done 3/3" << std::endl << std::endl;


	//-------------------------------//
	// Transverse dunes are created under unimodal wind, as well as medium to high sand supply
	//std::cout << "Transverse dunes res: 256, size 256" << std::endl;
	//DuneSediment dune1 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0),256,256);
	//for (int i = 0; i < 300; i++)
	//	dune1.SimulationStepMultiThreadAtomic();
	//dune1.ExportObj("transverseIndep.obj");
	//std::cout << "Done 1/2" << std::endl << std::endl;
	//Logger(" Low res dune done.");

	//time_t t1;
	//time(&t1);
	//logmsg = std::to_string(difftime(t1, start));
	//logmsg = "Time required = " + logmsg + " seconds";
	//Logger(logmsg);

	////-------------------------------//
	//Logger("--- Start of next ---");
	//std::cout << "Transverse dunes hyprid res: 256+512, size 256" << std::endl;
	//DuneSediment dune2 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0), 256, 256);
	//for (int i = 0; i < 290; i++)
	//	dune2.SimulationStepMultiThreadAtomic();
	////dune2.ExportObj("transverseHybridRes.obj");
	//DuneSediment dune2h = DuneSediment( dune2,2);
	//for (int i = 0; i < 10; i++)
	//	dune2h.SimulationStepMultiThreadAtomic();
	//dune2h.ExportObj("transverseHybridRes.obj");
	//std::cout << "Done 2/2" << std::endl << std::endl;
	//time_t t2;
	//time(&t2);
	//Logger("time hybrid dune done.");
	//logmsg = std::to_string(difftime(t2, t1));
	//logmsg = "Time required = " + logmsg + " seconds";
	//Logger(logmsg);

	////-------------------------------//
	
	//Logger("--- Start of next ---");
	//std::cout << "Transverse dunes hyprid res with back snaps: 256+512, size 256" << std::endl;


	//DuneSediment dune2 = DuneSediment(Box2D(Vector2(0), Vector2(256)), 1.0, 3.0, Vector2(3, 0), 256, 256);
	//for (int i = 0; i < backSteps; i++) {
	//	dune2.SimulationStepMultiThreadAtomic();
	//	duneQue.push(dune2);
	//	//to check that copy is pushed, not reference -> done!
	//}
	////std::cout << &(duneQue.front()) << std::endl << &(duneQue.back())<<std::endl; check for copy or ref
	//for (int i = backSteps; i < iter; i++) {
	//	dune2.SimulationStepMultiThreadAtomic();
	//	duneQue.pop();
	//	duneQue.push(dune2);
	//}
	////DuneSediment dune2h = DuneSediment(duneQue.front(), 2);
	////for (int i = 0; i < backSteps; i++)
	////	dune2h.SimulationStepMultiThreadAtomic();
	////dune2h.ExportObj("transverseHybridResBS.obj"); //BS as in back snaps


	////std::cout << "Done 2/2" << std::endl << std::endl;
	////time_t t2;
	////time(&t2);
	////Logger("time hybrid dune with backsnaps done.");
	////logmsg = std::to_string(difftime(t2, t1));
	////logmsg = "Time required = " + logmsg + " seconds";
	////Logger(logmsg);
	//dune2.ExportObj("transverseBase.obj");
	//DuneSediment dune2h = DuneSediment(duneQue.front(), 2);
	//for (int i = 0; i < backSteps; i++)
	//	dune2h.SimulationStepMultiThreadAtomic();
	//dune2h.ExportObj("transverseTH.obj");


	//std::cout << "Done 2/2" << std::endl << std::endl;
	//time_t t2;
	//time(&t2);
	//Logger("time hybrid (512) dune with backsnaps done.");
	//logmsg = std::to_string(difftime(t2, t1));
	//logmsg = "Time required = " + logmsg + " seconds";
	//Logger(logmsg);



	//-------------------------------//

	time(&t1);
	std::cout << "Transverse dunes res: 128, size 128" << std::endl;
	DuneSediment dune1 = DuneSediment(Box2D(Vector2(0), Vector2(128)), 1.0, 3.0, Vector2(3, 0), 128, 128);

	for (int i = 0; i < backSteps; i++) {
		dune1.SimulationStepMultiThreadAtomic();
		duneQue.push(dune1);
	}

	for (int i = backSteps; i < iter; i++) {
		dune1.SimulationStepMultiThreadAtomic();
		duneQue.pop();
		duneQue.push(dune1);
	}


	dune1.ExportObj("transverse_basic.obj");
	std::cout << "Done 1/4" << std::endl << std::endl;
	Logger(" Low res dune done.");

	time(&t2);
	logmsg = std::to_string(difftime(t2, t1));
	logmsg = "Time required = " + logmsg + " seconds";
	Logger(logmsg);
	
//-------------------------------//
	time(&t1);
	Logger("--- Start of next ---");
	std::cout << "Transverse dunes: Time Hybrid" << std::endl;
	DuneSediment duneTH = DuneSediment(duneQue.front(),2);
	//duneQue = std::queue<DuneSediment>();

	for (int i = 0; i < backSteps -qsteps ; i++) {
		duneTH.SimulationStepMultiThreadAtomic();
	}
	//DuneSediment duneTHsl = DuneSediment(duneTH);
	//duneQue.push(duneTH);
	for (int i=qsteps ; i < backSteps ; i++) {
		duneTH.SimulationStepMultiThreadAtomic();
	}


	duneTH.ExportObj("transverseTH.obj");
	std::cout << "Done 2/4" << std::endl << std::endl;
	Logger(" Time Hybrid dune done.");

	time(&t2);
	logmsg = std::to_string(difftime(t2, t1));
	logmsg = "Time required = " + logmsg + " seconds";
	Logger(logmsg);
	
//-------------------------------//
	//time_t t2;
	time(&t1);
	Logger("--- Start of next ---");
	std::cout << "Transverse dunes: Quad tree enrichment of basic" << std::endl;


	//DuneSediment duneQT = DuneSediment(duneQue.front(), &cellQueryFunction, 2);
	//DuneSediment duneQT = DuneSediment(duneQue.front(), &cellQueryFunction, 2);
	DuneSediment duneQT = DuneSediment(dune1, &cellQueryFunction, 2);

	std::cout << "Dune QT (NU) generated" << std::endl << std::endl;
	std::cout << "QT simulation in progress" << std::endl << std::endl;
	//for(int i=0; i<10; i++)
	duneQT.SimulationStepQTMultiThreadAtomic(2,false);

	std::cout << "Exporting obj file ..." << std::endl << std::endl;

	duneQT.ExportQTObj("transverseQT.obj");


	std::cout << "Done 3/4" << std::endl << std::endl;

	time(&t2);
	Logger("Non Uniform QT dune done.");
	logmsg = std::to_string(difftime(t2, t1));
	logmsg = "Time required = " + logmsg + " seconds";
	Logger(logmsg);



	//-------------------------------//
	time(&t1);
	Logger("--- Start of next ---");
	std::cout << "Transverse dunes: Quad tree enrichment of TH" << std::endl;


	//DuneSediment duneQT = DuneSediment(duneQue.front(), &cellQueryFunction, 2);
	DuneSediment duneQTonTH = DuneSediment(duneTH, &cellQueryFunction, 2);

	std::cout << "Dune QT (NU) generated" << std::endl << std::endl;
	std::cout << "QT simulation in progress" << std::endl << std::endl;
	//for(int i=0; i<10; i++)
	duneQTonTH.SimulationStepQTMultiThreadAtomic(2,false);

	std::cout << "Exporting obj file ..." << std::endl << std::endl;

	duneQTonTH.ExportQTObj("transverseQTonTH.obj");


	std::cout << "Done 4/4" << std::endl << std::endl;

	time(&t2);
	Logger("Non Uniform QT enrichment of TH done.");
	logmsg = std::to_string(difftime(t2, t1));
	logmsg = "Time required = " + logmsg + " seconds";
	Logger(logmsg);

	
	/////////////////////////////////////////////////////////////////////////////////
	
	time(&finish);
	logmsg = std::to_string( difftime(finish, start) );
	logmsg = "Total Time required = " + logmsg + " seconds";
	Logger(logmsg);
	Logger("--- Run Ends ----");

	return 0;
}
