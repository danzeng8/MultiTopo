// MultiTopo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdint.h>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include <iostream>
#include <queue>
#include <chrono>
#define BOOST_TYPEOF_EMULATION
#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <boost/heap/pairing_heap.hpp>
#include <typeinfo>
#include <stdlib.h>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include <map>
#include <stack>
#include <unordered_map>
#include <bitset>
#include "stats.h"
#include "options.h"
#include "procstatus.h"
#include "timer.h"
#include "ds.h"
#include "bbtree.h"
#include "prep.h"
#include <stdio.h>
#include "util.h"
#include "tbb/parallel_for.h"
#include "tbb/concurrent_vector.h"
#include <tiffio.h>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::heap;
using namespace boost::filesystem;

#pragma region 
typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty; // define edge weight property
typedef boost::property<boost::vertex_name_t, float> VertexWeightProperty; // define node weight property; note that: vertex_index_t is not mutable
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexWeightProperty, EdgeWeightProperty> grapht; // define all the graph properties
typedef boost::graph_traits<grapht>::adjacency_iterator AdjacencyIterator;
typedef boost::graph_traits<grapht >::vertex_descriptor vertex_t;

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, VertexWeightProperty, EdgeWeightProperty> graphL; // define all the graph properties
typedef boost::graph_traits<graphL>::adjacency_iterator AdjacencyIteratorL;
typedef boost::graph_traits<graphL>::vertex_descriptor vertex_tL;

#pragma endregion define graph property 2018Äê3ÔÂ23ÈÕ20:15:27

#include <fstream>


struct node {
	unsigned char type;
	unsigned char inFg;
	int64_t labelCost; //intensity costs: only needed for cuts and fills
	float floatCost;
	int64_t v = 0; int64_t e = 0;  int64_t f = 0;  int64_t c = 0; //cell complex costs: only needed for cuts and fills
	int index;
	bool valid = true;
	bool isArticulate = false;
	bool isArticulateFg = false;
	bool isNew = false;
	int compIndex = -1;
	int tin = 0;
	int low = 0;
	int64_t intensity;
	int overallCompIndexFg = -1; int overallCompIndexBg = -1;
	int level;
	int totalNodeIndex;
	float greatestDiff = 0;
	bool inConflictLower;
	bool inConflictUpper;
};

//#pragma endregion define heaps
struct VertexProps {
	float weight;
};

struct EdgeProps {
	bool strong;
};
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, VertexProps, EdgeProps> Graph;

struct coord {
	short x, y, z;
};

//Use 26-connectivity for foreground, 6-connectivity for background 
std::vector<std::vector<int>> structCube = {
	{-1, -1, -1},
	{-1, -1, 0},
	{-1, -1, 1},
	{-1, 0, -1},
	{-1, 0, 0},
	{-1,0, 1},
	{-1, 1, -1},
	{-1, 1, 0},
	{-1, 1, 1},
	{0, -1, -1},
	{0, -1, 0},
	{0, -1, 1},
	{0, 0, -1},
	{0, 0, 1},
	{0, 1, -1},
	{0, 1, 0},
	{0, 1,1},
	{1, -1, -1},
	{1, -1, 0},
	{1, -1, 1},
	{1, 0, -1},
	{1, 0, 0},
	{1,0, 1},
	{1, 1, -1},
	{1, 1, 0},
	{1, 1, 1} };

std::vector<std::vector<int>> struct18 = {
	{-1, -1, 0},
	{-1, 0, -1},
	{-1, 0, 0},
	{-1,0, 1},
	{-1, 1, 0},
	{0, -1, -1},
	{0, -1, 0},
	{0, -1, 1},
	{0, 0, -1},
	{0, 0, 1},
	{0, 1, -1},
	{0, 1, 0},
	{0, 1,1},
	{1, -1, 0},
	{1, 0, -1},
	{1, 0, 0},
	{1,0, 1},
	{1, 1, 0}
};

std::vector< std::vector<int>> structCross3D = {
	{-1, 0, 0},
	{0, -1, 0},
	{0, 0, -1},
	{0, 0, 1},
	{0, 1, 0},
	{1, 0, 0}
};

std::vector<std::vector<int>> structCubeFull = {
	{-1, -1, -1},
	{-1, -1, 0},
	{-1, -1, 1},
	{-1, 0, -1},
	{-1, 0, 0},
	{-1,0, 1},
	{-1, 1, -1},
	{-1, 1, 0},
	{-1, 1, 1},
	{0, -1, -1},
	{0, -1,0},
	{0, -1, 1},
	{0, 0, -1},
	{0, 0, 0},
	{0, 0, 1},
	{0, 1, -1},
	{0, 1,0},
	{0, 1, 1},
	{1, -1, -1},
	{1, -1, 0},
	{1, -1, 1},
	{1, 0, -1},
	{1, 0, 0},
	{1, 0, 1},
	{1, 1, -1},
	{1, 1, 0},
	{1, 1, 1}
};

std::vector<std::vector<int>> structSquare = {
	{-1, -1, 0},
	{-1, 0, 0},
	{-1, 1, 0},
	{0, -1, 0},
	{0, 1, 0},
	{1, -1, 0},
	{1, 0, 0},
	{1, 1, 0},
};


bool inStruct18[3][3][3] = {
	{{false,true,false},{true,true,true},{false,true,false}},
	{{true,true,true},{true,true,true},{true,true,true}},
	{{false,true,false},{true,true,true},{false,true,false}}
};

struct weightedCoord {
	short x, y, z;  float intensityDiff; //Records intensity difference of voxel compared to shape
	bool operator<(const weightedCoord& rhs) const
	{
		return intensityDiff < rhs.intensityDiff;
	}
};

std::vector< std::vector<int> > cubeFrontMask = {
	{0,0,0},
	{1,0,0},
	{0,1,0},
	{1,1,0},
	{1,0,1},
	{1,1,1},
	{0,1,1},
	{0,0,1}
};

typedef std::tuple<int, int, int, int, int> coor;

struct key_hash : public std::unary_function<coor, int>
{
	std::size_t operator()(const coor& k) const
	{
		return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
	}
};

typedef std::unordered_map<const coor, int, key_hash, std::equal_to<coor>> mapC;

struct Coordinate {
	int x, y, z;
	Coordinate(int x, int y, int z) : x(x), y(y), z(z) {}

	bool operator<(const Coordinate& coord) const {
		if (x < coord.x) return true;
		if (x > coord.x) return false;
		//x == coord.x
		if (y < coord.y) return true;
		if (y > coord.y) return false;
		//x == coord.x && y == coord.y
		if (z < coord.z) return true;
		if (z > coord.z) return false;
		//*this == coord
		return false;
	}

	bool operator==(const Coordinate& coord) const {
		if (x == coord.x && y == coord.y && z == coord.z)
			return true;
		return false;
	}

	inline bool isInRange(Coordinate coord, int range) const {
		if (pow(coord.x - this->x, 2) + pow(coord.y - this->y, 2) + pow(coord.z - this->z, 2) <= range * range)
			return true;
		return false;
	}
};

uint32_t unvisited = 2097152;

class Compare
{
public:
	bool operator() (node& a, node& b)
	{
		return abs(a.labelCost) > abs(b.labelCost);
	}
};

enum type {
	CORE = 0,
	N = 1,
	CUT = 2,
	FILL = 3,
	MIDT = 4,
	PINODE = 5,
	HYPERNODE = 6
};

int numDigits(int x)
{
	//Accomodate image volumes with up to 10^5 slices
	x = abs(x);
	return (x < 10 ? 1 :
		(x < 100 ? 2 :
			(x < 1000 ? 3 :
				(x < 10000 ? 4 : 5))));
}

int coordToIndex(int x, int y, int z, int width, int height, int depth) {
	return x + height * (y + depth * z);
}

void getCoordFromIndex(int index, int& x, int& y, int& z, int width, int height, int depth) {
	x = index - (height * (y + depth * z));
	y = ((index - x) / height) - (depth * z);
	z = (((index - x) / height) - y) / depth;
}


bool IsBitSet(uint32_t num, int bit)
{
	return 1 == ((num >> bit) & 1);
}

bool Visited(uint32_t num)
{
	return 1 == ((num >> 31) & 1);
}

//Extract k bits from position p
uint32_t bitExtracted(uint32_t number, int k, int p)
{
	uint32_t mask;
	mask = ((1 << k) - 1) << p;
	return number & mask;
}

uint32_t Label(uint32_t number) {
	return bitExtracted(number, 22, 0);
}

void nBitToX(uint32_t& number, int n, int x) {
	number ^= (-x ^ number) & (1UL << n);
}

void changeLabel(uint32_t& n1, uint32_t& n2) {
	for (int i = 0; i < 23; i++) {
		//Set ith bit of n2 to ith bit of n1
		nBitToX(n2, i, ((n1 >> i) & 1));
	}
}

void setVisitedFlag(uint32_t& n, int i) {
	nBitToX(n, 31, i);
}

void loadImages(std::vector<float*>& g_Image3D, int& numSlices, std::string dir, int& width, int& height,
	int& maxVal, int& minVal, vector<int>& maxVoxel, int inFileType) {

	//File format is .tif
	if (inFileType == 1) {
		TIFF* tif = TIFFOpen(dir.c_str(), "r");
		if (tif) {
			int dircount = 0;
			do {
				tdata_t buf;
				uint32 row;
				TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
				TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
				buf = _TIFFmalloc(TIFFScanlineSize(tif));
				float* imageF = new float[width * height];
				for (row = 0; row < height; row++) {
					TIFFReadScanline(tif, buf, row);
					float* imgF = (float*)buf;
					for (int i = 0; i < width; i++) {
						//cout << imgF[i] << endl;
						imageF[i + row * width] = imgF[i];

						if ((int)imageF[i + row * width] > maxVal) {
							maxVal = (int)imageF[i + row * width];
							maxVoxel = { i,(int)row,dircount };
							cout << " set to " << i << " " << row << " " << dircount << endl;
						}
					}
				}
				g_Image3D.push_back(imageF);
				_TIFFfree(buf);

				dircount++;
			} while (TIFFReadDirectory(tif));
			TIFFClose(tif);
			numSlices = g_Image3D.size();
			cout << width << " " << height << " " << numSlices << endl;
		}
		else {
			std::cout << "Error Reading Tif File  (not existent, not accessible or no TIFF file)" << std::endl;
		}
		return;
	}

	path p(dir);
	cout << "Loading files " << endl;
	int ct = 0;
	for (auto i = directory_iterator(p); i != directory_iterator(); i++)
	{
		if (!is_directory(i->path()))
		{
			stbi_uc* g_image = NULL;

			int bpp;
			g_image = stbi_load(i->path().string().c_str(), &width, &height, &bpp, 1);
			float* imageF = new float[width * height];
			for (int x = 0; x < width; x++) {
				for (int y = 0; y < height; y++) {
					imageF[x + y * width] = (float)g_image[x + y * width];
					if ((int)g_image[x + y * width] > maxVal) {
						maxVal = (int)g_image[x + y * width];
						maxVoxel = { x,y,ct };
					}
					minVal = min((int)g_image[x + y * width], minVal);
				}
			}
			ct++;
			g_Image3D.push_back(imageF);
		}
		else {
			continue;
		}
	}
	numSlices = g_Image3D.size();
}

void parseArguments(int argc, string& inFile, string& outFile, char** argv, vector<float>& shapes, float& vizEps, int& hypernodeSize, int& productThresh,
	int& globalSteinerTime, int& localSteinerTime, int& bbTime, bool& help, bool& shapeTopo,
	int& beamSize, int& oneConflict, float& epsilon, bool& propagate, int& geomCost, int& distMode,
	bool& outputBB, string& bbFile, bool& totalTopo, int& inFileType, int& outFileType, bool& createLevels, vector<int>& shapeIndices) {
	for (int i = 0; i < argc; i++) {
		if (((string)argv[i]).substr(0, 2) == "--") {
			string arg = (string)argv[i];
			arg.erase(0, 2);

			if (arg == "in") {
				inFile = (string)argv[i + 1]; //Input image slices directory
			}
			if (arg == "out") {
				outFile = (string)argv[i + 1]; //Output image slices directory
			}

			if (arg == "S") {
				std::stringstream ss(argv[i + 1]);
				std::string s;
				while (std::getline(ss, s, ' ')) {
					shapes.push_back(std::stof(s));
				}
			}

			if (arg == "indices") {
				std::stringstream ss(argv[i + 1]);
				std::string s;
				while (std::getline(ss, s, ' ')) {
					shapeIndices.push_back(std::stoi(s));
				}
			}

			if (arg == "shapeTopo") {
				shapeTopo = true;
			}

			if (arg == "help") {
				cout << "help menu: under construction" << endl;
			}

			if (arg == "beam") {
				beamSize = stoi(argv[i + 1]);
			}

			if (arg == "oneConflict") {
				oneConflict = stoi(argv[i + 1]);
			}
			if (arg == "productThresh") {
				productThresh = stoi(argv[i + 1]);
			}

			if (arg == "hypernodeSize") {
				hypernodeSize = stoi(argv[i + 1]);
			}

			if (arg == "epsilon") {
				epsilon = stof(argv[i + 1]);
			}

			if (arg == "propagate") {
				cout << "propagate mode " << endl;
				propagate = true;
			}

			if (arg == "geomCost") {
				geomCost = stoi(argv[i + 1]);
			}

			if (arg == "distMode") {
				distMode = stoi(argv[i + 1]);
			}

			if (arg == "outputBB") {
				outputBB = true;
				bbFile = argv[i + 1];

			}

			if (arg == "totalTopo") {
				totalTopo = true;
			}

			if (arg == "inFileType") {
				inFileType = stoi(argv[i + 1]);
			}

			if (arg == "outFileType") {
				outFileType = stoi(argv[i + 1]);
			}

			if (arg == "localSteinerTime") {
				localSteinerTime = stoi(argv[i + 1]);
			}

			if (arg == "globalSteinerTime") {
				globalSteinerTime = stoi(argv[i + 1]);
			}

			if (arg == "createLevels") {
				createLevels = true;
			}
		}
	}
}

//Flood fill component inside core, starting from above the core threshold C 
void floodFillCore(std::vector<float*>& g_Image3D, std::vector<float*>& distC, int distMode, std::vector<uint32_t*>& labels, int numSlices, int width, int height, float C, float S, int x, int y, int z, int& labelCtI, priority_queue<weightedCoord>& corePQ) {
	std::queue<coord> coordinates;
	coord c = { x,y,z };
	std::vector<std::vector<int>> mask = structCross3D;
	coordinates.push(c);
	uint32_t labelCt = (uint32_t)labelCtI;
	labels[z][x + y * width] = labelCt;
	cout << x << " " << y << " " << z << " labeled with " << labelCt << " intensity " << (int)g_Image3D[z][x + y * width] << endl;
	setVisitedFlag(labels[z][x + width * y], 1);
	while (!coordinates.empty()) {
		coord v = coordinates.front();
		setVisitedFlag(labels[v.z][v.x + width * v.y], 1);
		coordinates.pop();
		for (int i = 0; i < mask.size(); i++) { //Using 26-connectivity for foreground
			coord n = { v.x + mask[i][0], v.y + mask[i][1], v.z + mask[i][2] };
			if (n.x >= 0 && n.x < width && n.y >= 0 && n.y < height && n.z >= 0 && n.z < numSlices) { //Bounds check
				if (Label(labels[n.z][n.x + n.y * width]) == unvisited) {
					if ((float)g_Image3D[n.z][n.x + n.y * width] > C) { //neighboring voxel still in core
						coordinates.push(n);
						uint32_t label32 = labelCt;
						changeLabel(label32, labels[n.z][n.x + n.y * width]);
					}
					else {
						//cout << "neighbor gt? " << (float)g_Image3D[n.z][n.x + n.y * width] << " " << S << endl;
						if ((float)g_Image3D[n.z][n.x + n.y * width] > S) { //neighboring voxel in shape, to be expanded during core inflation

							if (Visited(labels[n.z][n.x + width * n.y]) == false) { //place in core priority queue if not in a priority queue already
								float value = ((float)g_Image3D[n.z][n.x + width * n.y] - S);
								if (distMode == 2) {
									value = distC[n.z][n.x + width * n.y];
								}
								weightedCoord wc = { n.x, n.y, n.z, value };
								corePQ.push(wc);
								setVisitedFlag(labels[n.z][n.x + width * n.y], 1);
							}
						}
					}
				}
			}
		}
	}
}

bool isPotentialNeighborToFront(int label, std::vector<node>& nodes) {
	if (label == unvisited) {
		return true;
	}
	if (label >= nodes.size()) {
		return true;
	}
	if (nodes[label].type != 1) {
		return true;
	}
	return false;
}

bool isType(int label, int qtype, std::vector<node>& nodes) {
	if (label == unvisited) {
		return false;
	}
	if (label > nodes.size()) {
		return false;
	}
	if ((int)(nodes[label].type) == qtype) {
		return true;
	}
	return false;

}

//Flood fill component inside complement of neighborhood, starting from below the neighborhood threshold N 
void floodFillNeighborhood(std::vector<float*>& g_Image3D, std::vector<float*>& distC, int distMode, std::vector<uint32_t*>& labels, int numSlices, int width, int height, float S, float N, int x, int y, int z, int& labelCtI, priority_queue< weightedCoord>& nPQ, std::vector<node>& nodes) {
	std::queue<coord> coordinates; //std::vector<bool *> & inPQ,
	coord c = { x,y,z };
	std::vector<std::vector<int>> mask = structCube;
	coordinates.push(c);
	uint32_t labelCt = (uint32_t)labelCtI;
	labels[z][x + y * width] = labelCt;
	setVisitedFlag(labels[z][x + width * y], 1);
	std::vector<std::vector<bool>> newPQFlag;
	for (int i = 0; i < numSlices; i++) {
		std::vector<bool> layer(width * height, false);
		newPQFlag.push_back(layer);
	}
	while (!coordinates.empty()) {
		coord v = coordinates.front();
		coordinates.pop();
		setVisitedFlag(labels[v.z][v.x + width * v.y], 1);
		for (int i = 0; i < mask.size(); i++) { //Using 26-connectivity for background
			coord n = { v.x + mask[i][0], v.y + mask[i][1], v.z + mask[i][2] };
			if (n.x >= 0 && n.x < width && n.y >= 0 && n.y < height && n.z >= 0 && n.z < numSlices) { //Bounds check

				if (Label(labels[n.z][n.x + n.y * width]) == unvisited) {
					if ((float)g_Image3D[n.z][n.x + width * n.y] <= N) { //neighboring voxel still in core


						coordinates.push(n);
						uint32_t label32 = labelCt;
						changeLabel(label32, labels[n.z][n.x + n.y * width]);
					}

					//else {

					//}

				}
				if ((float)g_Image3D[n.z][n.x + width * n.y] > N && (float)g_Image3D[n.z][n.x + width * n.y] <= S) { //neighboring voxel in shape, to be deflated to during neighborhood deflation

					if (!newPQFlag[n.z][n.x + width * n.y]) {
						//if (n.x == 22 && n.y == 56) {
							//cout << "push into next pq " << n.x << " " << n.y << " " << n.z << " " << (float)g_Image3D[n.z][n.x + n.y * width] << " " << S << endl;
						//}
						float value = abs((float)g_Image3D[n.z][n.x + width * n.y] - S);
						if (distMode == 2) {
							value = distC[n.z][n.x + width * n.y];
						}
						weightedCoord wc = { n.x, n.y, n.z, value }; //-abs((float)g_Image3D[n.z][n.x + width * n.y] - S) };
						nPQ.push(wc);
						newPQFlag[n.z][n.x + width * n.y] = true;
					}
				}
				//else {
					//if ((float)g_Image3D[n.z][n.x + width * n.y] != 0) {
					//	cout << (float)g_Image3D[n.z][n.x + width * n.y] << " N " << N << " S " << S << endl;
					//}
				//}
			}
		}
	}

	/**for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int z = 0; z < numSlices; z++) {
				//if (x == 22 && y == 56) {
				////	cout << "condering push into next pq " << x << " " << y << " " << (float)g_Image3D[z][x + y * width] << " " << isPotentialNeighborToFront(Label(labels[z][x + y * width]), nodes) << " " << S << endl;
				//}

				if (isPotentialNeighborToFront(Label(labels[z][x + y * width]), nodes)) {
					if ((float)g_Image3D[z][x + y * width] <= S && !newPQFlag[z][x + y * width]) {
						//if (x == 22 && y == 56) {
						//	cout << "layer 1 " << x << " " << y << " " << (float)g_Image3D[z][x + y * width] << " " << isPotentialNeighborToFront(Label(labels[z][x + y * width]), nodes) << " " << S << endl;
						//}
						for (int i = 0; i < mask.size(); i++) {
							int nx = x + mask[i][0]; int ny = y + mask[i][1]; int nz = z + mask[i][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
								//if (x == 22 && y == 56) {
								//	cout << "layer 2 " << nx << " " << ny << " " << Label(labels[nz][nx + ny * width]) << " " << isType(Label(labels[nz][nx + ny * width]), 1, nodes) << " " << isType(Label(labels[nz][nx + ny * width]), 0, nodes) << " " << isType(Label(labels[nz][nx + ny * width]), 2, nodes) << " " << S << " " << i << " " << mask.size() << endl;
								//}
								if (isType(Label(labels[nz][nx + ny * width]), 1, nodes)) {
									weightedCoord wc = { x, y, z , abs((float)g_Image3D[z][x + y * width] - S) };
									nPQ.push(wc);
									newPQFlag[z][x + y * width] = true;
									//if (x == 22 && y == 56) {
									//	cout << "pushed " << x << " " << y << " " << (float)g_Image3D[z][x + y * width] << " " << isPotentialNeighborToFront(Label(labels[z][x + y * width]), nodes) << " " << S << endl;
								//	}
									//setVisitedFlag(labels[nz][nx + ny * width], 1);
								}
							}
						}
					}
				}
			}
		}
	}**/
}

int neighborhoodToIndex(std::vector<std::vector<std::vector<int>>>& neighborhood) {
	int ct = 0;
	uint32_t index = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				nBitToX(index, ct, neighborhood[i][j][k]);
				ct += 1;
			}
		}
	}
	return index;
}

bool simple3DDictionary(const std::vector<uint32_t*>& labels, const std::vector<node>& nodes, int x, int y, int z, int numSlices, int width, int height, const std::vector<unsigned char>& simpleDictionary3D, int lType, bool inFg, int conn) {
	int minX = max(x - 1, 0); int minY = max(y - 1, 0); int minZ = max(z - 1, 0);
	int maxX = min(x + 2, width); int maxY = min(y + 2, height); int maxZ = min(z + 2, numSlices);
	//Correct version
	std::vector< std::vector< std::vector<int>>> cubeN(3, std::vector<std::vector<int>>(3, std::vector<int>(3, 0)));
	for (int i = minX; i < maxX; i++) {
		for (int j = minY; j < maxY; j++) {
			for (int k = minZ; k < maxZ; k++) {
				if (Label(labels[k][i + j * width]) != unvisited) {

					if (((int)nodes[Label(labels[k][i + j * width])].type) == lType) {
						cubeN[i - x + 1][j - y + 1][k - z + 1] = 1;

					}
					else {
						cubeN[i - x + 1][j - y + 1][k - z + 1] = 0;
					}
				}
				else {
					cubeN[i - x + 1][j - y + 1][k - z + 1] = 0;
				}
			}
		}
	}


	if (conn == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					cubeN[i][j][k] = 1 - cubeN[i][j][k];
				}
			}
		}
	}
	return ((int)simpleDictionary3D[neighborhoodToIndex(cubeN)]) == 49;

}


std::vector<int> getEulerNumbersFromBImg(const std::vector<uint8_t*>& bImg, int width, int height, int numSlices) {
	int v = 0; int e = 0; int f = 0; int c = 0;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				if (bImg[k][i + j * width] == 1) {
					v += 1;

					if (i + 1 < width) {
						//compare (i,j,k), (i+1, j, k)
						if ((bImg[k][(i + 1) + (width * j)]) == 1) {
							e += 1;
						}
					}

					if (j + 1 < height) {
						//compare (i,j,k), (i, j+1, k)
						if ((bImg[k][i + (width * (j + 1))]) == 1) {
							e += 1;
						}
					}

					if (k + 1 < numSlices) {
						//compare (i,j,k), (i, j+1, k)
						if ((bImg[k + 1][i + (width * (j))]) == 1) {
							e += 1;
						}
					}

					if (i + 1 < width && j + 1 < height) {
						//compare (i,j,k), (i, j+1, k)
						if ((bImg[k][(i + 1) + (width * (j))]) == 1 &&
							(bImg[k][(i)+(width * (j + 1))]) == 1 && (bImg[k][(i + 1) + (width * (j + 1))]) == 1) {
							f += 1;
						}
					}

					if (i + 1 < width && k + 1 < numSlices) {
						//compare (i,j,k), (i, j+1, k)
						if ((bImg[k][(i + 1) + (width * (j))]) == 1 &&
							(bImg[k + 1][(i)+(width * (j))]) == 1 && (bImg[k + 1][(i + 1) + (width * (j))]) == 1) {
							f += 1;

						}
					}

					//add up faces in yz plane
					if (j + 1 < height && k + 1 < numSlices) {
						//compare (i,j,k), (i, j+1, k)
						if (bImg[k][(i)+(width * (j + 1))] == 1 &&
							(bImg[k + 1][(i)+(width * (j))]) == 1 && (bImg[k + 1][(i)+(width * (j + 1))]) == 1
							) {
							f += 1;
						}
					}

					//Add up cubes
					if (i + 1 < width && j + 1 < height && k + 1 < numSlices) {
						bool hasCube = true;
						for (int o = 0; o < cubeFrontMask.size(); o++) {
							int coord[3] = { i + cubeFrontMask[o][0], j + cubeFrontMask[o][1], k + cubeFrontMask[o][2] };
							if (bImg[coord[2]][(coord[0]) + (width * (coord[1]))] == 0) {

								hasCube = false;
								break;
							}
						}
						if (hasCube) {
							c += 1;
						}

					}
				}
			}
		}
	}
	return { v,e,f,c };
}

int labelComponentsB(const std::vector<uint8_t*>& bImg, int conn, int width, int height, int numSlices) {
	int ct = 0;
	std::vector< std::vector<int > > mask;
	if (conn == 0) {
		mask = structCube;
	}
	else {
		mask = structCross3D;
	}

	std::vector< std::vector< std::vector<bool>>> visited(width, std::vector<std::vector<bool>>(height, std::vector<bool>(numSlices, false)));

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				//Unvisited foreground voxels
				if (((int)bImg[k][i + j * width]) == 1 && visited[i][j][k] == false) {
					ct += 1;
					std::queue<Coordinate> q;
					q.push(Coordinate(i, j, k));
					visited[i][j][k] = true;
					while (!q.empty()) {
						Coordinate qp = q.front();
						q.pop();
						for (int s = 0; s < mask.size(); s++) {
							Coordinate np(qp.x + mask[s][0], qp.y + mask[s][1], qp.z + mask[s][2]);
							if (np.x >= 0 && np.x < width && np.y >= 0 && np.y < height && np.z >= 0 && np.z < numSlices) {
								if (((int)bImg[np.z][np.x + (np.y * width)]) == 1 && visited[np.x][np.y][np.z] == false) {
									visited[np.x][np.y][np.z] = true;
									q.push(np);
								}
							}
						}
					}

				}
			}
		}
	}
	return ct;
}


int labelComponentsBBg(const std::vector<uint8_t*>& bImg, int conn, int width, int height, int numSlices, int& bgCompsAll) {
	int bgCt = 0;
	std::vector< std::vector<int > > mask;
	if (conn == 0) {
		mask = structCube;
	}
	else {
		mask = structCross3D;
	}

	std::vector< std::vector< std::vector<bool>>> visited(width, std::vector<std::vector<bool>>(height, std::vector<bool>(numSlices, false)));

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				//Unvisited foreground voxels
				if (((int)bImg[k][i + j * width]) == 0 && visited[i][j][k] == false) {
					bgCompsAll++;
					std::queue<Coordinate> q;
					q.push(Coordinate(i, j, k));
					visited[i][j][k] = true;
					bool isTrueCavity = true;
					if (i == 0 || i == width - 1 || j == 0 || j == height - 1 || k == 0 || k == numSlices - 1) {
						isTrueCavity = false;
					}
					while (!q.empty()) {
						Coordinate qp = q.front();
						q.pop();
						if (qp.x == 0 || qp.x == width - 1 || qp.y == 0 || qp.y == height - 1 || qp.z == 0 || qp.z == numSlices - 1) {
							isTrueCavity = false;
						}
						for (int s = 0; s < mask.size(); s++) {
							Coordinate np(qp.x + mask[s][0], qp.y + mask[s][1], qp.z + mask[s][2]);
							if (np.x >= 0 && np.x < width && np.y >= 0 && np.y < height && np.z >= 0 && np.z < numSlices) {
								if (((int)bImg[np.z][np.x + (np.y * width)]) == 0 && visited[np.x][np.y][np.z] == false) {
									visited[np.x][np.y][np.z] = true;
									q.push(np);
								}
							}
						}
					}
					if (isTrueCavity) {
						bgCt++;
					}
				}
			}
		}
	}
	return bgCt;
}



std::vector<int> getTopoFromBImg(const std::vector<uint8_t*> bImg, int fgconn, int width, int height, int numSlices) {
	int h0 = labelComponentsB(bImg, fgconn, width, height, numSlices);
	int bgCompsAtEdges = 0;
	int h2 = labelComponentsBBg(bImg, 1 - fgconn, width, height, numSlices, bgCompsAtEdges); // - 1
	std::vector<int> eulNums = getEulerNumbersFromBImg(bImg, width, height, numSlices);
	int h1 = h0 + h2 - (eulNums[0] - eulNums[1] + eulNums[2] - eulNums[3]);
	return { h0,h2,h1 };
}

void inflateCoreMulti(std::vector<uint32_t*>& labels, const std::vector<float*>& g_Image3D, const std::vector<float*>& distC, std::vector<node>& nodes, priority_queue<weightedCoord>& corePQ, priority_queue<weightedCoord>& nextPQ, int numSlices, int width,
	int height, float C, float S, float nextS, const std::vector<unsigned char>& simpleDictionary3D, int distMode) {
	priority_queue<weightedCoord> unsimpleVoxels;
	std::vector<std::vector<int>> mask = structCube;

	while (!corePQ.empty()) {
		weightedCoord wc = corePQ.top(); //Pop voxel on core boundary with closest intensity to shape
		corePQ.pop();
		auto sTime = std::chrono::high_resolution_clock::now();
		if (simple3DDictionary(labels, nodes, wc.x, wc.y, wc.z, numSlices, width, height, simpleDictionary3D, 0, true, 1)) {
			auto sTime2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = sTime2 - sTime;
			bool wasSet = false;
			for (int i = 0; i < mask.size(); i++) { //Since 3 x3 neighborhood is simple with respect to inflated core, can pick label of any neighbor to be current one,
				Coordinate n(wc.x + mask[i][0], wc.y + mask[i][1], wc.z + mask[i][2]);
				if (n.x >= 0 && n.y >= 0 && n.z >= 0 && n.x < width && n.y < height && n.z < numSlices) { //bounds check
					if (abs(mask[i][0]) + abs(mask[i][1]) + abs(mask[i][2]) == 1) {
						if (Label(labels[n.z][n.x + n.y * width]) != unvisited) {
							if (((int)nodes[Label(labels[n.z][n.x + n.y * width])].type) == 0) { //current voxel on boundary, which means one of neighbors is already in core
								if (Label(labels[wc.z][wc.x + wc.y * width]) == unvisited) {
									labels[wc.z][wc.x + wc.y * width] = labels[n.z][n.x + n.y * width];
								}
							}
						}
					}

					//if (S == 2) {
						//cout << "expanded to neighbor?" << (float)g_Image3D[n.z][n.x + n.y * width] << " " << S << " " << Visited(labels[n.z][n.x + n.y * width]) << endl;
					//}
					//Continue expanding to unvisited neighbors in correct intensity range
					if ((float)g_Image3D[n.z][n.x + n.y * width] > S && Visited(labels[n.z][n.x + n.y * width]) == false) {
						if (distMode == 2) {
							weightedCoord wnc = { n.x, n.y, n.z,  ((float)distC[n.z][n.x + n.y * width]) }; //distC[n.z][n.x + n.y * width]
							corePQ.push(wnc);
							setVisitedFlag(labels[n.z][n.x + n.y * width], 1);
						}
						else {
							weightedCoord wnc = { n.x, n.y, n.z,  ((float)g_Image3D[n.z][n.x + n.y * width] - S) };
							corePQ.push(wnc);
							setVisitedFlag(labels[n.z][n.x + n.y * width], 1);
						}


					}

				}
			}
		}
		else {
			//Pixel not simple, but may become simple later, so unflag in inPQ table to allow for a later visit
			auto sTime2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = sTime2 - sTime;
			setVisitedFlag(labels[wc.z][wc.x + wc.y * width], 0);
		}
	}

	/**for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int z = 0; z < numSlices; z++) {
				if (g_Image3D[z][x + y * width] > S && g_Image3D[z][x + y * width] <= C) {
					if (Label(labels[z][x + y * width]) == unvisited) {
						if (simple3DDictionary(labels, nodes, x, y, z, numSlices, width, height, simpleDictionary3D, 0, true, 1)) {
							cout << "remaining simple " << endl;
						}
					}
				}
			}
		}
	}**/

	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int z = 0; z < numSlices; z++) {
				if (Label(labels[z][x + y * width]) != unvisited) {
					if (nodes[Label(labels[z][x + y * width])].type == 0) {
						for (int i = 0; i < mask.size(); i++) {
							int nx = x + mask[i][0]; int ny = y + mask[i][1]; int nz = z + mask[i][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
								if (Label(labels[nz][nx + ny * width]) == unvisited) {
									if ((float)g_Image3D[nz][nx + ny * width] > nextS) {
										//g_Image3D
										//if (distMode == 1) {
										float value = ((float)g_Image3D[nz][nx + ny * width] - nextS);
										if (distMode == 2) {
											value = distC[nz][nx + ny * width];
										}

										weightedCoord wnc = { nx, ny, nz , ((float)distC[nz][nx + ny * width] - nextS) };
										nextPQ.push(wnc);
										setVisitedFlag(labels[nz][nx + ny * width], 1);

									}
								}
							}
						}
					}
				}
			}
		}
	}
	/**
	vector< uint8_t*> topoImg;
	for (int s = 0; s < numSlices; s++) {
		uint8_t* label8T = new uint8_t[width * height]();
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				uint32_t label = Label(labels[s][i + j * width]);//levelNewToOldIndices[l][
				//if (((int)origLvlNodes[l][label].type) == 0 || ((int)origLvlNodes[l][label].type) == 2 || ((int)origLvlNodes[l][label].type) == 3) {
				if (label != unvisited) {
					if (nodes[label].type == CORE) {
						label8T[i + j * width] = 1;
					}
				}
				else {
					label8T[i + j * width] = 0;
				}
			}
		}
		topoImg.push_back(label8T);
	}
	std::vector<int> topoNums = getTopoFromBImg(topoImg, 1, width, height, numSlices);
	cout << "topo nums " << topoNums[0] << " " << topoNums[1] << " " << topoNums[2] << endl;**/
}



void deflateNeighborhoodMulti(std::vector<uint32_t*>& labels, const std::vector<float*>& g_Image3D, const std::vector<float*>& priorityImg, std::vector<node>& nodes, priority_queue<weightedCoord>& nPQ, priority_queue<weightedCoord>& nextPQ, int numSlices, int width,
	int height, float N, float S, float nextS, const std::vector<unsigned char>& simpleDictionary3D, int distMode) {
	priority_queue<weightedCoord> unsimpleVoxels;
	std::vector<std::vector<int>> mask = structCube;
	std::vector<std::vector<bool>> newPQFlag;
	for (int i = 0; i < numSlices; i++) {
		std::vector<bool> layer(width * height, false);
		newPQFlag.push_back(layer);
	}
	while (!nPQ.empty()) {
		weightedCoord wc = nPQ.top(); //Pop voxel on core boundary with closest intensity to shape
		nPQ.pop();
		if (S == 136) {
			if (wc.x == 102 && wc.y == 51) {
				cout << "is coord simple? " << simple3DDictionary(labels, nodes, wc.x, wc.y, wc.z, numSlices, width, height, simpleDictionary3D, 1, false, 0) << endl;


			}
		}
		/**if (S == 190) {
			cout << wc.x << " " << wc.y << " " << wc.z << " " << (int)g_Image3D[wc.z][wc.x + wc.y * width] <<
				" " << simple3DDictionary(labels, nodes, wc.x, wc.y, wc.z, numSlices, width, height, simpleDictionary3D, 1, false, 0)
				<<
				endl;
		}**/
		//	cout << wc.x << " " << wc.y << " " << wc.z << " " << (int)g_Image3D[wc.z][wc.x + wc.y * width] <<
			//	" " << S << " " << simple3DDictionary(labels, nodes, wc.x, wc.y, wc.z, numSlices, width, height, simpleDictionary3D, 1, false, 0)
			//	<<
			//	endl;
		auto sTime = std::chrono::high_resolution_clock::now();
		if (simple3DDictionary(labels, nodes, wc.x, wc.y, wc.z, numSlices, width, height, simpleDictionary3D, 1, false, 0)) {
			auto sTime2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = sTime2 - sTime;
			bool wasSet = false;
			for (int i = 0; i < mask.size(); i++) { //Since 3 x3 neighborhood is simple with respect to inflated core, can pick label of any neighbor to be current one,
				Coordinate n(wc.x + mask[i][0], wc.y + mask[i][1], wc.z + mask[i][2]);
				if (n.x >= 0 && n.y >= 0 && n.z >= 0 && n.x < width && n.y < height && n.z < numSlices) { //bounds check
					if (Label(labels[n.z][n.x + n.y * width]) != unvisited) {
						if (((int)nodes[Label(labels[n.z][n.x + n.y * width])].type) == 1) { //current voxel on boundary, which means one of neighbors is already in core
							labels[wc.z][wc.x + wc.y * width] = labels[n.z][n.x + n.y * width];
						}
					}
					//Continue expanding to unvisited neighbors in correct intensity range
					if ((float)g_Image3D[n.z][n.x + n.y * width] <= S && Visited(labels[n.z][n.x + n.y * width]) == false
						&& Label(labels[n.z][n.x + n.y * width]) == unvisited
						) {//&& Label(labels[n.z][n.x + n.y * width]) == unvisited
						if (distMode == 2) {
							weightedCoord wnc = { n.x, n.y, n.z,  abs((float)priorityImg[n.z][n.x + width * n.y]) };
							nPQ.push(wnc);
							setVisitedFlag(labels[n.z][n.x + n.y * width], 1);
						}
						else {
							weightedCoord wnc = { n.x, n.y, n.z,  abs((float)g_Image3D[n.z][n.x + width * n.y] - S) };
							nPQ.push(wnc);
							setVisitedFlag(labels[n.z][n.x + n.y * width], 1);
						}

					}
					//else {
					//	cout << "didn't push " << n.x << " " << n.y << " " << n.z << " I; " << (float)g_Image3D[n.z][n.x + width * n.y] << " S: " << S << " is visited: " << Visited(labels[n.z][n.x + n.y * width])
					//	<< " label: " << Label(labels[n.z][n.x + n.y * width]) << " unvisited: " << unvisited << " are equal:  " << (Label(labels[n.z][n.x + n.y * width]) == unvisited) << endl;
					//}

				}
			}
			if ((Label(labels[wc.z][wc.x + wc.y * width]) == unvisited)) {
				cout << "current remains unlabeled" << endl;
			}
		}
		else {
			//Pixel not simple, but may become simple later, so unflag in inPQ table to allow for a later visit
			auto sTime2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = sTime2 - sTime;
			setVisitedFlag(labels[wc.z][wc.x + wc.y * width], 0);
			newPQFlag[wc.z][wc.x + wc.y * width] = false;
		}

	}
	cout << "finished deflating this round " << endl;
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int z = 0; z < numSlices; z++) {
				//if (nextS == 135) {
					//if (x == 102 && y == 51 && nextS == 136) {
						//cout << "next push? " << isPotentialNeighborToFront(Label(labels[z][x + y * width]), nodes)
							//<< " " << (float)g_Image3D[z][x + y * width] << " " << nextS << endl;
					//}
				//}
				if (isPotentialNeighborToFront(Label(labels[z][x + y * width]), nodes)) {
					if ((float)g_Image3D[z][x + y * width] <= nextS && !newPQFlag[z][x + y * width]) {
						for (int i = 0; i < mask.size(); i++) {
							int nx = x + mask[i][0]; int ny = y + mask[i][1]; int nz = z + mask[i][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {

								if (isType(Label(labels[nz][nx + ny * width]), 1, nodes)) {
									float value = abs((float)g_Image3D[z][x + y * width] - nextS);
									if (distMode == 2) {
										value = priorityImg[z][x + width * y];
									}
									weightedCoord wc = { x, y, z , value };
									nextPQ.push(wc);
									newPQFlag[z][x + y * width] = true;

								}
							}
						}
					}
				}
			}
		}
	}
}

float gradient(int x, int y, int z, const std::vector<float*>& g_Image3D, int width, int height, int numSlices) {
	int xA = min(x + 1, width - 1); int xB = max(x - 1, 0);
	float gX = (float)g_Image3D[z][xA + y * width] - (float)g_Image3D[z][xB + y * width];

	int yA = min(y + 1, height - 1); int yB = max(y - 1, 0);
	float gY = (float)g_Image3D[z][x + yA * width] - (float)g_Image3D[z][x + yB * width];

	int zA = min(z + 1, numSlices - 1); int zB = max(z - 1, 0);
	float gZ = (float)g_Image3D[zA][x + y * width] - (float)g_Image3D[zB][x + y * width];
	return sqrtf((gX * gX) + (gY * gY) + (gZ * gZ));
}

bool ccNeighbor(const std::vector<node>& nodes, int px, int py, int pz, int nx, int ny, int nz, const std::vector<uint32_t*>& labels, int width, int height, int numSlices, const std::vector<float*>& g_Image3D, float S) {
	//is a six neighbor
	bool c1 = abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1;
	if (c1) {
		return true;
	}

	if ((abs(nx - px) + abs(ny - py) + abs(nz - pz) == 2)) {
		if (pz == nz) {
			bool isNFillOrN1 = true;
			//if (Label(labels[pz][nx + (width * py)]) != unvisited) {
			if (((float)g_Image3D[pz][nx + (width * py)]) <= S) {
				isNFillOrN1 = false;
			}
			//}

			bool isNFillOrN2 = true;
			//if (Label(labels[pz][px + (width * ny)]) != unvisited) {
			if (((float)g_Image3D[pz][px + (width * ny)]) <= S) {
				isNFillOrN2 = false;
			}
			//}
			//2D diagonals are not fill or neighbor
			return isNFillOrN1 && isNFillOrN2;
		}

		if (px == nx) {
			bool isNFillOrN1 = true;
			//if (Label(labels[nz][px + (width * py)]) != unvisited) {
			if (((float)g_Image3D[nz][px + (width * py)]) <= S) {
				isNFillOrN1 = false;
			}
			//}

			bool isNFillOrN2 = true;
			//if (Label(labels[pz][px + (width * ny)]) != unvisited) {
			if (((float)g_Image3D[pz][px + (width * ny)]) <= S) {
				isNFillOrN2 = false;
			}
			//}
			//2D diagonals are not fill or neighbor
			return isNFillOrN1 && isNFillOrN2;
		}

		if (py == ny) {
			bool isNFillOrN1 = true;
			//if (Label(labels[nz][px + (width * py)]) != unvisited) {
			if (((float)g_Image3D[nz][px + (width * py)]) <= S) {
				isNFillOrN1 = false;
			}
			//}

			bool isNFillOrN2 = true;
			//if (Label(labels[pz][nx + (width * py)]) != unvisited) {
			if (((float)g_Image3D[pz][nx + (width * py)]) <= S) {
				isNFillOrN2 = false;
			}
			//}
			return isNFillOrN1 && isNFillOrN2;
		}
	}

	if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 3) {
		bool isNFillOrN = true;
		//if (Label(labels[nz][px + (width * ny)]) != unvisited) {
		if (((float)g_Image3D[nz][px + (width * ny)]) <= S) {
			isNFillOrN = false;
		}
		//}
		//if (Label(labels[nz][nx + (width * py)]) != unvisited) {
		if (((float)g_Image3D[nz][nx + (width * py)]) <= S) {
			isNFillOrN = false;
		}
		//}
		//if (Label(labels[pz][nx + (width * py)]) != unvisited) {
		if (((float)g_Image3D[pz][nx + (width * py)]) <= S) {
			isNFillOrN = false;
		}
		//}
		//if (Label(labels[nz][px + (width * py)]) != unvisited) {
		if (((float)g_Image3D[nz][px + (width * py)]) <= S) {
			isNFillOrN = false;
		}
		//}
		//if (Label(labels[pz][px + (width * ny)]) != unvisited) {
		if (((float)g_Image3D[pz][px + (width * ny)]) <= S) {
			isNFillOrN = false;
		}
		//}
		//if (Label(labels[pz][nx + (width * ny)]) != unvisited) {
		if (((float)g_Image3D[pz][nx + (width * ny)]) <= S) {
			isNFillOrN = false;
		}
		//}
		return isNFillOrN;
	}
	cout << "case not addressed (bug) at coordinate " << px << " " << py << " " << pz << " with neighbor " << nx << " " << ny << " " << nz << endl;
	//3D diagonals are not core
	return false;
}

void addEdge(Graph& G, map<std::vector<int>, int>& edgeWt, int i1, int i2, int px, int py, int pz, int nx, int ny, int nz) {
	if (edgeWt.find({ i1,i2 }) != edgeWt.end()) { //Indicate edge is 6-connected
		int wt = edgeWt[{i1, i2}];
		if (wt == 0 && (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1)) {
			edgeWt[{i1, i2}] = 1; edgeWt[{i2, i1}] = 1;
		}
	}
	else { // Edge doesnt exist
		add_edge(i1, i2, G);
		if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
			edgeWt[{i1, i2}] = 1; edgeWt[{i2, i1}] = 1;
		}
		else {
			edgeWt[{i1, i2}] = 0; edgeWt[{i2, i1}] = 0;
		}

	}
}

void identifyCutFromPixel3D(std::vector< std::vector<uint32_t*> >& levelLabels, int lvl, std::vector<float*>& g_Image3D, Graph& G, int x, int y, int z, int width, int height,
	int numSlices, int labelCtI, node& n, float S, std::vector < std::vector<node> >& lvlNodes, map<std::vector<int>, int>& edgeWt, int geomCost) {
	uint32_t labelCt = (uint32_t)labelCtI;
	std::vector<node> nodes = lvlNodes[lvl];
	std::vector<uint32_t*> labels = levelLabels[lvl];
	std::queue<tuple<int, int, int>> q; q.push({ x,y,z });
	n.labelCost = 0; n.floatCost = 0.0;  n.intensity = 0.0;
	while (!q.empty()) {
		tuple<int, int, int> pt = q.front(); q.pop();
		int px = get<0>(pt); int py = get<1>(pt); int pz = get<2>(pt);
		if (Label(labels[pz][px + (width * py)]) == unvisited) {
			//1; //
			if (geomCost == 0) {
				n.floatCost += (gradient(x, y, z, g_Image3D, width, height, numSlices));
			}
			else {
				n.floatCost += 1;
			}
			n.greatestDiff = max(n.greatestDiff, abs((float)g_Image3D[pz][px + (width * py)] - S));
			changeLabel(labelCt, labels[pz][px + (width * py)]);
			setVisitedFlag(labels[pz][px + (width * py)], 1);

			//flood fill to nearby neighbor depending on priority.
			for (int i = 0; i < structCube.size(); i++) {
				int nx = px + structCube[i][0]; int ny = py + structCube[i][1]; int nz = pz + structCube[i][2];
				if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
					if ((float)g_Image3D[nz][nx + (width * ny)] > S) {
						if (Label(labels[nz][nx + (width * ny)]) == unvisited) { //unvisited
							if (ccNeighbor(nodes, px, py, pz, nx, ny, nz, labels, width, height, numSlices, g_Image3D, S)) {
								q.push(std::make_tuple(nx, ny, nz));
							}
						}
					}

				}
			}

			for (int i = 0; i < structCube.size(); i++) {
				int xn = px + structCube[i][0]; int yn = py + structCube[i][1]; int zn = pz + structCube[i][2];
				if (xn >= 0 && xn < width && yn >= 0 && yn < height && zn >= 0 && zn < numSlices) {
					if (Label(labels[zn][xn + (width * yn)]) != unvisited && Label(labels[pz][px + (width * py)]) != Label(labels[zn][xn + (width * yn)])) { //neighboring voxel has a label, is potential edge
						if (ccNeighbor(nodes, px, py, pz, xn, yn, zn, labels, width, height, numSlices, g_Image3D, S)) {
							addEdge(G, edgeWt, labelCt, Label(labels[zn][xn + (width * yn)]), px, py, pz, xn, yn, zn);
						}
					}
				}
			}
		}
	}
}

bool ccNeighborFill6Conn(const std::vector<node>& nodes, int px, int py, int pz, int nx, int ny, int nz, const std::vector<uint32_t*>& labels, int width) {
	//is a six neighbor
	bool c1 = abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1;
	if (c1) {
		return true;
	}
	if ((abs(nx - px) + abs(ny - py) + abs(nz - pz) == 2)) {
		if (pz == nz) {

			bool isNFillOrN1 = true;
			if (Label(labels[pz][nx + (width * py)]) != unvisited) {
				if (((int)nodes[Label(labels[pz][nx + (width * py)])].type) == 1) {
					isNFillOrN1 = false;
				}
			}

			bool isNFillOrN2 = true;
			if (Label(labels[pz][px + (width * ny)]) != unvisited) {
				if (((int)nodes[Label(labels[pz][px + (width * ny)])].type) == 1) {
					isNFillOrN2 = false;
				}
			}

			//2D diagonals are not fill or neighbor
			return isNFillOrN1 && isNFillOrN2;
		}

		if (px == nx) {
			bool isNFillOrN1 = true;
			if (Label(labels[nz][px + (width * py)]) != unvisited) {
				if (((int)nodes[Label(labels[nz][px + (width * py)])].type) == 1) {
					isNFillOrN1 = false;
				}
			}

			bool isNFillOrN2 = true;
			if (Label(labels[pz][px + (width * ny)]) != unvisited) {
				if (((int)nodes[Label(labels[pz][px + (width * ny)])].type) == 1) {
					isNFillOrN2 = false;
				}
			}
			//2D diagonals are not fill or neighbor
			return isNFillOrN1 && isNFillOrN2;
		}

		if (py == ny) {
			bool isNFillOrN1 = true;
			if (Label(labels[nz][px + (width * py)]) != unvisited) {
				if (((int)nodes[Label(labels[nz][px + (width * py)])].type) == 1) {
					isNFillOrN1 = false;
				}
			}

			bool isNFillOrN2 = true;
			if (Label(labels[pz][nx + (width * py)]) != unvisited) {
				if (((int)nodes[Label(labels[pz][nx + (width * py)])].type) == 1) {
					isNFillOrN2 = false;
				}
			}
			return isNFillOrN1 && isNFillOrN2;
		}
	}

	if ((abs(nx - px) + abs(ny - py) + abs(nz - pz) == 3)) {
		bool isNFillOrN = true;
		if (Label(labels[nz][px + (width * ny)]) != unvisited) {
			if (((int)nodes[Label(labels[nz][px + (width * ny)])].type) == 1) {
				isNFillOrN = false;
			}
		}
		if (Label(labels[nz][nx + (width * py)]) != unvisited) {
			if (((int)nodes[Label(labels[nz][nx + (width * py)])].type) == 1) {
				isNFillOrN = false;
			}
		}
		if (Label(labels[pz][nx + (width * py)]) != unvisited) {
			if (((int)nodes[Label(labels[pz][nx + (width * py)])].type) == 1) {
				isNFillOrN = false;
			}
		}
		if (Label(labels[nz][px + (width * py)]) != unvisited) {
			if (((int)nodes[Label(labels[nz][px + (width * py)])].type) == 1) {
				isNFillOrN = false;
			}
		}
		if (Label(labels[pz][px + (width * ny)]) != unvisited) {
			if (((int)nodes[Label(labels[pz][px + (width * ny)])].type) == 1) {
				isNFillOrN = false;
			}
		}
		if (Label(labels[pz][nx + (width * ny)]) != unvisited) {
			if (((int)nodes[Label(labels[pz][nx + (width * ny)])].type) == 1) {
				isNFillOrN = false;
			}
		}
		return isNFillOrN;
	}
	cout << "case not addressed for fill connectivity, with coordinate " << px << " " << py << " " << pz << " and neighbor " << nx << " " << ny << " " << nz << endl;
	return false;
}

void identifyFillFromPixel3D(std::vector< std::vector<uint32_t*> >& levelLabels, int lvl, std::vector<float*>& g_Image3D, Graph& G, int x, int y, int z, int width, int height,
	int numSlices, int labelCtI, node& n, float S, std::vector< std::vector<node> >& lvlNodes, map<std::vector<int>, int>& edgeWt, int geomCost) {
	uint32_t labelCt = (uint32_t)labelCtI;
	std::vector<uint32_t*> labels = levelLabels[lvl];
	std::queue<tuple<int, int, int>> q; q.push(std::make_tuple(x, y, z));
	n.labelCost = 0; n.floatCost = 0.0;
	std::vector<float> diffs;
	std::vector<node> nodes = lvlNodes[lvl];
	while (!q.empty()) {
		tuple<int, int, int> pt = q.front(); q.pop();
		int px = std::get<0>(pt); int py = std::get<1>(pt); int pz = std::get<2>(pt);
		if (Label(labels[pz][px + (width * py)]) == unvisited) {
			if (geomCost == 0) {
				n.floatCost += (-gradient(x, y, z, g_Image3D, width, height, numSlices));//-1;//
			}
			else {
				n.floatCost -= 1;
			}
			n.greatestDiff = max(n.greatestDiff, abs((float)g_Image3D[pz][px + (width * py)] - S));
			changeLabel(labelCt, labels[pz][px + (width * py)]);
			setVisitedFlag(labels[pz][px + (width * py)], 1);

			//flood fill to nearby neighbor depending on priority.
			for (int i = 0; i < structCube.size(); i++) {
				int xn = px + structCube[i][0]; int yn = py + structCube[i][1]; int zn = pz + structCube[i][2];
				if (xn >= 0 && xn < width && yn >= 0 && yn < height && zn >= 0 && zn < numSlices) {
					if ((float)g_Image3D[zn][xn + (width * yn)] <= S) {
						if (Label(labels[zn][xn + (width * yn)]) == unvisited) { //unvisited
							//, g_Image3D, S
							if (ccNeighborFill6Conn(nodes, px, py, pz, xn, yn, zn, labels, width)) {
								q.push(std::make_tuple(xn, yn, zn));
							}
						}
					}

				}
			}

			//Add edge depending on how well connected
			for (int i = 0; i < structCube.size(); i++) {
				int xn = px + structCube[i][0]; int yn = py + structCube[i][1]; int zn = pz + structCube[i][2];
				if (xn >= 0 && xn < width && yn >= 0 && yn < height && zn >= 0 && zn < numSlices) {
					if (Label(labels[zn][xn + (width * yn)]) != unvisited && Label(labels[pz][px + (width * py)]) != Label(labels[zn][xn + (width * yn)])) { //neighboring voxel has a label, is potential edge
						if (ccNeighborFill6Conn(nodes, px, py, pz, xn, yn, zn, labels, width)) {
							addEdge(G, edgeWt, labelCt, Label(labels[zn][xn + (width * yn)]), px, py, pz, xn, yn, zn);
						}
					}
				}
			}
		}
	}
}

std::vector<int> getEulerNumbers(std::vector<node>& nodes, std::vector<uint32_t*>& labels, int width, int height, int numSlices) {
	int v = 0; int e = 0; int f = 0; int c = 0;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				//add up vertices
				int l = Label(labels[k][i + (width * j)]);
				if ((int)nodes[l].type == 1) { continue; }
				switch ((int)nodes[l].type) {
				case 0:
					v += 1;
					break;
				case 2:
					v += 1;
					nodes[l].v += 1; //add to generator count
					break;
				case 3:
					//Fill, so dont add to fg v count
					nodes[l].v += 1; //add to generator count
					break;
				default:
					break;

				}

				//add up edges in x direction
				if (i + 1 < width) {
					//compare (i,j,k), (i+1, j, k)
					int l2 = Label(labels[k][(i + 1) + (width * j)]);
					//both vts have to not be in neighborhood
					if (((int)nodes[l].type) != 1 && ((int)nodes[l2].type) != 1) {
						//0 is core, 1 is cut, 2 is fill. Priority is fill because 6-conn, so should take max
						int maxNode = l;
						if (((int)nodes[l2].type) > ((int)nodes[l].type)) {
							maxNode = l2;
						}
						switch ((int)nodes[maxNode].type) {
						case 0:
							e += 1;
							break;
						case 2:
							e += 1;
							nodes[maxNode].e += 1; //add to generator count
							break;
						case 3:
							//Fill, so dont add to fg e count
							nodes[maxNode].e += 1; //add to generator count
							break;
						default:
							break;

						}

					}
				}

				//add up edges in y direction
				if (j + 1 < height) {
					//compare (i,j,k), (i, j+1, k)
					int l2 = Label(labels[k][i + (width * (j + 1))]);
					//both vts have to not be in neighborhood
					if (((int)nodes[l].type) != 1 && ((int)nodes[l2].type) != 1) {
						//0 is core, 1 is cut, 2 is fill. Priority is fill because 6-conn, so should take max
						int maxNode = l;
						if (((int)nodes[l2].type) > ((int)nodes[l].type)) {
							maxNode = l2;
						}
						switch ((int)nodes[maxNode].type) {
						case 0:
							e += 1;
							break;
						case 2:
							e += 1;
							nodes[maxNode].e += 1; //add to generator count
							break;
						case 3:
							//Fill, so dont add to fg e count
							nodes[maxNode].e += 1; //add to generator count
							break;
						default:
							break;

						}

					}
				}

				//add up edges in z direction
				if (k + 1 < numSlices) {
					//compare (i,j,k), (i, j+1, k)
					int l2 = Label(labels[k + 1][i + (width * j)]);
					//both vts have to not be in neighborhood
					if (((int)nodes[l].type) != 1 && ((int)nodes[l2].type) != 1) {
						//0 is core, 1 is cut, 2 is fill. Priority is fill because 6-conn, so should take max
						int maxNode = l;
						if (((int)nodes[l2].type) > ((int)nodes[l].type)) {
							maxNode = l2;
						}
						switch ((int)nodes[maxNode].type) {
						case 0:
							e += 1;
							break;
						case 2:
							e += 1;
							nodes[maxNode].e += 1; //add to generator count
							break;
						case 3:
							//Fill, so dont add to fg e count
							nodes[maxNode].e += 1; //add to generator count
							break;
						default:
							break;

						}

					}
				}

				//add up faces in xy plane
				if (i + 1 < width && j + 1 < height) {
					//compare (i,j,k), (i, j+1, k)
					int l2 = Label(labels[k][(i + 1) + (width * (j))]);
					int l3 = Label(labels[k][(i)+(width * (j + 1))]);
					int l4 = Label(labels[k][(i + 1) + (width * (j + 1))]);
					//both vts have to not be in neighborhood
					if (((int)nodes[l].type) != 1 && ((int)nodes[l2].type) != 1 && ((int)nodes[l3].type) != 1 && ((int)nodes[l4].type) != 1) {
						//0 is core, 1 is cut, 2 is fill. Priority is fill because 6-conn, so should take max
						int maxNode = l;
						if (((int)nodes[l2].type) > ((int)nodes[maxNode].type)) {
							maxNode = l2;
						}
						if (((int)nodes[l3].type) > ((int)nodes[maxNode].type)) {
							maxNode = l3;
						}
						if (((int)nodes[l4].type) > ((int)nodes[maxNode].type)) {
							maxNode = l4;
						}
						switch ((int)nodes[maxNode].type) {
						case 0:
							f += 1;
							break;
						case 2:
							f += 1;
							nodes[maxNode].f += 1; //add to generator count
							break;
						case 3:
							//Fill, so dont add to fg e count
							nodes[maxNode].f += 1; //add to generator count

							break;
						default:
							break;

						}

					}
				}

				//add up faces in xz plane
				if (i + 1 < width && k + 1 < numSlices) {
					//compare (i,j,k), (i, j+1, k)
					int l2 = Label(labels[k][(i + 1) + (width * (j))]);
					int l3 = Label(labels[k + 1][(i)+(width * (j))]);
					int l4 = Label(labels[k + 1][(i + 1) + (width * (j))]);
					//both vts have to not be in neighborhood
					if (((int)nodes[l].type) != 1 && ((int)nodes[l2].type) != 1 && ((int)nodes[l3].type) != 1 && ((int)nodes[l4].type) != 1) {						  //0 is core, 1 is cut, 2 is fill. Priority is fill because 6-conn, so should take max
						int maxNode = l;
						if (((int)nodes[l2].type) > ((int)nodes[maxNode].type)) {
							maxNode = l2;
						}
						if (((int)nodes[l3].type) > ((int)nodes[maxNode].type)) {
							maxNode = l3;
						}
						if (((int)nodes[l4].type) > ((int)nodes[maxNode].type)) {
							maxNode = l4;
						}
						switch ((int)nodes[maxNode].type) {
						case 0:
							f += 1;
							break;
						case 2:
							f += 1;
							nodes[maxNode].f += 1; //add to generator count
							break;
						case 3:

							//Fill, so dont add to fg e count
							nodes[maxNode].f += 1; //add to generator count
							break;
						default:
							break;

						}

					}
				}

				//add up faces in yz plane
				if (j + 1 < height && k + 1 < numSlices) {
					//compare (i,j,k), (i, j+1, k)
					int l2 = Label(labels[k][(i)+(width * (j + 1))]);
					int l3 = Label(labels[k + 1][(i)+(width * (j))]);
					int l4 = Label(labels[k + 1][(i)+(width * (j + 1))]);
					//both vts have to not be in neighborhood
					if (((int)nodes[l].type) != 1 && ((int)nodes[l2].type) != 1 && ((int)nodes[l3].type) != 1 && ((int)nodes[l4].type) != 1) {						  //0 is core, 1 is cut, 2 is fill. Priority is fill because 6-conn, so should take max
						int maxNode = l;
						if (((int)nodes[l2].type) > ((int)nodes[maxNode].type)) {
							maxNode = l2;
						}
						if (((int)nodes[l3].type) > ((int)nodes[maxNode].type)) {
							maxNode = l3;
						}
						if (((int)nodes[l4].type) > ((int)nodes[maxNode].type)) {
							maxNode = l4;
						}
						switch ((int)nodes[maxNode].type) {
						case 0:
							f += 1;
							break;
						case 2:
							f += 1;
							nodes[maxNode].f += 1; //add to generator count
							break;
						case 3:
							//Fill, so dont add to fg e count
							nodes[maxNode].f += 1; //add to generator count

							break;
						default:
							break;

						}

					}
				}

				//Add up cubes
				if (i + 1 < width && j + 1 < height && k + 1 < numSlices) {
					bool hasCube = true;
					int maxNode = l;
					for (int o = 0; o < cubeFrontMask.size(); o++) {
						int coord[3] = { i + cubeFrontMask[o][0], j + cubeFrontMask[o][1], k + cubeFrontMask[o][2] };
						if (((int)nodes[Label(labels[coord[2]][(coord[0]) + (width * (coord[1]))])].type) == 1) {

							hasCube = false;
							break;
						}
						if (((int)nodes[Label(labels[coord[2]][(coord[0]) + (width * (coord[1]))])].type) >
							((int)nodes[maxNode].type)
							) {
							maxNode = Label(labels[coord[2]][(coord[0]) + (width * (coord[1]))]);
						}
					}
					if (hasCube) {
						switch ((int)nodes[maxNode].type) {
						case 0:
							c += 1;
							break;
						case 2:
							c += 1;
							nodes[maxNode].c += 1; //add to generator count
							break;
						case 3:
							//Fill, so dont add to fg e count
							nodes[maxNode].c += 1; //add to generator count
							break;
						default:
							break;
						}
					}
				}

			}
		}
	}
	return { v,e,f,c };
}

int labelComponentsTBg(const std::vector<float*>& g_Image3D, float t, int conn, int width, int height, int numSlices) {
	int ct = 0;
	std::vector< std::vector<int > > mask;
	if (conn == 0) {
		mask = structCube;
	}
	else {
		mask = structCross3D;
	}

	std::vector< std::vector< std::vector<bool>>> visited(width, std::vector<std::vector<bool>>(height, std::vector<bool>(numSlices, false)));

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				//Unvisited foreground voxels
				if (g_Image3D[k][i + j * width] <= t && visited[i][j][k] == false) {
					ct += 1;
					std::queue<Coordinate> q;
					q.push(Coordinate(i, j, k));
					visited[i][j][k] = true;
					while (!q.empty()) {
						Coordinate qp = q.front();
						q.pop();
						for (int s = 0; s < mask.size(); s++) {
							Coordinate np(qp.x + mask[s][0], qp.y + mask[s][1], qp.z + mask[s][2]);
							if (np.x >= 0 && np.x < width && np.y >= 0 && np.y < height && np.z >= 0 && np.z < numSlices) {
								if (g_Image3D[np.z][np.x + (np.y * width)] <= t && visited[np.x][np.y][np.z] == false) {
									visited[np.x][np.y][np.z] = true;
									q.push(np);
								}
							}
						}
					}

				}

			}
		}
	}
	return ct;
}


int labelComponentsT(const std::vector<float*>& g_Image3D, float t, int conn, int width, int height, int numSlices) {
	int ct = 0;
	std::vector< std::vector<int > > mask;
	if (conn == 0) {
		mask = structCube;
	}
	else {
		mask = structCross3D;
	}

	std::vector< std::vector< std::vector<bool>>> visited(width, std::vector<std::vector<bool>>(height, std::vector<bool>(numSlices, false)));

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				//Unvisited foreground voxels

				if (((float)g_Image3D[k][i + j * width]) > t && visited[i][j][k] == false) {
					ct += 1;
					std::queue<Coordinate> q;
					q.push(Coordinate(i, j, k));
					visited[i][j][k] = true;
					while (!q.empty()) {
						Coordinate qp = q.front();
						q.pop();
						for (int s = 0; s < mask.size(); s++) {
							Coordinate np(qp.x + mask[s][0], qp.y + mask[s][1], qp.z + mask[s][2]);
							if (np.x >= 0 && np.x < width && np.y >= 0 && np.y < height && np.z >= 0 && np.z < numSlices) {
								if (((float)g_Image3D[np.z][np.x + (np.y * width)]) > t && visited[np.x][np.y][np.z] == false) {
									visited[np.x][np.y][np.z] = true;
									q.push(np);
								}
							}
						}
					}

				}

			}
		}
	}
	return ct;
}

std::vector<int> getEulerNumbersFromGImg(const std::vector<float*>& g_Image3D, float t, int width, int height, int numSlices, bool complement) {
	int v = 0; int e = 0; int f = 0; int c = 0;

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				if (((float)g_Image3D[k][i + j * width] > t && !complement) || ((float)g_Image3D[k][i + j * width] <= t && complement)) {
					v += 1;

					if (i + 1 < width) {
						//compare (i,j,k), (i+1, j, k)
						if ((((float)g_Image3D[k][(i + 1) + (width * j)]) > t && !complement) || (((float)g_Image3D[k][(i + 1) + (width * j)]) <= t && complement)) {
							e += 1;
						}
					}

					if (j + 1 < height) {
						//compare (i,j,k), (i, j+1, k)
						if ((((float)g_Image3D[k][i + (width * (j + 1))]) > t && !complement) || (((float)g_Image3D[k][i + (width * (j + 1))]) <= t && complement)) {
							e += 1;
						}
					}

					if (k + 1 < numSlices) {
						//compare (i,j,k), (i, j+1, k)
						if ((((float)g_Image3D[k + 1][i + (width * (j))]) > t && !complement) || (((float)g_Image3D[k + 1][i + (width * (j))]) <= t && complement)) {
							e += 1;
						}
					}

					if (i + 1 < width && j + 1 < height) {
						//compare (i,j,k), (i, j+1, k)
						if (
							(((float)g_Image3D[k][(i + 1) + (width * (j))]) > t &&
								((float)g_Image3D[k][(i)+(width * (j + 1))]) > t && ((float)g_Image3D[k][(i + 1) + (width * (j + 1))]) > t && !complement) ||
							(((float)g_Image3D[k][(i + 1) + (width * (j))]) <= t &&
								((float)g_Image3D[k][(i)+(width * (j + 1))]) <= t && ((float)g_Image3D[k][(i + 1) + (width * (j + 1))]) <= t && complement)
							) {
							f += 1;
						}
					}

					if (i + 1 < width && k + 1 < numSlices) {
						//compare (i,j,k), (i, j+1, k)
						if (
							(((float)g_Image3D[k][(i + 1) + (width * (j))]) > t &&
								((float)g_Image3D[k + 1][(i)+(width * (j))]) > t && ((float)g_Image3D[k + 1][(i + 1) + (width * (j))]) > t && !complement)
							||
							(((float)g_Image3D[k][(i + 1) + (width * (j))]) <= t &&
								((float)g_Image3D[k + 1][(i)+(width * (j))]) <= t && ((float)g_Image3D[k + 1][(i + 1) + (width * (j))]) <= t && complement)
							) {
							f += 1;

						}
					}

					//add up faces in yz plane
					if (j + 1 < height && k + 1 < numSlices) {
						//compare (i,j,k), (i, j+1, k)
						if (
							((float)g_Image3D[k][(i)+(width * (j + 1))] > t &&
								((float)g_Image3D[k + 1][(i)+(width * (j))]) > t && ((float)g_Image3D[k + 1][(i)+(width * (j + 1))]) > t && !complement)
							||
							((float)g_Image3D[k][(i)+(width * (j + 1))] <= t &&
								((float)g_Image3D[k + 1][(i)+(width * (j))]) <= t && ((float)g_Image3D[k + 1][(i)+(width * (j + 1))]) <= t && complement)

							) {
							f += 1;
						}
					}

					//Add up cubes
					if (i + 1 < width && j + 1 < height && k + 1 < numSlices) {
						bool hasCube = true;
						for (int o = 0; o < cubeFrontMask.size(); o++) {
							int coord[3] = { i + cubeFrontMask[o][0], j + cubeFrontMask[o][1], k + cubeFrontMask[o][2] };
							if (((g_Image3D[coord[2]][(coord[0]) + (width * (coord[1]))] <= t) && !complement) || ((g_Image3D[coord[2]][(coord[0]) + (width * (coord[1]))] > t) && complement)) {

								hasCube = false;
								break;
							}
						}
						if (hasCube) {
							c += 1;
						}

					}
				}
			}
		}

	}
	return { v,e,f,c };
}

std::vector<int> getTopoFromGImg(const std::vector<float*>& g_Image3D, float t, int fgconn, int width, int height, int numSlices) {
	int h0 = labelComponentsT(g_Image3D, t, fgconn, width, height, numSlices);
	int h2 = labelComponentsTBg(g_Image3D, t, 1 - fgconn, width, height, numSlices) - 1;
	std::vector<int> eulNums = getEulerNumbersFromGImg(g_Image3D, t, width, height, numSlices, false);
	int h1 = h0 + h2 - (eulNums[0] - eulNums[1] + eulNums[2] - eulNums[3]);
	return { h0,h2,h1 };
}

struct Level {
	std::vector<node> nodes;
	std::vector<node> origNodes;
	double cost;
	int cfCt;
	int64_t wtSum;
	std::vector<std::vector<int> > newToOldComp;
	int h0;
	int h1;
	int h2;
	double geomCost;
};

struct Contradiction {
	int l1;
	int l2;
	int n1;
	int n2;
	std::vector<int> lowerIndices; //overall indices across levels
	std::vector<int> upperIndices;
};

struct State {
	std::vector<Level> levels;
	double cost;
	std::vector<Contradiction> contradictions;
	vector<vector<vector<int> > > levelNewToOldComps;
	std::vector< std::vector<int> > levelEulerNums;
	int64_t totalWtSum;
};

class doubleNode {
public:
	doubleNode(node n, int specialType, int side);
	int getType(); int getSide();
	std::vector<int> hypernodes; int localIndex; node origNode;
private:
	int type;
	int side;
};

std::vector< std::vector<int> > findGraphComponents(grapht& g, std::vector<node>& nodes, bool inFg) {

	int numNodes = nodes.size();
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(g, &nodeToComp[0]);
	std::vector< std::vector<int> > wholeComps(n);
	std::vector<bool> isCompIndex(n, false);
	for (int i = 0; i < nodeToComp.size(); i++) {
		if (nodes[i].valid) {
			wholeComps[nodeToComp[i]].push_back(i);
			if (!isCompIndex[nodeToComp[i]]) {
				isCompIndex[nodeToComp[i]] = true;
			}
		}
	}

	std::vector< std::vector<int> > comps;
	for (int i = 0; i < wholeComps.size(); i++) {
		std::vector<int> comp = wholeComps[i];
		if (comp.size() == 0) {
			continue;
		}
		bool addComp = true;
		for (int j = 0; j < comp.size(); j++) {
			if (nodes[comp[j]].valid) {
				if (((int)nodes[comp[j]].type) == 1 && inFg) {
					addComp = false;
				}
				if (((int)nodes[comp[j]].type) == 0 && !inFg) {
					addComp = false;
				}
			}
			else {
				addComp = false;
			}
		}
		if (addComp) {
			comps.push_back(comp);
		}
	}
	return comps;
}

void pruneIsolatedNonTerminalGroups(grapht& fgG, grapht& bgG, std::vector<node>& nodes, std::vector<bool>& isAffectedComp, std::vector<int>& affectedComponents, map<int, int>& nodeToComp, bool trackComponents) {
	std::vector< std::vector<int> > fgComps = findGraphComponents(fgG, nodes, true); std::vector< std::vector<int> > bgComps = findGraphComponents(bgG, nodes, false);
	for (int i = 0; i < fgComps.size(); i++) {
		bool hasCore = false;
		std::vector<int> comp = fgComps[i];
		for (int j = 0; j < comp.size(); j++) {
			if (((int)nodes[comp[j]].type) == 0) {
				hasCore = true;
				break;
			}
		}
		if (!hasCore) {
			for (int j = 0; j < comp.size(); j++) {
				nodes[comp[j]].type = 1;
				nodes[comp[j]].inFg = 0;
				if (trackComponents) {
					isAffectedComp[nodeToComp[comp[j]]] = true; affectedComponents.push_back(nodeToComp[comp[j]]);
				}
			}
		}
	}
	for (int i = 0; i < bgComps.size(); i++) {
		bool hasN = false;
		std::vector<int> comp = bgComps[i];
		for (int j = 0; j < comp.size(); j++) {
			if (((int)nodes[comp[j]].type) == 1) {
				hasN = true;
				break;
			}
		}
		if (!hasN) {
			for (int j = 0; j < comp.size(); j++) {
				nodes[comp[j]].type = 0;
				nodes[comp[j]].inFg = 1;
				if (trackComponents) {
					isAffectedComp[nodeToComp[comp[j]]] = true; affectedComponents.push_back(nodeToComp[comp[j]]);
				}
			}
		}
	}
}

void updateEdges(Graph& G, grapht& g, std::vector<node>& nodes, std::map< std::vector<int>, int >& edgeWtMap, bool inFg) {
	g = grapht();
	for (int i = 0; i < nodes.size(); i++) {
		add_vertex(g);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (inFg) {
			if ((int)nodes[v1].type == 1 || (int)nodes[v2].type == 1) {
				continue;
			}
			if (edgeWtMap[{v1, v2}] == 1) {
				add_edge(v1, v2, 0, g);
			}
		}
		else {
			if ((int)nodes[v1].type == 0 || (int)nodes[v2].type == 0) {
				continue;
			}
			add_edge(v1, v2, 0, g);
		}
	}
}

std::vector<int> findCriticalArticulationNodes(grapht& g, std::vector<node>& nodes) {
	std::vector<vertex_t> art_points;
	boost::articulation_points(g, std::back_inserter(art_points));
	std::vector<int> isArticulation(nodes.size(), 0);
	for (int i = 0; i < art_points.size(); i++) {
		isArticulation[(int)art_points[i]] = 1;
	}

	std::vector<int> isArt(nodes.size(), 0);

	for (int i = 0; i < art_points.size(); i++) {
		int av = (int)art_points[i];
		std::vector<int> adjNodes;
		auto neighbours = adjacent_vertices(av, g);
		std::vector<bool> visited(nodes.size(), false);
		visited[av] = true;
		int termCt = 0;
		for (auto u : make_iterator_range(neighbours)) {
			//Each neighbour of articulation point forms potential new component
			if (!visited[u]) {
				bool hasTerm = false;
				queue<int> q;
				q.push(u); visited[u] = true;
				while (!q.empty()) {
					int p = q.front();
					q.pop();
					if (((int)nodes[p].type) == 0 || ((int)nodes[p].type) == 1) {
						hasTerm = true;
					}
					auto neighboursp = adjacent_vertices(p, g);
					for (auto up : make_iterator_range(neighboursp)) {
						if (!visited[up]) {
							q.push(up);
							visited[up] = true;
						}
					}
				}
				if (hasTerm) {
					termCt += 1;
				}
			}
		}
		if (termCt > 1) {
			isArt[av] = termCt;
		}
	}
	return isArt;
}

int assignArticulationNodes(grapht& fgA, grapht& bgA, std::vector<node>& nodes, std::vector<bool>& isAffectedComp, std::vector<int>& affectedComponents, map<int, int>& nodeToComp, bool trackComponents, Graph& G) {
	std::vector<int> isArtFg = findCriticalArticulationNodes(fgA, nodes); std::vector<int> isArtBg = findCriticalArticulationNodes(bgA, nodes);

	int artNodes = 0;
	int64_t maxIRange = 0;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].valid) {
			if ((int)nodes[i].type == 2 || (int)nodes[i].type == 3) {
				maxIRange += abs(nodes[i].labelCost);
			}
		}
	}

	maxIRange *= 2.0;
	maxIRange = max((int64_t)1, maxIRange);
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].valid) {
			if ((int)nodes[i].type == 2 || (int)nodes[i].type == 3) {
				if (isArtFg[i] > 1 || isArtBg[i] > 1) {
					artNodes += 1;
					if (isArtFg[i] > 1 && isArtBg[i] == 0) {
						if ((int)nodes[i].type == 3) {
							//If fill, assign all neighbouring cuts to be fg
							auto neighbours = adjacent_vertices(i, G);
							for (auto vd : make_iterator_range(neighbours)) {
								if ((int)nodes[vd].type == 2) {
									nodes[vd].type = 0;
									nodes[vd].inFg = 1;
									if (trackComponents) {
										isAffectedComp[nodeToComp[vd]] = true; affectedComponents.push_back(nodeToComp[vd]);
									}
								}
							}
						}

						//assign as fg terminal
						nodes[i].type = 0;
						nodes[i].inFg = 1;
						if (trackComponents) {
							isAffectedComp[nodeToComp[i]] = true; affectedComponents.push_back(nodeToComp[i]);
						}
					}
					if (isArtFg[i] == 0 && isArtBg[i] > 1) {
						if ((int)nodes[i].type == 2) {
							auto neighbours = adjacent_vertices(i, G);
							for (auto vd : make_iterator_range(neighbours)) {
								if ((int)nodes[vd].type == 3) {
									nodes[vd].type = 1;
									nodes[vd].inFg = 0;
									if (trackComponents) {
										isAffectedComp[nodeToComp[vd]] = true; affectedComponents.push_back(nodeToComp[vd]);
									}
								}
							}
						}
						nodes[i].type = 1;
						nodes[i].inFg = 0;
						//If cut, assign all neighbouring fills to be bg
						if (trackComponents) {
							isAffectedComp[nodeToComp[i]] = true; affectedComponents.push_back(nodeToComp[i]);
						}
					}
					if (isArtFg[i] > 1 && isArtBg[i] > 1) {
						int64_t lexCost = (maxIRange * (isArtFg[i] - isArtBg[i])) + nodes[i].labelCost;
						if (lexCost > 0) { //separates more fg than bg comps
							if ((int)nodes[i].type == 3) {
								//If fill, assign all neighbouring cuts to be fg
								auto neighbours = adjacent_vertices(i, G);
								for (auto vd : make_iterator_range(neighbours)) {
									if ((int)nodes[vd].type == 2) {
										nodes[vd].type = 0;
										nodes[vd].inFg = 1;
										if (trackComponents) {
											isAffectedComp[nodeToComp[vd]] = true; affectedComponents.push_back(nodeToComp[vd]);
										}
									}
								}
							}
							nodes[i].type = 0;
							nodes[i].inFg = 1;
							//If fill, assign all neighbouring cuts to be fg
							if (trackComponents) {
								isAffectedComp[nodeToComp[i]] = true; affectedComponents.push_back(nodeToComp[i]);
							}
						}
						else {
							if ((int)nodes[i].type == 2) {
								auto neighbours = adjacent_vertices(i, G);
								for (auto vd : make_iterator_range(neighbours)) {
									if ((int)nodes[vd].type == 3) {
										nodes[vd].type = 1;
										nodes[vd].inFg = 0;
										if (trackComponents) {
											isAffectedComp[nodeToComp[vd]] = true; affectedComponents.push_back(nodeToComp[vd]);
										}
									}
								}
							}
							nodes[i].type = 1;
							nodes[i].inFg = 0;
							//If cut, assign all neighbouring fills to be bg
							if (trackComponents) {
								isAffectedComp[nodeToComp[i]] = true; affectedComponents.push_back(nodeToComp[i]);
							}
						}
					}
				}
			}
		}
	}
	return artNodes;
}

void removeCAndNEdges(Graph& g, std::vector<node>& nodes) {
	Graph gN;
	for (int i = 0; i < nodes.size(); i++) {
		add_vertex(gN);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if ((int)nodes[v1].type == 0 && (int)nodes[v2].type == 1) {
			continue;
		}
		if ((int)nodes[v1].type == 1 && (int)nodes[v2].type == 0) {
			continue;
		}

		add_edge(v1, v2, gN);
	}
	g = gN;
}

void preprocessGraph(Graph& G, std::vector<node>& nodes, std::map< std::vector<int>, int>& edgeWtMap, std::vector<bool>& isAffectedComp,
	std::vector<int>& affectedComponents, map<int, int>& nodeToComp, bool trackComponents) {
	int numArt = 100000;
	grapht fgGA, bgGA;
	for (int i = 0; i < nodes.size(); i++) {
		add_vertex(fgGA); add_vertex(bgGA);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if ((((int)nodes[v1].type) == 0 || ((int)nodes[v1].type) == 2 || ((int)nodes[v1].type) == 3) &&
			(((int)nodes[v2].type) == 0 || ((int)nodes[v2].type) == 2 || ((int)nodes[v2].type) == 3) &&
			(nodes[v1].valid) && (nodes[v2].valid) && edgeWtMap[{v1, v2}] == 1) {
			//If src vtx or tgt vtx is core, cut, or fill, add to fg graph

			add_edge(v1, v2, 0, fgGA);
		}

		if ((((int)nodes[v1].type) == 1 || ((int)nodes[v1].type) == 2 || ((int)nodes[v1].type) == 3) &&
			(((int)nodes[v2].type) == 1 || ((int)nodes[v2].type) == 2 || ((int)nodes[v2].type) == 3) &&
			(nodes[v1].valid) && (nodes[v2].valid)) {
			//If src vtx or tgt vtx is core, cut, or fill, add to fg graph 
			add_edge(v1, v2, 0, bgGA);
		}
	}

	while (numArt > 0) {
		pruneIsolatedNonTerminalGroups(fgGA, bgGA, nodes, isAffectedComp, affectedComponents, nodeToComp, trackComponents);

		updateEdges(G, fgGA, nodes, edgeWtMap, true); updateEdges(G, bgGA, nodes, edgeWtMap, false);

		numArt = assignArticulationNodes(fgGA, bgGA, nodes, isAffectedComp, affectedComponents, nodeToComp, trackComponents, G);

		updateEdges(G, fgGA, nodes, edgeWtMap, true); updateEdges(G, bgGA, nodes, edgeWtMap, false);
	}
	removeCAndNEdges(G, nodes);
}

class componentGraph {
public:
	componentGraph() {
		numNodes = 0;
		nodeAssignment = 0;
	};
	componentGraph(int graphSize) {
		numNodes = graphSize;
		for (int i = 0; i < graphSize; i++) {
			add_vertex(g);
		}
		nodeAssignment = 0;
	};
	void clearLocalVertex(int localIndex);
	std::vector<std::vector<int> > getGlobalComponentsExcluding(std::vector<int>& excludeIndices);
	Graph g;
	int numNodes; int nodeAssignment;
	std::vector<int> globalNodes;
};


void componentGraph::clearLocalVertex(int localIndex) {
	clear_vertex(localIndex, g);
}

std::vector<std::vector<int> > componentGraph::getGlobalComponentsExcluding(std::vector<int>& excludeIndices) {
	std::vector< std::vector<int> > components;
	std::vector<int> nodeToComp(numNodes);

	if (numNodes == 1) {
		return { {globalNodes[0]} };
	}
	std::vector<bool> isExcluded(numNodes, false);
	for (int i = 0; i < excludeIndices.size(); i++) {
		isExcluded[excludeIndices[i]] = true;
	}
	if (num_edges(g) == 0) {
		for (int i = 0; i < numNodes; i++) {
			if (!isExcluded[i]) {
				components.push_back({ globalNodes[i] });
			}
		}
		return components;
	}

	int n = (int)boost::connected_components(g, &nodeToComp[0]);

	int numComps = 0; std::vector<int> isCompIndex(n, -1);

	for (int i = 0; i < numNodes; i++) {
		if (!isExcluded[i]) {
			if (isCompIndex[nodeToComp[i]] == -1) {

				isCompIndex[nodeToComp[i]] = numComps;
				std::vector<int> newComp = { globalNodes[i] };
				components.push_back(newComp);
				numComps += 1;
			}
			else {
				components[isCompIndex[nodeToComp[i]]].push_back(globalNodes[i]);
			}
		}
	}
	return components;
}

class doubleGraph {
public:
	doubleGraph(std::vector<node>& nodes, std::vector<int>& clusterLocalNodesGlobalIndices, std::vector< std::vector<int> >& clusterLocalEdgesGlobalIndices, std::vector< std::vector<int> >& complexEdges,
		map<int, int>& globalToLocal, map< std::vector<int>, int>& edgeWts);
	doubleGraph(std::vector<node>& nodes, int numNodes, Graph& g, map< std::vector<int>, int>& edgeWts);
	doubleGraph() {};
	std::vector<std::vector<int>> findTerminalComponents(bool inFg); void findNTComponentsConstraintEdges(bool inFg);
	Graph getBoostGraph();
	std::vector< std::vector<int> > getFgCompsComplex(); std::vector< std::vector<int> > getBgCompsComplex(); void assignEdgesToFgAndBgComps(); int getSize();
	std::vector<int> cores; std::vector<int> ns; std::vector<int> midTs;
	std::vector< std::vector<int> > findRemainingCFComponents(bool inFg, std::vector<int>& excludeIndices);
	int getPiIndex(); std::vector<doubleNode> doubleNodes; Graph complexG; Graph fgG; Graph bgG;
	map<int, int> doubleNodeFgCompMapping; map<int, int> doubleNodeBgCompMapping; std::vector<componentGraph> fgCompGraphs; std::vector<componentGraph> bgCompGraphs;
	std::vector< std::vector<int> > doubleEdges;
	map<int, std::vector<int> > nodeToDoubleMapping;
private:

	std::vector< std::vector<int> > doubleEdgesComplex; //Graph doubleG;
	int numDoubleNodes; int piIndex; std::vector< std::vector<int> > fgCompsComplex; std::vector< std::vector<int> > bgCompsComplex;

};

std::vector< std::vector<int> > doubleGraph::getFgCompsComplex() {
	return fgCompsComplex;
}

std::vector< std::vector<int> > doubleGraph::getBgCompsComplex() {
	return bgCompsComplex;
}

int doubleGraph::getPiIndex() {
	return piIndex;
}

int doubleGraph::getSize() {
	return numDoubleNodes;
}

std::vector<std::vector<int>>  doubleGraph::findTerminalComponents(bool inFg) {
	if (inFg) {
		Graph fgG;
		for (int i = 0; i < doubleNodes.size(); i++) {
			add_vertex(fgG);
		}
		for (int i = 0; i < doubleEdges.size(); i++) {
			if (doubleNodes[doubleEdges[i][0]].getSide() == 1 && doubleNodes[doubleEdges[i][1]].getSide() == 1) {
				add_edge(doubleEdges[i][0], doubleEdges[i][1], fgG);
			}
		}
		std::vector< std::vector<int> > components;
		std::vector<int> nodeToComp(numDoubleNodes);
		int n = (int)boost::connected_components(fgG, &nodeToComp[0]);
		int numComps = 0; std::vector<int> isCompIndex(n, -1);
		for (int i = 0; i < numDoubleNodes; i++) {
			if (((int)doubleNodes[i].getSide()) == 1) {
				if (isCompIndex[nodeToComp[i]] == -1) {

					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i);
				}

			}
		}

		std::vector<std::vector<int>> terminalComponents;
		for (int i = 0; i < components.size(); i++) {
			for (int j = 0; j < components[i].size(); j++) {
				if (((int)doubleNodes[components[i][j]].getType()) == CORE) {
					terminalComponents.push_back(components[i]);
					break;
				}
			}
		}
		int coreCt = 0;
		for (int i = 0; i < doubleNodes.size(); i++) {
			if ((int)doubleNodes[i].getType() == CORE) {
				coreCt++;
			}
		}
		return terminalComponents;
	}
	else {
		Graph bgG;
		for (int i = 0; i < doubleNodes.size(); i++) {
			add_vertex(bgG);
		}
		for (int i = 0; i < doubleEdges.size(); i++) {
			if (doubleNodes[doubleEdges[i][0]].getSide() == -1 && doubleNodes[doubleEdges[i][1]].getSide() == -1) {
				add_edge(doubleEdges[i][0], doubleEdges[i][1], bgG);
			}
		}
		std::vector< std::vector<int> > components;
		std::vector<int> nodeToComp(numDoubleNodes);
		int n = (int)boost::connected_components(bgG, &nodeToComp[0]);
		int numComps = 0; std::vector<int> isCompIndex(n, -1);
		for (int i = 0; i < numDoubleNodes; i++) {
			if (((int)doubleNodes[i].getSide()) == -1) {
				if (isCompIndex[nodeToComp[i]] == -1) {
					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i);
				}
			}
		}
		std::vector<std::vector<int>> terminalComponents;
		for (int i = 0; i < components.size(); i++) {
			for (int j = 0; j < components[i].size(); j++) {
				if (((int)doubleNodes[components[i][j]].getType()) == N) {
					terminalComponents.push_back(components[i]);
					break;
				}
			}
		}
		return terminalComponents;
	}
}

void doubleGraph::findNTComponentsConstraintEdges(bool inFg) {
	if (inFg) {
		for (int i = 0; i < doubleNodes.size(); i++) {
			add_vertex(fgG);
		}
		for (int i = 0; i < doubleEdgesComplex.size(); i++) {

			if ((doubleNodes[doubleEdgesComplex[i][0]].getSide() == 1 && doubleNodes[doubleEdgesComplex[i][0]].getType() != CORE) &&
				(doubleNodes[doubleEdgesComplex[i][1]].getSide() == 1 && doubleNodes[doubleEdgesComplex[i][1]].getType() != CORE)) {
				add_edge(doubleEdgesComplex[i][0], doubleEdgesComplex[i][1], fgG);
			}
		}
		std::vector< std::vector<int> > components;
		std::vector<int> nodeToComp(numDoubleNodes);
		int n = (int)boost::connected_components(fgG, &nodeToComp[0]);
		int numComps = 0; std::vector<int> isCompIndex(n, -1);
		for (int i = 0; i < numDoubleNodes; i++) {
			if (((int)doubleNodes[i].getSide()) == 1 && doubleNodes[i].getType() != CORE) {
				if (isCompIndex[nodeToComp[i]] == -1) {

					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1; doubleNodeFgCompMapping[i] = numComps - 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i); doubleNodeFgCompMapping[i] = isCompIndex[nodeToComp[i]];
				}
			}
		}
		for (int i = 0; i < components.size(); i++) {
			for (int j = 0; j < components[i].size(); j++) {
				if (((int)doubleNodes[components[i][j]].getSide()) == 1 && doubleNodes[components[i][j]].getType() != CORE) {
					fgCompsComplex.push_back(components[i]); fgCompGraphs.push_back(componentGraph(components[i].size()));
					break;
				}
			}
		}
	}
	else {
		for (int i = 0; i < doubleNodes.size(); i++) {
			add_vertex(bgG);
		}
		for (int i = 0; i < doubleEdgesComplex.size(); i++) {
			if ((doubleNodes[doubleEdgesComplex[i][0]].getSide() == -1 && doubleNodes[doubleEdgesComplex[i][0]].getType() != N) &&
				(doubleNodes[doubleEdgesComplex[i][1]].getSide() == -1 && doubleNodes[doubleEdgesComplex[i][1]].getType() != N)) {
				add_edge(doubleEdgesComplex[i][0], doubleEdgesComplex[i][1], bgG);
			}
		}
		std::vector< std::vector<int> > components;
		std::vector<int> nodeToComp(numDoubleNodes);
		int n = (int)boost::connected_components(bgG, &nodeToComp[0]);
		int numComps = 0; std::vector<int> isCompIndex(n, -1);
		for (int i = 0; i < numDoubleNodes; i++) {
			if (((int)doubleNodes[i].getSide()) == -1 && doubleNodes[i].getType() != N) {
				if (isCompIndex[nodeToComp[i]] == -1) {
					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1; doubleNodeBgCompMapping[i] = numComps - 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i); doubleNodeBgCompMapping[i] = isCompIndex[nodeToComp[i]];
				}
			}
		}
		for (int i = 0; i < components.size(); i++) {
			for (int j = 0; j < components[i].size(); j++) {
				if (((int)doubleNodes[components[i][j]].getSide()) == -1 && doubleNodes[components[i][j]].getType() != N) {
					bgCompsComplex.push_back(components[i]); bgCompGraphs.push_back(componentGraph(components[i].size()));
					break;
				}
			}
		}
	}
}

std::vector< std::vector<int> > doubleGraph::findRemainingCFComponents(bool inFg, std::vector<int>& excludeIndices) {
	if (inFg) {
		Graph fgRemaining = fgG; std::vector<bool> exclude(numDoubleNodes, false);
		for (int i = 0; i < excludeIndices.size(); i++) {
			clear_vertex(excludeIndices[i], fgRemaining); exclude[excludeIndices[i]] = true;
		}

		std::vector< std::vector<int> > components;
		std::vector<int> nodeToComp(numDoubleNodes); std::vector< std::vector<int> > fgCompsComplexR;
		int n = (int)boost::connected_components(fgRemaining, &nodeToComp[0]);
		int numComps = 0; std::vector<int> isCompIndex(n, -1);
		for (int i = 0; i < numDoubleNodes; i++) {
			if (((int)doubleNodes[i].getSide()) == 1 && doubleNodes[i].getType() != CORE && !exclude[i]) {
				if (isCompIndex[nodeToComp[i]] == -1) {

					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i);
				}
			}
		}
		for (int i = 0; i < components.size(); i++) {
			for (int j = 0; j < components[i].size(); j++) {
				if (((int)doubleNodes[components[i][j]].getSide()) == 1 && doubleNodes[components[i][j]].getType() != CORE && !exclude[components[i][j]]) {
					fgCompsComplexR.push_back(components[i]);
					break;
				}
			}
		}
		return fgCompsComplexR;
	}
	else {
		Graph bgGRemaining = bgG; std::vector<bool> exclude(numDoubleNodes, false);
		for (int i = 0; i < excludeIndices.size(); i++) {
			clear_vertex(excludeIndices[i], bgGRemaining); exclude[excludeIndices[i]] = true;
		}
		std::vector< std::vector<int> > components;
		std::vector<int> nodeToComp(numDoubleNodes); std::vector< std::vector<int> > bgCompsComplexR;
		int n = (int)boost::connected_components(bgGRemaining, &nodeToComp[0]);
		int numComps = 0; std::vector<int> isCompIndex(n, -1);
		for (int i = 0; i < numDoubleNodes; i++) {
			if (((int)doubleNodes[i].getSide()) == -1 && doubleNodes[i].getType() != N && !exclude[i]) {
				if (isCompIndex[nodeToComp[i]] == -1) {
					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i);
				}
			}
		}
		for (int i = 0; i < components.size(); i++) {
			for (int j = 0; j < components[i].size(); j++) {
				if (((int)doubleNodes[components[i][j]].getSide()) == -1 && doubleNodes[components[i][j]].getType() != N && !exclude[components[i][j]]) {
					bgCompsComplexR.push_back(components[i]);
					break;
				}
			}
		}

		return bgCompsComplexR;
	}

}

void doubleGraph::assignEdgesToFgAndBgComps() {
	//Assign local node index to double node
	//Within local graph object, create vector of global ints
	for (int i = 0; i < doubleNodes.size(); i++) {
		if (doubleNodes[i].getSide() == 1 && doubleNodes[i].getType() != CORE) {
			doubleNodes[i].localIndex = fgCompGraphs[doubleNodeFgCompMapping[i]].nodeAssignment;
			fgCompGraphs[doubleNodeFgCompMapping[i]].nodeAssignment = fgCompGraphs[doubleNodeFgCompMapping[i]].nodeAssignment + 1;
			fgCompGraphs[doubleNodeFgCompMapping[i]].globalNodes.push_back(i);
		}
		else {
			if (doubleNodes[i].getSide() == -1 && doubleNodes[i].getType() != N) {
				doubleNodes[i].localIndex = bgCompGraphs[doubleNodeBgCompMapping[i]].nodeAssignment;
				bgCompGraphs[doubleNodeBgCompMapping[i]].nodeAssignment = bgCompGraphs[doubleNodeBgCompMapping[i]].nodeAssignment + 1;
				bgCompGraphs[doubleNodeBgCompMapping[i]].globalNodes.push_back(i);
			}
		}
	}

	for (int i = 0; i < doubleEdgesComplex.size(); i++) {
		if (doubleNodes[doubleEdgesComplex[i][0]].getSide() == 1 && doubleNodes[doubleEdgesComplex[i][1]].getSide() == 1) {
			if (doubleNodes[doubleEdgesComplex[i][0]].getType() != CORE && doubleNodes[doubleEdgesComplex[i][1]].getType() != CORE) { //Must be cut or fill
				add_edge(doubleNodes[doubleEdgesComplex[i][0]].localIndex, doubleNodes[doubleEdgesComplex[i][1]].localIndex, fgCompGraphs[doubleNodeFgCompMapping[i]].g);
			}
		}
		else {
			if (doubleNodes[doubleEdgesComplex[i][0]].getSide() == -1 && doubleNodes[doubleEdgesComplex[i][1]].getSide() == -1) {
				if (doubleNodes[doubleEdgesComplex[i][0]].getType() != N && doubleNodes[doubleEdgesComplex[i][1]].getType() != N) { //Must be cut or fill
					add_edge(doubleNodes[doubleEdgesComplex[i][0]].localIndex, doubleNodes[doubleEdgesComplex[i][1]].localIndex, bgCompGraphs[doubleNodeBgCompMapping[i]].g);
				}
			}
		}
	}
}

doubleGraph::doubleGraph(std::vector<node>& nodes, std::vector<int>& clusterLocalNodesGlobalIndices, std::vector< std::vector<int> >& clusterLocalEdgesGlobalIndices,
	std::vector< std::vector<int> >& complexEdges, map<int, int>& globalToLocal, map< std::vector<int>, int>& edgeWts) {
	//Add nodes to double graph
	//Construct boost graph object to find components: could be useful later on as well
	for (int i = 0; i < clusterLocalNodesGlobalIndices.size(); i++) {
		if (((int)nodes[clusterLocalNodesGlobalIndices[i]].type) == CORE) {
			doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], CORE, 1));
			add_vertex(complexG); cores.push_back(((int)doubleNodes.size()) - 1);
			nodeToDoubleMapping[globalToLocal[clusterLocalNodesGlobalIndices[i]]] = { ((int)doubleNodes.size()) - 1 };
		}
		else {
			if (((int)nodes[clusterLocalNodesGlobalIndices[i]].type) == N) {
				doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], N, -1));
				add_vertex(complexG); ns.push_back(((int)doubleNodes.size()) - 1);
				nodeToDoubleMapping[globalToLocal[clusterLocalNodesGlobalIndices[i]]] = { ((int)doubleNodes.size()) - 1 };
			}
			else {
				if (((int)nodes[clusterLocalNodesGlobalIndices[i]].type) == CUT) {
					doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], CUT, 1));
					add_vertex(complexG);
					doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], MIDT, 0));
					add_vertex(complexG); midTs.push_back(doubleNodes.size() - 1);
					doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], CUT, -1));
					add_vertex(complexG);
					nodeToDoubleMapping[globalToLocal[clusterLocalNodesGlobalIndices[i]]] = { ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 };
					//Add edges between node on fg side, pi node, and node on bg side
					doubleEdges.push_back({ ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2 });
					doubleEdges.push_back({ ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 });
					//Do same for complex edges
					doubleEdgesComplex.push_back({ ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2 }); doubleEdgesComplex.push_back({ ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 });
					add_edge(((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2, complexG); add_edge(((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1, complexG);
				}
				else {
					if (((int)nodes[clusterLocalNodesGlobalIndices[i]].type) == FILL) {
						doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], FILL, 1));
						add_vertex(complexG);
						doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], MIDT, 0));
						add_vertex(complexG); midTs.push_back(doubleNodes.size() - 1);
						doubleNodes.push_back(doubleNode(nodes[clusterLocalNodesGlobalIndices[i]], FILL, -1));
						add_vertex(complexG);
						nodeToDoubleMapping[globalToLocal[clusterLocalNodesGlobalIndices[i]]] = { ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 };
						//Add edges between node on fg side, pi node, and node on bg side
						doubleEdges.push_back({ ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2 });
						doubleEdges.push_back({ ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 });
						//Do same for complex edges
						doubleEdgesComplex.push_back({ ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2 });
						doubleEdgesComplex.push_back({ ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 });
						add_edge(((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2, complexG);
						add_edge(((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1, complexG);
					}
				}
			}
		}
	}
	node n; n.type = PINODE; doubleNodes.push_back(doubleNode(n, PINODE, 0)); piIndex = doubleNodes.size() - 1; numDoubleNodes = doubleNodes.size(); add_vertex(complexG); //add_vertex(doubleG); 

	//Construct double graph edges for all edges from original graph
	for (int i = 0; i < clusterLocalEdgesGlobalIndices.size(); i++) {
		int v1 = clusterLocalEdgesGlobalIndices[i][0]; int v2 = clusterLocalEdgesGlobalIndices[i][1];
		int local1 = globalToLocal[v1]; int local2 = globalToLocal[v2];

		if (nodeToDoubleMapping[local1].size() == 1 && nodeToDoubleMapping[local2].size() > 1) {

			if (((int)nodes[v1].type) == 0) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
				}
			}
			else {
				if (((int)nodes[v1].type) == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][2] });
				}
			}
		}
		if (nodeToDoubleMapping[local2].size() == 1 && nodeToDoubleMapping[local1].size() > 1) {

			if (((int)nodes[v2].type) == CORE) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[local2][0], nodeToDoubleMapping[local1][0] });
				}
			}
			else {
				if (((int)nodes[v2].type) == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[local2][0], nodeToDoubleMapping[local1][2] });
				}
			}
		}
		if (nodeToDoubleMapping[local1].size() > 1 && nodeToDoubleMapping[local2].size() > 1) {

			if (edgeWts[{v1, v2}] == 1) {
				doubleEdges.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
			}
			doubleEdges.push_back({ nodeToDoubleMapping[local1][2], nodeToDoubleMapping[local2][2] });
		}
		if (nodeToDoubleMapping[local1].size() == 1 && nodeToDoubleMapping[local2].size() == 1) {
			if (((int)nodes[v1].type) == CORE && ((int)nodes[v2].type) == CORE) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
				}
			}
			else {
				if (((int)nodes[v1].type) == N && ((int)nodes[v2].type) == N) {
					doubleEdges.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
				}
			}
		}
	}
	//Construct complex double graph edges for all edges from original graph
	for (int i = 0; i < complexEdges.size(); i++) {

		int v1 = complexEdges[i][0]; int v2 = complexEdges[i][1];
		int local1 = globalToLocal[v1]; int local2 = globalToLocal[v2];
		if (nodeToDoubleMapping[local1].size() == 1 && nodeToDoubleMapping[local2].size() > 1) {
			if (((int)nodes[v1].type) == 0) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdgesComplex.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
					add_edge(nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0], complexG);
				}
			}
			else {
				if (((int)nodes[v1].type) == 1) {
					doubleEdgesComplex.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][2] });
					add_edge(nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][2], complexG);
				}
			}
		}
		if (nodeToDoubleMapping[local2].size() == 1 && nodeToDoubleMapping[local1].size() > 1) {
			if (((int)nodes[v2].type) == 0) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdgesComplex.push_back({ nodeToDoubleMapping[local2][0], nodeToDoubleMapping[local1][0] });
					add_edge(nodeToDoubleMapping[local2][0], nodeToDoubleMapping[local1][0], complexG);
				}
			}
			else {
				if (((int)nodes[v2].type) == 1) {
					doubleEdgesComplex.push_back({ nodeToDoubleMapping[local2][0], nodeToDoubleMapping[local1][2] });
					add_edge(nodeToDoubleMapping[local2][0], nodeToDoubleMapping[local1][2], complexG);
				}
			}
		}
		if (nodeToDoubleMapping[local1].size() > 1 && nodeToDoubleMapping[local2].size() > 1) {
			doubleEdgesComplex.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
			add_edge(nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0], complexG);
			doubleEdgesComplex.push_back({ nodeToDoubleMapping[local1][2], nodeToDoubleMapping[local2][2] });
			add_edge(nodeToDoubleMapping[local1][2], nodeToDoubleMapping[local2][2], complexG);
		}

		if (nodeToDoubleMapping[local1].size() == 1 && nodeToDoubleMapping[local2].size() == 1) {
			if (((int)nodes[v1].type) == CORE && ((int)nodes[v2].type) == CORE) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdgesComplex.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
					add_edge(nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0], complexG);
				}
			}
			else {
				if (((int)nodes[v1].type) == N && ((int)nodes[v2].type) == N) {
					doubleEdgesComplex.push_back({ nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0] });
					add_edge(nodeToDoubleMapping[local1][0], nodeToDoubleMapping[local2][0], complexG);
				}
			}
		}

	}
	int coreCompEs = 0;
	//Find terminal components
	std::vector<std::vector<int> > coreComps = findTerminalComponents(true);
	for (int i = 0; i < coreComps.size(); i++) {
		for (int j = 0; j < coreComps[i].size(); j++) {
			if (doubleNodes[coreComps[i][j]].getType() == CORE) {
				doubleEdges.push_back({ coreComps[i][j], piIndex }); //add_edge(coreComps[i][j], piIndex, doubleG);
				coreCompEs += 1;
				break;
			}
		}
	}
	int nCompEs = 0;
	std::vector<std::vector<int> > nComps = findTerminalComponents(false);
	for (int i = 0; i < nComps.size(); i++) {
		for (int j = 0; j < nComps[i].size(); j++) {
			if (doubleNodes[nComps[i][j]].getType() == N) {
				doubleEdges.push_back({ nComps[i][j], piIndex }); //add_edge(nComps[i][j], piIndex, doubleG);
				nCompEs += 1;
				break;
			}
		}
	}
	findNTComponentsConstraintEdges(true); findNTComponentsConstraintEdges(false);

	assignEdgesToFgAndBgComps(); //Will make it faster to build up hypernodes of connected components excluding particular nodes
}

doubleGraph::doubleGraph(std::vector<node>& nodes, int numNodes, Graph& g, map< std::vector<int>, int>& edgeWts) {
	//Add nodes to double graph
	//Construct boost graph object to find components: could be useful later on as well
	int numCores = 0;
	for (int i = 0; i < numNodes; i++) {
		if (((int)nodes[i].type) == CORE) {
			doubleNodes.push_back(doubleNode(nodes[i], CORE, 1));
			cores.push_back(((int)doubleNodes.size()) - 1);
			nodeToDoubleMapping[i] = { ((int)doubleNodes.size()) - 1 };
			numCores += 1;
		}
		else {
			if (((int)nodes[i].type) == N) {
				doubleNodes.push_back(doubleNode(nodes[i], N, -1));
				ns.push_back(((int)doubleNodes.size()) - 1);
				nodeToDoubleMapping[i] = { ((int)doubleNodes.size()) - 1 };
			}
			else {
				if (((int)nodes[i].type) == CUT) {
					doubleNodes.push_back(doubleNode(nodes[i], CUT, 1));
					doubleNodes.push_back(doubleNode(nodes[i], MIDT, 0));
					midTs.push_back(doubleNodes.size() - 1);
					doubleNodes.push_back(doubleNode(nodes[i], CUT, -1)); //add_vertex(doubleG); 
					nodeToDoubleMapping[i] = { ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 };
					//Add edges between node on fg side, pi node, and node on bg side
					doubleEdges.push_back({ ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2 });
					doubleEdges.push_back({ ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 });
				}
				else {
					if (((int)nodes[i].type) == FILL) {
						doubleNodes.push_back(doubleNode(nodes[i], FILL, 1));
						doubleNodes.push_back(doubleNode(nodes[i], MIDT, 0)); midTs.push_back(doubleNodes.size() - 1);
						doubleNodes.push_back(doubleNode(nodes[i], FILL, -1));
						nodeToDoubleMapping[i] = { ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 };
						//Add edges between node on fg side, pi node, and node on bg side
						doubleEdges.push_back({ ((int)doubleNodes.size()) - 3, ((int)doubleNodes.size()) - 2 });
						doubleEdges.push_back({ ((int)doubleNodes.size()) - 2, ((int)doubleNodes.size()) - 1 });
					}
				}
			}
		}
	}
	node n; n.type = PINODE; doubleNodes.push_back(doubleNode(n, PINODE, 0)); piIndex = doubleNodes.size() - 1; numDoubleNodes = doubleNodes.size();


	//Construct double graph edges for all edges from original graph
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		//cout << "edge " << ((int)nodes[v1].type) << " " << ((int)nodes[v2].type) << endl;
		if (nodeToDoubleMapping[v1].size() == 1 && nodeToDoubleMapping[v2].size() > 1) {

			if (((int)nodes[v1].type) == CORE) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[v1][0], nodeToDoubleMapping[v2][0] });
				}
			}
			else {
				if (((int)nodes[v1].type) == N) {
					doubleEdges.push_back({ nodeToDoubleMapping[v1][0], nodeToDoubleMapping[v2][2] });
				}
			}
		}
		if (nodeToDoubleMapping[v2].size() == 1 && nodeToDoubleMapping[v1].size() > 1) {

			if (((int)nodes[v2].type) == CORE) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[v2][0], nodeToDoubleMapping[v1][0] });
				}
			}
			else {
				if (((int)nodes[v2].type) == N) {
					doubleEdges.push_back({ nodeToDoubleMapping[v2][0], nodeToDoubleMapping[v1][2] });
				}
			}
		}
		if (nodeToDoubleMapping[v1].size() > 1 && nodeToDoubleMapping[v2].size() > 1) {
			if (edgeWts[{v1, v2}] == 1) {
				doubleEdges.push_back({ nodeToDoubleMapping[v1][0], nodeToDoubleMapping[v2][0] });
			}
			doubleEdges.push_back({ nodeToDoubleMapping[v1][2], nodeToDoubleMapping[v2][2] });
		}
		if (nodeToDoubleMapping[v1].size() == 1 && nodeToDoubleMapping[v2].size() == 1) {
			if (((int)nodes[v1].type) == CORE && ((int)nodes[v2].type) == CORE) {
				if (edgeWts[{v1, v2}] == 1) {
					doubleEdges.push_back({ nodeToDoubleMapping[v1][0], nodeToDoubleMapping[v2][0] });
				}
			}
			else {
				if (((int)nodes[v1].type) == N && ((int)nodes[v2].type) == N) {
					doubleEdges.push_back({ nodeToDoubleMapping[v1][0], nodeToDoubleMapping[v2][0] });
				}
			}
		}
	}

	int coreCompEs = 0;
	//Find terminal components
	std::vector<std::vector<int> > coreComps = findTerminalComponents(true);
	for (int i = 0; i < coreComps.size(); i++) {
		for (int j = 0; j < coreComps[i].size(); j++) {
			if (doubleNodes[coreComps[i][j]].getType() == CORE) {
				doubleEdges.push_back({ coreComps[i][j], piIndex });
				coreCompEs += 1;
				break;
			}
		}

	}
	int nCompEs = 0;
	std::vector<std::vector<int> > nComps = findTerminalComponents(false);
	for (int i = 0; i < nComps.size(); i++) {
		for (int j = 0; j < nComps[i].size(); j++) {
			if (doubleNodes[nComps[i][j]].getType() == N) {
				doubleEdges.push_back({ nComps[i][j], piIndex });

				nCompEs += 1;
				break;
			}
		}
	}

}

class hyperNode {
public:
	hyperNode(std::vector<int>& subnodes, int specialType, int fgSide);
	int getType(); void setWeight(int64_t wt);
	int getSide(); void assignHyperNodeWt(std::vector<doubleNode>& doubleNodes, int64_t wtSum);
	std::vector<int> doubleSubnodes; double getWeight(); int decIndex;
private:
	int type; int side; int64_t weight;
};

hyperNode::hyperNode(std::vector<int>& subnodes, int specialType, int fgSide) {
	doubleSubnodes = subnodes; type = specialType; side = fgSide; decIndex = -1;
}

int hyperNode::getType() {
	return type;
}

int hyperNode::getSide() {
	return side;
}

double hyperNode::getWeight() {
	return weight;
}

void hyperNode::setWeight(int64_t wt) {
	weight = wt;
}



doubleNode::doubleNode(node n, int nodeType, int nodeSide) {
	type = nodeType;
	side = nodeSide;
	origNode = n;
}

int doubleNode::getType() {
	return type;
}

int doubleNode::getSide() {
	return side;
}

void hyperNode::assignHyperNodeWt(std::vector<doubleNode>& doubleNodes, int64_t wtSum) {

	if ((type == CORE) || (type == N) || (type == MIDT) || (type == PINODE)) {
		weight = INFINITY;
	}
	else {
		if (side == 1) {
			weight = 0;

			for (int i = 0; i < doubleSubnodes.size(); i++) {
				weight += (-(int64_t)doubleNodes[doubleSubnodes[i]].origNode.labelCost + wtSum);
			}
		}
		else {//Side must be  -1 (AKA background)
			weight = 0;
			for (int i = 0; i < doubleSubnodes.size(); i++) {
				weight += (((int64_t)doubleNodes[doubleSubnodes[i]].origNode.labelCost) + wtSum);
			}
		}
	}

}

class hyperGraph {
public:
	hyperGraph(doubleGraph& doubleGIn, int nodeSize, int64_t wtSum); //Constructor local stage
	hyperGraph(std::vector<node>& nodes, doubleGraph& doubleGIn, map< std::vector<int>, int>& edgeWts, Graph& G, tbb::concurrent_vector< hyperNode >& globalHypernodes, int64_t wtSum); //Constructor for global stage
	void constructHyperEdges();
	std::vector<hyperNode> hyperNodes;
	std::vector<std::vector<int> > hyperEdges;
	doubleGraph doubleG; int numHypernodes; int numTerminals;
	std::vector<int> coreIndices; std::vector<int> nIndices; int piIndex;
};


void hyperGraph::constructHyperEdges() {
	map<std::vector<int>, bool> edgeExists; Graph G;
	for (int i = 0; i < hyperNodes.size(); i++) {
		add_vertex(G);
	}
	for (int i = 0; i < doubleG.doubleEdges.size(); i++) {
		int v1 = doubleG.doubleEdges[i][0]; int v2 = doubleG.doubleEdges[i][1];
		std::vector<int> v1Hypernodes = doubleG.doubleNodes[v1].hypernodes; std::vector<int> v2Hypernodes = doubleG.doubleNodes[v2].hypernodes;

		for (int j = 0; j < v1Hypernodes.size(); j++) {
			for (int k = 0; k < v2Hypernodes.size(); k++) {
				if (v1Hypernodes[j] != v2Hypernodes[k]) {
					if (hyperNodes[v1Hypernodes[j]].getType() == HYPERNODE && hyperNodes[v2Hypernodes[k]].getType() == HYPERNODE) {
						continue;
					}

					std::vector<int> hEdge = { v1Hypernodes[j], v2Hypernodes[k] }; sort(hEdge.begin(), hEdge.end());
					if (edgeExists.find(hEdge) == edgeExists.end()) {

						hyperEdges.push_back(hEdge);
						edgeExists[hEdge] = true; edgeExists[{hEdge[1], hEdge[0]}] = true;
						add_edge(hEdge[0], hEdge[1], G);

					}

				}
			}
		}
	}

	edgeExists.clear();
	return;
}

hyperGraph::hyperGraph(doubleGraph& doubleGIn, int nodeSize, int64_t wtSum) {
	doubleG = doubleGIn; //Initialize hypergraph to have its first nodes be terminal nodes: cores, neighborhoods, middle terminals, pi node
	numHypernodes = 0; numTerminals = 0;
	for (int i = 0; i < doubleG.cores.size(); i++) {
		std::vector<int> coreV = { doubleG.cores[i] };
		hyperNodes.push_back(hyperNode(coreV, CORE, 1));
		numHypernodes += 1;
		numTerminals += 1;
		hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		coreIndices.push_back(numHypernodes - 1);
		doubleG.doubleNodes[doubleG.cores[i]].hypernodes.push_back(hyperNodes.size() - 1);
	}
	for (int i = 0; i < doubleG.ns.size(); i++) {
		std::vector<int> nV = { doubleG.ns[i] };
		hyperNodes.push_back(hyperNode(nV, N, -1));
		numHypernodes += 1; numTerminals += 1;
		hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		nIndices.push_back(numHypernodes - 1);
		doubleG.doubleNodes[doubleG.ns[i]].hypernodes.push_back(hyperNodes.size() - 1);

	}
	for (int i = 0; i < doubleG.midTs.size(); i++) {
		std::vector<int> midTV = { doubleG.midTs[i] };
		hyperNodes.push_back(hyperNode(midTV, MIDT, 0));
		numHypernodes += 1; numTerminals += 1;
		hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		doubleG.doubleNodes[doubleG.midTs[i]].hypernodes.push_back(hyperNodes.size() - 1);
	}
	std::vector<int> piV = { doubleG.getPiIndex() };
	hyperNodes.push_back(hyperNode(piV, PINODE, 0));
	numHypernodes += 1;
	numTerminals += 1;
	hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
	piIndex = numHypernodes - 1;
	doubleG.doubleNodes[doubleG.getPiIndex()].hypernodes.push_back(hyperNodes.size() - 1);
	//Add individual cuts on fg side and individual fills on bg side
	for (int i = 0; i < doubleG.doubleNodes.size(); i++) {
		if (doubleG.doubleNodes[i].getType() == CUT) {
			if (doubleG.doubleNodes[i].getSide() == 1) {
				std::vector<int> cutV = { i };
				hyperNodes.push_back(hyperNode(cutV, HYPERNODE, 1));
				numHypernodes += 1; hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
				doubleG.doubleNodes[i].hypernodes.push_back(hyperNodes.size() - 1);
			}
		}

		if (doubleG.doubleNodes[i].getType() == FILL) {
			if (doubleG.doubleNodes[i].getSide() == -1) {
				std::vector<int> fillV = { i };
				hyperNodes.push_back(hyperNode(fillV, HYPERNODE, -1));
				numHypernodes += 1;
				hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
				doubleG.doubleNodes[i].hypernodes.push_back(hyperNodes.size() - 1);
			}
		}
	}
	int hs0 = hyperNodes.size();
	map< std::vector<int>, bool> subsetExists;
	std::vector< std::vector<int> > fgComplexEdges = doubleG.getFgCompsComplex();
	std::vector< std::vector<int> > bgComplexEdges = doubleG.getBgCompsComplex();
	for (int i = 0; i < fgComplexEdges.size(); i++) {
		std::vector<int> fgComplexComp = fgComplexEdges[i];
		if (fgComplexComp.size() == 1 && doubleG.doubleNodes[fgComplexComp[0]].getType() == CUT) {
			continue;
		}
		else {
			std::vector<int> setFills;
			for (int j = 0; j < fgComplexComp.size(); j++) {
				if (doubleG.doubleNodes[fgComplexComp[j]].getType() == FILL) { setFills.push_back(fgComplexComp[j]); }
			}
			sort(setFills.begin(), setFills.end()); subsetExists[setFills] = true; //Record that this subset of fills exists
			std::vector<int> allCutsAndFills = setFills;
			for (int j = 0; j < setFills.size(); j++) {
				auto neighboursFg = adjacent_vertices(setFills[j], doubleG.complexG);
				for (auto vd : make_iterator_range(neighboursFg)) {
					if (doubleG.doubleNodes[vd].getType() == CUT) {
						allCutsAndFills.push_back(vd);
					}
				}
			}
			sort(allCutsAndFills.begin(), allCutsAndFills.end());
			allCutsAndFills.erase(unique(allCutsAndFills.begin(), allCutsAndFills.end()), allCutsAndFills.end()); //Sort and delete duplicates
			hyperNodes.push_back(hyperNode(allCutsAndFills, HYPERNODE, 1));
			numHypernodes += 1;
			hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
			for (int j = 0; j < allCutsAndFills.size(); j++) {

				doubleG.doubleNodes[allCutsAndFills[j]].hypernodes.push_back(hyperNodes.size() - 1);
			}
		}
	}

	for (int i = 0; i < bgComplexEdges.size(); i++) { //k = -1 (all components on fg side)
		std::vector<int> bgComplexComp = bgComplexEdges[i];
		if (bgComplexComp.size() == 1 && doubleG.doubleNodes[bgComplexComp[0]].getType() == FILL) {
			continue;
		}
		else {
			std::vector<int> setCuts;
			for (int j = 0; j < bgComplexComp.size(); j++) {
				if (doubleG.doubleNodes[bgComplexComp[j]].getType() == CUT) {
					setCuts.push_back(bgComplexComp[j]);
				}
			}
			sort(setCuts.begin(), setCuts.end()); subsetExists[setCuts] = true; //Record that this subset of cuts exists
			std::vector<int> allCutsAndFills = setCuts;
			for (int j = 0; j < setCuts.size(); j++) {
				auto neighboursBg = adjacent_vertices(setCuts[j], doubleG.complexG);
				for (auto vd : make_iterator_range(neighboursBg)) {
					if (doubleG.doubleNodes[vd].getType() == FILL) {
						allCutsAndFills.push_back(vd);
					}
				}
			}
			sort(allCutsAndFills.begin(), allCutsAndFills.end());
			allCutsAndFills.erase(unique(allCutsAndFills.begin(), allCutsAndFills.end()), allCutsAndFills.end()); //Sort and delete duplicates
			hyperNodes.push_back(hyperNode(allCutsAndFills, HYPERNODE, -1));
			numHypernodes += 1;
			hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
			for (int j = 0; j < allCutsAndFills.size(); j++) {

				doubleG.doubleNodes[allCutsAndFills[j]].hypernodes.push_back(hyperNodes.size() - 1);
			}
		}
	}
	int hs1 = hyperNodes.size() - hs0;
	//Construct size 0 nodes
	for (int i = 0; i < doubleG.doubleNodes.size(); i++) {
		if (doubleG.doubleNodes[i].getSide() == 1) {
			if (doubleG.doubleNodes[i].getType() == CUT) {
				//Find bg complement of individual cut
				int excludeIndex = i + 2;

				std::vector<int> excludeIndices = { excludeIndex };
				std::vector< std::vector<int> >  bgComplementComps = doubleG.findRemainingCFComponents(false, excludeIndices);

				for (int q = 0; q < bgComplementComps.size(); q++) {
					std::vector<int> currentBgComp = bgComplementComps[q];
					if (currentBgComp.size() == 0) { continue; }
					if (currentBgComp.size() == 1 && doubleG.doubleNodes[currentBgComp[0]].getType() == FILL) { continue; }
					else {
						std::vector<int> setCuts;
						for (int k = 0; k < currentBgComp.size(); k++) {
							if (doubleG.doubleNodes[currentBgComp[k]].getType() == CUT) {
								setCuts.push_back(currentBgComp[k]);
							}
						}
						if (setCuts.size() == 0) { continue; }
						sort(setCuts.begin(), setCuts.end());
						if (subsetExists.find(setCuts) == subsetExists.end()) { //Add node if the set of e-linked cuts in the component does not exist yet

							subsetExists[setCuts] = true;
							std::vector<int> allCutsAndFillsC = setCuts;
							for (int l = 0; l < setCuts.size(); l++) {
								auto neighboursBg = adjacent_vertices(setCuts[l], doubleG.complexG);
								for (auto vd : make_iterator_range(neighboursBg)) {
									if (doubleG.doubleNodes[vd].getType() == FILL) {
										allCutsAndFillsC.push_back(vd);
									}
								}
							}
							sort(allCutsAndFillsC.begin(), allCutsAndFillsC.end());
							allCutsAndFillsC.erase(unique(allCutsAndFillsC.begin(), allCutsAndFillsC.end()), allCutsAndFillsC.end()); //Sort and delete duplicates
							hyperNodes.push_back(hyperNode(allCutsAndFillsC, HYPERNODE, -1));
							numHypernodes += 1;
							hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
							for (int l = 0; l < allCutsAndFillsC.size(); l++) {
								doubleG.doubleNodes[allCutsAndFillsC[l]].hypernodes.push_back(hyperNodes.size() - 1);
							}
						}
					}
				}
			}
		}
		else {
			if (doubleG.doubleNodes[i].getSide() == -1) {
				if (doubleG.doubleNodes[i].getType() == FILL) {
					//Find fg complement of individual fill
					int excludeIndex = i - 2;

					std::vector<int> excludeIndices = { excludeIndex };
					std::vector< std::vector<int> > fgComplementComps = doubleG.findRemainingCFComponents(true, excludeIndices);
					for (int q = 0; q < fgComplementComps.size(); q++) {
						std::vector<int> currentFgComp = fgComplementComps[q];
						if (currentFgComp.size() == 0) { continue; }
						if (currentFgComp.size() == 1 && doubleG.doubleNodes[currentFgComp[0]].getType() == CUT) { continue; }
						else {
							std::vector<int> setFills;
							for (int k = 0; k < currentFgComp.size(); k++) {
								if (doubleG.doubleNodes[currentFgComp[k]].getType() == FILL) {
									setFills.push_back(currentFgComp[k]);
								}
							}
							if (setFills.size() == 0) { continue; }
							sort(setFills.begin(), setFills.end());
							if (subsetExists.find(setFills) == subsetExists.end()) { //Add node if the set of e-linked fills in the component does not exist yet
								subsetExists[setFills] = true;
								std::vector<int> allCutsAndFillsC = setFills;
								for (int l = 0; l < setFills.size(); l++) {
									auto neighboursFg = adjacent_vertices(setFills[l], doubleG.complexG);
									for (auto vd : make_iterator_range(neighboursFg)) {
										if (doubleG.doubleNodes[vd].getType() == CUT) {
											allCutsAndFillsC.push_back(vd);
										}
									}
								}
								sort(allCutsAndFillsC.begin(), allCutsAndFillsC.end());
								allCutsAndFillsC.erase(unique(allCutsAndFillsC.begin(), allCutsAndFillsC.end()), allCutsAndFillsC.end()); //Sort and delete duplicates
								hyperNodes.push_back(hyperNode(allCutsAndFillsC, HYPERNODE, 1));
								numHypernodes += 1;
								hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
								for (int l = 0; l < allCutsAndFillsC.size(); l++) {
									doubleG.doubleNodes[allCutsAndFillsC[l]].hypernodes.push_back(hyperNodes.size() - 1);
								}
							}
						}
					}

				}
			}
		}
	}

	//Recursively build up hypernodes up to the maximum hypernode size
	std::vector< std::vector<int> > nodesLastRoundFg; std::vector< std::vector<int> > nodesLastRoundBg; int fgCt = 0; int bgCt = 0;
	for (int i = 0; i < doubleG.getSize(); i++) {
		if (doubleG.doubleNodes[i].getSide() == 1) {
			fgCt += 1;
			if (doubleG.doubleNodes[i].getType() == FILL) {
				nodesLastRoundFg.push_back({ i });
			}
		}
		else {
			if (doubleG.doubleNodes[i].getSide() == -1) {
				bgCt += 1;
				if (doubleG.doubleNodes[i].getType() == CUT) {
					nodesLastRoundBg.push_back({ i });
				}
			}
		}
	}
	int maxNodeSize = min(nodeSize, max(fgCt, bgCt)); //If user accidentally specifies an infinitely large node size, this limits the amount of node building iterations

	//cout << "hyp size " << nodeSize << " " << maxNodeSize << endl;
//Construct nodes up to user-specified maximum size
	for (int i = 0; i < maxNodeSize; i++) {
		//Construct FG nodes and its BG complement
		std::vector< std::vector<int> > currentHyperNodesFg = nodesLastRoundFg;
		std::vector< std::vector<int> > newNodesFg;
		map< std::vector<int>, bool> newNodeExistsFg;
		for (int j = 0; j < currentHyperNodesFg.size(); j++) {
			std::vector<int> eLinkedFills = currentHyperNodesFg[j];
			std::vector<int> allCutsAndFills = eLinkedFills;
			std::sort(eLinkedFills.begin(), eLinkedFills.end());
			if (subsetExists.find(eLinkedFills) == subsetExists.end()) {
				for (int k = 0; k < eLinkedFills.size(); k++) {
					auto neighbourCuts = adjacent_vertices(eLinkedFills[k], doubleG.complexG);
					for (auto vd : make_iterator_range(neighbourCuts)) {
						if (doubleG.doubleNodes[vd].getType() == CUT) {
							allCutsAndFills.push_back(vd);
							auto neighbourFills = adjacent_vertices(vd, doubleG.complexG);
							for (auto vf : make_iterator_range(neighbourFills)) {
								if (doubleG.doubleNodes[vf].getType() == FILL) {
									if (std::find(eLinkedFills.begin(), eLinkedFills.end(), vf) == eLinkedFills.end()) {
										std::vector<int> newNode = eLinkedFills;
										newNode.push_back(vf);
										std::sort(newNode.begin(), newNode.end());
										if (newNodeExistsFg.find(newNode) == newNodeExistsFg.end()) {
											newNodesFg.push_back(newNode);
											newNodeExistsFg[newNode] = true;
										}
									}
								}
							}
						}
					}
				}
				std::sort(allCutsAndFills.begin(), allCutsAndFills.end());
				allCutsAndFills.erase(unique(allCutsAndFills.begin(), allCutsAndFills.end()), allCutsAndFills.end());
				hyperNodes.push_back(hyperNode(allCutsAndFills, HYPERNODE, 1));
				numHypernodes += 1;
				hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
				subsetExists[eLinkedFills] = true;
				std::vector<int> excludeIndices;
				for (int k = 0; k < allCutsAndFills.size(); k++) {
					doubleG.doubleNodes[allCutsAndFills[k]].hypernodes.push_back(hyperNodes.size() - 1);
					//for each cut and fill simulate removing from fg graph, find complement of component
					int excludeIndex = allCutsAndFills[k] + 2;
					excludeIndices.push_back(excludeIndex);

				}
				if (excludeIndices.size() > 0) {
					std::vector< std::vector<int> > bgComplementComps = doubleG.findRemainingCFComponents(false, excludeIndices);

					for (int q = 0; q < bgComplementComps.size(); q++) {
						std::vector<int> currentBgComp = bgComplementComps[q];
						if (currentBgComp.size() == 0) { continue; }
						if (currentBgComp.size() == 1 && doubleG.doubleNodes[currentBgComp[0]].getType() == FILL) { continue; }
						else {
							std::vector<int> setCuts;
							for (int k = 0; k < currentBgComp.size(); k++) {
								if (doubleG.doubleNodes[currentBgComp[k]].getType() == CUT) {
									setCuts.push_back(currentBgComp[k]);
								}
							}
							if (setCuts.size() == 0) { continue; }
							sort(setCuts.begin(), setCuts.end());
							if (subsetExists.find(setCuts) == subsetExists.end()) { //Add node if the set of e-linked cuts in the component does not exist yet

								subsetExists[setCuts] = true;
								std::vector<int> allCutsAndFillsC = setCuts;
								for (int l = 0; l < setCuts.size(); l++) {
									auto neighboursBg = adjacent_vertices(setCuts[l], doubleG.complexG);
									for (auto vd : make_iterator_range(neighboursBg)) {
										if (doubleG.doubleNodes[vd].getType() == FILL) {
											allCutsAndFillsC.push_back(vd);
										}
									}
								}
								sort(allCutsAndFillsC.begin(), allCutsAndFillsC.end());
								allCutsAndFillsC.erase(unique(allCutsAndFillsC.begin(), allCutsAndFillsC.end()), allCutsAndFillsC.end()); //Sort and delete duplicates
								hyperNodes.push_back(hyperNode(allCutsAndFillsC, HYPERNODE, -1));
								numHypernodes += 1;
								hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
								for (int l = 0; l < allCutsAndFillsC.size(); l++) {
									doubleG.doubleNodes[allCutsAndFillsC[l]].hypernodes.push_back(hyperNodes.size() - 1);
								}
							}
						}
					}
				}
			}
		}
		nodesLastRoundFg = newNodesFg;

		//Construct BG nodes and its FG complement
		std::vector< std::vector<int> > currentHyperNodesBg = nodesLastRoundBg;
		std::vector< std::vector<int> > newNodesBg;
		map< std::vector<int>, bool> newNodeExistsBg;
		for (int j = 0; j < currentHyperNodesBg.size(); j++) {
			std::vector<int> eLinkedCuts = currentHyperNodesBg[j];
			std::vector<int> allCutsAndFills = eLinkedCuts;
			std::sort(eLinkedCuts.begin(), eLinkedCuts.end());
			if (subsetExists.find(eLinkedCuts) == subsetExists.end()) {
				for (int k = 0; k < eLinkedCuts.size(); k++) {
					auto neighbourFills = adjacent_vertices(eLinkedCuts[k], doubleG.complexG);
					for (auto vd : make_iterator_range(neighbourFills)) {
						if (doubleG.doubleNodes[vd].getType() == FILL) {
							allCutsAndFills.push_back(vd);
							auto neighbourCuts = adjacent_vertices(vd, doubleG.complexG);
							for (auto vf : make_iterator_range(neighbourCuts)) {
								if (doubleG.doubleNodes[vf].getType() == CUT) {
									if (std::find(eLinkedCuts.begin(), eLinkedCuts.end(), vf) == eLinkedCuts.end()) {
										std::vector<int> newNode = eLinkedCuts; newNode.push_back(vf);
										std::sort(newNode.begin(), newNode.end());
										if (newNodeExistsBg.find(newNode) == newNodeExistsBg.end()) {
											newNodesBg.push_back(newNode);
											newNodeExistsBg[newNode] = true;
										}
									}
								}
							}
						}
					}
				}

				subsetExists[eLinkedCuts] = true;
				std::sort(allCutsAndFills.begin(), allCutsAndFills.end());
				allCutsAndFills.erase(unique(allCutsAndFills.begin(), allCutsAndFills.end()), allCutsAndFills.end());
				hyperNodes.push_back(hyperNode(allCutsAndFills, HYPERNODE, -1));
				numHypernodes += 1;
				hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);

				std::vector<int> excludeIndices;
				for (int k = 0; k < allCutsAndFills.size(); k++) {

					doubleG.doubleNodes[allCutsAndFills[k]].hypernodes.push_back(hyperNodes.size() - 1);
					//for each cut and fill simulate removing from bg graph, find complement of component
					int excludeIndex = allCutsAndFills[k] - 2;
					excludeIndices.push_back(excludeIndex);
				}

				if (excludeIndices.size() > 0) {
					std::vector< std::vector<int> > fgComplementComps = doubleG.findRemainingCFComponents(true, excludeIndices);
					for (int q = 0; q < fgComplementComps.size(); q++) {
						std::vector<int> currentFgComp = fgComplementComps[q];
						if (currentFgComp.size() == 0) { continue; }
						if (currentFgComp.size() == 1 && doubleG.doubleNodes[currentFgComp[0]].getType() == CUT) { continue; }
						else {
							std::vector<int> setFills;
							for (int k = 0; k < currentFgComp.size(); k++) {
								if (doubleG.doubleNodes[currentFgComp[k]].getType() == FILL) {
									setFills.push_back(currentFgComp[k]);
								}
							}
							if (setFills.size() == 0) { continue; }
							sort(setFills.begin(), setFills.end());
							if (subsetExists.find(setFills) == subsetExists.end()) { //Add node if the set of e-linked fills in the component does not exist yet
								subsetExists[setFills] = true;
								std::vector<int> allCutsAndFillsC = setFills;
								for (int l = 0; l < setFills.size(); l++) {
									auto neighboursFg = adjacent_vertices(setFills[l], doubleG.complexG);
									for (auto vd : make_iterator_range(neighboursFg)) {
										if (doubleG.doubleNodes[vd].getType() == CUT) {
											allCutsAndFillsC.push_back(vd);
										}
									}
								}
								sort(allCutsAndFillsC.begin(), allCutsAndFillsC.end());
								allCutsAndFillsC.erase(unique(allCutsAndFillsC.begin(), allCutsAndFillsC.end()), allCutsAndFillsC.end()); //Sort and delete duplicates
								hyperNodes.push_back(hyperNode(allCutsAndFillsC, HYPERNODE, 1));
								numHypernodes += 1; hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
								for (int l = 0; l < allCutsAndFillsC.size(); l++) {
									doubleG.doubleNodes[allCutsAndFillsC[l]].hypernodes.push_back(hyperNodes.size() - 1);
								}
							}
						}
					}
				}
			}
		}
		nodesLastRoundBg = newNodesBg;
	}

	//Construct edges of hypergraph
	constructHyperEdges();
	//cout << "hyperg size " << hyperNodes.size() << " " << hyperEdges.size() << " " << wtSum << endl;
}

hyperGraph::hyperGraph(std::vector<node>& nodes, doubleGraph& doubleGIn, map< std::vector<int>, int>& edgeWts, Graph& G, tbb::concurrent_vector< hyperNode >& globalHypernodes,
	int64_t wtSum) {
	doubleG = doubleGIn; //Initialize hypergraph to have its first nodes be terminal nodes: cores, neighborhoods, middle terminals, pi node
	numHypernodes = 0; numTerminals = 0;
	for (int i = 0; i < doubleG.cores.size(); i++) {
		std::vector<int> coreV = { doubleG.cores[i] };
		hyperNodes.push_back(hyperNode(coreV, CORE, 1));
		numHypernodes += 1;
		numTerminals += 1;
		hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		coreIndices.push_back(numHypernodes - 1);
		doubleG.doubleNodes[doubleG.cores[i]].hypernodes.push_back(hyperNodes.size() - 1);
	}
	for (int i = 0; i < doubleG.ns.size(); i++) {
		std::vector<int> nV = { doubleG.ns[i] };
		hyperNodes.push_back(hyperNode(nV, N, -1));
		numHypernodes += 1; numTerminals += 1;
		hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		nIndices.push_back(numHypernodes - 1);
		doubleG.doubleNodes[doubleG.ns[i]].hypernodes.push_back(hyperNodes.size() - 1);
	}
	for (int i = 0; i < doubleG.midTs.size(); i++) {
		std::vector<int> midTV = { doubleG.midTs[i] };
		hyperNodes.push_back(hyperNode(midTV, MIDT, 0));
		numHypernodes += 1; numTerminals += 1;
		hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		doubleG.doubleNodes[doubleG.midTs[i]].hypernodes.push_back(hyperNodes.size() - 1);
	}
	std::vector<int> piV = { doubleG.getPiIndex() };
	hyperNodes.push_back(hyperNode(piV, PINODE, 0));
	numHypernodes += 1;
	numTerminals += 1;
	hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
	piIndex = numHypernodes - 1;
	doubleG.doubleNodes[doubleG.getPiIndex()].hypernodes.push_back(hyperNodes.size() - 1);
	//Create hypernodes 
	for (int i = 0; i < globalHypernodes.size(); i++) {
		std::vector<int> subnodes = globalHypernodes[i].doubleSubnodes;
		std::vector<int> doubleSubnodes;
		int side = globalHypernodes[i].getSide(); //All subnodes must be cuts or fills
		if (side == 1) {
			for (int j = 0; j < subnodes.size(); j++) {
				doubleSubnodes.push_back(doubleG.nodeToDoubleMapping[subnodes[j]][0]);
				doubleG.doubleNodes[doubleG.nodeToDoubleMapping[subnodes[j]][0]].hypernodes.push_back(hyperNodes.size());
			}
			hyperNodes.push_back(hyperNode(doubleSubnodes, HYPERNODE, 1));
			numHypernodes += 1;
			hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		}
		else { // side must be -1
			for (int j = 0; j < subnodes.size(); j++) {
				doubleSubnodes.push_back(doubleG.nodeToDoubleMapping[subnodes[j]][2]);
				doubleG.doubleNodes[doubleG.nodeToDoubleMapping[subnodes[j]][2]].hypernodes.push_back(hyperNodes.size());
			}
			hyperNodes.push_back(hyperNode(doubleSubnodes, HYPERNODE, -1));
			numHypernodes += 1;
			hyperNodes[numHypernodes - 1].assignHyperNodeWt(doubleG.doubleNodes, wtSum);
		}
	}

	constructHyperEdges();
}

//Find 2 * sum of node weights, to ensure hypergraph nodes always have positive weights
int64_t getWtSum(std::vector<node>& nodes) {
	int64_t nodeSum = 0;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].valid) {
			if (((int)nodes[i].type) == 2 || ((int)nodes[i].type) == 3) {
				nodeSum += (2 * abs((int64_t)nodes[i].labelCost));
			}
		}
	}
	cout.precision(14);
	nodeSum = std::min(nodeSum, (int64_t)(WMAX / (nodes.size() * 2)));
	return nodeSum;
}

std::vector< std::vector<int>> getCfComponents(Graph& g, std::vector<node>& nodes, int numNodes) {
	std::vector< std::vector<int> > typeComponents;
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(g, &nodeToComp[0]);
	int numComps = 0; std::vector<int> isCompIndex(n, -1);
	for (int i = 0; i < numNodes; i++) {

		if (((int)nodes[i].type) == 2 || ((int)nodes[i].type) == 3) {
			if (isCompIndex[nodeToComp[i]] == -1) {

				isCompIndex[nodeToComp[i]] = numComps;
				std::vector<int> newComp = { i };
				typeComponents.push_back(newComp);
				numComps += 1;
			}
			else {
				typeComponents[isCompIndex[nodeToComp[i]]].push_back(i);
			}

		}
	}
	return typeComponents;
}

void getLocalGraphClusters(Graph& G, std::vector<node>& nodes, int numNodes, Graph& fgG, Graph& bgG, std::vector< std::vector< std::vector<int> > >& localGraphEdges, std::vector< std::vector<int> >& localNodesGlobalIndex, map< std::vector<int>, int>& edgeWt) {
	Graph cutFillGraph;
	for (int i = 0; i < numNodes; i++) {
		add_vertex(cutFillGraph); add_vertex(fgG); add_vertex(bgG);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 2) || (nodes[v2].type == 3))) {
			add_edge(v1, v2, cutFillGraph);
		}
		if (((nodes[v1].type == 0) || (nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 0) || (nodes[v2].type == 2) || (nodes[v2].type == 3)) && edgeWt[{v1, v2}] == 1) {
			add_edge(v1, v2, fgG);
		}
		if (((nodes[v1].type == 1) || (nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 1) || (nodes[v2].type == 2) || (nodes[v2].type == 3))) {
			add_edge(v1, v2, bgG);
		}

	}
	std::vector< std::vector<int>> cfComponents = getCfComponents(cutFillGraph, nodes, numNodes);
	std::vector< std::vector<int> > globalToLocalNodes(numNodes, std::vector<int>(0));

	for (int i = 0; i < cfComponents.size(); i++) {
		std::vector<int> cfWithTerminals = cfComponents[i];
		map<int, bool> termAdded;
		for (int j = 0; j < cfComponents[i].size(); j++) {
			globalToLocalNodes[cfComponents[i][j]].push_back(i); //Map each node to the cluster which contains it

			auto neighbours = adjacent_vertices(cfComponents[i][j], G);
			for (auto vd : make_iterator_range(neighbours)) {
				if (((int)nodes[vd].type) == 0 || ((int)nodes[vd].type) == 1) {
					if (termAdded.find(vd) == termAdded.end()) {
						cfWithTerminals.push_back(vd);
						termAdded[vd] = true;
					}
				}
			}
		}
		localNodesGlobalIndex.push_back(cfWithTerminals);
		std::vector< std::vector<int> > localEdges; localGraphEdges.push_back(localEdges);
	}

	//Iterate through edges, find subgraph for each local cut-fill cluster in one sweep
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((int)nodes[v1].type == 2) || ((int)nodes[v1].type == 3)) {
			localGraphEdges[globalToLocalNodes[v1][0]].push_back({ v1, v2 });
			continue;
		}
		else {
			if (((int)nodes[v2].type == 2) || ((int)nodes[v2].type == 3)) {
				localGraphEdges[globalToLocalNodes[v2][0]].push_back({ v1, v2 });
				continue;
			}

		}
	}
}

void getCoreG(Graph& G, std::vector<node>& nodes, int numNodes, Graph& coreG, map<std::vector<int>, int>& edgeWts) {
	for (int i = 0; i < numNodes; i++) {
		add_vertex(coreG);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if ((int)nodes[v1].type == CORE && (int)nodes[v2].type == CORE) {
			if (edgeWts[{v1, v2}] == 1) {
				add_edge(v1, v2, coreG);
			}
		}
	}
}

void getNG(Graph& G, std::vector<node>& nodes, int numNodes, Graph& nG) {
	for (int i = 0; i < numNodes; i++) {
		add_vertex(nG);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if ((int)nodes[v1].type == N && (int)nodes[v2].type == N) {
			add_edge(v1, v2, nG);
		}
	}
}


int findComponents(Graph& G, std::vector<node>& nodes, int numNodes, bool inFg) {
	std::vector< std::vector<int> > components;
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(G, &nodeToComp[0]);
	int numComps = 0; std::vector<int> isCompIndex(n, -1);
	//If in fg, core. Else if in bg, neighborhood
	int type = 0;
	if (!inFg) {
		type = 1;
	}
	for (int i = 0; i < numNodes; i++) {
		if (nodes[i].valid) {
			if (((int)nodes[i].type) == type || ((int)nodes[i].type) == 2 || ((int)nodes[i].type) == 3) {
				if (isCompIndex[nodeToComp[i]] == -1) {

					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i);
				}

			}
		}
	}
	int compsWithTerms = 0;
	for (int i = 0; i < components.size(); i++) {
		bool containsTerm = false;
		for (int j = 0; j < components[i].size(); j++) {
			if (((int)nodes[components[i][j]].type) == type) {
				containsTerm = true;
				break;
			}
		}
		if (containsTerm) {
			compsWithTerms += 1;
		}
	}
	return compsWithTerms;
}

std::vector< std::vector<int> > getComponents(Graph& G, std::vector<node>& nodes, int inFg) {
	std::vector< std::vector<int> > components;
	int numNodes = nodes.size();
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(G, &nodeToComp[0]);
	int numComps = 0; std::vector<int> isCompIndex(n, -1);
	//If in fg, core. Else if in bg, neighborhood
	for (int i = 0; i < numNodes; i++) {
		if (nodes[i].valid) {
			if (((int)nodes[i].inFg) == inFg) {
				if (isCompIndex[nodeToComp[i]] == -1) {

					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { i };
					components.push_back(newComp);
					numComps += 1;
				}
				else {
					components[isCompIndex[nodeToComp[i]]].push_back(i);
				}

			}
		}
	}
	return components;
}

void getCompToLocalIndex(std::vector<int>& clusterLocalNodesGlobalIndices, map<int, int>& localToGlobal, map<int, int>& globalToLocal, map< std::vector<int>, int>& edgeWts) {
	for (int i = 0; i < clusterLocalNodesGlobalIndices.size(); i++) {
		globalToLocal[clusterLocalNodesGlobalIndices[i]] = i; localToGlobal[i] = clusterLocalNodesGlobalIndices[i];
	}
}

void connectTerminalsExternally(hyperGraph& hypG, Graph& coreGIn, Graph& nGIn, int numNodes, std::vector<int>& clusterLocalNodesGlobalIndices, std::vector< std::vector<int> >& clusterLocalEdgesGlobalIndices, map<std::vector<int>, bool>& edgeExists) {
	Graph coreG = coreGIn;
	Graph nG = nGIn;
	for (int i = 0; i < clusterLocalEdgesGlobalIndices.size(); i++) {
		remove_edge(clusterLocalEdgesGlobalIndices[i][0], clusterLocalEdgesGlobalIndices[i][1], coreG);
		remove_edge(clusterLocalEdgesGlobalIndices[i][0], clusterLocalEdgesGlobalIndices[i][1], nG);
	}
	map<int, int> nodeToHyper;
	for (int i = 0; i < hypG.coreIndices.size(); i++) {
		nodeToHyper[hypG.doubleG.doubleNodes[hypG.hyperNodes[hypG.coreIndices[i]].doubleSubnodes[0]].origNode.index] = hypG.coreIndices[i];
	}
	for (int i = 0; i < hypG.nIndices.size(); i++) {
		nodeToHyper[hypG.doubleG.doubleNodes[hypG.hyperNodes[hypG.nIndices[i]].doubleSubnodes[0]].origNode.index] = hypG.nIndices[i];
	}

	std::vector<int> nodeToCompFg(numNodes);
	int nFg = (int)boost::connected_components(coreG, &nodeToCompFg[0]);
	int numCompsFg = 0; std::vector<int> isCompIndexFg(nFg, -1);
	std::vector< std::vector<int> > coreComponents;
	for (int i = 0; i < hypG.coreIndices.size(); i++) {
		int coreNIndex = hypG.doubleG.doubleNodes[hypG.hyperNodes[hypG.coreIndices[i]].doubleSubnodes[0]].origNode.index;
		if (isCompIndexFg[nodeToCompFg[coreNIndex]] == -1) {
			isCompIndexFg[nodeToCompFg[coreNIndex]] = numCompsFg;
			std::vector<int> newComp = { nodeToHyper[coreNIndex] };
			coreComponents.push_back(newComp);
			numCompsFg += 1;
		}
		else {
			coreComponents[isCompIndexFg[nodeToCompFg[coreNIndex]]].push_back(nodeToHyper[coreNIndex]);
		}
	}

	std::vector<int> nodeToCompBg(numNodes);
	int nBg = (int)boost::connected_components(nG, &nodeToCompBg[0]);
	int numCompsBg = 0;
	std::vector<int> isCompIndexBg(nBg, -1);
	std::vector< std::vector<int> > nComponents;
	for (int i = 0; i < hypG.nIndices.size(); i++) {
		int nNIndex = hypG.doubleG.doubleNodes[hypG.hyperNodes[hypG.nIndices[i]].doubleSubnodes[0]].origNode.index;
		if (isCompIndexBg[nodeToCompBg[nNIndex]] == -1) {
			isCompIndexBg[nodeToCompBg[nNIndex]] = numCompsBg;
			std::vector<int> newComp = { nodeToHyper[nNIndex] };
			nComponents.push_back(newComp);
			numCompsBg += 1;
		}
		else {
			nComponents[isCompIndexBg[nodeToCompBg[nNIndex]]].push_back(nodeToHyper[nNIndex]);
		}
	}



	for (int i = 0; i < coreComponents.size(); i++) {
		for (int j = 0; j < coreComponents[i].size() - 1; j++) {
			std::vector<int> e = { coreComponents[i][j], coreComponents[i][j + 1] };
			if (edgeExists.find(e) == edgeExists.end()) {
				hypG.hyperEdges.push_back(e);
				edgeExists[e] = true; edgeExists[{e[1], e[0]}] = true;
			}
		}
	}

	for (int i = 0; i < nComponents.size(); i++) {
		for (int j = 0; j < nComponents[i].size() - 1; j++) {
			std::vector<int> e = { nComponents[i][j], nComponents[i][j + 1] };
			if (edgeExists.find(e) == edgeExists.end()) {
				hypG.hyperEdges.push_back(e);
				edgeExists[e] = true; edgeExists[{e[1], e[0]}] = true;
			}
		}
	}


}

//Check if terminals connected outside this cluster
std::vector< std::vector<int> > getTermsConnectedOutsideCluster(Graph& inG, int type, std::vector<node>& nodes, int numNodes, std::vector<int>& clusterGlobalIndices, hyperGraph& hyp) {
	std::vector<bool> inCluster(numNodes, false); Graph g = inG; std::vector<bool> termInComp(numNodes, false); map<int, int> nodeToHyper;
	if (type == CORE) {
		for (int i = 0; i < hyp.coreIndices.size(); i++) {
			//For each term index in hypergraph, check if in same component in original graph
			termInComp[hyp.doubleG.doubleNodes[hyp.hyperNodes[hyp.coreIndices[i]].doubleSubnodes[0]].origNode.index] = true;
			nodeToHyper[hyp.doubleG.doubleNodes[hyp.hyperNodes[hyp.coreIndices[i]].doubleSubnodes[0]].origNode.index] = hyp.coreIndices[i];
		}
		typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
		edge_iter ei, ei_end;
		for (tie(ei, ei_end) = edges(inG); ei != ei_end; ++ei) {
			int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
			if ((int)nodes[v1].type == CORE && (int)nodes[v2].type == CORE) {
				remove_edge(v1, v2, g);
			}
		}

	}
	else { //type must be neighborhood
		for (int i = 0; i < hyp.nIndices.size(); i++) {
			//For each term index in hypergraph, check if in same component in original graph
			termInComp[hyp.doubleG.doubleNodes[hyp.hyperNodes[hyp.nIndices[i]].doubleSubnodes[0]].origNode.index] = true;
			nodeToHyper[hyp.doubleG.doubleNodes[hyp.hyperNodes[hyp.nIndices[i]].doubleSubnodes[0]].origNode.index] = hyp.nIndices[i];
		}
		typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
		edge_iter ei, ei_end;
		for (tie(ei, ei_end) = edges(inG); ei != ei_end; ++ei) {
			int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
			if ((int)nodes[v1].type == N && (int)nodes[v2].type == N) {
				remove_edge(v1, v2, g);
			}
		}
	}

	for (int i = 0; i < clusterGlobalIndices.size(); i++) {
		if (nodes[clusterGlobalIndices[i]].type == CORE || nodes[clusterGlobalIndices[i]].type == N) {
			continue;
		}
		inCluster[clusterGlobalIndices[i]] = true;
		clear_vertex(clusterGlobalIndices[i], g);
	}

	std::vector< std::vector<int> > termComponents;
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(g, &nodeToComp[0]);
	int numComps = 0; std::vector<int> isCompIndex(n, -1);
	for (int i = 0; i < numNodes; i++) {
		if (termInComp[i]) {
			if (isCompIndex[nodeToComp[i]] == -1) {

				isCompIndex[nodeToComp[i]] = numComps;
				std::vector<int> newComp = { nodeToHyper[i] };
				termComponents.push_back(newComp);
				numComps += 1;
			}
			else {
				termComponents[isCompIndex[nodeToComp[i]]].push_back(nodeToHyper[i]);
			}

		}
	}
	return termComponents;
}

//Check if terminals connected outside this cluster
std::vector< std::vector<int> > getTermsConnectedOutsideCluster(Graph& inG, int type, std::vector<node>& nodes, int numNodes, std::vector<int>& clusterGlobalIndices,
	hyperGraph& hyp, map<int, int>& nodeToCompIn) {
	std::vector<bool> inCluster(numNodes, false); Graph g = inG; std::vector<bool> termInComp(numNodes, false); map<int, int> nodeToHyper;


	for (int i = 0; i < clusterGlobalIndices.size(); i++) {
		if (nodes[clusterGlobalIndices[i]].type == CORE || nodes[clusterGlobalIndices[i]].type == N) {
			if (nodeToCompIn.find(clusterGlobalIndices[i]) == nodeToCompIn.end()) {
				continue;
			}
		}

		inCluster[clusterGlobalIndices[i]] = true;
		clear_vertex(clusterGlobalIndices[i], g);
	}

	for (int i = 0; i < hyp.hyperEdges.size(); i++) {
		if ((hyp.hyperNodes[hyp.hyperEdges[i][0]].getType() == CORE && hyp.hyperNodes[hyp.hyperEdges[i][1]].getType() == CORE)
			||
			(hyp.hyperNodes[hyp.hyperEdges[i][0]].getType() == N && hyp.hyperNodes[hyp.hyperEdges[i][1]].getType() == N)
			) {
			int v1 = hyp.doubleG.doubleNodes[hyp.hyperNodes[hyp.hyperEdges[i][0]].doubleSubnodes[0]].origNode.index;
			int v2 = hyp.doubleG.doubleNodes[hyp.hyperNodes[hyp.hyperEdges[i][1]].doubleSubnodes[0]].origNode.index;
		}
	}


	std::vector< std::vector<int> > termComponents;
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(g, &nodeToComp[0]);
	int numComps = 0; std::vector<int> isCompIndex(n, -1);
	for (int i = 0; i < numNodes; i++) {
		if (termInComp[i]) {
			if (!inCluster[i]) {
				if (isCompIndex[nodeToComp[i]] == -1) {

					isCompIndex[nodeToComp[i]] = numComps;
					std::vector<int> newComp = { nodeToHyper[i] };
					termComponents.push_back(newComp);
					numComps += 1;
				}
				else {
					termComponents[isCompIndex[nodeToComp[i]]].push_back(nodeToHyper[i]);
				}
			}
		}
	}
	return termComponents;
}

std::vector<std::vector<std::vector<int>>> getPartitions(const std::vector<int>& elements) {
	std::vector<std::vector<std::vector<int>>> fList;

	std::vector<std::vector<int>> lists;
	std::vector<int> indexes(elements.size(), 0);
	lists.emplace_back(std::vector<int>());
	lists[0].insert(lists[0].end(), elements.begin(), elements.end());

	int counter = -1;

	for (;;) {
		counter += 1;
		fList.emplace_back(lists);

		int i, index;
		bool obreak = false;
		for (i = indexes.size() - 1;; --i) {
			if (i <= 0) {
				obreak = true;
				break;
			}
			index = indexes[i];
			lists[index].erase(lists[index].begin() + lists[index].size() - 1);
			if (lists[index].size() > 0)
				break;
			lists.erase(lists.begin() + index);
		}
		if (obreak) break;

		++index;
		if (index >= lists.size())
			lists.emplace_back(std::vector<int>());
		for (; i < indexes.size(); ++i) {
			indexes[i] = index;
			lists[index].emplace_back(elements[i]);
			index = 0;
		}

	}
	return fList;
}

void findGenVoxels(const std::vector<uint32_t*>& labels, std::vector <Coordinate>& voxels, int labelIndex, int width, int height, int numSlices,
	int x, int y, int z
) {
	map<Coordinate, bool> visited;
	queue<Coordinate> q;
	Coordinate cp(x, y, z);
	q.push(cp);
	while (!q.empty()) {
		Coordinate c = q.front(); q.pop();
		visited[c] = true;
		voxels.push_back(c);
		for (int i = 0; i < structCube.size(); i++) {
			int nx = c.x + structCube[i][0]; int ny = c.y + structCube[i][1]; int nz = c.z + structCube[i][2];
			if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
				if (Label(labels[nz][nx + ny * width]) == labelIndex) {
					Coordinate cn(nx, ny, nz);
					if (visited.find(cn) == visited.end()) {
						q.push(cn);
						visited[cn] = true;

					}
				}
			}
		}
	}

}

struct genVoxels {
	std::vector<Coordinate> coords;
	node n;

};

struct CompareGen {
	bool operator()(genVoxels const& p1, genVoxels const& p2)
	{
		return abs(p1.n.labelCost) < abs(p2.n.labelCost);
	}
};

void dfsOverall(grapht& g, int v, std::vector<bool>& visited, int& timer, int p, int& overallIndex, std::vector<node>& nodes, bool inFg) {
	visited[v] = true;
	if (inFg) {
		nodes[v].overallCompIndexFg = overallIndex;
	}
	else {
		nodes[v].overallCompIndexBg = overallIndex;
	}
	//nodes[v].tin = timer; nodes[v].low = timer;
	timer += 1;

	std::queue<int> q;
	q.push(v);
	while (!q.empty()) {
		int n = q.front();
		q.pop();
		int numNs = 0;
		auto neighbours = adjacent_vertices(n, g);
		/**if (nodes.size() == 8091) {
			cout << n << " of type " << (int)nodes[n].type << endl;
			int ct = 0;
			typedef boost::graph_traits<grapht>::edge_iterator edge_iter;
			edge_iter ei, ei_end;
			for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
				int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
				if ( ((int)n) == v1 || ((int)n) == v2) {
					cout << v1 << " " << v2 << endl;
					ct++;
				}

			}
			cout << "num edges " << ct << endl;
		}**/
		for (auto u : make_iterator_range(neighbours)) {
			if (!visited[u]) {
				if (nodes[u].valid) {
					q.push(u);
					visited[u] = true;
					if (inFg) {
						nodes[u].overallCompIndexFg = overallIndex;
					}
					else {
						nodes[u].overallCompIndexBg = overallIndex;
					}
				}

			}
			numNs++;
		}
	}
}

void findOverallComponents(grapht& g, std::vector<node>& nodes, int fg, int genIndex) {
	int nV = nodes.size();
	std::vector<bool> visited(nV, false); int overallIndex = 0;

	if (fg == 1) {



		for (int i = 0; i < nodes.size(); i++) {
			//			cout << "node " << i << " of " << nodes.size() << endl;
			if (!nodes[i].valid) { continue; }
			if (i == genIndex) {
				visited[i] = true;
				nodes[i].overallCompIndexFg = overallIndex;
				overallIndex++;
				continue;
			}
			if (((int)nodes[i].inFg) == 1 || ((int)nodes[i].type == FILL)) {
				if (!visited[i]) {
					int timer = 0;
					dfsOverall(g, i, visited, timer, -1, overallIndex, nodes, true);
					overallIndex += 1;
				}
			}
			//cout << "finished node " << i << endl;
		}
		//cout << "after dfs 1" << endl;
		for (int i = 0; i < nodes.size(); i++) {
			if (!nodes[i].valid) { continue; }
			if (((int)nodes[i].inFg) == 1 || ((int)nodes[i].type == FILL)) {
				if (nodes[i].overallCompIndexFg == -1) {
					nodes[i].overallCompIndexFg = overallIndex;
					overallIndex += 1;
				}
			}
		}
		//cout << "overall comp index fg" << endl;
	}
	else {
		//cout << "before dfs 1" << endl;
		for (int i = 0; i < nodes.size(); i++) {
			if (!nodes[i].valid) { continue; }
			if (((int)nodes[i].inFg) == 0 || ((int)nodes[i].type == CUT)) {
				if (!visited[i]) {
					int timer = 0;
					dfsOverall(g, i, visited, timer, -1, overallIndex, nodes, false);
					overallIndex += 1;
				}
			}
		}
		//cout << "after dfs 1" << endl;
		for (int i = 0; i < nodes.size(); i++) {
			if (!nodes[i].valid) { continue; }
			if (i == genIndex) {
				visited[i] = true;
				nodes[i].overallCompIndexBg = overallIndex;
				overallIndex++;
				continue;
			}
			if (((int)nodes[i].inFg) == 0 || ((int)nodes[i].type == CUT)) {
				if (nodes[i].overallCompIndexBg == -1) {
					nodes[i].overallCompIndexBg = overallIndex;
					overallIndex += 1;
				}
			}
		}
	}
	//Find component indices for isolated comps
}

void dfs(Graph& g, int v, std::vector<bool>& visited, int& timer, int p, int& overallIndex, std::vector<node>& nodes) {

	visited[v] = true;
	nodes[v].compIndex = overallIndex;
	nodes[v].tin = timer; nodes[v].low = timer;
	timer += 1;
	int children = 0;
	auto neighbours = adjacent_vertices(v, g);
	for (auto u : make_iterator_range(neighbours)) {
		if (u != p) {
			if (visited[u]) {
				nodes[v].low = min(nodes[v].low, nodes[u].tin);
			}
			else {
				dfs(g, u, visited, timer, v, overallIndex, nodes);
				nodes[v].low = min(nodes[v].low, nodes[u].low);
				if (nodes[u].low >= nodes[v].tin) {
					if (nodes[u].tin > 1) {
						nodes[u].isNew = true;
					}
					if (p != -1) {
						nodes[v].isArticulate = true;
					}
				}
				children += 1;
			}
		}
	}
	if (p == -1 && children > 1) {
		nodes[v].isArticulate = true;
	}
}

void dfs(grapht& g, int v, std::vector<bool>& visited, int& timer, int p, int& overallIndex, std::vector<node>& nodes) {

	visited[v] = true;
	nodes[v].compIndex = overallIndex;
	nodes[v].tin = timer; nodes[v].low = timer;
	timer += 1;
	int children = 0;
	auto neighbours = adjacent_vertices(v, g);
	for (auto u : make_iterator_range(neighbours)) {
		if (u != p) {
			if (visited[u]) {
				nodes[v].low = min(nodes[v].low, nodes[u].tin);
			}
			else {
				dfs(g, u, visited, timer, v, overallIndex, nodes);
				nodes[v].low = min(nodes[v].low, nodes[u].low);
				if (nodes[u].low >= nodes[v].tin) {
					if (nodes[u].tin > 1) {
						nodes[u].isNew = true;
					}
					if (p != -1) {
						nodes[v].isArticulate = true;
					}
				}
				children += 1;
			}
		}
	}
	if (p == -1 && children > 1) {
		nodes[v].isArticulate = true;
	}
}

bool compareByTime(const node& a, const node& b)
{
	return a.tin < b.tin;
}

map<int, map<int, int> > findTermComponents(grapht& g, std::vector<node>& nodes, int type) {
	int nV = nodes.size();
	std::vector<bool> visited(nV, false); int overallIndex = 0;
	std::vector<bool> visited1(nV, false); int overallIndex1 = 0;
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].type) == type) {
			if (!visited[i]) {
				int timer = 0; int timer1 = 0;
				dfs(g, i, visited, timer, -1, overallIndex, nodes);

				overallIndex += 1;
			}
		}
	}
	//Find component indices for isolated comps
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].type) == type) {
			if (nodes[i].compIndex == -1) {
				nodes[i].compIndex = overallIndex;
				overallIndex += 1;
			}
		}
	}
	map<int, map<int, int> > nodeConnectivity;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].isArticulate) {
			auto neighbourItr = adjacent_vertices(i, g);
			std::vector<node> neighbours;
			map<int, int> componentMapping; componentMapping[0] = 0;
			for (auto u : make_iterator_range(neighbourItr)) {
				neighbours.push_back(nodes[u]);
				componentMapping[u] = 0;
			}
			std::sort(neighbours.begin(), neighbours.end(), compareByTime);
			int compIndex = 0;
			map<int, int> nodeToComp;
			for (int j = 0; j < neighbours.size(); j++) {
				if (neighbours[j].low < nodes[i].tin) {
					componentMapping[neighbours[j].tin] = componentMapping[neighbours[j].low]; //Use time of ancestor
				}
				else {
					if (neighbours[j].isNew) { //reporting node
						compIndex += 1;
						componentMapping[neighbours[j].tin] = compIndex;
					}
					else { //Inherit component of nearest left neighbor
						if (j - 1 > 0) {
							componentMapping[neighbours[j].tin] = componentMapping[neighbours[j - 1].tin];
						}
					}
				}
				nodeToComp[neighbours[j].index] = componentMapping[neighbours[j].tin];
			}
			nodeConnectivity[i] = nodeToComp;
		}
	}
	return nodeConnectivity;
}

map<int, map<int, int> > findTermComponents(Graph& g, std::vector<node>& nodes, int type) {
	int nV = nodes.size();
	std::vector<bool> visited(nV, false); int overallIndex = 0;
	std::vector<bool> visited1(nV, false); int overallIndex1 = 0;
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].type) == type) {
			if (!visited[i]) {
				int timer = 0; int timer1 = 0;
				dfs(g, i, visited, timer, -1, overallIndex, nodes);

				overallIndex += 1;
			}
		}
	}
	//Find component indices for isolated comps
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].type) == type) {
			if (nodes[i].compIndex == -1) {
				nodes[i].compIndex = overallIndex;
				overallIndex += 1;
			}
		}
	}
	map<int, map<int, int> > nodeConnectivity;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].isArticulate) {
			auto neighbourItr = adjacent_vertices(i, g);
			std::vector<node> neighbours;
			map<int, int> componentMapping; componentMapping[0] = 0;
			for (auto u : make_iterator_range(neighbourItr)) {
				neighbours.push_back(nodes[u]);
				componentMapping[u] = 0;
			}
			std::sort(neighbours.begin(), neighbours.end(), compareByTime);
			int compIndex = 0;
			map<int, int> nodeToComp;
			for (int j = 0; j < neighbours.size(); j++) {
				if (neighbours[j].low < nodes[i].tin) {
					componentMapping[neighbours[j].tin] = componentMapping[neighbours[j].low]; //Use time of ancestor
				}
				else {
					if (neighbours[j].isNew) { //reporting node
						compIndex += 1;
						componentMapping[neighbours[j].tin] = compIndex;
					}
					else { //Inherit component of nearest left neighbor
						if (j - 1 > 0) {
							componentMapping[neighbours[j].tin] = componentMapping[neighbours[j - 1].tin];
						}
					}
				}
				nodeToComp[neighbours[j].index] = componentMapping[neighbours[j].tin];
			}
			nodeConnectivity[i] = nodeToComp;
		}
	}
	return nodeConnectivity;
}

int getParent(map<int, int>& parentComp, int node) {
	int parent = parentComp[node];
	while (parent != parentComp[parent]) {
		parent = parentComp[parent];

	}
	return parent;

}

bool compareByIntensity(const weightedCoord& a, const weightedCoord& b)
{
	return a.intensityDiff > b.intensityDiff;
}

int minigraphComps(const grapht& miniG, const std::vector<bool>& miniGValid, int numNodes) {
	std::vector<int> nodeToComp(numNodes);
	int n = (int)boost::connected_components(miniG, &nodeToComp[0]);
	int numComps = 0;
	std::vector<bool> isCompIndex(n, false);
	for (int i = 0; i < nodeToComp.size(); i++) {
		// cout << i << " " << nodeToComp[i] << endl;
		if (miniGValid[i]) {
			if (!isCompIndex[nodeToComp[i]]) {
				isCompIndex[nodeToComp[i]] = true;
				numComps += 1;
			}
		}
	}
	return numComps;
}

std::vector< std::vector< std::vector<int> > >  adjFaces6Conn = {
   {{1, 0, 0}, {1, 0, 1}, {0, 0, 1}},
   {{1, 0, 0}, {1, 0, -1}, {0, 0, -1}},
   {{-1, 0, 0}, {-1, 0, -1}, {0, 0, -1}},
   {{-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}},
   {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}},
   {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}},
   {{-1, 0, 0}, {-1, -1, 0}, {0, -1, 0}},
   {{1, 0, 0}, {1, -1, 0}, {0, -1, 0}},
   {{0, 0, 1}, {0, 1, 1}, {0, 1, 0}},
   {{0, 0, -1}, {0, 1, -1}, {0, 1, 0}},
   {{0, 0, -1}, {0, -1, -1}, {0, -1, 0}},
   {{0, 0, 1}, {0, -1, 1}, {0, -1, 0}}
};

std::vector< std::vector< std::vector<int> > > adjCubes6Conn = { {{0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1,
   0}, {1, 1, 1}}, {{0, 0, -1}, {0, 1, -1}, {0, 1, 0}, {1, 0, -1}, {1,
	0, 0}, {1, 1, -1}, {1, 1, 0}}, {{0, -1, 0}, {0, -1, 1}, {0, 0,
   1}, {1, -1, 0}, {1, -1, 1}, {1, 0, 0}, {1, 0,
   1}}, {{0, -1, -1}, {0, -1, 0}, {0, 0, -1}, {1, -1, -1}, {1, -1,
   0}, {1, 0, -1}, {1, 0, 0}}, {{-1, 0, 0}, {-1, 0, 1}, {-1, 1,
   0}, {-1, 1, 1}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}}, {{-1,
   0, -1}, {-1, 0, 0}, {-1, 1, -1}, {-1, 1, 0}, {0, 0, -1}, {0,
   1, -1}, {0, 1, 0}}, {{-1, -1, 0}, {-1, -1, 1}, {-1, 0, 0}, {-1, 0,
   1}, {0, -1, 0}, {0, -1, 1}, {0, 0, 1}}, {{-1, -1, -1}, {-1, -1,
   0}, {-1, 0, -1}, {-1, 0, 0}, {0, -1, -1}, {0, -1, 0}, {0, 0, -1}} };

bool simple3DLabel(const std::vector<uint32_t*>& labels, const std::vector<node>& nodes, int x, int y, int z, int numSlices, int width, int height, const std::vector<unsigned char>& simpleDictionary3D, int labelIndex, bool inFg, int& conn) {
	int minX = max(x - 1, 0); int minY = max(y - 1, 0); int minZ = max(z - 1, 0);
	int maxX = min(x + 2, width); int maxY = min(y + 2, height); int maxZ = min(z + 2, numSlices);
	//Check core simplicity in 26-neighborhood
	std::vector< std::vector< std::vector<int>>> cubeN(3, std::vector<std::vector<int>>(3, std::vector<int>(3, 0)));
	for (int i = minX; i < maxX; i++) {
		for (int j = minY; j < maxY; j++) {
			for (int k = minZ; k < maxZ; k++) {
				if (Label(labels[k][i + j * width]) != unvisited) {

					if (((int)Label(labels[k][i + j * width])) == labelIndex ||
						((int)nodes[((int)Label(labels[k][i + j * width]))].type) == 0 ||
						((int)nodes[((int)Label(labels[k][i + j * width]))].type) == 2
						) {
						cubeN[i - x + 1][j - y + 1][k - z + 1] = 1;
					}
					else {
						cubeN[i - x + 1][j - y + 1][k - z + 1] = 0;
					}
				}
				else {
					cubeN[i - x + 1][j - y + 1][k - z + 1] = 0;
				}
			}
		}
	}

	if (conn == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					cubeN[i][j][k] = 1 - cubeN[i][j][k];
				}
			}
		}
	}
	return ((int)simpleDictionary3D[neighborhoodToIndex(cubeN)]) == 49;
}

bool simple3DLabelBg(const std::vector<uint32_t*>& labels, const std::vector<node>& nodes, int x, int y, int z, int numSlices, int width, int height, const std::vector<unsigned char>& simpleDictionary3D, int labelIndex, bool inFg, int& conn) {
	int minX = max(x - 1, 0); int minY = max(y - 1, 0); int minZ = max(z - 1, 0);
	int maxX = min(x + 2, width); int maxY = min(y + 2, height); int maxZ = min(z + 2, numSlices);
	//Check core simplicity in 26-neighborhood
	//Correct version
	std::vector< std::vector< std::vector<int>>> cubeN(3, std::vector<std::vector<int>>(3, std::vector<int>(3, 0)));
	for (int i = minX; i < maxX; i++) {
		for (int j = minY; j < maxY; j++) {
			for (int k = minZ; k < maxZ; k++) {
				if (Label(labels[k][i + j * width]) != unvisited) {

					if (((int)Label(labels[k][i + j * width])) == labelIndex ||
						((int)nodes[((int)Label(labels[k][i + j * width]))].type) == 1 ||
						((int)nodes[((int)Label(labels[k][i + j * width]))].type) == 3
						) {
						cubeN[i - x + 1][j - y + 1][k - z + 1] = 1;

					}
					else {
						cubeN[i - x + 1][j - y + 1][k - z + 1] = 0;
					}
				}
				else {
					cubeN[i - x + 1][j - y + 1][k - z + 1] = 0;
				}
			}
		}
	}


	if (conn == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					cubeN[i][j][k] = 1 - cubeN[i][j][k];
				}
			}
		}
	}
	return ((int)simpleDictionary3D[neighborhoodToIndex(cubeN)]) == 49;
}

void dfsOverallComps(grapht& g, int v, std::vector<bool>& visited, int& timer, int p, int& overallIndex, std::vector<node>& nodes, bool inFg) {
	visited[v] = true;
	if (inFg) {
		nodes[v].overallCompIndexFg = overallIndex;
	}
	else {
		nodes[v].overallCompIndexBg = overallIndex;
	}
	nodes[v].tin = timer; nodes[v].low = timer;
	timer += 1;
	int children = 0;
	auto neighbours = adjacent_vertices(v, g);
	for (auto u : make_iterator_range(neighbours)) {
		if (u != p) {
			if (visited[u]) {
				nodes[v].low = min(nodes[v].low, nodes[u].tin);
			}
			else {
				dfsOverallComps(g, u, visited, timer, v, overallIndex, nodes, inFg);
				nodes[v].low = min(nodes[v].low, nodes[u].low);
				if (nodes[u].low >= nodes[v].tin) {
					if (nodes[u].tin > 1) {
						nodes[u].isNew = true;
					}
					if (p != -1) {
						nodes[v].isArticulateFg = true;
					}
				}
				children += 1;
			}
		}
	}
	if (p == -1 && children > 1) {
		nodes[v].isArticulateFg = true;
	}

}

void findArticulationPoints(grapht& g, std::vector<node>& nodes, int fg, int genIndex) {
	int nV = nodes.size();
	std::vector<bool> visited(nV, false); int overallIndex = 0;


	if (fg == 1) {

		for (int i = 0; i < nodes.size(); i++) {
			//			cout << "node " << i << " of " << nodes.size() << endl;
			if (!nodes[i].valid) { continue; }
			if (((int)nodes[i].inFg) == 1 || ((int)nodes[i].type == FILL)) {
				if (!visited[i]) {
					int timer = 0;
					dfsOverallComps(g, i, visited, timer, -1, overallIndex, nodes, true);
					overallIndex += 1;
				}
			}
			//cout << "finished node " << i << endl;
		}
		//cout << "after dfs 1" << endl;
		for (int i = 0; i < nodes.size(); i++) {
			if (!nodes[i].valid) { continue; }
			if (((int)nodes[i].inFg) == 1 || ((int)nodes[i].type == FILL)) {
				if (nodes[i].overallCompIndexFg == -1) {
					nodes[i].overallCompIndexFg = overallIndex;
					overallIndex += 1;
				}
			}
		}
	}
	else {
		//cout << "before dfs 1" << endl;
		for (int i = 0; i < nodes.size(); i++) {
			if (!nodes[i].valid) { continue; }
			if (((int)nodes[i].inFg) == 0 || ((int)nodes[i].type == CUT)) {
				if (!visited[i]) {
					int timer = 0;
					dfsOverallComps(g, i, visited, timer, -1, overallIndex, nodes, false);
					overallIndex += 1;
				}
			}
		}
		//cout << "after dfs 1" << endl;
		for (int i = 0; i < nodes.size(); i++) {
			if (!nodes[i].valid) { continue; }
			if (((int)nodes[i].inFg) == 0 || ((int)nodes[i].type == CUT)) {
				if (nodes[i].overallCompIndexBg == -1) {
					nodes[i].overallCompIndexBg = overallIndex;
					overallIndex += 1;
				}
			}
		}
	}
	//Find component indices for isolated comps


}



void getAllCombinations(const std::vector< std::vector<int> >& comboVecs, int index, std::vector<int> combo, std::vector< std::vector<int> >& combinations) {
	if (index >= comboVecs.size()) {
		combinations.push_back(combo);
		return;
	}
	for (int i = 0; i < comboVecs[index].size(); i++) {
		std::vector<int> newCombo = combo;
		newCombo.push_back(comboVecs[index][i]);
		getAllCombinations(comboVecs, index + 1, newCombo, combinations);

	}
}

void setSteinerParams(int timelimit, int bbtime) {
	params.timelimit = timelimit;
	params.memlimit = 10000000000;
	params.cutoff = -1;
	params.cutoffopt = false;
	params.absgap = 0;
	params.bbinfofreq = 100;
	params.branchtype = 0;
	params.nodeselect = 0;
	params.daiterations = 10;
	params.perturbedheur = true;
	params.nodelimit = -1;
	params.daeager = 1.25;
	params.dasat = -1;
	params.daguide = true;
	params.seed = 0;

	params.rootonly = false;
	params.heuronly = false;
	params.initheur = true;
	params.initprep = true;
	params.redrootonly = false;
	params.bigM = false;
	params.semiBigM = false;
	params.eagerroot = false;

	params.heureps = -1;
	params.heurroots = 1;
	params.heurbb = true;
	params.heurbbtime = timelimit;
	params.heursupportG = true;

	// reduction tests
	params.d1 = true;
	params.d2 = true;
	params.ma = true;
	params.ms = true;
	params.ss = true;
	params.lc = true;
	params.nr = true;
	params.boundbased = true;

	int p = 10;
	params.precision = 1;

	for (int i = 0; i < p; i++) {
		params.precision *= 10;
	}
	params.type = "nwstp";
}

Inst loadNWSTP(hyperGraph& hyperG) {
	Inst inst;
	inst.offset = 0; //double w, prize;
	inst.r = -1;
	inst.isInt = true;
	inst.resizeNodes(hyperG.numHypernodes); inst.resizeEdges(2 * hyperG.hyperEdges.size()); inst.t = hyperG.numTerminals; std::vector<weight_t> nw(hyperG.numHypernodes, 0.0);
	int ij = 0;
	for (int i = 0; i < hyperG.hyperEdges.size(); i++) { //Add edges
		inst.newArc(hyperG.hyperEdges[i][0], hyperG.hyperEdges[i][1], ij, ij + 1, (weight_t)0); ij++;
		inst.newArc(hyperG.hyperEdges[i][1], hyperG.hyperEdges[i][0], ij, ij - 1, (weight_t)0); ij++;
	}
	int assignedT = 0;
	for (int i = 0; i < hyperG.numHypernodes; i++) {
		if (hyperG.hyperNodes[i].getType() == CORE || hyperG.hyperNodes[i].getType() == N || hyperG.hyperNodes[i].getType() == MIDT || hyperG.hyperNodes[i].getType() == PINODE) {
			inst.T[i] = true;
			inst.f1[i] = 1;
			inst.p[i] = WMAX;
		}
	}

	for (int i = 0; i < inst.n; i++) {
		nw[i] = (weight_t)hyperG.hyperNodes[i].getWeight();
		if (nw[i] == INFINITY) {
			nw[i] = 1;
		}
		for (int ij : inst.din[i]) {
			inst.c[ij] = inst.c[ij] + nw[i];
		}
	}

	for (int i = 0; i < inst.n; i++) {
		if (inst.T[i]) {
			inst.r = i;
			break;
		}
	}

	inst.offset += nw[inst.r];

	for (int i = 0; i < inst.n; i++) {
		inst.din[i].shrink_to_fit();
		inst.dout[i].shrink_to_fit();
	}

	inst.t = 0;
	for (int i = 0; i < inst.n; i++) {
		if (inst.p[i] > 0 || inst.f1[i]) inst.t++;
		inst.T[i] = (inst.p[i] != 0);

	}
	stats.initial = inst.countInstSize();

	// associates the anti-parallel arc to each arc if it exists (-1 otherwise)
	int antiparallelArcs = 0;
	if (inst.isAsym) {
		for (int ij = 0; ij < inst.m; ij++)
			inst.opposite[ij] = -1;
		for (int ij = 0; ij < inst.m; ij++) {
			if (inst.opposite[ij] != -1) continue;
			const int i = inst.tail[ij];
			const int j = inst.head[ij];
			for (int jk : inst.dout[j]) {
				const int k = inst.head[jk];
				if (k == i) {
					inst.opposite[ij] = jk;
					inst.opposite[jk] = ij;
					break;
				}
			}
		}
		for (int ij = 0; ij < inst.m; ij++) {
			if (inst.opposite[ij] == -1) continue;
			antiparallelArcs++;
		}
	}
	else {
		antiparallelArcs = inst.m;
	}

	// compute the ratio of bidirected edges / arcs that have an antiparallel counterpart
	stats.bidirect = (double)antiparallelArcs / (2 * inst.m - antiparallelArcs);
	if (params.bigM) {
		if (inst.r == -1) {
			inst = inst.createRootedBigMCopy();
		}
		else {
			params.bigM = false;
		}
	}
	return inst;
}

void solveSteinerTree(hyperGraph& hyperG, std::vector<node>& nodes, std::vector< std::vector<int> >& slnEdges, int steinerTime, int bbTime) {
	setSteinerParams(steinerTime, bbTime);
	double bestKnown = -1;
	ProcStatus::setMemLimit(params.memlimit);
	params.timelimit = steinerTime;
	srand(params.seed);
	Timer tLoad(true);
	Inst nwstp = loadNWSTP(hyperG);
	if (params.bigM) {
		if (nwstp.r == -1) {
			nwstp = nwstp.createRootedBigMCopy();
		}
		else {
			params.bigM = false;
		}
	}
	BBTree bbtree(nwstp);
	if (params.initprep) {
		Timer tPrep(true);

		bbtree.initPrep();
		stats.prep = nwstp.countInstSize();
		stats.preptime = tPrep.elapsed().getSeconds();
	}
	if (params.cutoff > 0.0)
		bbtree.setCutUp(params.cutoff);
	if (params.cutoffopt)
		bbtree.setCutUp(bestKnown);

	bbtree.setBestKnown(bestKnown);
	// solve
	if (params.heuronly) {
		bbtree.initHeur();
	}
	else if (params.rootonly) {

		bbtree.initHeur();

		if (params.timelimit >= 0)
			bbtree.setTimeLim(max(0.0, params.timelimit));
		if (params.nodelimit >= 0)
			bbtree.setNodeLim(params.nodelimit);
		bbtree.processRoots();
	}

	else {
		bbtree.initHeur();
		if (params.timelimit >= 0) {
			bbtree.setTimeLim(max(0.0, params.timelimit));
		}
		if (params.nodelimit >= 0)
			bbtree.setNodeLim(params.nodelimit);
		bbtree.solve();
	}
	stats.bbnodes = bbtree.getNnodes();
	stats.isInt = nwstp.isInt;
	stats.isAsym = nwstp.isAsym;

	// timing
	stats.bbtime = bbtree.getTime();
	stats.timeBest = bbtree.getTimeBest();
	stats.roottime = bbtree.getRootTime();
	stats.heurtime = bbtree.getHeurTime();
	stats.heurbbtime = bbtree.getHeurBBTime();

	// bounds
	stats.ub = format(bbtree.getUB(), nwstp);
	stats.lb = format(bbtree.getLB(), nwstp);
	stats.gap = gapP(stats.lb, stats.ub);

	// root
	stats.rootub = format(bbtree.getRootUB(), nwstp);
	stats.rootlb = format(bbtree.getRootLB(), nwstp);
	stats.rootgap = gapP(stats.rootlb, stats.rootub);
	stats.roots = bbtree.getNroots();
	stats.oroots = bbtree.getNrootsOpen();
	stats.proots = bbtree.getNrootsProcessed();

	if (bbtree.getState() == BBTree::State::BB_MEMLIMIT)
		stats.memout = 1;

	// get running time and backmapped solution
	stats.time = Timer::total.elapsed().getSeconds();
	auto S = bbtree.getInc1();
	weight_t ub = S.obj;
	weight_t lb = bbtree.getLB();

	// check if value from bound file matches computed optimum value
	const bool match = (abs(format(ub, nwstp) - bestKnown) < 1e-5);
	// validates solution
	stats.valid = S.validate();
	auto inst = bbtree.getInst1(); auto sln = bbtree.getInc1();
	int slnSize = 0;
	for (int ij = 0; ij < inst.m; ij++) {
		if (sln.arcs[ij]) {
			slnSize += 1;
			slnEdges.push_back({ hyperG.hyperEdges[ij / 2][0] , hyperG.hyperEdges[ij / 2][1] });
		}
	}

}

int64_t findNodeCost(int index, int side, std::vector<node>& nodes, int numNodes, Graph& G, Graph& fgG, Graph& bgG, int h0, int h2, map<std::vector<int>, int>& edgeWt) {
	std::vector<node> nodesTemp = nodes; int64_t c = 1000000000;
	Graph fgGTemp = fgG; Graph bgGTemp = bgG;
	if (side == 1) {


		nodesTemp[index].type = CORE; nodesTemp[index].inFg = 1;
		if (nodes[index].type == FILL) {
			auto neighboursT = adjacent_vertices(index, G);
			for (auto u : make_iterator_range(neighboursT)) {
				if (nodes[u].type == CUT) {
					nodesTemp[u].type = CORE;
					clear_vertex(u, bgGTemp);
				}
			}
		}
		int newH0 = findComponents(fgGTemp, nodesTemp, numNodes, true);
		clear_vertex(index, bgGTemp);
		int newH2 = findComponents(bgGTemp, nodesTemp, numNodes, false);
		int64_t labelCost = -nodes[index].labelCost;
		int deltaE = nodes[index].v - nodes[index].e + nodes[index].f - nodes[index].c;
		if (nodes[index].type == CUT) {
			labelCost = 0;
		}
		int64_t cost = (int64_t)(((newH0 - h0) + (newH2 - h2)) * c + labelCost);
		return cost;
	}
	else {
		nodesTemp[index].type = N;
		nodesTemp[index].inFg = 0;
		if (nodes[index].type == CUT) {
			auto neighboursT = adjacent_vertices(index, G);
			for (auto u : make_iterator_range(neighboursT)) {
				if (nodes[u].type == FILL) {
					nodesTemp[u].type = N;
					clear_vertex(u, fgGTemp);
				}
			}
		}
		clear_vertex(index, fgGTemp);
		int newH0 = findComponents(fgGTemp, nodesTemp, numNodes, true);
		int newH2 = findComponents(bgGTemp, nodesTemp, numNodes, false);
		int64_t labelCost = nodes[index].labelCost;
		int deltaE = nodes[index].v - nodes[index].e + nodes[index].f - nodes[index].c;
		if (nodes[index].type == FILL) {
			labelCost = 0;
		}
		int64_t cost = (int64_t)(((newH0 - h0) + (newH2 - h2)) * c + labelCost);
		return cost;

	}
}

int findLocalNodeToFix(std::vector<node>& nodes, int numNodes, map<int, int>& globalToLocal, std::vector<bool>& inFg, std::vector<bool>& inBg, hyperGraph& hypG,
	std::vector< std::vector<int> >& slnEdges, Graph& G, Graph& fgG, Graph& bgG, int fgComps, int bgComps, int& termSide,
	map<std::vector<int>, int>& edgeWt, std::vector<int>& clusterLocalNodesGlobalIndices
) {
	int violatingNode = -1;
	int64_t lowestCost = WMAX;
	for (int i = 0; i < slnEdges.size(); i++) {
		int v1 = slnEdges[i][0];
		switch (hypG.hyperNodes[v1].getType()) {
		case CORE:
			inFg[globalToLocal[hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[0]].origNode.index]] = true;
			break;
		case N:
			inBg[globalToLocal[hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[0]].origNode.index]] = true;
			break;
		case HYPERNODE:
			if (hypG.hyperNodes[v1].getSide() == 1) {
				for (int j = 0; j < hypG.hyperNodes[v1].doubleSubnodes.size(); j++) {
					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[j]].origNode.index;
					inFg[globalToLocal[globalI]] = true;
					if (inBg[globalToLocal[globalI]]) {
						int64_t cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost1;
							termSide = 1;
						}
						int64_t cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost0;
							termSide = 0;
						}
					}
				}
			}
			else { //Must be background side hypernode
				for (int j = 0; j < hypG.hyperNodes[v1].doubleSubnodes.size(); j++) {

					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[j]].origNode.index;
					inBg[globalToLocal[globalI]] = true;
					if (inFg[globalToLocal[globalI]]) { // fgCompNT, bgCompNT,
						int64_t cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost1;
							termSide = 1;
						}
						int64_t cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost0;
							termSide = 0;
						}
					}
				}
			}
			break;
		default:
			break;
		}
		int v2 = slnEdges[i][1];
		switch (hypG.hyperNodes[v2].getType()) {
		case CORE:
			inFg[globalToLocal[hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[0]].origNode.index]] = true;
			break;
		case N:
			inBg[globalToLocal[hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[0]].origNode.index]] = true;
			break;
		case HYPERNODE:
			if (hypG.hyperNodes[v2].getSide() == 1) {
				for (int j = 0; j < hypG.hyperNodes[v2].doubleSubnodes.size(); j++) {
					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[j]].origNode.index;
					inFg[globalToLocal[globalI]] = true;
					if (inBg[globalToLocal[globalI]]) {
						int64_t cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost1;
							termSide = 1;
						}
						int64_t cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost0;
							termSide = 0;
						}
					}
				}
			}
			else { //Must be background side hypernode
				for (int j = 0; j < hypG.hyperNodes[v2].doubleSubnodes.size(); j++) {
					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[j]].origNode.index;
					inBg[globalToLocal[globalI]] = true;
					if (inFg[globalToLocal[globalI]]) {
						int64_t cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost1;
							termSide = 1;
						}
						int64_t cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalToLocal[globalI];
							lowestCost = cost0;
							termSide = 0;
						}
					}
				}
			}
			break;
		default:
			break;
		}
	}
	return violatingNode;
}

void getGlobalToLocalCFGraphs(Graph& localG, std::vector< std::vector<int> >& clusterLocalEdgesGlobalIndices, std::vector<int>& clusterLocalNodesGlobalIndices, map<int, int>& globalToLocal, std::vector<node>& nodes) {
	for (int i = 0; i < clusterLocalNodesGlobalIndices.size(); i++) {
		add_vertex(localG);
	}

	for (int i = 0; i < clusterLocalEdgesGlobalIndices.size(); i++) {
		int v1 = clusterLocalEdgesGlobalIndices[i][0]; int v2 = clusterLocalEdgesGlobalIndices[i][1];
		if ((nodes[v1].type == CUT || nodes[v1].type == FILL) && (nodes[v2].type == CUT || nodes[v2].type == FILL)) {
			add_edge(globalToLocal[v1], globalToLocal[v2], localG);
		}
	}
}

void updateGraphs(Graph& G, Graph& fgG, Graph& bgG, std::vector<node>& nodes, map<std::vector<int>, int>& edgeWts) {
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((nodes[v1].type == 0) || (nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 0) || (nodes[v2].type == 2) || (nodes[v2].type == 3)) && edgeWts[{v1, v2}] == 1) {
			add_edge(v1, v2, fgG);
		}
		if (((nodes[v1].type == 1) || (nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 1) || (nodes[v2].type == 2) || (nodes[v2].type == 3))) {
			add_edge(v1, v2, bgG);
		}
	}
}

void getSteinerGroupingsForDecision(std::vector<bool>& steinerDecision, Graph& localG, std::vector<int>& clusterLocalNodesGlobalIndices, std::vector<node>& nodes,
	std::vector< std::vector<int> >& fgComponents, std::vector< std::vector<int> >& bgComponents) {
	std::vector<bool> cfInFg(clusterLocalNodesGlobalIndices.size(), false); std::vector<bool> cfInBg(clusterLocalNodesGlobalIndices.size(), false);
	Graph fgClusterGraph = localG; Graph bgClusterGraph = localG;
	for (int d = 0; d < steinerDecision.size(); d++) {
		if (nodes[clusterLocalNodesGlobalIndices[d]].type == CORE || nodes[clusterLocalNodesGlobalIndices[d]].type == N) {
			clear_vertex(d, fgClusterGraph);
			clear_vertex(d, bgClusterGraph);
			continue;
		}
		else {
			if (steinerDecision[d]) {
				cfInFg[d] = true; clear_vertex(d, bgClusterGraph);
			}
			else {
				cfInBg[d] = true; clear_vertex(d, fgClusterGraph);
			}
		}
	}


	std::vector<int> nodeToCompFg(clusterLocalNodesGlobalIndices.size());
	int nFg = (int)boost::connected_components(fgClusterGraph, &nodeToCompFg[0]);
	int numCompsFg = 0;
	std::vector<int> isCompIndexFg(nFg, -1);
	for (int i = 0; i < clusterLocalNodesGlobalIndices.size(); i++) {
		if (cfInFg[i]) {
			if (isCompIndexFg[nodeToCompFg[i]] == -1) {
				isCompIndexFg[nodeToCompFg[i]] = numCompsFg;
				std::vector<int> newComp = { i };
				fgComponents.push_back(newComp);
				numCompsFg += 1;
			}
			else {
				fgComponents[isCompIndexFg[nodeToCompFg[i]]].push_back(i);
			}

		}
	}

	std::vector<int> nodeToCompBg(clusterLocalNodesGlobalIndices.size());
	int nBg = (int)boost::connected_components(bgClusterGraph, &nodeToCompBg[0]);
	int numCompsBg = 0;
	std::vector<int> isCompIndexBg(nBg, -1);
	for (int i = 0; i < clusterLocalNodesGlobalIndices.size(); i++) {
		if (cfInBg[i]) {
			if (isCompIndexBg[nodeToCompBg[i]] == -1) {
				isCompIndexBg[nodeToCompBg[i]] = numCompsBg;
				std::vector<int> newComp = { i };
				bgComponents.push_back(newComp);
				numCompsBg += 1;
			}
			else {
				bgComponents[isCompIndexBg[nodeToCompBg[i]]].push_back(i);
			}

		}
	}
}

void solveLocalGraphs(Graph& G, std::vector<node>& nodes, int numNodes, map< std::vector<int>, int>& edgeWts, int hypernodeSize, int64_t wtSum, int productThresh, tbb::concurrent_vector<hyperNode>& globalHypernodes,
	int localSteinerTime,
	int bbTime, int& maxLocalGraphNodes, int& maxLocalGraphEdges) {
	std::vector<node> newNodes = nodes;
	Graph fgG, bgG;
	std::vector< std::vector< std::vector<int> > > localEdgesGlobalIndex;
	std::vector< std::vector<int> > localNodesGlobalIndex;
	getLocalGraphClusters(G, nodes, numNodes, fgG, bgG, localEdgesGlobalIndex, localNodesGlobalIndex, edgeWts);
	int numClusters = localNodesGlobalIndex.size();
	//These two graphs used to connect terminals together in local graph that must be connected outside graph
	Graph coreG;
	getCoreG(G, nodes, numNodes, coreG, edgeWts);
	Graph nG;
	getNG(G, nodes, numNodes, nG);

	int decisionIndex = 0;
	//Find global components first, so don't need to re-find for each run of pruneLocalImpossibleComplexEdges
	int fgComps = findComponents(fgG, nodes, numNodes, true); int bgComps = findComponents(bgG, nodes, numNodes, false); std::vector<hyperGraph> localHyperGraphs;
	//cout << fgComps << " " << bgComps << endl;
	//for (int i = 0; i < localNodesGlobalIndex.size(); i++) {
	tbb::parallel_for(tbb::blocked_range<int>(0, localNodesGlobalIndex.size()),
		[&](tbb::blocked_range<int> r)
	{
		for (int i = r.begin(); i < r.end(); ++i)
		{
			std::vector<bool> steinerDecisions;
			std::vector< std::vector<int> > clusterLocalEdgesGlobalIndices = localEdgesGlobalIndex[i];
			std::vector<int> clusterLocalNodesGlobalIndices = localNodesGlobalIndex[i];

			std::vector< std::vector<int> > complexEdges = clusterLocalEdgesGlobalIndices;
			map<int, int> localToGlobal; map<int, int> globalToLocal;
			getCompToLocalIndex(clusterLocalNodesGlobalIndices, localToGlobal, globalToLocal, edgeWts);
			//cout << "construct double g " << endl;
			doubleGraph localDoubleG(nodes, clusterLocalNodesGlobalIndices, clusterLocalEdgesGlobalIndices, complexEdges, globalToLocal, edgeWts);
			//cout << "done double g " << endl;
			hyperGraph localHyperG(localDoubleG, hypernodeSize, wtSum);
			//cout << "done hyper g" << endl;
			maxLocalGraphNodes = max((int)localHyperG.numHypernodes, maxLocalGraphNodes);
			maxLocalGraphEdges = max((int)localHyperG.hyperEdges.size(), maxLocalGraphEdges);
			//If two nodes have to be connected outside cluster, create edge in hypergraph
			map< std::vector<int>, bool> edgeExists;
			for (int i = 0; i < localHyperG.hyperEdges.size(); i++) {
				edgeExists[{localHyperG.hyperEdges[i][0], localHyperG.hyperEdges[i][1]}] = true; edgeExists[{localHyperG.hyperEdges[i][1], localHyperG.hyperEdges[i][0]}] = true;
			}
			connectTerminalsExternally(localHyperG, coreG, nG, numNodes, clusterLocalNodesGlobalIndices, clusterLocalEdgesGlobalIndices, edgeExists);
			;
			//Solve hypergraph for every combination of terminals
			std::vector< std::vector<int> > connectedCores = getTermsConnectedOutsideCluster(fgG, CORE, nodes, numNodes, clusterLocalNodesGlobalIndices, localHyperG);
			std::vector< std::vector<int> > connectedNs = getTermsConnectedOutsideCluster(bgG, N, nodes, numNodes, clusterLocalNodesGlobalIndices, localHyperG);
			int numCombinations = 1;
			for (int j = 0; j < connectedCores.size(); j++) { numCombinations *= connectedCores[j].size(); }
			for (int j = 0; j < connectedNs.size(); j++) { numCombinations *= connectedNs[j].size(); }
			if (numCombinations > productThresh) {
				//Defer nodes to global graph
				for (int j = 0; j < localHyperG.hyperNodes.size(); j++) {
					if (localHyperG.hyperNodes[j].getType() == HYPERNODE) {
						std::vector<int> hyperNodeGlobalIndices;
						for (int s = 0; s < localHyperG.hyperNodes[j].doubleSubnodes.size(); s++) {
							hyperNodeGlobalIndices.push_back(localDoubleG.doubleNodes[localHyperG.hyperNodes[j].doubleSubnodes[s]].origNode.index);
						}
						std::sort(hyperNodeGlobalIndices.begin(), hyperNodeGlobalIndices.end());
						globalHypernodes.push_back(hyperNode(hyperNodeGlobalIndices, HYPERNODE, localHyperG.hyperNodes[j].getSide()));
					}
				}
			}
			else {
				std::vector< std::vector<std::vector<std::vector<int> > > > cCompPartitions;
				for (int j = 0; j < connectedCores.size(); j++) {
					std::vector<std::vector<std::vector<int> > > corePartitions = getPartitions(connectedCores[j]);
					cCompPartitions.push_back(corePartitions);
				}
				std::vector< std::vector<std::vector<std::vector<int> > > > nCompPartitions;
				for (int j = 0; j < connectedNs.size(); j++) {
					std::vector<std::vector<std::vector<int> > > nPartitions = getPartitions(connectedNs[j]);
					nCompPartitions.push_back(nPartitions);
				}
				std::vector< std::vector<int> > cCompIndices;
				for (int j = 0; j < cCompPartitions.size(); j++) {
					std::vector<int> compPartIndices;
					for (int k = 0; k < cCompPartitions[j].size(); k++) {
						compPartIndices.push_back(k);
					}
					cCompIndices.push_back(compPartIndices);
				}
				std::vector< std::vector<int> > nCompIndices;
				for (int j = 0; j < nCompPartitions.size(); j++) {
					std::vector<int> compPartIndices;
					for (int k = 0; k < nCompPartitions[j].size(); k++) {
						compPartIndices.push_back(k);
					}
					nCompIndices.push_back(compPartIndices);
				}

				std::vector<int> coreCombo; std::vector< std::vector<int> > coreCombinations;
				getAllCombinations(cCompIndices, 0, coreCombo, coreCombinations);
				std::vector<int> nCombo; std::vector< std::vector<int> > nCombinations;
				getAllCombinations(nCompIndices, 0, nCombo, nCombinations);
				int totalCombinations = coreCombinations.size() * nCombinations.size(); std::vector< std::vector<bool> > steinerDecisions;
				std::vector<int> combIndices; for (int index = 0; index < totalCombinations; index++) { combIndices.push_back(index); }
				for (int itr = 0; itr < combIndices.size(); itr++) {
					//tbb::parallel_for(tbb::blocked_range<int>(0, combIndices.size()),
						//[&](tbb::blocked_range<int> range)
						//{
						//for (int itr = range.begin(); itr < range.end(); ++itr)
							//{
					int c = itr / ((int)nCombinations.size());
					int n = itr % ((int)nCombinations.size());
					std::vector<int> cPartitionSelection = coreCombinations[c]; map< std::vector<int>, bool > edgeExistsTemp = edgeExists; //map< std::vector<int>, int > edgeWtsTemp = edgeWts;
					hyperGraph partitionHyperG = localHyperG; std::vector<node> nodesTemp = nodes;
					int fgCompsTemp = fgComps;
					int bgCompsTemp = bgComps;
					Graph tempG = G;
					Graph fgGTemp = fgG;
					Graph bgGTemp = bgG;
					std::vector< std::vector<int> > termEdges;
					std::vector<int> connTermIndices;
					for (int p = 0; p < cPartitionSelection.size(); p++) {
						std::vector<std::vector<int> > componentPartition = cCompPartitions[p][cPartitionSelection[p]];
						for (int g = 0; g < componentPartition.size(); g++) {
							std::vector<int> group = componentPartition[g];
							for (int gi = 0; gi < group.size() - 1; gi++) {
								if (edgeExistsTemp.find({ group[gi], group[gi + 1] }) == edgeExistsTemp.end()) {
									partitionHyperG.hyperEdges.push_back({ group[gi], group[gi + 1] }); termEdges.push_back({ group[gi], group[gi + 1] }); connTermIndices.push_back(group[gi]); connTermIndices.push_back(group[gi + 1]);
									edgeExistsTemp[{ group[gi + 1], group[gi] }] = true; edgeExistsTemp[{ group[gi], group[gi + 1] }] = true;

								}
							}
						}
					}
					std::vector<int> nPartitionSelection = nCombinations[n];
					for (int p = 0; p < nPartitionSelection.size(); p++) {
						std::vector<std::vector<int> > componentPartition = nCompPartitions[p][nPartitionSelection[p]];
						for (int g = 0; g < componentPartition.size(); g++) {
							std::vector<int> group = componentPartition[g];
							for (int gi = 0; gi < group.size() - 1; gi++) {
								if (edgeExistsTemp.find({ group[gi], group[gi + 1] }) == edgeExistsTemp.end()) {
									partitionHyperG.hyperEdges.push_back({ group[gi], group[gi + 1] }); termEdges.push_back({ group[gi], group[gi + 1] }); termEdges.push_back({ group[gi], group[gi + 1] }); connTermIndices.push_back(group[gi]); connTermIndices.push_back(group[gi + 1]);
									edgeExistsTemp[{ group[gi + 1], group[gi] }] = true; edgeExistsTemp[{ group[gi], group[gi + 1] }] = true;
								}
							}
						}
					}

					bool violating = true;
					while (violating) {
						std::vector< std::vector<int> > slnEdges; hyperGraph currentHyperG = partitionHyperG;
						//cout << "before steiner solve " << endl;
						solveSteinerTree(currentHyperG, nodesTemp, slnEdges, localSteinerTime, bbTime);
						//cout << "after solve " << endl;
						violating = false;
						std::vector<bool> inFg(clusterLocalNodesGlobalIndices.size(), false);
						std::vector<bool> inBg(clusterLocalNodesGlobalIndices.size(), false);
						int termSide;
						int newTerminal = findLocalNodeToFix(nodesTemp, numNodes, globalToLocal, inFg, inBg, currentHyperG, slnEdges, G, fgGTemp, bgGTemp, fgCompsTemp, bgCompsTemp, termSide, edgeWts, clusterLocalNodesGlobalIndices);

						if (newTerminal != -1) {
							newTerminal = clusterLocalNodesGlobalIndices[newTerminal];
							violating = true;
							if (termSide == 1) { //Terminal set to the foreground
								if (nodesTemp[newTerminal].type == FILL) { //If type fill sent to FG, adjacent cuts must be sent to FG as well
									auto neighboursT = adjacent_vertices(newTerminal, G);
									for (auto u : make_iterator_range(neighboursT)) {
										if (nodesTemp[u].type == CUT) {

											nodesTemp[u].type = CORE;
										}
									}
								}
								nodesTemp[newTerminal].type = CORE; nodesTemp[newTerminal].inFg = 1;
							}
							else {
								if (nodesTemp[newTerminal].type == CUT) { //If type cut sent to BG, adjacent fills must be sent to BG as well
									auto neighboursT = adjacent_vertices(newTerminal, G);
									for (auto u : make_iterator_range(neighboursT)) {
										if (nodesTemp[u].type == FILL) {
											nodesTemp[u].type = N;
										}
									}
								}
								nodesTemp[newTerminal].type = N; nodesTemp[newTerminal].inFg = 0;
							}
							removeCAndNEdges(tempG, nodesTemp);
							updateGraphs(tempG, fgGTemp, bgGTemp, nodesTemp, edgeWts);
							fgCompsTemp = findComponents(fgGTemp, nodesTemp, numNodes, true);
							int bgCompsTemp = findComponents(bgGTemp, nodesTemp, numNodes, false);
							doubleGraph localDoubleGTemp(nodesTemp, clusterLocalNodesGlobalIndices, clusterLocalEdgesGlobalIndices, complexEdges, globalToLocal, edgeWts);
							partitionHyperG = hyperGraph(localDoubleGTemp, hypernodeSize, wtSum);
							//Add edges between terminals using indices of updated partition graph
							std::sort(connTermIndices.begin(), connTermIndices.end());
							connTermIndices.erase(unique(connTermIndices.begin(), connTermIndices.end()), connTermIndices.end());
							map<int, int> connTermIndicesMapping;
							for (int j = 0; j < connTermIndices.size(); j++) {
								connTermIndicesMapping[localHyperG.doubleG.doubleNodes[localHyperG.hyperNodes[connTermIndices[j]].doubleSubnodes[0]].origNode.index] = -1;
							}
							for (int j = 0; j < partitionHyperG.coreIndices.size(); j++) {
								int localI = localDoubleGTemp.doubleNodes[partitionHyperG.hyperNodes[partitionHyperG.coreIndices[j]].doubleSubnodes[0]].origNode.index;
								if (connTermIndicesMapping.find(localI) != connTermIndicesMapping.end()) {
									connTermIndicesMapping[localI] = partitionHyperG.coreIndices[j];
								}
							}
							for (int j = 0; j < partitionHyperG.nIndices.size(); j++) {
								int localI = localDoubleGTemp.doubleNodes[partitionHyperG.hyperNodes[partitionHyperG.nIndices[j]].doubleSubnodes[0]].origNode.index;
								if (connTermIndicesMapping.find(localI) != connTermIndicesMapping.end()) {
									connTermIndicesMapping[localI] = partitionHyperG.nIndices[j];
								}
							}
							for (int j = 0; j < termEdges.size(); j++) {
								partitionHyperG.hyperEdges.push_back({
								connTermIndicesMapping[localHyperG.doubleG.doubleNodes[localHyperG.hyperNodes[termEdges[j][0]].doubleSubnodes[0]].origNode.index],
									connTermIndicesMapping[localHyperG.doubleG.doubleNodes[localHyperG.hyperNodes[termEdges[j][1]].doubleSubnodes[0]].origNode.index]
									});
							}

						}
						else {
							steinerDecisions.push_back(inFg);
						}
					}
					//}

				//});
				}
				//if a node always makes the same decisions across all steiner trees, assign it to fg or bg
				int newTerms = 0; int aFg = 0; int aBg = 0;
				for (int j = 0; j < clusterLocalNodesGlobalIndices.size(); j++) {

					if (nodes[clusterLocalNodesGlobalIndices[j]].type == CORE || nodes[clusterLocalNodesGlobalIndices[j]].type == N) { continue; };
					bool allFg = true;
					bool allBg = true;
					for (int k = 0; k < steinerDecisions.size(); k++) {
						if (steinerDecisions[k][j]) {
							allBg = false;
						}
						else {
							allFg = false;
						}
					}

					if (allFg && allBg) { cout << "all fg and all bg: bug! " << endl; }
					else {
						if (allFg) { //Assign to be core
							aFg += 1;
							if (nodes[clusterLocalNodesGlobalIndices[j]].type == FILL) { //If type fill sent to FG, adjacent cuts must be sent to FG as well
								auto neighboursT = adjacent_vertices(clusterLocalNodesGlobalIndices[j], G);
								for (auto u : make_iterator_range(neighboursT)) {
									if (newNodes[u].type == CUT) {
										newNodes[u].type = CORE;
										newNodes[u].inFg = 1;
									}
								}
							}
							newNodes[clusterLocalNodesGlobalIndices[j]].type = CORE;
							newNodes[clusterLocalNodesGlobalIndices[j]].inFg = 1;
						}
						else {
							if (allBg) { //Assign to be N
								aBg += 1;
								if (nodes[clusterLocalNodesGlobalIndices[j]].type == CUT) { //If type fill sent to FG, adjacent cuts must be sent to FG as well
									auto neighboursT = adjacent_vertices(clusterLocalNodesGlobalIndices[j], G);
									for (auto u : make_iterator_range(neighboursT)) {
										if (newNodes[u].type == FILL) {
											newNodes[u].type = N;
											newNodes[u].inFg = 0;
										}
									}
								}
								newNodes[clusterLocalNodesGlobalIndices[j]].type = N;
								newNodes[clusterLocalNodesGlobalIndices[j]].inFg = 0;
							}
						}
					}
				}
				Graph localG;
				getGlobalToLocalCFGraphs(localG, complexEdges, clusterLocalNodesGlobalIndices, globalToLocal, newNodes);
				//Cluster nodes which make same decisions together
				int pushed = 0;
				map< std::vector<int>, bool> hypAdded;
				for (int j = 0; j < steinerDecisions.size(); j++) {
					decisionIndex += 1;
					std::vector<bool> hasPos(clusterLocalNodesGlobalIndices.size(), false);
					std::vector<bool> hasNeg(clusterLocalNodesGlobalIndices.size(), false);

					std::vector< std::vector<int> > cfFgComps;
					std::vector< std::vector<int> > cfBgComps;
					getSteinerGroupingsForDecision(steinerDecisions[j], localG, clusterLocalNodesGlobalIndices, newNodes, cfFgComps, cfBgComps);
					for (int cfComp = 0; cfComp < cfFgComps.size(); cfComp++) {
						std::vector<int> fgCFCompLocalIndices = cfFgComps[cfComp]; std::vector<int> fgCFComp;
						bool hasIndex = false;
						for (int k = 0; k < fgCFCompLocalIndices.size(); k++) {
							fgCFComp.push_back(clusterLocalNodesGlobalIndices[fgCFCompLocalIndices[k]]);
						}


						std::sort(fgCFComp.begin(), fgCFComp.end());
						fgCFComp.erase(unique(fgCFComp.begin(), fgCFComp.end()), fgCFComp.end());
						std::vector<int> nodeWithSign = fgCFComp;
						nodeWithSign.push_back(1);
						for (int k = 0; k < fgCFComp.size(); k++) {
							hasPos[globalToLocal[fgCFComp[k]]] = true;
						}

						if (hypAdded.find(nodeWithSign) == hypAdded.end()) {
							globalHypernodes.push_back(hyperNode(fgCFComp, HYPERNODE, 1));
							pushed += 1;
							hypAdded[nodeWithSign] = true;
						}
					}
					for (int cfComp = 0; cfComp < cfBgComps.size(); cfComp++) {
						std::vector<int> bgCFCompLocalIndices = cfBgComps[cfComp];
						std::vector<int> bgCFComp;
						bool hasIndex = false;
						for (int k = 0; k < bgCFCompLocalIndices.size(); k++) {
							bgCFComp.push_back(clusterLocalNodesGlobalIndices[bgCFCompLocalIndices[k]]);
						}
						std::sort(bgCFComp.begin(), bgCFComp.end());
						bgCFComp.erase(unique(bgCFComp.begin(), bgCFComp.end()), bgCFComp.end());
						for (int k = 0; k < bgCFComp.size(); k++) {
							hasNeg[globalToLocal[bgCFComp[k]]] = true;
						}
						std::vector<int> nodeWithSign = bgCFComp;
						nodeWithSign.push_back(-1);
						if (hypAdded.find(nodeWithSign) == hypAdded.end()) {
							globalHypernodes.push_back(hyperNode(bgCFComp, HYPERNODE, -1));
							hypAdded[nodeWithSign] = true;
							pushed += 1;
						}
					}
					for (int k = 0; k < clusterLocalNodesGlobalIndices.size(); k++) {
						if (newNodes[clusterLocalNodesGlobalIndices[k]].type == CUT || newNodes[clusterLocalNodesGlobalIndices[k]].type == FILL) {
							if (!hasPos[k] && !hasNeg[k]) {
								cout << "doesn't have positive or negative node " << k << " " << clusterLocalNodesGlobalIndices[k] << endl;
							}

						}
					}
				}

			}
		}
	});
	//}
	nodes = newNodes;

	//cout << "done local solve " << endl;
}

std::vector< std::vector<int>> getTypeComponents(grapht& g, std::vector<node>& nodes, int nodesSize, int type) {
	std::vector< std::vector<int> > typeComponents;
	std::vector<int> nodeToComp(nodesSize);
	int n = (int)boost::connected_components(g, &nodeToComp[0]);
	int numComps = 0; std::vector<int> isCompIndex(n, -1);
	for (int i = 0; i < nodesSize; i++) {

		if (((int)nodes[i].type) == type) {
			if (isCompIndex[nodeToComp[i]] == -1) {

				isCompIndex[nodeToComp[i]] = numComps;
				std::vector<int> newComp = { i };
				typeComponents.push_back(newComp);
				numComps += 1;
			}
			else {
				typeComponents[isCompIndex[nodeToComp[i]]].push_back(i);
			}

		}
	}
	return typeComponents;
}

std::vector<std::vector<int> > mergeAdjacentTerminals(Graph& G, std::vector<node>& nodes, map< std::vector<int>, int>& edgeWt, map<int, int>& oldToNew) {
	std::vector<node> newNodes;
	grapht coreG; grapht neighborhoodG;
	int nodesSize = nodes.size();
	for (int i = 0; i < nodesSize; i++) {
		add_vertex(coreG); add_vertex(neighborhoodG);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((int)nodes[v1].type == 0) && nodes[v1].valid && ((int)nodes[v2].type == 0) && nodes[v2].valid && edgeWt[{v1, v2}] == 1) {
			add_edge(v1, v2, 0, coreG);
		}
		if (((int)nodes[v1].type == 1) && nodes[v1].valid && ((int)nodes[v2].type == 1) && nodes[v2].valid) {
			add_edge(v1, v2, 0, neighborhoodG);
		}
	}
	std::vector< std::vector<int>> coreComponents = getTypeComponents(coreG, nodes, nodesSize, 0); std::vector< std::vector<int>> neighborhoodComponents = getTypeComponents(neighborhoodG, nodes, nodesSize, 1);
	std::vector<std::vector<int> > newToOldComp;
	for (int i = 0; i < coreComponents.size(); i++) {
		newToOldComp.push_back(coreComponents[i]);
		node n; n.type = 0; n.index = newNodes.size(); n.inFg = 1; newNodes.push_back(n);
		for (int j = 0; j < coreComponents[i].size(); j++) {
			oldToNew[coreComponents[i][j]] = newNodes.size() - 1;
		}
	}
	for (int i = 0; i < neighborhoodComponents.size(); i++) {
		newToOldComp.push_back(neighborhoodComponents[i]);
		node n; n.type = 1; n.index = newNodes.size(); n.inFg = 0; newNodes.push_back(n);
		for (int j = 0; j < neighborhoodComponents[i].size(); j++) {
			oldToNew[neighborhoodComponents[i][j]] = newNodes.size() - 1;
		}
	}
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].type) == 2) {
			newToOldComp.push_back({ i });
			node n = nodes[i]; n.index = newNodes.size();  n.type = 2; n.inFg = 1; n.labelCost = nodes[i].labelCost; n.intensity = nodes[i].intensity;
			newNodes.push_back(n);
			oldToNew[i] = newNodes.size() - 1;
		}
		else {
			if (((int)nodes[i].type) == 3) {
				newToOldComp.push_back({ i });
				//node n = nodes[i];
				node n = nodes[i]; n.index = newNodes.size();
				n.type = 3; n.inFg = 0; n.index = newNodes.size(); n.labelCost = nodes[i].labelCost; n.intensity = nodes[i].intensity;
				newNodes.push_back(n);
				oldToNew[i] = newNodes.size() - 1;
			}
		}
	}
	map< std::vector<int>, int > newEdgeWt;
	Graph newG;
	for (int i = 0; i < newNodes.size(); i++) {
		add_vertex(newG);
	}
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (oldToNew[v1] == oldToNew[v2]) { continue; }
		if (newEdgeWt.find({ oldToNew[v1], oldToNew[v2] }) == newEdgeWt.end()) {
			add_edge(oldToNew[v1], oldToNew[v2], newG);
			newEdgeWt[{ oldToNew[v1], oldToNew[v2] }] = edgeWt[{v1, v2}]; newEdgeWt[{ oldToNew[v2], oldToNew[v1] }] = edgeWt[{v1, v2}];
		}

		else {
			if (newEdgeWt[{ oldToNew[v1], oldToNew[v2] }] == 0) {
				if (edgeWt[{v1, v2}] == 1) {
					newEdgeWt[{ oldToNew[v1], oldToNew[v2] }] = edgeWt[{v1, v2}]; newEdgeWt[{ oldToNew[v2], oldToNew[v1] }] = edgeWt[{v1, v2}];
				}
			}
		}
	}
	G = newG;
	edgeWt = newEdgeWt;
	nodes = newNodes;
	return newToOldComp;
}

void getFgBgGraphs(Graph& G, std::vector<node>& nodes, int numNodes, Graph& fgG, Graph& bgG, map< std::vector<int>, int>& edgeWts) {
	Graph cutFillGraph;
	for (int i = 0; i < numNodes; i++) {
		add_vertex(fgG); add_vertex(bgG);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((nodes[v1].type == 0) || (nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 0) || (nodes[v2].type == 2) || (nodes[v2].type == 3)) && edgeWts[{v1, v2}] == 1) {
			add_edge(v1, v2, fgG);
		}
		if (((nodes[v1].type == 1) || (nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 1) || (nodes[v2].type == 2) || (nodes[v2].type == 3))) {
			add_edge(v1, v2, bgG);
		}

	}
}

void getLocalGraphClustersMapping(Graph& G, std::vector<node>& nodes, int numNodes, std::vector< std::vector< std::vector<int> > >& localGraphEdges, std::vector< std::vector<int> >& localNodesGlobalIndex, map<int, int>& nodeToComp) {
	Graph cutFillGraph;
	for (int i = 0; i < numNodes; i++) {
		add_vertex(cutFillGraph);
	}
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((nodes[v1].type == 2) || (nodes[v1].type == 3)) && ((nodes[v2].type == 2) || (nodes[v2].type == 3))) {
			add_edge(v1, v2, cutFillGraph);
		}

	}
	std::vector< std::vector<int>> cfComponents = getCfComponents(cutFillGraph, nodes, numNodes);
	std::vector< std::vector<int> > globalToLocalNodes(numNodes, std::vector<int>(0));

	for (int i = 0; i < cfComponents.size(); i++) {
		std::vector<int> cfWithTerminals = cfComponents[i];
		map<int, bool> termAdded;
		for (int j = 0; j < cfComponents[i].size(); j++) {
			nodeToComp[cfComponents[i][j]] = i;
			globalToLocalNodes[cfComponents[i][j]].push_back(i); //Map each node to the cluster which contains it

			auto neighbours = adjacent_vertices(cfComponents[i][j], G);
			for (auto vd : make_iterator_range(neighbours)) {
				if (((int)nodes[vd].type) == 0 || ((int)nodes[vd].type) == 1) {
					if (termAdded.find(vd) == termAdded.end()) {
						cfWithTerminals.push_back(vd);
						termAdded[vd] = true;
					}
				}
			}
		}
		localNodesGlobalIndex.push_back(cfWithTerminals);
		std::vector< std::vector<int> > localEdges; localGraphEdges.push_back(localEdges);
	}

	//Iterate through edges, find subgraph for each local cut-fill cluster in one sweep
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (((int)nodes[v1].type == 2) || ((int)nodes[v1].type == 3)) {
			localGraphEdges[globalToLocalNodes[v1][0]].push_back({ v1, v2 });
			continue;
		}
		else {
			if (((int)nodes[v2].type == 2) || ((int)nodes[v2].type == 3)) {
				localGraphEdges[globalToLocalNodes[v2][0]].push_back({ v1, v2 });
				continue;
			}

		}
	}
}

int findGlobalNodeToFix(std::vector<node>& nodes, int numNodes, std::vector<bool>& inFg, std::vector<bool>& inBg, hyperGraph& hypG,
	std::vector< std::vector<int> >& slnEdges, Graph& G, Graph& fgG, Graph& bgG, int fgComps, int bgComps, int& termSide,
	map<std::vector<int>, int>& edgeWt, std::vector<bool>& nodeToFix, int64_t& cost
) {
	int violatingNode = -1; int64_t lowestCost = WMAX;
	for (int i = 0; i < slnEdges.size(); i++) {
		int v1 = slnEdges[i][0];
		switch (hypG.hyperNodes[v1].getType()) {
		case CORE:
			inFg[hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[0]].origNode.index] = true;
			break;
		case N:
			inBg[hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[0]].origNode.index] = true;
			break;
		case HYPERNODE:
			if (hypG.hyperNodes[v1].getSide() == 1) {
				for (int j = 0; j < hypG.hyperNodes[v1].doubleSubnodes.size(); j++) {
					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[j]].origNode.index;
					inFg[globalI] = true;
					if (inBg[globalI]) {
						if (nodeToFix[globalI]) { continue; }
						double cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalI; lowestCost = cost1; termSide = 1;
						}
						double cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalI; lowestCost = cost0; termSide = 0;
						}
					}
				}
			}
			else { //Must be background side hypernode
				for (int j = 0; j < hypG.hyperNodes[v1].doubleSubnodes.size(); j++) {

					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v1].doubleSubnodes[j]].origNode.index;
					inBg[globalI] = true;
					if (inFg[globalI]) {
						if (nodeToFix[globalI]) { continue; }
						double cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalI; lowestCost = cost1; termSide = 1;
						}
						double cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalI; lowestCost = cost0; termSide = 0;
						}
					}
				}
			}
			break;
		default:
			break;
		}
		int v2 = slnEdges[i][1];
		switch (hypG.hyperNodes[v2].getType()) {
		case CORE:
			inFg[hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[0]].origNode.index] = true;
			break;
		case N:
			inBg[hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[0]].origNode.index] = true;
			break;
		case HYPERNODE:
			if (hypG.hyperNodes[v2].getSide() == 1) {
				for (int j = 0; j < hypG.hyperNodes[v2].doubleSubnodes.size(); j++) {
					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[j]].origNode.index;
					inFg[globalI] = true;
					if (inBg[globalI]) {
						if (nodeToFix[globalI]) { continue; }
						double cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalI; lowestCost = cost1; termSide = 1;
						}
						double cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalI; lowestCost = cost0; termSide = 0;
						}
					}
				}
			}
			else { //Must be background side hypernode
				for (int j = 0; j < hypG.hyperNodes[v2].doubleSubnodes.size(); j++) {
					int globalI = hypG.doubleG.doubleNodes[hypG.hyperNodes[v2].doubleSubnodes[j]].origNode.index;
					inBg[globalI] = true;
					if (inFg[globalI]) {
						if (nodeToFix[globalI]) { continue; }
						double cost1 = findNodeCost(globalI, 1, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost1 < lowestCost) {
							violatingNode = globalI; lowestCost = cost1; termSide = 1;
						}
						double cost0 = findNodeCost(globalI, 0, nodes, numNodes, G, fgG, bgG, fgComps, bgComps, edgeWt);
						if (cost0 < lowestCost) {
							violatingNode = globalI; lowestCost = cost0; termSide = 0;
						}
					}
				}
			}
			break;
		default:
			break;
		}
	}
	cost = lowestCost;
	return violatingNode;
}

void solveComponentGraph(Graph& G, std::vector<node>& nodes, std::vector<node>& newNodes, int numNodes, map< std::vector<int>, int>& edgeWts, int hypernodeSize, int64_t wtSum, int productThresh, tbb::concurrent_vector< hyperNode >& globalHypernodes,
	std::vector< std::vector<int> >& clusterLocalEdgesGlobalIndices, std::vector<int>& clusterLocalNodesGlobalIndices, int fgComps, int bgComps, int localSteinerTime, Graph& fgG, Graph& bgG, int compIndex, map<int, int>& nodeToComp
	, int bbTime) {
	newNodes = nodes;

	//These two graphs used to connect terminals together in local graph that must be connected outside graph
	Graph coreG; getCoreG(G, nodes, numNodes, coreG, edgeWts);
	Graph nG; getNG(G, nodes, numNodes, nG);

	std::vector<bool> steinerDecisions;
	std::vector< std::vector<int> > complexEdges = clusterLocalEdgesGlobalIndices;
	map<int, int> localToGlobal; map<int, int> globalToLocal;
	getCompToLocalIndex(clusterLocalNodesGlobalIndices, localToGlobal, globalToLocal, edgeWts);

	doubleGraph localDoubleG(nodes, clusterLocalNodesGlobalIndices, clusterLocalEdgesGlobalIndices, complexEdges, globalToLocal, edgeWts);

	hyperGraph localHyperG(localDoubleG, hypernodeSize, wtSum);
	map< std::vector<int>, bool> edgeExists;
	for (int i = 0; i < localHyperG.hyperEdges.size(); i++) {
		edgeExists[{localHyperG.hyperEdges[i][0], localHyperG.hyperEdges[i][1]}] = true;
		edgeExists[{localHyperG.hyperEdges[i][1], localHyperG.hyperEdges[i][0]}] = true;
	}
	connectTerminalsExternally(localHyperG, coreG, nG, numNodes, clusterLocalNodesGlobalIndices, clusterLocalEdgesGlobalIndices, edgeExists);
	//Solve hypergraph for every combination of terminals

	std::vector< std::vector<int> > connectedCores = getTermsConnectedOutsideCluster(fgG, CORE, nodes, numNodes, clusterLocalNodesGlobalIndices, localHyperG, nodeToComp);
	std::vector< std::vector<int> > connectedNs = getTermsConnectedOutsideCluster(bgG, N, nodes, numNodes, clusterLocalNodesGlobalIndices, localHyperG, nodeToComp);
	int numCombinations = 1;
	for (int j = 0; j < connectedCores.size(); j++) { numCombinations *= connectedCores[j].size(); } for (int j = 0; j < connectedNs.size(); j++) { numCombinations *= connectedNs[j].size(); }
	if (numCombinations > productThresh) {
		for (int j = 0; j < localHyperG.hyperNodes.size(); j++) {
			if (localHyperG.hyperNodes[j].getType() == HYPERNODE) {
				std::vector<int> hyperNodeGlobalIndices;
				for (int s = 0; s < localHyperG.hyperNodes[j].doubleSubnodes.size(); s++) {
					hyperNodeGlobalIndices.push_back(localDoubleG.doubleNodes[localHyperG.hyperNodes[j].doubleSubnodes[s]].origNode.index);

				}
				globalHypernodes.push_back(hyperNode(hyperNodeGlobalIndices, HYPERNODE, localHyperG.hyperNodes[j].getSide()));
			}
		}
	}
	else {
		std::vector< std::vector<std::vector<std::vector<int> > > > cCompPartitions;
		for (int j = 0; j < connectedCores.size(); j++) {
			std::vector<std::vector<std::vector<int> > > corePartitions = getPartitions(connectedCores[j]);
			cCompPartitions.push_back(corePartitions);
		}
		std::vector< std::vector<std::vector<std::vector<int> > > > nCompPartitions;
		for (int j = 0; j < connectedNs.size(); j++) {
			std::vector<std::vector<std::vector<int> > > nPartitions = getPartitions(connectedNs[j]);
			nCompPartitions.push_back(nPartitions);
		}
		std::vector< std::vector<int> > cCompIndices;
		for (int j = 0; j < cCompPartitions.size(); j++) {
			std::vector<int> compPartIndices;
			for (int k = 0; k < cCompPartitions[j].size(); k++) {
				compPartIndices.push_back(k);
			}
			cCompIndices.push_back(compPartIndices);
		}
		std::vector< std::vector<int> > nCompIndices;
		for (int j = 0; j < nCompPartitions.size(); j++) {
			std::vector<int> compPartIndices;
			for (int k = 0; k < nCompPartitions[j].size(); k++) {
				compPartIndices.push_back(k);
			}
			nCompIndices.push_back(compPartIndices);
		}

		std::vector<int> coreCombo; std::vector< std::vector<int> > coreCombinations;
		getAllCombinations(cCompIndices, 0, coreCombo, coreCombinations);
		std::vector<int> nCombo; std::vector< std::vector<int> > nCombinations;
		getAllCombinations(nCompIndices, 0, nCombo, nCombinations);

		int totalCombinations = coreCombinations.size() * nCombinations.size(); std::vector< std::vector<bool> > steinerDecisions;
		std::vector<int> combIndices; for (int index = 0; index < totalCombinations; index++) { combIndices.push_back(index); }
		tbb::parallel_for(tbb::blocked_range<int>(0, combIndices.size()),
			[&](tbb::blocked_range<int> range)
		{
			for (int itr = range.begin(); itr < range.end(); ++itr)
			{
				int c = itr / ((int)nCombinations.size());
				int n = itr % ((int)nCombinations.size());
				std::vector<int> cPartitionSelection = coreCombinations[c];
				map< std::vector<int>, bool > edgeExistsTemp = edgeExists;
				hyperGraph partitionHyperG = localHyperG;
				std::vector<node> nodesTemp = nodes;
				int fgCompsTemp = fgComps;
				int bgCompsTemp = bgComps;
				Graph tempG = G;
				Graph fgGTemp = fgG;
				Graph bgGTemp = bgG;
				std::vector< std::vector<int> > termEdges;
				std::vector<int> connTermIndices;
				for (int p = 0; p < cPartitionSelection.size(); p++) {
					std::vector<std::vector<int> > componentPartition = cCompPartitions[p][cPartitionSelection[p]];
					for (int g = 0; g < componentPartition.size(); g++) {
						std::vector<int> group = componentPartition[g];
						for (int gi = 0; gi < group.size() - 1; gi++) {
							if (edgeExistsTemp.find({ group[gi], group[gi + 1] }) == edgeExistsTemp.end()) {
								partitionHyperG.hyperEdges.push_back({ group[gi], group[gi + 1] }); termEdges.push_back({ group[gi], group[gi + 1] }); connTermIndices.push_back(group[gi]); connTermIndices.push_back(group[gi + 1]);
								edgeExistsTemp[{ group[gi + 1], group[gi] }] = true; edgeExistsTemp[{ group[gi], group[gi + 1] }] = true;
							}
						}
					}
				}
				std::vector<int> nPartitionSelection = nCombinations[n];
				for (int p = 0; p < nPartitionSelection.size(); p++) {
					std::vector<std::vector<int> > componentPartition = nCompPartitions[p][nPartitionSelection[p]];
					for (int g = 0; g < componentPartition.size(); g++) {
						std::vector<int> group = componentPartition[g];
						for (int gi = 0; gi < group.size() - 1; gi++) {
							if (edgeExistsTemp.find({ group[gi], group[gi + 1] }) == edgeExistsTemp.end()) {
								partitionHyperG.hyperEdges.push_back({ group[gi], group[gi + 1] }); termEdges.push_back({ group[gi], group[gi + 1] }); termEdges.push_back({ group[gi], group[gi + 1] }); connTermIndices.push_back(group[gi]); connTermIndices.push_back(group[gi + 1]);
								edgeExistsTemp[{ group[gi + 1], group[gi] }] = true; edgeExistsTemp[{ group[gi], group[gi + 1] }] = true;
							}
						}
					}
				}
				bool violating = true;
				while (violating) {
					std::vector< std::vector<int> > slnEdges;
					hyperGraph currentHyperG = partitionHyperG;
					solveSteinerTree(currentHyperG, nodesTemp, slnEdges, localSteinerTime, bbTime);

					violating = false; std::vector<bool> inFg(clusterLocalNodesGlobalIndices.size(), false);
					std::vector<bool> inBg(clusterLocalNodesGlobalIndices.size(), false); int termSide;
					int newTerminal = findLocalNodeToFix(nodesTemp, numNodes, globalToLocal, inFg, inBg, currentHyperG, slnEdges, G, fgGTemp, bgGTemp, fgCompsTemp, bgCompsTemp, termSide, edgeWts, clusterLocalNodesGlobalIndices);

					if (newTerminal != -1) {
						newTerminal = clusterLocalNodesGlobalIndices[newTerminal];
						violating = true;
						if (termSide == 1) { //Terminal set to the foreground

							if (nodesTemp[newTerminal].type == FILL) { //If type fill sent to FG, adjacent cuts must be sent to FG as well
								auto neighboursT = adjacent_vertices(newTerminal, G);
								for (auto u : make_iterator_range(neighboursT)) {
									if (nodesTemp[u].type == CUT) {
										nodesTemp[u].type = CORE;
									}
								}
							}
							nodesTemp[newTerminal].type = CORE; nodesTemp[newTerminal].inFg = 1;
						}
						else {

							if (nodesTemp[newTerminal].type == CUT) { //If type cut sent to BG, adjacent fills must be sent to BG as well
								auto neighboursT = adjacent_vertices(newTerminal, G);
								for (auto u : make_iterator_range(neighboursT)) {
									if (nodesTemp[u].type == FILL) {
										nodesTemp[u].type = N;
									}
								}
							}
							nodesTemp[newTerminal].type = N; nodesTemp[newTerminal].inFg = 0;
						}
						removeCAndNEdges(tempG, nodesTemp); updateGraphs(tempG, fgGTemp, bgGTemp, nodesTemp, edgeWts);
						fgCompsTemp = findComponents(fgGTemp, nodesTemp, numNodes, true);
						int bgCompsTemp = findComponents(bgGTemp, nodesTemp, numNodes, false);
						doubleGraph localDoubleGTemp(nodesTemp, clusterLocalNodesGlobalIndices, clusterLocalEdgesGlobalIndices, complexEdges, globalToLocal, edgeWts);

						partitionHyperG = hyperGraph(localDoubleGTemp, hypernodeSize, wtSum);
						//Add edges between terminals using indices of updated partition graph
						std::sort(connTermIndices.begin(), connTermIndices.end());
						connTermIndices.erase(unique(connTermIndices.begin(), connTermIndices.end()), connTermIndices.end()); map<int, int> connTermIndicesMapping;
						for (int j = 0; j < connTermIndices.size(); j++) {
							connTermIndicesMapping[localHyperG.doubleG.doubleNodes[localHyperG.hyperNodes[connTermIndices[j]].doubleSubnodes[0]].origNode.index] = -1;
						}
						for (int j = 0; j < partitionHyperG.coreIndices.size(); j++) {
							int localI = localDoubleGTemp.doubleNodes[partitionHyperG.hyperNodes[partitionHyperG.coreIndices[j]].doubleSubnodes[0]].origNode.index;
							if (connTermIndicesMapping.find(localI) != connTermIndicesMapping.end()) {
								connTermIndicesMapping[localI] = partitionHyperG.coreIndices[j];
							}
						}
						for (int j = 0; j < partitionHyperG.nIndices.size(); j++) {
							int localI = localDoubleGTemp.doubleNodes[partitionHyperG.hyperNodes[partitionHyperG.nIndices[j]].doubleSubnodes[0]].origNode.index;
							if (connTermIndicesMapping.find(localI) != connTermIndicesMapping.end()) {
								connTermIndicesMapping[localI] = partitionHyperG.nIndices[j];
							}
						}
						for (int j = 0; j < termEdges.size(); j++) {
							partitionHyperG.hyperEdges.push_back({
							connTermIndicesMapping[localHyperG.doubleG.doubleNodes[localHyperG.hyperNodes[termEdges[j][0]].doubleSubnodes[0]].origNode.index],
								connTermIndicesMapping[localHyperG.doubleG.doubleNodes[localHyperG.hyperNodes[termEdges[j][1]].doubleSubnodes[0]].origNode.index]
								});
						}
					}
					else {
						steinerDecisions.push_back(inFg);
					}
				}
			}
		});
		//if a node always makes the same decisions across all steiner trees, assign it to fg or bg
		int newTerms = 0; int aFg = 0; int aBg = 0;
		for (int j = 0; j < clusterLocalNodesGlobalIndices.size(); j++) {
			if (nodes[clusterLocalNodesGlobalIndices[j]].type == CORE || nodes[clusterLocalNodesGlobalIndices[j]].type == N) { continue; };
			bool allFg = true;
			bool allBg = true;
			for (int k = 0; k < steinerDecisions.size(); k++) {
				if (steinerDecisions[k][j]) { allBg = false; }
				else { allFg = false; }
			}

			if (allFg && allBg) { cout << "all fg and all bg: bug! " << endl; }
			if (allFg) { //Assign to be core
				aFg += 1;
				if (nodes[clusterLocalNodesGlobalIndices[j]].type == FILL) { //If type fill sent to FG, adjacent cuts must be sent to FG as well
					auto neighboursT = adjacent_vertices(clusterLocalNodesGlobalIndices[j], G);
					for (auto u : make_iterator_range(neighboursT)) {
						if (newNodes[u].type == CUT) {
							newNodes[u].type = CORE; newNodes[u].inFg = 1;
							newTerms += 1;
						}
					}
				}
				newTerms += 1;
				newNodes[clusterLocalNodesGlobalIndices[j]].type = CORE; newNodes[clusterLocalNodesGlobalIndices[j]].inFg = 1;
			}
			if (allBg) { //Assign to be N
				aBg += 1;
				if (nodes[clusterLocalNodesGlobalIndices[j]].type == CUT) { //If type fill sent to FG, adjacent cuts must be sent to FG as well
					auto neighboursT = adjacent_vertices(clusterLocalNodesGlobalIndices[j], G);
					for (auto u : make_iterator_range(neighboursT)) {
						if (newNodes[u].type == FILL) {
							newTerms += 1;
							newNodes[u].type = N; newNodes[u].inFg = 0;
						}
					}
				}
				newTerms += 1;
				newNodes[clusterLocalNodesGlobalIndices[j]].type = N; newNodes[clusterLocalNodesGlobalIndices[j]].inFg = 0;
			}
		}
		Graph localG;
		getGlobalToLocalCFGraphs(localG, complexEdges, clusterLocalNodesGlobalIndices, globalToLocal, newNodes);
		//Cluster nodes which make same decisions together


		map< std::vector<int>, bool> hypAdded;
		for (int j = 0; j < steinerDecisions.size(); j++) {
			std::vector< std::vector<int> > cfFgComps; std::vector< std::vector<int> > cfBgComps;
			getSteinerGroupingsForDecision(steinerDecisions[j], localG, clusterLocalNodesGlobalIndices, newNodes, cfFgComps, cfBgComps);
			for (int cfComp = 0; cfComp < cfFgComps.size(); cfComp++) {
				std::vector<int> fgCFCompLocalIndices = cfFgComps[cfComp]; std::vector<int> fgCFComp;
				for (int k = 0; k < fgCFCompLocalIndices.size(); k++) {
					fgCFComp.push_back(clusterLocalNodesGlobalIndices[fgCFCompLocalIndices[k]]);

				}

				std::sort(fgCFComp.begin(), fgCFComp.end());
				fgCFComp.erase(unique(fgCFComp.begin(), fgCFComp.end()), fgCFComp.end());
				std::vector<int> nodeWithSign = fgCFComp; nodeWithSign.push_back(1);
				if (hypAdded.find(nodeWithSign) == hypAdded.end()) {
					globalHypernodes.push_back(hyperNode(fgCFComp, HYPERNODE, 1));
					hypAdded[nodeWithSign] = true;
				}
			}
			for (int cfComp = 0; cfComp < cfBgComps.size(); cfComp++) {
				std::vector<int> bgCFCompLocalIndices = cfBgComps[cfComp];
				std::vector<int> bgCFComp;
				for (int k = 0; k < bgCFCompLocalIndices.size(); k++) {
					bgCFComp.push_back(clusterLocalNodesGlobalIndices[bgCFCompLocalIndices[k]]);

				}
				std::sort(bgCFComp.begin(), bgCFComp.end());
				bgCFComp.erase(unique(bgCFComp.begin(), bgCFComp.end()), bgCFComp.end());
				std::vector<int> nodeWithSign = bgCFComp;
				nodeWithSign.push_back(-1);
				if (hypAdded.find(nodeWithSign) == hypAdded.end()) {
					globalHypernodes.push_back(hyperNode(bgCFComp, HYPERNODE, -1));
					hypAdded[nodeWithSign] = true;
				}
			}
		}
	}
	nodes = newNodes;
}

// std::vector< std::vector<int> >& newToOldComp,
void solveGlobalGraph(std::vector<node>& nodes, int numNodes, Graph& G, Graph& origG, tbb::concurrent_vector< hyperNode >& globalHypernodes, int64_t wtSum, map< std::vector<int>, int>& edgeWts,
	int hypernodeSize, int productThresh, int globalSteinerTime, int localSteinerTime,
	int nodesToFix, int bbTime) {

	doubleGraph globalDoubleG(nodes, numNodes, G, edgeWts);
	hyperGraph hyperG(nodes, globalDoubleG, edgeWts, G, globalHypernodes, wtSum);
	hyperGraph currentHyperG = hyperG; Graph fgG; Graph bgG;

	getFgBgGraphs(G, nodes, numNodes, fgG, bgG, edgeWts);
	std::vector< std::vector< std::vector<int> > > localGraphEdges;
	std::vector< std::vector<int> > localNodesGlobalIndex; map<int, int> nodeToComp;
	getLocalGraphClustersMapping(G, nodes, numNodes, localGraphEdges, localNodesGlobalIndex, nodeToComp);
	int fgComps = findComponents(fgG, nodes, numNodes, true);
	int bgComps = findComponents(bgG, nodes, numNodes, false);
	bool conflicting = true; //Iterate until no more nodes which are doubly labeled
	int numItr = 0;
	//cout << "Hypergraph size " << currentHyperG.hyperNodes.size() << " " << currentHyperG.hyperEdges.size() << endl;
	while (conflicting) {
		//cout << "Global steiner tree iteration  " << numItr + 1 << " Graph nodes: " << nodes.size() << " Double graph nodes: " << currentHyperG.doubleG.doubleNodes.size() << " Double graph edges: " << currentHyperG.doubleG.doubleEdges.size() << " Hypergraph nodes: " << currentHyperG.numHypernodes << " Hypergraph edges: " << currentHyperG.hyperEdges.size() << endl;
		std::vector< std::vector<int> > slnEdges;

		solveSteinerTree(currentHyperG, nodes, slnEdges, globalSteinerTime, bbTime);

		conflicting = false; std::vector<bool> inFg(numNodes, false);
		std::vector<bool> inBg(numNodes, false); int termSide;
		std::vector<bool> nodeFixed(numNodes, false);
		int64_t termCost;
		int newTerminal = findGlobalNodeToFix(nodes, numNodes, inFg, inBg, currentHyperG, slnEdges, G, fgG, bgG, fgComps, bgComps, termSide, edgeWts, nodeFixed, termCost);

		if (newTerminal != -1) {
			nodeFixed[newTerminal] = true;
			conflicting = true;
			std::vector<bool> isAffectedComp(localNodesGlobalIndex.size(), false); std::vector<int> affectedComponents; //After a new terminal is upon a conflicting node, these variables represents the components of cuts and fills which are affected
			//Set new terminal and potentially adjacent nodes to be terminals to satisfy cell complex constraint 
			if (termSide == 1) {

				if (((int)nodes[newTerminal].type) == FILL) {
					auto neighbours = adjacent_vertices(newTerminal, G);
					for (auto u : make_iterator_range(neighbours)) {
						if (((int)nodes[u].type) == CUT) {
							nodes[u].type = CORE; isAffectedComp[nodeToComp[u]] = true; affectedComponents.push_back(nodeToComp[u]); nodeFixed[u] = true;
						}
					}
				}
				nodes[newTerminal].type = CORE; isAffectedComp[nodeToComp[newTerminal]] = true; affectedComponents.push_back(nodeToComp[newTerminal]);
			}
			else {
				if (((int)nodes[newTerminal].type) == CUT) {
					auto neighbours = adjacent_vertices(newTerminal, G);
					for (auto u : make_iterator_range(neighbours)) {
						if (((int)nodes[u].type) == FILL) {
							nodes[u].type = N; isAffectedComp[nodeToComp[u]] = true; affectedComponents.push_back(nodeToComp[u]);  nodeFixed[u] = true;
						}
					}
				}
				nodes[newTerminal].type = N; isAffectedComp[nodeToComp[newTerminal]] = true; affectedComponents.push_back(nodeToComp[newTerminal]);
			}

			for (int n = 1; n < nodesToFix; n++) {
				int nextT = findGlobalNodeToFix(nodes, numNodes, inFg, inBg, currentHyperG, slnEdges, G, fgG, bgG, fgComps, bgComps, termSide, edgeWts, nodeFixed, termCost);
				if (nextT != -1 && termCost < 0) {
					if (termSide == 1) {

						if (((int)nodes[nextT].type) == FILL) {
							auto neighbours = adjacent_vertices(nextT, G);
							for (auto u : make_iterator_range(neighbours)) {
								if (((int)nodes[u].type) == CUT) {
									nodes[u].type = CORE; isAffectedComp[nodeToComp[u]] = true; affectedComponents.push_back(nodeToComp[u]);  nodeFixed[u] = true;
								}
							}
						}
						nodes[nextT].type = CORE; isAffectedComp[nodeToComp[nextT]] = true; affectedComponents.push_back(nodeToComp[nextT]); nodeFixed[nextT] = true;
					}
					else {
						if (((int)nodes[nextT].type) == CUT) {
							auto neighbours = adjacent_vertices(nextT, G);
							for (auto u : make_iterator_range(neighbours)) {
								if (((int)nodes[u].type) == FILL) {
									nodes[u].type = N; isAffectedComp[nodeToComp[u]] = true; affectedComponents.push_back(nodeToComp[u]);  nodeFixed[u] = true;
								}
							}
						}
						nodes[nextT].type = N; isAffectedComp[nodeToComp[nextT]] = true; affectedComponents.push_back(nodeToComp[nextT]);  nodeFixed[nextT] = true;
					}
				}
				else {
					break;
				}
			}
			//Preprocess graph globally

			preprocessGraph(G, nodes, edgeWts, isAffectedComp, affectedComponents, nodeToComp, true);
			tbb::concurrent_vector< hyperNode > newHypernodes; //Also need to make newHypernodeSides thread safe
																										  //Find hypernodes for graph in next iteration 
			for (int i = 0; i < currentHyperG.hyperNodes.size(); i++) {
				if (currentHyperG.hyperNodes[i].getType() == HYPERNODE) {
					//For all components not affected by new terminals, propagate hypernodes to next iteration
					if (!isAffectedComp[nodeToComp[currentHyperG.doubleG.doubleNodes[currentHyperG.hyperNodes[i].doubleSubnodes[0]].origNode.index]]) {
						std::vector<int> hypernode;
						for (int j = 0; j < currentHyperG.hyperNodes[i].doubleSubnodes.size(); j++) {
							hypernode.push_back(currentHyperG.doubleG.doubleNodes[currentHyperG.hyperNodes[i].doubleSubnodes[j]].origNode.index);
						}
						newHypernodes.push_back(hyperNode(hypernode, HYPERNODE, currentHyperG.hyperNodes[i].getSide()));
					}
				}
			}
			std::sort(affectedComponents.begin(), affectedComponents.end()); affectedComponents.erase(unique(affectedComponents.begin(), affectedComponents.end()), affectedComponents.end()); //Find unique components which are affected by new terminals
			fgG = Graph();
			bgG = Graph();
			getFgBgGraphs(G, nodes, numNodes, fgG, bgG, edgeWts);
			fgComps = findComponents(fgG, nodes, numNodes, true); bgComps = findComponents(bgG, nodes, numNodes, false);
			for (int i = 0; i < affectedComponents.size(); i++) {
				int compIndex = affectedComponents[i];
				//If all nodes in component are now terminals, continue
				bool hasCF = false;
				for (int j = 0; j < localNodesGlobalIndex[compIndex].size(); j++) {
					if ((int)nodes[localNodesGlobalIndex[compIndex][j]].type == CUT || (int)nodes[localNodesGlobalIndex[compIndex][j]].type == FILL) {
						hasCF = true;
						break;
					}
				}
				if (hasCF) {
					std::vector<node> newNodes = nodes;
					solveComponentGraph(G, nodes, newNodes, numNodes, edgeWts, hypernodeSize, wtSum, productThresh, newHypernodes,
						localGraphEdges[compIndex], localNodesGlobalIndex[compIndex], fgComps, bgComps, localSteinerTime, fgG, bgG, compIndex, nodeToComp, bbTime);
				}
			}


			removeCAndNEdges(G, nodes);
			/**map<int, int> oldToNew2;
			std::vector< std::vector<int> > newToOldComp2 = mergeAdjacentTerminals(G, nodes, edgeWts, oldToNew2);
			std::vector< std::vector<int> > newToOldTemp;
			numNodes = nodes.size();
			for (int i = 0; i < newToOldComp2.size(); i++) {

				std::vector<int> combinedNewToOld;
				for (int j = 0; j < newToOldComp2[i].size(); j++) {
					int oldIndex = newToOldComp2[i][j];
					for (int k = 0; k < newToOldComp[oldIndex].size(); k++) {
						combinedNewToOld.push_back(newToOldComp[oldIndex][k]);
					}
				}
				newToOldTemp.push_back(combinedNewToOld);
			}
			newToOldComp = newToOldTemp;**/
			tbb::concurrent_vector<hyperNode> globalNodesTemp;
			for (int i = 0; i < newHypernodes.size(); i++) {
				std::vector<int> subnodes;
				for (int j = 0; j < newHypernodes[i].doubleSubnodes.size(); j++) {
					//subnodes.push_back(oldToNew2[newHypernodes[i].doubleSubnodes[j]]);
					subnodes.push_back(newHypernodes[i].doubleSubnodes[j]);
				}
				std::sort(subnodes.begin(), subnodes.end());
				subnodes.erase(unique(subnodes.begin(), subnodes.end()), subnodes.end());

				globalNodesTemp.push_back(hyperNode(subnodes, HYPERNODE, newHypernodes[i].getSide()));
			}
			newHypernodes = globalNodesTemp;
			removeCAndNEdges(G, nodes); numNodes = nodes.size();
			doubleGraph globalDoubleGN(nodes, numNodes, G, edgeWts);
			hyperGraph hyperGN(nodes, globalDoubleGN, edgeWts, G, newHypernodes, wtSum);
			currentHyperG = hyperGN;
			fgG = Graph();
			bgG = Graph();
			getFgBgGraphs(G, nodes, numNodes, fgG, bgG, edgeWts);
			localGraphEdges.clear();
			localNodesGlobalIndex.clear();
			nodeToComp.clear();
			getLocalGraphClustersMapping(G, nodes, numNodes, localGraphEdges, localNodesGlobalIndex, nodeToComp);
			fgComps = findComponents(fgG, nodes, numNodes, true);
			bgComps = findComponents(bgG, nodes, numNodes, false);
		}
		else {
			//No more conflicting nodes
			for (int i = 0; i < nodes.size(); i++) {
				nodes[i].inFg = inFg[i];
			}
		}
		numItr += 1;
	}
	//cout << "Finished global steiner tree stage " << endl;
}

void findContradictions(State& state, vector< vector<int> >& overallToLvlNodes, Graph& intersectionGraph, vector<vector<node> >& origLvlNodes) {

	for (int i = 1; i < state.levels.size(); i++) {
		std::vector<node> nodes = state.levels[i].origNodes; // state.levels[i].nodes;

		for (int j = 0; j < nodes.size(); j++) {

			//cout << "lvl " << i << " " << j << " " << (int)nodes[j].inFg << " " << (int)origLvlNodes[i][j].type << " " << nodes.size() << endl;
			//&& ((int)origLvlNodes[i][j].type == FILL || (int)origLvlNodes[i][j].type == CUT)
			if ((int)nodes[j].inFg == 1 && nodes[j].valid)
			{
				auto neighbours = adjacent_vertices(nodes[j].totalNodeIndex, intersectionGraph);
				for (auto u : make_iterator_range(neighbours)) {

					int nLvl = overallToLvlNodes[(int)u][0];
					int nIndex = overallToLvlNodes[(int)u][1];
					if (nLvl < i) { //cut or fill at lower level is in BG
						if ((int)state.levels[nLvl].origNodes[nIndex].inFg == 0) {
							// &&
							//((int)origLvlNodes[nLvl][nIndex].type == CUT
							//	|| (int)origLvlNodes[nLvl][nIndex].type == FILL)
							Contradiction c;
							c.l1 = nLvl;
							c.n1 = nIndex;
							c.l2 = i;
							c.n2 = j;
							//find upper and lower contradictions for this node
							//find lower
							vector<bool> visited(overallToLvlNodes.size(), false);
							queue<int> q;
							q.push((int)u); //us is lower node
							visited[(int)u] = true;
							visited[nodes[j].totalNodeIndex] = true;
							while (!q.empty()) { //search in lower levels
								int top = q.front();
								int cLvl = overallToLvlNodes[(int)top][0];
								int cIndex = overallToLvlNodes[(int)top][1];
								if ((int)state.levels[cLvl].origNodes[cIndex].inFg == 0) {
									c.lowerIndices.push_back((int)top); //only push if current label 0
								}
								q.pop();
								auto neighboursLower = adjacent_vertices(top, intersectionGraph);
								for (auto a : make_iterator_range(neighboursLower)) {
									int nLvlLower = overallToLvlNodes[(int)a][0];
									int nIndexLower = overallToLvlNodes[(int)a][1];
									if (nLvlLower < cLvl) { //cut or fill at lower level is in BG
										if (((int)origLvlNodes[nLvlLower][nIndexLower].type == CUT
											|| (int)origLvlNodes[nLvlLower][nIndexLower].type == FILL)) {
											if (!visited[(int)a]) {
												q.push((int)a);
												visited[(int)a] = true;
											}
										}
									}
								}
							}

							visited = vector<bool>(overallToLvlNodes.size(), false);
							q.push(nodes[j].totalNodeIndex); //search in upper levels
							visited[(int)u] = true;
							visited[nodes[j].totalNodeIndex] = true;
							while (!q.empty()) { //search in upper levels
								int top = q.front();
								int cLvl = overallToLvlNodes[(int)top][0];
								int cIndex = overallToLvlNodes[(int)top][1];
								//cout << cLvl << " " << cIndex << " " << (int)state.levels[cLvl].nodes[cIndex].type << " " << (int)state.levels[cLvl].nodes[cIndex].inFg << endl;
								if ((int)state.levels[cLvl].origNodes[cIndex].inFg == 1) {
									c.upperIndices.push_back((int)top); //only push if current label 1
								}
								q.pop();
								auto neighboursUpper = adjacent_vertices(top, intersectionGraph);
								for (auto a : make_iterator_range(neighboursUpper)) {
									int nLvlUpper = overallToLvlNodes[(int)a][0];
									int nIndexUpper = overallToLvlNodes[(int)a][1];
									if (nLvlUpper > cLvl) { //cut or fill at lower level is in BG
										if (((int)origLvlNodes[nLvlUpper][nIndexUpper].type == CUT
											|| (int)origLvlNodes[nLvlUpper][nIndexUpper].type == FILL)) {
											if (!visited[(int)a]) {
												q.push((int)a);
												visited[(int)a] = true;
											}
										}
									}
								}
							}

							state.contradictions.push_back(c);
						}
					}
				}
			}
		}
		//cout << "finished level " << endl;
	}
	//cout << "end fn " << endl;
}

struct CompareState {
	bool operator()(State const& s1, State const& s2) {
		return s1.cost > s2.cost;
	}
};

void solveLevel(State& state, int l, vector<Graph>& levelGraphsIn, vector<map<vector<int>, int>>& levelEdgeWtsIn,
	int hypernodeSize, int productThresh, int localSteinerTime, int bbTime, int globalSteinerTime, vector<node>& origNodes
) {
	//Assign labeling for isolated cut-fill components and critical articulation points
	std::vector<bool> fillerB; std::vector<int> fillerI; map<int, int> fillerMap; //These variables are filler variables here; this function takes more arguments during the global steiner tree stage.
	state.levels[l].nodes = state.levels[l].origNodes;
	map<int, int> oldToNew;
	vector<Graph> levelGraphs = levelGraphsIn;
	vector<map<vector<int>, int> > levelEdgeWts = levelEdgeWtsIn;
	//cout << "before preproces " << endl;
	preprocessGraph(levelGraphs[l], state.levels[l].nodes, levelEdgeWts[l], fillerB, fillerI, fillerMap, false);
	std::vector< std::vector<int> > newToOldComp = mergeAdjacentTerminals(levelGraphs[l], state.levels[l].nodes, levelEdgeWts[l], oldToNew);
	tbb::concurrent_vector< hyperNode > globalHypernodes;
	//Represents nodes after local stage: some have been assigned to core and neighborhood
	//Local steiner tree stage on clusters of cuts and fills
	//cout << "Original graph size: Nodes: " << nodes.size() << ", Edges: " << num_edges(G) << endl;
	int maxLocalGraphNodes = 0; int maxLocalGraphEdges = 0;
	//cout << "steiner time " << globalSteinerTime << " " << localSteinerTime << " " << bbTime << endl;
	//cout << "before solve local " << endl;
	solveLocalGraphs(levelGraphs[l], state.levels[l].nodes, state.levels[l].nodes.size(), levelEdgeWts[l], hypernodeSize,
		state.levels[l].wtSum, productThresh, globalHypernodes, localSteinerTime, bbTime,
		maxLocalGraphNodes, maxLocalGraphEdges
	);
	//cout << "Max size of local hypergraph: Nodes: " << maxLocalGraphNodes << " Edges: " << maxLocalGraphEdges <<
		//" wtSum " << level.wtSum <<
		//endl;

	removeCAndNEdges(levelGraphs[l], state.levels[l].nodes); map<int, int> oldToNew2;

	std::vector< std::vector<int> > newToOldComp2 = mergeAdjacentTerminals(levelGraphs[l], state.levels[l].nodes, levelEdgeWts[l], oldToNew2);
	std::vector< std::vector<int> > newToOldTemp;
	for (int i = 0; i < newToOldComp2.size(); i++) {
		std::vector<int> combinedNewToOld;
		for (int j = 0; j < newToOldComp2[i].size(); j++) {
			int oldIndex = newToOldComp2[i][j];
			for (int k = 0; k < newToOldComp[oldIndex].size(); k++) {
				combinedNewToOld.push_back(newToOldComp[oldIndex][k]);
			}
		}
		newToOldTemp.push_back(combinedNewToOld);
	}

	newToOldComp = newToOldTemp;
	//level.newToOldComp = newToOldComp;
	state.levelNewToOldComps[l] = newToOldComp;

	for (int i = 0; i < globalHypernodes.size(); i++) {
		std::vector<int> subnodes;
		for (int j = 0; j < globalHypernodes[i].doubleSubnodes.size(); j++) {
			subnodes.push_back(oldToNew2[globalHypernodes[i].doubleSubnodes[j]]);
		}
		std::sort(subnodes.begin(), subnodes.end());
		subnodes.erase(unique(subnodes.begin(), subnodes.end()), subnodes.end());
		globalHypernodes[i] = hyperNode(subnodes, HYPERNODE, globalHypernodes[i].getSide());
	}

	//Global steiner tree stage
	int nodesToFix = 1; int numNodes = state.levels[l].nodes.size();
	Graph origGraph = levelGraphs[l];
	//cout << "before solve global " << endl;
	solveGlobalGraph(state.levels[l].nodes, numNodes, levelGraphs[l], origGraph, globalHypernodes, state.levels[l].wtSum, levelEdgeWts[l], hypernodeSize, productThresh, globalSteinerTime, localSteinerTime,
		nodesToFix, bbTime); //newToOldComp,

	for (int i = 0; i < state.levels[l].nodes.size(); i++) {
		vector<int> subnodes = state.levelNewToOldComps[l][i]; //; level.newToOldComp[i];
		for (int j = 0; j < subnodes.size(); j++) {
			state.levels[l].origNodes[subnodes[j]].inFg = state.levels[l].nodes[i].inFg;
		}
	}
	vector<int> currentEulerNums = state.levelEulerNums[l];

	double geometryCost = 0;
	int numFills = 0; int numCuts = 0;
	int oFills = 0; int oCuts = 0;
	for (int j = 0; j < origNodes.size(); j++) {
		if ((int)origNodes[j].type == FILL) {
			oFills++;
		}
		if ((int)origNodes[j].type == CUT) {
			oCuts++;
		}
		if ((int)state.levels[l].origNodes[j].type == FILL) {
			numFills++;
		}
		if ((int)state.levels[l].origNodes[j].type == CUT) {
			numCuts++;
		}
		/**if (((int)origNodes[j].type == FILL && (int)state.levels[l].origNodes[j].inFg == 1) ||
			((int)origNodes[j].type == CUT && (int)state.levels[l].origNodes[j].inFg == 0)) {
			geometryCost += abs(state.levels[l].origNodes[j].intensity);
		}**/

		if (((int)origNodes[j].type == FILL && (int)state.levels[l].origNodes[j].inFg == 1)) {
			currentEulerNums[0] += origNodes[j].v;
			currentEulerNums[1] += origNodes[j].e;
			currentEulerNums[2] += origNodes[j].f;
			currentEulerNums[3] += origNodes[j].c;
			geometryCost += abs(state.levels[l].origNodes[j].floatCost);
		}

		if (((int)origNodes[j].type == CUT && (int)state.levels[l].origNodes[j].inFg == 0)) {
			currentEulerNums[0] -= origNodes[j].v;
			currentEulerNums[1] -= origNodes[j].e;
			currentEulerNums[2] -= origNodes[j].f;
			currentEulerNums[3] -= origNodes[j].c;
			geometryCost += abs(state.levels[l].origNodes[j].floatCost);
		}
	}

	Graph fgG = Graph();
	Graph bgG = Graph();
	for (int i = 0; i < state.levels[l].nodes.size(); i++) {
		add_vertex(fgG);
		add_vertex(bgG);
	}//
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(levelGraphs[l]); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if ((int)state.levels[l].nodes[v1].inFg == 1 && (int)state.levels[l].nodes[v2].inFg == 1) {
			if (levelEdgeWts[l][{v1, v2}] == 1) {
				add_edge(v1, v2, fgG);
			}
		}
		/**if ((int)state.levels[l].nodes[v1].inFg == 1) {
			if ((int)state.levels[l].nodes[v1].type == FILL) {
				if ((int)state.levels[l].nodes[v2].type == CUT) {
					if ((int)state.levels[l].nodes[v2].inFg == 0) {
						cout << " contradiction" << endl;
					}
				}
			}
		}
		if ((int)state.levels[l].nodes[v2].inFg == 1) {
			if ((int)state.levels[l].nodes[v2].type == FILL) {
				if ((int)state.levels[l].nodes[v1].type == CUT) {
					if ((int)state.levels[l].nodes[v1].inFg == 0) {
						cout << " contradiction1" << endl;
					}
				}
			}
		}**/

		if ((int)state.levels[l].nodes[v1].inFg == 0 && (int)state.levels[l].nodes[v2].inFg == 0) {
			add_edge(v1, v2, bgG);
		}


	}

	//find lexicographical cost for this level
	vector<vector<int>> fgComps = getComponents(fgG, state.levels[l].nodes, 1);
	vector<vector<int>> bgComps = getComponents(bgG, state.levels[l].nodes, 0);

	state.levels[l].h0 = fgComps.size();
	state.levels[l].h2 = bgComps.size() - 1;
	int h1 = fgComps.size() + (bgComps.size() - 1) - (currentEulerNums[0] - currentEulerNums[1] + currentEulerNums[2] - currentEulerNums[3]);
	state.levels[l].h1 = h1;
	state.levels[l].geomCost = geometryCost; //state.levels[l].wtSum
	state.levels[l].cost = (fgComps.size() + state.levels[l].h2 + h1) * state.totalWtSum + geometryCost;


	//cout << "resolved level " << l << " has " << fgComps.size() << " fg components " << bgComps.size() << " bg components  " << 
		//geometryCost << " geometry cost " << numCuts << " " << numFills << " " << oCuts << " " << oFills << endl;

}

State constrainLowerState(State& stateIn, int contraIndex, vector< vector<int> >& overallToLvlNodes, vector<vector<node>>& lvlOrigNodesIn,
	Graph& intersectionG, vector<Graph>& levelGraphs, vector<map<vector<int>, int>>& levelEdgeWts, int hypernodeSize, int productThresh, int localSteinerTime,
	int bbTime, int globalSteinerTime, vector< vector<vector<int>> >& levelNewToOldComps, vector<vector<node> >& origLvlNodes, int& levelSolve) {
	State state = stateIn;
	Contradiction c = state.contradictions[contraIndex];
	//std::vector< std::vector<node>> origLvlNodes = origLvlNodesIn;
	vector<int> lvls;
	for (int i = 0; i < c.lowerIndices.size(); i++) {

		int lvl = overallToLvlNodes[c.lowerIndices[i]][0];
		//cout << "Contra index " << contraIndex << " Lower index " << lvl << " " << c.lowerIndices[i] << " " << (int)state.levels[lvl].origNodes[overallToLvlNodes[c.lowerIndices[i]][1]].type << endl;

		//int index = state.levels[lvl].origNodes[overallToLvlNodes[c.lowerIndices[i]][1]];
		//cout << "lower index " << c.lowerIndices[i] << " lvl " << (int)state.levels[lvl].nodes[index].type << " " << (int)state.levels[lvl].nodes[index].inFg << endl;
		//cout << "Lower index " << lvl << " " << c.lowerIndices[i] << " " << (int)state.levels[lvl].origNodes[overallToLvlNodes[c.lowerIndices[i]][1]].type << endl;
		if (origLvlNodes[lvl][overallToLvlNodes[c.lowerIndices[i]][1]].type == FILL) {
			auto neighboursp = adjacent_vertices(overallToLvlNodes[c.lowerIndices[i]][1], levelGraphs[lvl]);
			for (auto up : make_iterator_range(neighboursp)) {
				if (origLvlNodes[lvl][(int)up].type == CUT) {
					state.levels[lvl].origNodes[(int)up].inFg = 1;
					state.levels[lvl].origNodes[(int)up].type = CORE;
					//cout << "neig set to 1 " << endl;
				}
			}
		}
		state.levels[lvl].origNodes[overallToLvlNodes[c.lowerIndices[i]][1]].inFg = 1;
		state.levels[lvl].origNodes[overallToLvlNodes[c.lowerIndices[i]][1]].type = CORE;
		lvls.push_back(lvl);

	}
	std::sort(lvls.begin(), lvls.end());
	lvls.erase(unique(lvls.begin(), lvls.end()), lvls.end());
	for (int i = 0; i < lvls.size(); i++) {
		solveLevel(state, lvls[i], levelGraphs, levelEdgeWts, hypernodeSize, productThresh,
			localSteinerTime, bbTime, globalSteinerTime, origLvlNodes[lvls[i]]
		);
		levelSolve++;
	}

	state.contradictions.clear();
	state.cost = 0;
	for (int i = 1; i < state.levels.size() - 1; i++) {
		state.cost += state.levels[i].cost;
	}
	findContradictions(state, overallToLvlNodes, intersectionG, origLvlNodes);
	return state;
}

State constrainUpperState(State& stateIn, int contraIndex, vector< vector<int> >& overallToLvlNodes, vector<vector<node>>& lvlOrigNodes,
	Graph& intersectionG, vector<Graph>& levelGraphs, vector<map<vector<int>, int>>& levelEdgeWts, int hypernodeSize, int productThresh, int localSteinerTime,
	int bbTime, int globalSteinerTime, vector< vector<vector<int>> >& levelNewToOldComps, std::vector< std::vector<node>>& origLvlNodes, int& levelSolve) {
	State state = stateIn;
	Contradiction c = state.contradictions[contraIndex];
	vector<int> lvls;
	for (int i = 0; i < c.upperIndices.size(); i++) {
		int lvl = overallToLvlNodes[c.upperIndices[i]][0];
		//cout << "Upper index " << lvl << " " << c.upperIndices[i] << " " << (int)state.levels[lvl].origNodes[overallToLvlNodes[c.upperIndices[i]][1]].type << endl;
		if (origLvlNodes[lvl][overallToLvlNodes[c.upperIndices[i]][1]].type == CUT) {
			auto neighboursp = adjacent_vertices(overallToLvlNodes[c.upperIndices[i]][1], levelGraphs[lvl]);
			for (auto up : make_iterator_range(neighboursp)) {
				if (origLvlNodes[lvl][(int)up].type == FILL) {
					state.levels[lvl].origNodes[(int)up].inFg = 0;
					state.levels[lvl].origNodes[(int)up].type = N;
					//cout << "neig set to 0 " << endl;
				}
			}
		}
		state.levels[lvl].origNodes[overallToLvlNodes[c.upperIndices[i]][1]].inFg = 0;
		state.levels[lvl].origNodes[overallToLvlNodes[c.upperIndices[i]][1]].type = N;
		lvls.push_back(lvl);

	}
	std::sort(lvls.begin(), lvls.end());
	lvls.erase(unique(lvls.begin(), lvls.end()), lvls.end());
	for (int i = 0; i < lvls.size(); i++) {
		solveLevel(state, lvls[i], levelGraphs, levelEdgeWts, hypernodeSize, productThresh,
			localSteinerTime, bbTime, globalSteinerTime, origLvlNodes[lvls[i]]
		);
		levelSolve++;
	}
	state.contradictions.clear();
	state.cost = 0;
	for (int i = 1; i < state.levels.size() - 1; i++) {
		state.cost += state.levels[i].cost;
	}
	findContradictions(state, overallToLvlNodes, intersectionG, origLvlNodes);
	return state;
}

std::vector<int> getStateEncoding(State& state) {
	vector<int> types;
	for (int i = 0; i < state.levels.size(); i++) {
		types.push_back(-i - 1);
		for (int j = 0; j < state.levels[i].nodes.size(); j++) {
			for (int k = 0; k < state.levelNewToOldComps[i][j].size(); k++) {
				types.push_back((int)state.levels[i].origNodes[state.levelNewToOldComps[i][j][k]].type); //state encodes types
			}
		}
	}
	return types;
}




void reindexGraph(std::vector<node>& nodes, Graph& G, map<int, int>& oldToNewIndex, map<int, int>& newToOldIndex, map< std::vector<int>, int>& edgeWt) {
	int index = 0; std::vector<node> newNodes; Graph newG; map< std::vector<int>, int> newEdgeWt;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].valid) {
			node n = nodes[i];
			n.index = newNodes.size(); n.intensity = nodes[i].intensity;
			oldToNewIndex[i] = newNodes.size(); newToOldIndex[newNodes.size()] = i;
			newNodes.push_back(n);
			add_vertex(newG);
		}
	}

	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
		if (nodes[v1].valid && nodes[v2].valid) {
			add_edge(oldToNewIndex[v1], oldToNewIndex[v2], newG);
			newEdgeWt[{oldToNewIndex[v1], oldToNewIndex[v2]}] = edgeWt[{v1, v2}];
			newEdgeWt[{oldToNewIndex[v2], oldToNewIndex[v1]}] = newEdgeWt[{oldToNewIndex[v1], oldToNewIndex[v2]}];
		}
	}
	nodes = newNodes;
	G = newG;
	edgeWt = newEdgeWt;
}

vector<int> getNums(State state) {
	int numCore = 0;
	int numN = 0;
	int numCuts = 0;
	int numFills = 0;
	for (int i = 0; i < state.levels.size(); i++) {
		for (int j = 0; j < state.levels[i].origNodes.size(); j++) {
			if (state.levels[i].origNodes[j].type == CORE) {
				numCore++;
			}
			if (state.levels[i].origNodes[j].type == N) {
				numN++;
			}
			if (state.levels[i].origNodes[j].type == FILL) {
				numFills++;
			}
			if (state.levels[i].origNodes[j].type == CUT) {
				numCuts++;
			}
		}
	}
	return { numCore, numN, numCuts, numFills };
}

int getContradictionScore(State stateIn, int cindex, int type, vector<Graph>& levelGraphs, vector<map<vector<int>, int>>& levelEdgeWts) {
	if (type == 0) {
		int diff = -abs((int)stateIn.contradictions[cindex].lowerIndices.size() - (int)stateIn.contradictions[cindex].upperIndices.size());
		return diff;
	}
	else {
		State state = stateIn;
		typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
		int topoTotal = 0;
		for (int l = 0; l < state.levels.size(); l++) {
			Graph fgG = Graph();
			Graph bgG = Graph();
			for (int i = 0; i < state.levels[l].nodes.size(); i++) {
				add_vertex(fgG);
				add_vertex(bgG);
			}


			edge_iter ei, ei_end;
			for (tie(ei, ei_end) = edges(levelGraphs[l]); ei != ei_end; ++ei) {
				int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
				if ((int)state.levels[l].nodes[v1].inFg == 1 && (int)state.levels[l].nodes[v1].inFg == 1) {
					if (levelEdgeWts[l][{v1, v2}] == 1) {
						add_edge(v1, v2, fgG);
					}
				}

				if ((int)state.levels[l].nodes[v1].inFg == 0 && (int)state.levels[l].nodes[v1].inFg == 0) {
					add_edge(v1, v2, bgG);
				}
			}
			//find lexicographical cost for this level
			int h0 = findComponents(fgG, state.levels[l].nodes, state.levels[l].nodes.size(), 1);
			int h2 = findComponents(bgG, state.levels[l].nodes, state.levels[l].nodes.size(), 0);
			//int h1 = h0 + h2 - (eulNums[0] - eulNums[1] + eulNums[2] - eulNums[3]);
			return 0;
		}
	}
}

map<int, map<int, int> > findComponents(Graph& g, std::vector<node>& nodes, int fg) {
	int nV = nodes.size();
	std::vector<bool> visited(nV, false); int overallIndex = 0;
	std::vector<bool> visited1(nV, false); int overallIndex1 = 0;
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].inFg) == fg) {
			if (!visited[i]) {
				int timer = 0; int timer1 = 0;
				dfs(g, i, visited, timer, -1, overallIndex, nodes);

				overallIndex += 1;
			}
		}
	}
	//Find component indices for isolated comps
	for (int i = 0; i < nodes.size(); i++) {
		if (((int)nodes[i].inFg) == fg) {
			if (nodes[i].compIndex == -1) {
				nodes[i].compIndex = overallIndex;
				overallIndex += 1;
			}
		}
	}
	map<int, map<int, int> > nodeConnectivity;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].isArticulate) {
			auto neighbourItr = adjacent_vertices(i, g);
			std::vector<node> neighbours;
			map<int, int> componentMapping; componentMapping[0] = 0;
			for (auto u : make_iterator_range(neighbourItr)) {
				neighbours.push_back(nodes[u]);
				componentMapping[u] = 0;
			}
			std::sort(neighbours.begin(), neighbours.end(), compareByTime);
			int compIndex = 0;
			map<int, int> nodeToComp;
			for (int j = 0; j < neighbours.size(); j++) {
				if (neighbours[j].low < nodes[i].tin) {
					componentMapping[neighbours[j].tin] = componentMapping[neighbours[j].low]; //Use time of ancestor
				}
				else {
					if (neighbours[j].isNew) { //reporting node
						compIndex += 1;
						componentMapping[neighbours[j].tin] = compIndex;
					}
					else { //Inherit component of nearest left neighbor
						if (j - 1 > 0) {
							componentMapping[neighbours[j].tin] = componentMapping[neighbours[j - 1].tin];
						}
					}
				}
				nodeToComp[neighbours[j].index] = componentMapping[neighbours[j].tin];
			}
			nodeConnectivity[i] = nodeToComp;
		}
	}

	return nodeConnectivity;
}

bool violatesComplex(std::vector<node>& nodes, int v1, int v2, map<std::vector<int>, int>& edgeWt, int fgConn) {

	if (((((int)nodes[v1].inFg) == 1) && (((int)nodes[v1].type) == 3)) // v2 is fill in fg
		&& ((((int)nodes[v2].inFg) == 0) && ((int)nodes[v2].type) == 2) //v1 is cut in bg
		) {
		return true;
	}

	if ((((int)nodes[v2].inFg == 1) && ((int)nodes[v2].type) == 3) // v2 is fill in fg
		&& (((int)nodes[v1].inFg == 0) && ((int)nodes[v1].type) == 2) //v1 is cut in bg
		) {
		return true;
	}

	return false;
}

void findNodeToSwap(Graph& G, Graph& fgG, Graph bgG, std::vector<node>& nodes, std::vector<node>& nodesIn, int& nodeToSwap, map< std::vector<int>, int>& edgeWt, int fgConn,
	int origFgComps, int origBgComps, int& newFgComps, int& newBgComps, int costType) {
	int bgConn = 1 - fgConn;

	nodeToSwap = -1;
	std::vector< tuple<int, int, int, int, int, int, int, int> > nodeCosts; //std::vector< tuple<int, int, int, int, int, int, int, int> > nodeCosts;
	newFgComps = origFgComps;
	newBgComps = origBgComps;
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i].isArticulate = false; nodes[i].compIndex = -1;
		nodes[i].isNew = false;
		nodes[i].tin = 0;
		nodes[i].low = 0;
		nodes[i].overallCompIndexFg = -1; nodes[i].overallCompIndexBg = -1;
	}
	map<int, map<int, int> > fgLocalArtConnectivity = findComponents(fgG, nodes, 1);
	map<int, map<int, int> > bgLocalArtConnectivity = findComponents(bgG, nodes, 0);

	int artCt = 0;
	for (int i = 0; i < nodes.size(); i++) {
		if (!nodes[i].valid) { continue; }

		if (((int)nodes[i].inFg) == 1 && ((((int)nodesIn[i].type) == 2))) {
			int changeH0 = 0;
			if (nodes[i].isArticulate) {
				artCt += 1;
				map<int, int> nodeToComp = fgLocalArtConnectivity[i];
				std::vector<int> uniqueComps;
				for (std::map<int, int>::iterator iter = nodeToComp.begin(); iter != nodeToComp.end(); ++iter)
				{
					uniqueComps.push_back(iter->second);
				}
				std::sort(uniqueComps.begin(), uniqueComps.end());
				uniqueComps.erase(unique(uniqueComps.begin(), uniqueComps.end()), uniqueComps.end());
				changeH0 = uniqueComps.size() - 1;
			}
			std::vector<int> nBgComps;
			auto neighbours = adjacent_vertices(i, G); //If connect to BG, how many components would it connect?
			int adjFg = 0;
			for (auto u : make_iterator_range(neighbours)) {
				if (nodes[u].valid) {
					if (((int)nodes[u].inFg) == 0 && nodes[u].compIndex != -1) {
						nBgComps.push_back(nodes[u].compIndex);
					}
					if (((int)nodes[u].inFg) == 1) {
						adjFg += 1;
					}
				}
			}
			if (adjFg == 0 && changeH0 == 0) {
				changeH0 = -1;
			}
			std::sort(nBgComps.begin(), nBgComps.end());
			nBgComps.erase(unique(nBgComps.begin(), nBgComps.end()), nBgComps.end());
			int changeH2 = -((int)nBgComps.size() - 1); //Would connect together this many background components

			nodes[i].inFg = 0;
			auto neighbours1 = adjacent_vertices(i, G);
			bool violatesCellComplex = false;
			for (auto vd : make_iterator_range(neighbours1)) {
				if (nodes[vd].valid) {
					if (violatesComplex(nodes, i, vd, edgeWt, bgConn)) {
						violatesCellComplex = true;
						break;

					}
				}
			}
			if (violatesCellComplex) {
				nodes[i].inFg = 1;
				continue;
			}

			int eulerNum = nodes[i].v - nodes[i].e + nodes[i].f - nodes[i].c;
			int deltaH1 = (changeH0)+(changeH2)-(-eulerNum);
			cout << changeH0 << " " << changeH2 << " " << deltaH1 << endl;
			if (costType == 0) {
				nodeCosts.push_back(std::make_tuple((changeH0)+(changeH2)+deltaH1, nodes[i].intensity, i, origFgComps + changeH0, origBgComps + changeH2, changeH0, changeH2, deltaH1));
			}
			if (costType == 1) {
				nodeCosts.push_back(std::make_tuple((changeH0)+(changeH2), nodes[i].intensity, i, origFgComps + changeH0, origBgComps + changeH2, changeH0, changeH2, deltaH1));

			}
			nodes[i].inFg = 1;
		}
		else {
			if (((int)nodes[i].inFg) == 0 && (((int)nodesIn[i].type) == 3)) {
				int changeH2 = 0;
				if (nodes[i].isArticulate) {
					artCt += 1;
					map<int, int> nodeToComp = bgLocalArtConnectivity[i];
					std::vector<int> uniqueComps;
					for (std::map<int, int>::iterator iter = nodeToComp.begin(); iter != nodeToComp.end(); ++iter)
					{
						uniqueComps.push_back(iter->second);
					}
					std::sort(uniqueComps.begin(), uniqueComps.end());
					uniqueComps.erase(unique(uniqueComps.begin(), uniqueComps.end()), uniqueComps.end());
					changeH2 = uniqueComps.size() - 1;
				}
				std::vector<int> nFgComps;
				auto neighbours = adjacent_vertices(i, G); //If connect to BG, how many components would it connect?
				int adjBg = 0;
				for (auto u : make_iterator_range(neighbours)) {
					if (nodes[u].valid) {
						if (((int)nodes[u].inFg) == 1 && nodes[u].compIndex != -1) {
							nFgComps.push_back(nodes[u].compIndex);
						}
						if (((int)nodes[u].inFg) == 0) {
							adjBg += 1;
						}
					}
				}
				std::sort(nFgComps.begin(), nFgComps.end());
				nFgComps.erase(unique(nFgComps.begin(), nFgComps.end()), nFgComps.end());
				int changeH0 = -((int)nFgComps.size() - 1); //Would connect together this many background components
				if (adjBg == 0 && changeH2 == 0) {
					changeH2 = -1;
				}

				nodes[i].inFg = 1;
				auto neighbours1 = adjacent_vertices(i, G);
				bool violatesCellComplex = false;
				for (auto vd : make_iterator_range(neighbours1)) {
					if (nodes[vd].valid) {
						if (violatesComplex(nodes, i, vd, edgeWt, bgConn)) {
							violatesCellComplex = true;
							break;

						}
					}
				}
				if (violatesCellComplex) {
					nodes[i].inFg = 0;
					continue;
				}

				int eulerNum = nodes[i].v - nodes[i].e + nodes[i].f - nodes[i].c;
				int deltaH1 = (changeH0)+(changeH2)-(eulerNum);
				if (costType == 0) {
					nodeCosts.push_back(std::make_tuple(changeH0 + changeH2 + deltaH1, -nodes[i].intensity, i, origFgComps + changeH0, origBgComps + changeH2, changeH0, changeH2, deltaH1));
				}
				if (costType == 1) {
					nodeCosts.push_back(std::make_tuple((changeH0)+(changeH2), nodes[i].intensity, i, origFgComps + changeH0, origBgComps + changeH2, changeH0, changeH2, deltaH1));

				}

				nodes[i].inFg = 0;
			}
		}
	}
	if (nodeCosts.size() > 0) {
		std::sort(nodeCosts.begin(), nodeCosts.end());

		if (get<0>(nodeCosts[0]) < 0) {
			nodeToSwap = get<2>(nodeCosts[0]);
			newFgComps = get<3>(nodeCosts[0]);
			newBgComps = get<4>(nodeCosts[0]);
		}
		else {
			if (get<0>(nodeCosts[0]) == 0) {
				if (get<1>(nodeCosts[0]) < 0) {
					nodeToSwap = get<2>(nodeCosts[0]);
					newFgComps = get<3>(nodeCosts[0]);
					newBgComps = get<4>(nodeCosts[0]);
				}
			}
		}
	}
}

int addToSubgraph(Graph& G, Graph& subG, std::vector<node>& nodes, node n, int fg, int conn, map< std::vector<int>, int>& edgeWt) {
	auto neighbours = adjacent_vertices(n.index, G);
	int nEdgesAdded = 0;
	for (auto vd : make_iterator_range(neighbours)) {

		if (((int)nodes[vd].inFg) == fg) {
			if (nodes[vd].valid) {
				if (conn == 0 || (conn == 1 && edgeWt[{(int)n.index, (int)vd}] == 1)) { //If fg connectivity is cube or fg connectivity is structCross3D and edge is strong
					//add_edge(n.index, vd, 0, subG);
					add_edge(n.index, vd, subG);
					nEdgesAdded += 1;
				}
			}
		}
	}

	return nEdgesAdded;
}

int getNumComponents(Graph& g, std::vector<node>& nodes, int fg) {
	std::vector<int> nodeToComp(nodes.size());
	int n = (int)boost::connected_components(g, &nodeToComp[0]);
	int numComps = 0;
	std::vector<bool> isCompIndex(nodes.size(), false);
	for (int i = 0; i < nodeToComp.size(); i++) {
		if (((int)nodes[i].inFg) == fg && nodes[i].valid) {
			if (!isCompIndex[nodeToComp[i]]) {
				isCompIndex[nodeToComp[i]] = true;
				numComps += 1;
			}
		}
	}
	return numComps;
}

void swapLabelsGreedy(Graph& G, Graph& fgG, Graph& bgG, std::vector<node>& nodes, map< std::vector<int>, int>& edgeWt, int fgConn, int costType) {
	std::vector<node> newNodes = nodes;
	std::vector<node> origNodes = nodes;
	int bgConn = 1 - fgConn;
	int fgComps = getNumComponents(fgG, newNodes, 1); int bgComps = getNumComponents(bgG, newNodes, 0);
	int nodeToSwap = -1;
	int newFgComps = -1; int newBgComps = -1;
	for (int i = 0; i < nodes.size(); i++) {
		if ((int)nodes[i].type == CUT || (int)nodes[i].type == FILL) {

			nodes[i].intensity = nodes[i].labelCost;
			origNodes[i].intensity = nodes[i].labelCost;
			newNodes[i].intensity = nodes[i].labelCost;
		}
	}
	findNodeToSwap(G, fgG, bgG, origNodes, newNodes, nodeToSwap, edgeWt, fgConn, fgComps, bgComps, newFgComps, newBgComps, costType);
	int itr = 0;
	while (nodeToSwap != -1) {
		itr += 1;
		if (((int)newNodes[nodeToSwap].inFg) == 1) {
			newNodes[nodeToSwap].inFg = 0;
			origNodes[nodeToSwap].inFg = 0;
			fgComps = newFgComps; bgComps = newBgComps;
			clear_vertex(nodeToSwap, fgG);
			addToSubgraph(G, bgG, newNodes, newNodes[nodeToSwap], 0, bgConn, edgeWt);
		}
		else {
			if (((int)newNodes[nodeToSwap].inFg) == 0) {
				newNodes[nodeToSwap].inFg = 1;
				origNodes[nodeToSwap].inFg = 1;
				fgComps = newFgComps; bgComps = newBgComps;
				clear_vertex(nodeToSwap, bgG);
				addToSubgraph(G, fgG, newNodes, newNodes[nodeToSwap], 1, fgConn, edgeWt);
			}
		}
		for (int i = 0; i < newNodes.size(); i++) {
			origNodes[i].inFg = newNodes[i].inFg;
		}
		findNodeToSwap(G, fgG, bgG, origNodes, newNodes, nodeToSwap, edgeWt, fgConn, fgComps, bgComps, newFgComps, newBgComps, costType);
	}
	for (int i = 0; i < newNodes.size(); i++) {
		nodes[i].inFg = newNodes[i].inFg;
	}
}

float euclideanDistance(vector<float>& v1, vector<float>& v2) {
	return sqrtf((v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1]) + (v1[2] - v2[2]) * (v1[2] - v2[2]));
}

void simplifyGenerator(Graph& G, std::vector<uint32_t*>& labels, const std::vector<float*>& g_Image3D,
	map<int, map<int, int> >& fgLocalArtConnectivity,
	map<int, map<int, int> >& bgLocalArtConnectivity, int numSlices, int width, int height, int fgConn,
	std::vector<node>& nodes, int nodeIndex, map<std::vector<int>, int>& edgeWt, float S,
	const std::vector<unsigned char>& simpleDictionary3D, const std::vector<Coordinate>& nodeVoxels, int genMode, map<int, int>& parentCompFg, map<int, int>& parentCompBg, bool& generatorChanged
	, grapht& fgGWithFills, grapht& bgGWithCuts, grapht& coreG, grapht& nG, float& minigraphT, float& addETime, float& otherTime, int costType
) {
	//grapht & fgG, grapht & bgG, 
	auto t1 = std::chrono::high_resolution_clock::now();
	std::vector< std::vector<int> > fgMask, bgMask;
	if (fgConn == 0) {
		fgMask = structCube;
		bgMask = structCross3D;
	}
	else {
		fgMask = structCross3D;
		bgMask = structCube;
	}
	auto t_start = std::chrono::high_resolution_clock::now();
	if (((int)nodes[nodeIndex].type) == 3) { //Fill

		auto neighbours = adjacent_vertices(nodeIndex, G);
		std::vector<int> compIndices; int nodeAdjFg = 0;
		for (auto u : make_iterator_range(neighbours)) {

			if (((int)nodes[u].inFg) == 1) {
				if (edgeWt[{nodeIndex, (int)u}] == 1) {
					if (((int)nodes[u].type) == 0 || ((int)nodes[u].type) == 2) {
						if (nodes[u].valid) {
							if (nodes[u].overallCompIndexFg != -1) {
								compIndices.push_back(nodes[u].overallCompIndexFg);
								nodeAdjFg += 1;
							}
						}
					}
				}
			}
		}
		//Delete duplicate connected components
		std::sort(compIndices.begin(), compIndices.end());
		compIndices.erase(unique(compIndices.begin(), compIndices.end()), compIndices.end());

		//Fill not connected to any fg comps
		int vx, vy, vz, nx, ny, nz;

		map<Coordinate, int> vtxToGNode;
		int  newVtxIndex = 1;
		map<int, int> fgCompToMini;
		if (compIndices.size() > 0) {
			for (int i = 0; i < compIndices.size(); i++) {
				fgCompToMini[compIndices[i]] = i;
			}
			newVtxIndex = compIndices.size();
		}
		grapht miniFgG = grapht();
		std::vector<weightedCoord> P;
		//Create vertex indices
		std::vector<bool> miniFgValid(newVtxIndex + nodeVoxels.size(), false);

		//std::vector<node> minigraphNodes;
		for (int i = 0; i < nodeVoxels.size(); i++) {
			vx = nodeVoxels[i].x; vy = nodeVoxels[i].y; vz = nodeVoxels[i].z;
			vtxToGNode[Coordinate(vx, vy, vz)] = newVtxIndex;
			miniFgValid[newVtxIndex] = true;
			newVtxIndex += 1;


			weightedCoord wc = { vx, vy, vz, abs((float)g_Image3D[vz][vx + width * vy] - S) };
			P.push_back(wc);

		}


		for (int i = 0; i < compIndices.size(); i++) {
			miniFgValid[fgCompToMini[compIndices[i]]] = true;
		}
		for (int i = 0; i < newVtxIndex; i++) {
			add_vertex(miniFgG);
		}
		//Construct internal and external FG minigraph edges
		for (int i = 0; i < nodeVoxels.size(); i++) {
			vx = nodeVoxels[i].x; vy = nodeVoxels[i].y; vz = nodeVoxels[i].z;
			Coordinate vC(vx, vy, vz);
			for (int j = 0; j < fgMask.size(); j++) {
				nx = vx + fgMask[j][0]; ny = vy + fgMask[j][1]; nz = vz + fgMask[j][2];
				Coordinate vN(nx, ny, nz);
				if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
					if (Label(labels[nz][nx + (width * (ny))]) == Label(labels[vz][vx + (width * (vy))])) { //edge within voxel
						if (!boost::edge(vtxToGNode[vC], vtxToGNode[vN], miniFgG).second) {
							add_edge(vtxToGNode[vC], vtxToGNode[vN], miniFgG);
						}
					}
					else {
						if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0
							|| ((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 2
							) {
							if (!boost::edge(vtxToGNode[vC], fgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexFg], miniFgG).second) {
								add_edge(vtxToGNode[vC], fgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexFg], miniFgG);
							}
						}
					}
				}
			}
		}

		int h0 = 1; int h2;
		map<int, int> parentComp;
		map<int, int> nodeToComp;
		auto neighboursN = adjacent_vertices(nodeIndex, G);
		std::vector<int> uniqueComps;
		for (auto u : make_iterator_range(neighboursN)) {
			if (((int)nodes[u].type) == N) {
				uniqueComps.push_back(getParent(parentCompBg, nodes[u].compIndex));
			}
		}
		std::sort(uniqueComps.begin(), uniqueComps.end());
		uniqueComps.erase(unique(uniqueComps.begin(), uniqueComps.end()), uniqueComps.end());

		h2 = uniqueComps.size();
		//Next step sort voxels by intensity and go thru them
		std::sort(P.begin(), P.end(), compareByIntensity);
		int voxelsDeleted = 0;
		std::vector<vertex_t> art_points;
		//if (P.size() > 100) {
			//cout << "P size " << P.size() << endl;
		//}

		std::vector<bool> checkVoxel(num_vertices(miniFgG), true);
		for (int i = 0; i < P.size(); i++) {
			weightedCoord wc = P[i];

			if (Label(labels[wc.z][wc.x + (width * (wc.y))]) == nodeIndex) {

				//auto t1 = std::chrono::high_resolution_clock::now();
				Coordinate wcC(wc.x, wc.y, wc.z);
				if (!checkVoxel[vtxToGNode[wcC]]) {
					//cout << "checked voxel already " << endl;
					continue;
				}
				int changeH2;
				std::vector<int> nBgNodes;
				for (int j = 0; j < bgMask.size(); j++) {
					nx = wc.x + bgMask[j][0]; ny = wc.y + bgMask[j][1]; nz = wc.z + bgMask[j][2];
					if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
						if (Label(labels[nz][nx + (width * (ny))]) != nodeIndex) { //edge within voxel
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1
								) {
								nBgNodes.push_back(getParent(parentCompBg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex));
							}
						}
					}
				}
				//Number of neighboring bg components
				std::sort(nBgNodes.begin(), nBgNodes.end());
				nBgNodes.erase(unique(nBgNodes.begin(), nBgNodes.end()), nBgNodes.end());

				changeH2 = 1 - nBgNodes.size();
				int h0N; int changeH0;

				if (P.size() == 1 && nodeAdjFg == 0) {
					changeH0 = -1;
				}
				else {

					clear_vertex(vtxToGNode[wcC], miniFgG);
					miniFgValid[vtxToGNode[wcC]] = false;

					auto t1 = std::chrono::high_resolution_clock::now();
					/**if (P.size() > 100) {
						int nSame = 0;
						for (int j = 0; j < bgMask.size(); j++) {
							nx = wc.x + bgMask[j][0]; ny = wc.y + bgMask[j][1]; nz = wc.z + bgMask[j][2];

							if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
								if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) { //edge within voxel
									nSame++;
								}
							}
						}
						cout << "call to mini " << nSame << endl;

					}**/
					h0N = minigraphComps(miniFgG, miniFgValid, newVtxIndex);
					auto t2 = std::chrono::high_resolution_clock::now();
					std::chrono::duration<double> elapsed = t2 - t1;
					minigraphT += elapsed.count();
					changeH0 = h0N - h0;
				}

				int dV = 1, dE = 0, dF = 0, dC = 0;
				for (int j = 0; j < fgMask.size(); j++) {
					nx = wc.x + fgMask[j][0]; ny = wc.y + fgMask[j][1]; nz = wc.z + fgMask[j][2];
					if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
						if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0 ||
							((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3 ||
							((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 2
							) {
							dE += 1;
						}
					}
				}
				for (int j = 0; j < adjFaces6Conn.size(); j++) {
					std::vector< std::vector<int> > adjFace = adjFaces6Conn[j];
					bool hasN = false;
					for (int k = 0; k < adjFace.size(); k++) {
						nx = wc.x + adjFace[k][0]; ny = wc.y + adjFace[k][1]; nz = wc.z + adjFace[k][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1) {
								hasN = true;
							}
						}
					}
					if (!hasN) {
						dF += 1;
					}
				}
				for (int j = 0; j < adjCubes6Conn.size(); j++) {
					std::vector< std::vector<int> > adjCube = adjCubes6Conn[j];
					bool hasN = false;
					for (int k = 0; k < adjCube.size(); k++) {
						nx = wc.x + adjCube[k][0]; ny = wc.y + adjCube[k][1]; nz = wc.z + adjCube[k][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1) {
								hasN = true;
							}
						}
					}
					if (!hasN) {
						dC += 1;
					}
				}
				int dEuler = -(dV - dE + dF - dC);
				int changeH1 = changeH0 + changeH2 - dEuler;

				bool reduce = false;
				if (costType == 0) {
					if (changeH0 <= 0 && changeH2 <= 0 && changeH1 <= 0 && changeH0 + changeH2 + changeH1 < 0) {
						reduce = true;
					}
				}
				if (costType == 1) {
					if (changeH0 <= 0 && changeH2 <= 0 && changeH0 + changeH2 < 0) {
						reduce = true;
					}

				}

				if (reduce) {

					generatorChanged = true;
					h2 = h2 + changeH2;
					for (int j = 0; j < bgMask.size(); j++) {
						nx = wc.x + bgMask[j][0]; ny = wc.y + bgMask[j][1]; nz = wc.z + bgMask[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1) {
								labels[wc.z][wc.x + (width * (wc.y))] = labels[nz][nx + (width * (ny))];
								break;
							}
						}
					}

					voxelsDeleted += 1;


					//voxel has now been relabeled to a background component, create potential new edges for background component

					for (int j = 0; j < bgMask.size(); j++) {
						nx = wc.x + bgMask[j][0]; ny = wc.y + bgMask[j][1]; nz = wc.z + bgMask[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[wc.z][wc.x + (width * (wc.y))]) != Label(labels[nz][nx + (width * (ny))])) {
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1 ||
									((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3
									) {
									if (ccNeighborFill6Conn(nodes, wc.x, wc.y, wc.z, nx, ny, nz, labels, width)) { //, g_Image3D, S
										if (edgeWt.find({ (int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))]) }) == edgeWt.end()) {
											add_edge(Label(labels[nz][nx + (width * (ny))]), Label(labels[wc.z][wc.x + (width * (wc.y))]), G);
											if (abs(nx - wc.x) + abs(ny - wc.y) + abs(nz - wc.z) == 1) {
												if (
													(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CORE || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CUT || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == FILL)
													&&
													(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

													) {
													add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

												}
												if (
													(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CORE)
													&&
													(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

													) {
													add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), coreG);

												}
											}
											if (
												(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == N || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CUT || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == FILL)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

												) {
												add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), bgGWithCuts);

											}
											if (
												(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == N)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

												) {
												add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), nG);

											}
										}
									}

									if (abs(nx - wc.x) + abs(ny - wc.y) + abs(nz - wc.z) == 1) {
										if (edgeWt[{(int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))])}] == 0) {
											if (
												(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CORE || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CUT || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == FILL)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

												) {
												add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

											}
										}
										edgeWt[{(int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))])}] = 1;
										edgeWt[{(int)Label(labels[nz][nx + (width * (ny))]), (int)Label(labels[wc.z][wc.x + (width * (wc.y))])}] = 1;


									}
									else {
										if (ccNeighborFill6Conn(nodes, wc.x, wc.y, wc.z, nx, ny, nz, labels, width)) { //, g_Image3D, S
											if (edgeWt.find({ (int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))]) }) == edgeWt.end()) {
												edgeWt[{(int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))])}] = 0;
												edgeWt[{(int)Label(labels[nz][nx + (width * (ny))]), (int)Label(labels[wc.z][wc.x + (width * (wc.y))])}] = 0;
											}
										}
									}
								}
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1 //|| ((int)nodes[Label(labels[nz][nx + (width*(ny))])].type) == 3
									) {
									if (getParent(parentCompBg, nodes[Label(labels[wc.z][wc.x + (width * (wc.y))])].compIndex)
										!= getParent(parentCompBg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex)
										) {
										parentCompBg[getParent(parentCompBg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex)] = getParent(parentCompBg, nodes[Label(labels[wc.z][wc.x + (width * (wc.y))])].compIndex);

									}
								}
							}

						}


					}

					queue<Coordinate> q;
					//Iteratively remove simple voxels
					for (int j = 0; j < structCube.size(); j++) {
						nx = wc.x + structCube[j][0]; ny = wc.y + structCube[j][1]; nz = wc.z + structCube[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
								q.push(Coordinate(nx, ny, nz));
							}
						}
					}
					map<Coordinate, int> genVisited;
					for (int j = 0; j < P.size(); j++) {

						if (Label(labels[P[j].z][P[j].x + (width * (P[j].y))]) == nodeIndex) {
							Coordinate pNow(P[j].x, P[j].y, P[j].z);
							genVisited[pNow] = 0;
						}
					}




					Coordinate wCoord(wc.x, wc.y, wc.z);
					genVisited[wCoord] = 1;
					while (!q.empty()) {
						Coordinate p = q.front();
						int px = p.x; int py = p.y; int pz = p.z;
						q.pop();
						genVisited[p] = 1;
						if (simple3DLabel(labels, nodes, px, py, pz, numSlices, width, height, simpleDictionary3D, nodeIndex, true, fgConn)) {
							if ((int)Label(labels[pz][px + (width * (py))]) == nodeIndex) {
								voxelsDeleted += 1;
								for (int j = 0; j < structCube.size(); j++) {
									nx = px + structCube[j][0]; ny = py + structCube[j][1]; nz = pz + structCube[j][2];
									if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
										if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
											if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1) {
												labels[pz][px + (width * (py))] = labels[nz][nx + (width * (ny))];
												clear_vertex(vtxToGNode[p], miniFgG);
												miniFgValid[vtxToGNode[p]] = false;
												for (int a = 0; a < structCube.size(); a++) {
													int nxi = px + structCube[a][0];
													int nyi = py + structCube[a][1];
													int nzi = pz + structCube[a][2];
													if (nxi >= 0 && nyi >= 0 && nzi >= 0 && nxi < width && nyi < height && nzi < numSlices) {
														if (Label(labels[pz][px + (width * (py))]) != Label(labels[nzi][nxi + (width * (nyi))])) {
															if (((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 1 ||
																((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 2 ||
																((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 3
																) {
																if (ccNeighborFill6Conn(nodes, px, py, pz, nxi, nyi, nzi, labels, width)) {
																	if (edgeWt.find({ (int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))]) }) == edgeWt.end()) {
																		add_edge(Label(labels[nzi][nxi + (width * (nyi))]), Label(labels[pz][px + (width * (py))]), G);
																		if (abs(nxi - px) + abs(nyi - py) + abs(nzi - pz) == 1) {
																			if (
																				(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																				&&
																				(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CORE || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CUT || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == FILL)

																				) {
																				add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), fgGWithFills);

																			}
																			if (
																				(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
																				&&
																				(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CORE)

																				) {
																				add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), coreG);

																			}
																		}
																		if (
																			(((int)nodes[Label(labels[pz][px + py * width])].type) == N)
																			&&
																			(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == N)

																			) {
																			add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), nG);

																		}

																	}
																	if (abs(nxi - px) + abs(nyi - py) + abs(nzi - pz) == 1) {
																		edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))])}] = 1;
																		edgeWt[{(int)Label(labels[nzi][nxi + (width * (nyi))]), (int)Label(labels[pz][px + (width * (py))])}] = 1;
																		if (edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nz][nx + (width * (ny))])}] == 0) {
																			if (
																				(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																				&&
																				(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CORE || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CUT || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == FILL)

																				) {
																				add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), fgGWithFills);

																			}
																		}
																		edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))])}] = 1;
																		edgeWt[{(int)Label(labels[nzi][nxi + (width * (nyi))]), (int)Label(labels[pz][px + (width * (py))])}] = 1;



																	}
																	else {
																		if (edgeWt.find({ (int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))]) }) == edgeWt.end()) {
																			edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))])}] = 0;
																			edgeWt[{(int)Label(labels[nzi][nxi + (width * (nyi))]), (int)Label(labels[pz][px + (width * (py))])}] = 0;
																		}
																	}
																}
															}

														}

														if (((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == N) {
															if (getParent(parentCompBg, nodes[Label(labels[pz][px + (width * (py))])].compIndex)
																!= getParent(parentCompBg, nodes[Label(labels[nzi][nxi + (width * (nyi))])].compIndex)
																) {
																parentCompBg[getParent(parentCompBg, nodes[Label(labels[nzi][nxi + (width * (nyi))])].compIndex)] = getParent(parentCompBg, nodes[Label(labels[pz][px + (width * (py))])].compIndex);

															}
														}
													}

												}
											}
										}
										if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
											Coordinate np(nx, ny, nz);
											if (genVisited[np] == 0) {
												q.push(np);
												genVisited[np] = 1;
											}
										}
									}
								}
							}
						}
						else {
							genVisited[p] = 0;
						}

					}

				}
				else {
					//if h0, h1, h2 total does not decrease, add back to graph
					Coordinate wCoord(wc.x, wc.y, wc.z);
					miniFgValid[vtxToGNode[wCoord]] = true;
					//checkVoxel[vtxToGNode[wCoord]] = false;
					for (int j = 0; j < fgMask.size(); j++) {
						nx = wc.x + fgMask[j][0]; ny = wc.y + fgMask[j][1]; nz = wc.z + fgMask[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
								if (!boost::edge(vtxToGNode[{nx, ny, nz}], vtxToGNode[{wc.x, wc.y, wc.z}], miniFgG).second) {
									add_edge(vtxToGNode[Coordinate(nx, ny, nz)], vtxToGNode[wCoord], miniFgG);
								}
							}
							else {
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0
									|| ((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 2
									) {
									if (!boost::edge(vtxToGNode[{wc.x, wc.y, wc.z}], fgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexFg], miniFgG).second) {
										add_edge(vtxToGNode[wCoord], fgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexFg], miniFgG);
									}
								}

							}
						}
					}

					/**

					map<Coordinate, int> genVisited;
					std::vector<Coordinate> changedVoxels;
					for (int j = 0; j < P.size(); j++) {

						if (Label(labels[P[j].z][P[j].x + (width * (P[j].y))]) == nodeIndex) {
							Coordinate pNow(P[j].x, P[j].y, P[j].z);
							genVisited[pNow] = 0;

						}
					}
					queue<Coordinate> q;
					genVisited[wCoord] = 1;
					for (int j = 0; j < structCube.size(); j++) {
						nx = wCoord.x + structCube[j][0]; ny = wCoord.y + structCube[j][1]; nz = wCoord.z + structCube[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (abs(nx - wCoord.x) + abs(ny - wCoord.y) + abs(nz - wCoord.z) == 1) {
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1) {
									labels[wCoord.z][wCoord.x + (width * (wCoord.y))] = labels[nz][nx + (width * (ny))];

								}
							}
						}
					}
					//Iteratively remove simple voxels
					for (int j = 0; j < structCube.size(); j++) {
						nx = wCoord.x + structCube[j][0]; ny = wCoord.y + structCube[j][1]; nz = wCoord.z + structCube[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
								q.push(Coordinate(nx, ny, nz));
							}
						}
					}

					changedVoxels.push_back(wCoord);
					checkVoxel[vtxToGNode[wCoord]] = false;

					while (!q.empty()) {
						Coordinate p = q.front();
						int px = p.x; int py = p.y; int pz = p.z;
						q.pop();
						genVisited[p] = 1;
						if (simple3DLabel(labels, nodes, px, py, pz, numSlices, width, height, simpleDictionary3D, nodeIndex, true, fgConn)) {
							if ((int)Label(labels[pz][px + (width * (py))]) == nodeIndex) {
								changedVoxels.push_back(p);
								checkVoxel[vtxToGNode[p]] = false;
								//cout << "check voxel update 1 " << vtxToGNode[wCoord] << " " << vtxToGNode[p] << endl;
								for (int j = 0; j < structCube.size(); j++) {
									nx = px + structCube[j][0]; ny = py + structCube[j][1]; nz = pz + structCube[j][2];
									if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
										if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
											if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1) {
												labels[pz][px + (width * (py))] = labels[nz][nx + (width * (ny))];

											}
										}
										if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
											Coordinate np(nx, ny, nz);
											if (genVisited[np] == 0) {
												q.push(np);
												genVisited[np] = 1;
											}
										}
									}
								}
							}
						}
						else {
							genVisited[p] = 0;
						}

					}
					for (int c = 0; c < changedVoxels.size(); c++) {
						Coordinate cv = changedVoxels[c];
						labels[cv.z][cv.x + (width * (cv.y))] = nodeIndex;
					}**/
				}

			}
		}


		if (voxelsDeleted >= P.size()) {
			nodes[nodeIndex].valid = false;
			auto neighbours = adjacent_vertices(nodeIndex, G);
			for (auto u : make_iterator_range(neighbours)) {
				edgeWt.erase({ nodeIndex, (int)u }); edgeWt.erase({ (int)u , nodeIndex });
			}
			clear_vertex(nodeIndex, G);
			for (int j = 0; j < P.size(); j++) {
				int px = P[j].x; int py = P[j].y; int pz = P[j].z;
				for (int k = 0; k < structCube.size(); k++) {
					int nx = px + structCube[k][0]; int ny = py + structCube[k][1]; int nz = pz + structCube[k][2];
					if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
						if (Label(labels[pz][px + py * width]) != Label(labels[nz][nx + ny * width])) {
							if (((int)nodes[Label(labels[nz][nx + ny * width])].type) != CORE) {
								if (edgeWt.find({ (int)Label(labels[pz][px + py * width]) , (int)Label(labels[nz][nx + ny * width]) }) == edgeWt.end()) {
									if (ccNeighborFill6Conn(nodes, px, py, pz, nx, ny, nz, labels, width)) {
										add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), G);
										if (
											(((int)nodes[Label(labels[pz][px + py * width])].type) == N)
											&&
											(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

											) {
											add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), nG);

										}

										if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
											edgeWt[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
											edgeWt[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
											if (
												(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

												) {
												add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

											}
											if (
												(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

												) {
												add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), coreG);

											}
										}
										else {
											edgeWt[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 0;
											edgeWt[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 0;
										}
									}
								}
								else {
									if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
										edgeWt[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
										edgeWt[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
									}
								}
							}
						}
					}
				}


			}
		}

		else {
			if (voxelsDeleted > 0) {
				auto t1 = std::chrono::high_resolution_clock::now();

				map< std::vector<int>, int> edgeWtTemp = edgeWt;
				auto neighbours = adjacent_vertices(nodeIndex, G);
				for (auto u : make_iterator_range(neighbours)) {
					edgeWtTemp.erase({ nodeIndex, (int)u }); edgeWtTemp.erase({ (int)u , nodeIndex });
				}
				clear_vertex(nodeIndex, G);
				map<Coordinate, bool> visited;
				for (int j = 0; j < P.size(); j++) {
					visited[Coordinate(P[j].x, P[j].y, P[j].z)] = false;
				}
				int currentLabelIndex = nodeIndex;

				for (int j = 0; j < P.size(); j++) {
					Coordinate pNow(P[j].x, P[j].y, P[j].z);
					if (!visited[pNow]) {
						if (Label(labels[pNow.z][pNow.x + (width * (pNow.y))]) == nodeIndex) {


							queue<Coordinate> q;
							q.push(pNow);
							if (currentLabelIndex != nodeIndex) {
								add_vertex(G);
								node n; n.labelCost = 0; n.floatCost = 0.0;  n.type = 3; n.inFg = 0; n.index = currentLabelIndex;
								nodes.push_back(n);
							}
							else {
								//nodes[currentLabelIndex] = n;
								nodes[currentLabelIndex].labelCost = 0;
								nodes[currentLabelIndex].floatCost = 0.0;  nodes[currentLabelIndex].type = 3; nodes[currentLabelIndex].inFg = 0; nodes[currentLabelIndex].index = currentLabelIndex;
							}
							std::vector<float> diffs;
							while (!q.empty()) {
								Coordinate p = q.front();
								int px = p.x; int py = p.y; int pz = p.z;
								q.pop();
								visited[p] = true;

								nodes[currentLabelIndex].floatCost += (-gradient(px, py, pz, g_Image3D, width, height, numSlices));


								uint32_t label32 = currentLabelIndex;
								changeLabel(label32, labels[pz][px + py * width]);
								for (int k = 0; k < structCube.size(); k++) {
									nx = px + structCube[k][0]; ny = py + structCube[k][1]; nz = pz + structCube[k][2];
									if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
										if (Label(labels[pz][px + py * width]) != Label(labels[nz][nx + ny * width])) {
											if (nodeIndex != Label(labels[nz][nx + ny * width])) {
												if (ccNeighborFill6Conn(nodes, px, py, pz, nx, ny, nz, labels, width)) { //, g_Image3D, S
													if (edgeWtTemp.find({ (int)Label(labels[pz][px + py * width]),(int)Label(labels[nz][nx + ny * width]) }) == edgeWtTemp.end()) {
														add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), G);
														if (
															(((int)nodes[Label(labels[pz][px + py * width])].type) == N || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
															&&
															(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

															) {
															add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), bgGWithCuts);

														}
														if (
															(((int)nodes[Label(labels[pz][px + py * width])].type) == N)
															&&
															(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

															) {
															add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), nG);
														}
														if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
															edgeWtTemp[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
															edgeWtTemp[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
															if (
																(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																&&
																(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

																) {
																add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

															}

															if (
																(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
																&&
																(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

																) {
																add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), coreG);

															}
														}
														else {
															edgeWtTemp[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 0;
															edgeWtTemp[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 0;

														}
													}
													else {
														if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
															edgeWtTemp[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
															edgeWtTemp[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
														}
													}
												}
											}
										}
										if (Label(labels[nz][nx + ny * width]) == nodeIndex) {
											Coordinate np(nx, ny, nz);
											if (!visited[np]) {
												if (ccNeighborFill6Conn(nodes, px, py, pz, nx, ny, nz, labels, width)) {
													q.push(np);
													visited[np] = true;
												}
											}
										}
									}
								}
							}
							//If empty, do lexicographical ordering if necessary
							if (currentLabelIndex == nodeIndex) {
								currentLabelIndex = nodes.size();
							}
							else {
								currentLabelIndex += 1;
							}
						}
					}
				}
				edgeWt = edgeWtTemp;

				auto t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = t2 - t1;
				addETime += elapsed.count();
			}
		}

	}
	//node type is cut
	if (((int)nodes[nodeIndex].type) == 2) {
		int bgConn = 1 - fgConn;
		auto neighbours = adjacent_vertices(nodeIndex, G);
		std::vector<int> compIndices; int adjNodesBg = 0;
		for (auto u : make_iterator_range(neighbours)) {
			if (((int)nodes[u].inFg) == 0) {
				if (edgeWt[{nodeIndex, (int)u}] == 1 || (edgeWt[{nodeIndex, (int)u}] == 0 &&
					bgConn == 0
					)) {
					if (nodes[u].valid) {
						if (((int)nodes[u].type) == 1 || ((int)nodes[u].type) == 3) {
							if (nodes[u].overallCompIndexBg != -1) {
								compIndices.push_back(nodes[u].overallCompIndexBg);
								adjNodesBg += 1;
							}
						}
					}
				}
			}

		}




		//Delete duplicate connected components
		std::sort(compIndices.begin(), compIndices.end());
		compIndices.erase(unique(compIndices.begin(), compIndices.end()), compIndices.end());
		//Cut not connected to any bg comps
		int vx, vy, vz, nx, ny, nz;
		map<Coordinate, int> vtxToGNode;
		int  newVtxIndex = 1;
		map<int, int> bgCompToMini;
		if (compIndices.size() > 0) {
			for (int i = 0; i < compIndices.size(); i++) {
				bgCompToMini[compIndices[i]] = i;
			}
			newVtxIndex = compIndices.size();

		}
		grapht miniBgG = grapht();
		std::vector<weightedCoord> P;
		//Create vertex indices
		std::vector<bool> miniBgValid(newVtxIndex + nodeVoxels.size(), false);
		for (int i = 0; i < nodeVoxels.size(); i++) {
			vx = nodeVoxels[i].x; vy = nodeVoxels[i].y; vz = nodeVoxels[i].z;
			vtxToGNode[Coordinate(vx, vy, vz)] = newVtxIndex;
			miniBgValid[newVtxIndex] = true;
			newVtxIndex += 1;


			weightedCoord wc = { vx, vy, vz, abs((float)g_Image3D[vz][vx + width * vy] - S) };
			P.push_back(wc);

		}

		for (int i = 0; i < compIndices.size(); i++) {
			miniBgValid[bgCompToMini[compIndices[i]]] = true;
		}
		for (int i = 0; i < newVtxIndex; i++) {
			add_vertex(miniBgG);

		}
		//Construct internal and external BG minigraph edges
		for (int i = 0; i < nodeVoxels.size(); i++) {
			vx = nodeVoxels[i].x; vy = nodeVoxels[i].y; vz = nodeVoxels[i].z;
			Coordinate vC(vx, vy, vz);
			for (int j = 0; j < bgMask.size(); j++) {
				nx = vx + bgMask[j][0]; ny = vy + bgMask[j][1]; nz = vz + bgMask[j][2];
				Coordinate vN(nx, ny, nz);
				if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
					if (Label(labels[nz][nx + (width * (ny))]) == Label(labels[vz][vx + (width * (vy))])) { //edge within voxel
						if (!boost::edge(vtxToGNode[vC], vtxToGNode[vN], miniBgG).second) {
							add_edge(vtxToGNode[vC], vtxToGNode[vN], miniBgG);
						}
					}
					else {
						if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1
							|| ((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3
							) {
							if (!boost::edge(vtxToGNode[vC], bgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexBg], miniBgG).second) {
								add_edge(vtxToGNode[vC], bgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexBg], miniBgG);
							}
						}
					}
				}
			}
		}
		int h2 = 1;
		int h0;
		map<int, int> nodeToComp;
		auto neighboursN = adjacent_vertices(nodeIndex, G);
		std::vector<int> uniqueComps;
		for (auto u : make_iterator_range(neighboursN)) {
			if (nodes[u].valid) {
				if (((int)nodes[u].type) == CORE) {
					uniqueComps.push_back(getParent(parentCompFg, nodes[u].compIndex));
				}
			}
		}
		std::sort(uniqueComps.begin(), uniqueComps.end());
		uniqueComps.erase(unique(uniqueComps.begin(), uniqueComps.end()), uniqueComps.end());

		h0 = uniqueComps.size();
		//Next step sort voxels by intensity and go thru them
		std::sort(P.begin(), P.end(), compareByIntensity);
		int voxelsDeleted = 0;
		auto t21 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = t21 - t1;
		otherTime += elapsed1.count();
		std::vector<bool> checkVoxel(num_vertices(miniBgG), true);

		for (int i = 0; i < P.size(); i++) {
			weightedCoord wc = P[i];
			Coordinate wcC(wc.x, wc.y, wc.z);
			if (!checkVoxel[vtxToGNode[wcC]]) {
				continue;
			}

			if (Label(labels[wc.z][wc.x + (width * (wc.y))]) == nodeIndex) {
				Coordinate wcC(wc.x, wc.y, wc.z);

				int changeH0;
				std::vector<int> nFgNodes;
				for (int j = 0; j < fgMask.size(); j++) {
					nx = wc.x + fgMask[j][0]; ny = wc.y + fgMask[j][1]; nz = wc.z + fgMask[j][2];
					if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
						if (Label(labels[nz][nx + (width * (ny))]) != nodeIndex) { //edge within voxel
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0
								) {
								nFgNodes.push_back(getParent(parentCompFg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex));
							}
						}
					}
				}
				//Number of neighboring fg components
				std::sort(nFgNodes.begin(), nFgNodes.end());
				nFgNodes.erase(unique(nFgNodes.begin(), nFgNodes.end()), nFgNodes.end());

				changeH0 = 1 - nFgNodes.size();

				int h2N; int changeH2;

				clear_vertex(vtxToGNode[wcC], miniBgG);
				miniBgValid[vtxToGNode[wcC]] = false;
				auto t1 = std::chrono::high_resolution_clock::now();

				h2N = minigraphComps(miniBgG, miniBgValid, newVtxIndex);
				auto t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = t2 - t1;
				minigraphT += elapsed.count();
				changeH2 = h2N - h2;

				int dV = 1, dE = 0, dF = 0, dC = 0;
				for (int j = 0; j < fgMask.size(); j++) {
					nx = wc.x + fgMask[j][0]; ny = wc.y + fgMask[j][1]; nz = wc.z + fgMask[j][2];
					if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
						if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0
							) {


							dE += 1;
						}
					}
				}
				for (int j = 0; j < adjFaces6Conn.size(); j++) {
					std::vector< std::vector<int> > adjFace = adjFaces6Conn[j];
					bool hasN = false;
					for (int k = 0; k < adjFace.size(); k++) {
						nx = wc.x + adjFace[k][0]; ny = wc.y + adjFace[k][1]; nz = wc.z + adjFace[k][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1 ||
								((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3 ||
								((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 2 //is now core, cut has higher priority
								) {
								hasN = true;
							}
						}
					}
					if (!hasN) {
						dF += 1;
					}
				}
				for (int j = 0; j < adjCubes6Conn.size(); j++) {
					std::vector< std::vector<int> > adjCube = adjCubes6Conn[j];
					bool hasN = false;
					for (int k = 0; k < adjCube.size(); k++) {
						nx = wc.x + adjCube[k][0]; ny = wc.y + adjCube[k][1]; nz = wc.z + adjCube[k][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1 ||
								((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3 ||
								((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 2
								) {
								hasN = true;
							}
						}
					}
					if (!hasN) {
						dC += 1;
					}
				}

				int dEuler = (dV - dE + dF - dC);
				int changeH1 = changeH0 + changeH2 - dEuler;
				bool reduce = false;
				if (costType == 0) {
					if (changeH0 <= 0 && changeH2 <= 0 && changeH1 <= 0 && changeH0 + changeH2 + changeH1 < 0) {
						reduce = true;
					}
				}
				if (costType == 1) {
					if (changeH0 <= 0 && changeH2 <= 0 && changeH0 + changeH2 < 0) {
						reduce = true;
					}

				}
				if (reduce) {
					generatorChanged = true;
					h0 = h0 + changeH0;
					bool relabeled = false;
					for (int j = 0; j < fgMask.size(); j++) {
						nx = wc.x + fgMask[j][0]; ny = wc.y + fgMask[j][1]; nz = wc.z + fgMask[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) != nodeIndex) { //edge within voxel
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0) {
									relabeled = true;
									labels[wc.z][wc.x + (width * (wc.y))] = labels[nz][nx + (width * (ny))];
								}

							}
						}
					}

					voxelsDeleted += 1;
					//voxel has now been relabeled to a fg component, create potential new edges for fg component
					for (int j = 0; j < structCross3D.size(); j++) {
						nx = wc.x + structCross3D[j][0]; ny = wc.y + structCross3D[j][1]; nz = wc.z + structCross3D[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[wc.z][wc.x + (width * (wc.y))]) != Label(labels[nz][nx + (width * (ny))])) {
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0 ||
									((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 2 ||
									((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3
									) {

									if (ccNeighbor(nodes, wc.x, wc.y, wc.z, nx, ny, nz, labels, width, height, numSlices, g_Image3D, S)) {
										if (edgeWt.find({ (int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))]) }) == edgeWt.end()) {
											add_edge(Label(labels[nz][nx + (width * (ny))]), Label(labels[wc.z][wc.x + (width * (wc.y))]), G);
											if (abs(nx - wc.x) + abs(ny - wc.y) + abs(nz - wc.z) == 1) {
												if (
													(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CORE || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CUT || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == FILL)
													&&
													(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

													) {
													add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

												}
												if (
													(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CORE)
													&&
													(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

													) {
													add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), coreG);

												}
											}
											if (
												(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == N || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CUT || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == FILL)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

												) {
												add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), bgGWithCuts);

											}
											if (
												(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == N)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

												) {
												add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), nG);
											}
										}
									}

									if (abs(nx - wc.x) + abs(ny - wc.y) + abs(nz - wc.z) == 1) {
										edgeWt[{(int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))])}] = 1;
										edgeWt[{(int)Label(labels[nz][nx + (width * (ny))]), (int)Label(labels[wc.z][wc.x + (width * (wc.y))])}] = 1;
										if (edgeWt[{(int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))])}] == 0) {
											if (
												(((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CORE || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == CUT || ((int)nodes[Label(labels[wc.z][wc.x + wc.y * width])].type) == FILL)
												&&
												(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

												) {
												add_edge(Label(labels[wc.z][wc.x + wc.y * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

											}
										}
									}
									else {
										if (ccNeighbor(nodes, wc.x, wc.y, wc.z, nx, ny, nz, labels, width, height, numSlices, g_Image3D, S)) {
											if (edgeWt.find({ (int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))]) }) == edgeWt.end()) {
												edgeWt[{(int)Label(labels[wc.z][wc.x + (width * (wc.y))]), (int)Label(labels[nz][nx + (width * (ny))])}] = 0;
												edgeWt[{(int)Label(labels[nz][nx + (width * (ny))]), (int)Label(labels[wc.z][wc.x + (width * (wc.y))])}] = 0;
											}
										}
									}
								}

								if (abs(nx - wc.x) + abs(ny - wc.y) + abs(nz - wc.z) == 1) {
									if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0) {
										if (getParent(parentCompFg, nodes[Label(labels[wc.z][wc.x + (width * (wc.y))])].compIndex)
											!= getParent(parentCompFg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex)
											) {
											parentCompFg[getParent(parentCompFg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex)] = getParent(parentCompFg, nodes[Label(labels[wc.z][wc.x + (width * (wc.y))])].compIndex);
										}
									}
								}
							}

						}


					}

					queue<Coordinate> q;
					//Iteratively remove simple voxels
					for (int j = 0; j < structCube.size(); j++) {
						nx = wc.x + structCube[j][0]; ny = wc.y + structCube[j][1]; nz = wc.z + structCube[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
								q.push(Coordinate(nx, ny, nz));
							}
						}
					}
					map<Coordinate, int> genVisited;
					for (int j = 0; j < P.size(); j++) {

						if (Label(labels[P[j].z][P[j].x + (width * (P[j].y))]) == nodeIndex) {
							Coordinate pNow(P[j].x, P[j].y, P[j].z);
							genVisited[pNow] = 0;
						}
					}
					Coordinate wCoord(wc.x, wc.y, wc.z);
					genVisited[wCoord] = 1;
					while (!q.empty()) {
						Coordinate p = q.front();
						int px = p.x; int py = p.y; int pz = p.z;
						q.pop();
						genVisited[{px, py, pz}] = 1;

						if (simple3DLabelBg(labels, nodes, px, py, pz, numSlices, width,
							height, simpleDictionary3D, nodeIndex, true, bgConn)) {
							if (Label(labels[pz][px + (width * (py))]) == nodeIndex) {
								voxelsDeleted += 1;
								for (int j = 0; j < structCube.size(); j++) {
									nx = px + structCube[j][0]; ny = py + structCube[j][1]; nz = pz + structCube[j][2];
									if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
										if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
											if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0) {
												labels[pz][px + (width * (py))] = labels[nz][nx + (width * (ny))];

												clear_vertex(vtxToGNode[p], miniBgG);
												miniBgValid[vtxToGNode[p]] = false;

												for (int a = 0; a < structCross3D.size(); a++) {
													int nxi = px + structCross3D[a][0];
													int nyi = py + structCross3D[a][1];
													int nzi = pz + structCross3D[a][2];
													if (nxi >= 0 && nyi >= 0 && nzi >= 0 && nxi < width && nyi < height && nzi < numSlices) {
														if (Label(labels[pz][px + (width * (py))]) != Label(labels[nzi][nxi + (width * (nyi))])) {
															if (((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 0 ||
																((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 2 ||
																((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 3
																) {
																if (ccNeighbor(nodes, px, py, pz, nxi, nyi, nzi, labels, width, height, numSlices, g_Image3D, S)) {
																	if (edgeWt.find({ (int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))]) }) == edgeWt.end()) {
																		add_edge(Label(labels[nzi][nxi + (width * (nyi))]), Label(labels[pz][px + (width * (py))]), G);
																		if (abs(nxi - px) + abs(nyi - py) + abs(nzi - pz) == 1) {
																			if (
																				(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																				&&
																				(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CORE || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CUT || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == FILL)

																				) {
																				add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), fgGWithFills);

																			}
																			if (
																				(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
																				&&
																				(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CORE)

																				) {
																				add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), coreG);

																			}
																		}
																		if (
																			(((int)nodes[Label(labels[pz][px + py * width])].type) == N || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																			&&
																			(((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == N || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == CUT || ((int)nodes[Label(labels[nzi][nxi + nyi * width])].type) == FILL)

																			) {
																			add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), bgGWithCuts);

																		}
																		if (
																			(((int)nodes[Label(labels[pz][px + py * width])].type) == N)
																			&&
																			(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

																			) {
																			add_edge(Label(labels[pz][px + py * width]), Label(labels[nzi][nxi + nyi * width]), nG);

																		}

																		if (abs(nxi - px) + abs(nyi - py) + abs(nzi - pz) == 1) {
																			edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))])}] = 1;
																			edgeWt[{(int)Label(labels[nzi][nxi + (width * (nyi))]), (int)Label(labels[pz][px + (width * (py))])}] = 1;
																		}
																		else {
																			if (edgeWt.find({ (int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))]) }) == edgeWt.end()) {
																				edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))])}] = 0;
																				edgeWt[{(int)Label(labels[nzi][nxi + (width * (nyi))]), (int)Label(labels[pz][px + (width * (py))])}] = 0;
																			}
																		}
																	}
																	else {
																		if (abs(nxi - px) + abs(nyi - py) + abs(nzi - pz) == 1) {
																			edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nzi][nxi + (width * (nyi))])}] = 1;
																			edgeWt[{(int)Label(labels[nzi][nxi + (width * (nyi))]), (int)Label(labels[pz][px + (width * (py))])}] = 1;
																			if (edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nz][nx + (width * (ny))])}] == 0) {
																				if (
																					(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																					&&
																					(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

																					) {
																					add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);
																				}
																			}
																		}
																	}
																}


															}

															if (abs(nxi - px) + abs(nyi - py) + abs(nzi - pz) == 1) {
																if (((int)nodes[Label(labels[nzi][nxi + (width * (nyi))])].type) == 0) {
																	if (getParent(parentCompFg, nodes[Label(labels[pz][px + (width * (py))])].compIndex)
																		!= getParent(parentCompFg, nodes[Label(labels[nzi][nxi + (width * (nyi))])].compIndex)
																		) {
																		parentCompFg[getParent(parentCompFg, nodes[Label(labels[nzi][nxi + (width * (nyi))])].compIndex)] = getParent(parentCompFg, nodes[Label(labels[pz][px + (width * (py))])].compIndex);
																	}

																}
															}
														}
													}
												}

											}
										}
										if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
											Coordinate np(nx, ny, nz);
											if (genVisited[np] == 0) {
												q.push(np);
												genVisited[np] = 1;
											}
										}
									}
								}
							}
						}
						else {
							genVisited[p] = 0;
						}

					}

				}
				else {
					//if h0, h1, h2 total does not decrease, add back to graph
					Coordinate wCoord(wc.x, wc.y, wc.z);
					miniBgValid[vtxToGNode[wCoord]] = true;
					for (int j = 0; j < bgMask.size(); j++) {
						nx = wc.x + bgMask[j][0]; ny = wc.y + bgMask[j][1]; nz = wc.z + bgMask[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
								add_edge(vtxToGNode[Coordinate(nx, ny, nz)], vtxToGNode[wCoord], miniBgG);
							}
							else {
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 1
									|| ((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 3
									) {
									if (!boost::edge(vtxToGNode[wCoord], bgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexBg], miniBgG).second) {
										add_edge(vtxToGNode[wCoord], bgCompToMini[nodes[Label(labels[nz][nx + (width * (ny))])].overallCompIndexBg], miniBgG);
									}
								}

							}
						}
					}


					/**map<Coordinate, int> genVisited;
					std::vector<Coordinate> changedVoxels;
					for (int j = 0; j < P.size(); j++) {

						if (Label(labels[P[j].z][P[j].x + (width * (P[j].y))]) == nodeIndex) {
							Coordinate pNow(P[j].x, P[j].y, P[j].z);
							genVisited[pNow] = 0;

						}
					}
					queue<Coordinate> q;
					genVisited[wCoord] = 1;
					for (int j = 0; j < structCube.size(); j++) {
						nx = wCoord.x + structCube[j][0]; ny = wCoord.y + structCube[j][1]; nz = wCoord.z + structCube[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (abs(nx - wCoord.x) + abs(ny - wCoord.y) + abs(nz - wCoord.z) == 1) {
								if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0) {
									labels[wCoord.z][wCoord.x + (width * (wCoord.y))] = labels[nz][nx + (width * (ny))];

								}
							}
						}
					}
					//Iteratively remove simple voxels
					for (int j = 0; j < structCube.size(); j++) {
						nx = wCoord.x + structCube[j][0]; ny = wCoord.y + structCube[j][1]; nz = wCoord.z + structCube[j][2];
						if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
							if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
								q.push(Coordinate(nx, ny, nz));
							}
						}
					}

					changedVoxels.push_back(wCoord);
					checkVoxel[vtxToGNode[wCoord]] = false;

					while (!q.empty()) {
						Coordinate p = q.front();
						int px = p.x; int py = p.y; int pz = p.z;
						q.pop();
						genVisited[p] = 1;
						if (simple3DLabelBg(labels, nodes, px, py, pz, numSlices, width, height, simpleDictionary3D, nodeIndex, true, bgConn)) {
							if ((int)Label(labels[pz][px + (width * (py))]) == nodeIndex) {
								changedVoxels.push_back(p);
								checkVoxel[vtxToGNode[p]] = false;
								//cout << "check voxel update 1 " << vtxToGNode[wCoord] << " " << vtxToGNode[p] << endl;
								for (int j = 0; j < structCube.size(); j++) {
									nx = px + structCube[j][0]; ny = py + structCube[j][1]; nz = pz + structCube[j][2];
									if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
										if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {


											if (((int)nodes[Label(labels[nz][nx + (width * (ny))])].type) == 0) {
												labels[pz][px + (width * (py))] = labels[nz][nx + (width * (ny))];

											}
										}
										if (Label(labels[nz][nx + (width * (ny))]) == nodeIndex) {
											Coordinate np(nx, ny, nz);
											if (genVisited[np] == 0) {
												q.push(np);
												genVisited[np] = 1;
											}
										}
									}
								}
							}
						}
						else {
							genVisited[p] = 0;
						}

					}
					for (int c = 0; c < changedVoxels.size(); c++) {
						Coordinate cv = changedVoxels[c];
						labels[cv.z][cv.x + (width * (cv.y))] = nodeIndex;
					}**/
				}
			}
		}

		if (voxelsDeleted >= P.size()) {
			nodes[nodeIndex].valid = false;
			auto neighbours = adjacent_vertices(nodeIndex, G);
			for (auto u : make_iterator_range(neighbours)) {
				edgeWt.erase({ nodeIndex, (int)u }); edgeWt.erase({ (int)u , nodeIndex });
			}
			clear_vertex(nodeIndex, G);
			for (int j = 0; j < P.size(); j++) {
				int px = P[j].x; int py = P[j].y; int pz = P[j].z;
				for (int k = 0; k < structCross3D.size(); k++) {
					int nx = px + structCross3D[k][0]; int ny = py + structCross3D[k][1]; int nz = pz + structCross3D[k][2];
					if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
						if (Label(labels[pz][px + py * width]) != Label(labels[nz][nx + ny * width])) {
							if (((int)nodes[Label(labels[nz][nx + ny * width])].type) != N) {
								if (edgeWt.find({ (int)Label(labels[pz][px + py * width]) , (int)Label(labels[nz][nx + ny * width]) }) == edgeWt.end()) {
									add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), G);
									if (
										(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
										&&
										(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

										) {
										add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

									}
									if (
										(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
										&&
										(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

										) {
										add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), coreG);

									}

									if (
										(((int)nodes[Label(labels[pz][px + py * width])].type) == N || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
										&&
										(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

										) {
										add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), bgGWithCuts);

									}
									if (
										(((int)nodes[Label(labels[pz][px + py * width])].type) == N)
										&&
										(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

										) {
										add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), nG);

									}
								}
								edgeWt[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
								edgeWt[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
								if (edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nz][nx + (width * (ny))])}] == 0) {
									if (
										(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
										&&
										(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

										) {
										add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);
										if (
											(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
											&&
											(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

											) {
											add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), coreG);

										}
									}
								}
							}
						}
					}
				}
			}

		}

		else {

			if (voxelsDeleted > 0) {
				map< std::vector<int>, int> edgeWtTemp = edgeWt;
				auto neighbours = adjacent_vertices(nodeIndex, G);
				for (auto u : make_iterator_range(neighbours)) {
					edgeWtTemp.erase({ nodeIndex, (int)u }); edgeWtTemp.erase({ (int)u , nodeIndex });
				}
				clear_vertex(nodeIndex, G);
				map<Coordinate, bool> visited;
				for (int j = 0; j < P.size(); j++) {
					visited[Coordinate(P[j].x, P[j].y, P[j].z)] = false;
				}
				int currentLabelIndex = nodeIndex;
				auto t1 = std::chrono::high_resolution_clock::now();
				for (int j = 0; j < P.size(); j++) {
					Coordinate pNow(P[j].x, P[j].y, P[j].z);
					if (!visited[pNow]) {
						if (Label(labels[P[j].z][P[j].x + (width * (P[j].y))]) == nodeIndex) {
							node n; n.labelCost = 0; n.floatCost = 0.0;  n.type = 2; n.inFg = 1; n.index = currentLabelIndex; n.intensity = 0.0;
							queue<Coordinate> q;
							q.push(pNow);
							if (currentLabelIndex != nodeIndex) {
								add_vertex(G);// add_vertex(fgG);
								nodes.push_back(n);

							}
							else {
								nodes[currentLabelIndex] = n;
							}
							std::vector<float> diffs;
							while (!q.empty()) {
								Coordinate p = q.front();
								q.pop();
								int px = p.x; int py = p.y; int pz = p.z;
								visited[p] = true;

								nodes[currentLabelIndex].floatCost += (gradient(px, py, pz, g_Image3D, width, height, numSlices));


								uint32_t label32 = currentLabelIndex;
								changeLabel(label32, labels[pz][px + py * width]);
								for (int k = 0; k < structCube.size(); k++) {
									nx = px + structCube[k][0]; ny = py + structCube[k][1]; nz = pz + structCube[k][2];
									if (nx >= 0 && ny >= 0 && nz >= 0 && nx < width && ny < height && nz < numSlices) {
										if (Label(labels[pz][px + py * width]) != Label(labels[nz][nx + ny * width])) {
											if (nodeIndex != Label(labels[nz][nx + ny * width])) {
												if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {

													if (((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE) {
														if (getParent(parentCompFg, nodes[Label(labels[pz][px + (width * (py))])].compIndex)
															!= getParent(parentCompFg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex)
															) {
															parentCompFg[getParent(parentCompFg, nodes[Label(labels[nz][nx + (width * (ny))])].compIndex)] = getParent(parentCompFg, nodes[Label(labels[pz][px + (width * (py))])].compIndex);


														}
													}

												}
												if (ccNeighbor(nodes, px, py, pz, nx, ny, nz, labels, width, height, numSlices, g_Image3D, S)) {
													if (edgeWtTemp.find({ (int)Label(labels[pz][px + py * width]),(int)Label(labels[nz][nx + ny * width]) }) == edgeWtTemp.end()) {
														add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), G);
														if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
															if (
																(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																&&
																(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

																) {
																add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);

															}
															if (
																(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
																&&
																(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

																) {
																add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), coreG);

															}
														}

														if (
															(((int)nodes[Label(labels[pz][px + py * width])].type) == N || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
															&&
															(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

															) {
															add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), bgGWithCuts);

														}
														if (
															(((int)nodes[Label(labels[pz][px + py * width])].type) == N)
															&&
															(((int)nodes[Label(labels[nz][nx + ny * width])].type) == N)

															) {
															add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), nG);

														}
														//}
														if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
															edgeWtTemp[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
															edgeWtTemp[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
														}
														else {
															edgeWtTemp[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 0;
															edgeWtTemp[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 0;

														}
													}
													else {
														if (edgeWt[{(int)Label(labels[pz][px + (width * (py))]), (int)Label(labels[nz][nx + (width * (ny))])}] == 0) {
															if (
																(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE || ((int)nodes[Label(labels[pz][px + py * width])].type) == CUT || ((int)nodes[Label(labels[pz][px + py * width])].type) == FILL)
																&&
																(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == CUT || ((int)nodes[Label(labels[nz][nx + ny * width])].type) == FILL)

																) {
																add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), fgGWithFills);
																if (
																	(((int)nodes[Label(labels[pz][px + py * width])].type) == CORE)
																	&&
																	(((int)nodes[Label(labels[nz][nx + ny * width])].type) == CORE)

																	) {
																	add_edge(Label(labels[pz][px + py * width]), Label(labels[nz][nx + ny * width]), coreG);

																}
															}
														}

														if (abs(nx - px) + abs(ny - py) + abs(nz - pz) == 1) {
															edgeWtTemp[{(int)Label(labels[pz][px + py * width]), (int)Label(labels[nz][nx + ny * width])}] = 1;
															edgeWtTemp[{(int)Label(labels[nz][nx + ny * width]), (int)Label(labels[pz][px + py * width])}] = 1;
														}
													}
												}
												//}


											}
										}
										if (Label(labels[nz][nx + ny * width]) == nodeIndex) {
											Coordinate np(nx, ny, nz);
											if (!visited[np]) {
												if (ccNeighbor(nodes, px, py, pz, nx, ny, nz, labels, width, height, numSlices, g_Image3D, S)) {
													q.push(np);
													visited[np] = true;
												}
											}
										}
									}
								}
							}
							if (currentLabelIndex == nodeIndex) {
								currentLabelIndex = nodes.size();
							}
							else {
								currentLabelIndex += 1;
							}
						}
					}
				}
				edgeWt = edgeWtTemp;
				auto t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = t2 - t1;
				addETime += elapsed.count();

			}
		}

	}
	auto t_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedt = t_end - t_start;
	otherTime += elapsedt.count();

}

void simplifyGenerators(std::vector<uint32_t*>& labels, const std::vector<float*>& g_Image3D, int numSlices, int width, int height, float S,
	std::vector<node>& nodes, Graph& G, grapht& fgG, grapht& fgGWithFills, grapht& bgGWithCuts, grapht& coreG, grapht& bgG, grapht& nG, map< std::vector<int>, int>& edgeWt,
	int fgConn, const std::vector<unsigned char>& simpleDictionary3D, int genMode, bool rootMorphoSingle, int costType) {
	auto before = std::chrono::high_resolution_clock::now();
	cout << "Begin DFS to find articulation nodes " << endl;

	map<int, map<int, int> > fgLocalArtConnectivity, bgLocalArtConnectivity;
	auto after = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = after - before;
	cout << "DFS time: " << elapsed.count() << endl;
	int n = nodes.size();
	int avgGenSize = 0;
	int numGens = 0;
	std::vector<bool> simplified(n, false); map<int, int> parentCompFg; map<int, int> parentCompBg;
	//std::vector<bool> nextToNeighbor(nodes.size(), false);
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].type == CORE) {
			parentCompFg[nodes[i].compIndex] = nodes[i].compIndex;

		}
		if (nodes[i].type == N) {
			parentCompBg[nodes[i].compIndex] = nodes[i].compIndex;
			if (rootMorphoSingle) {
				auto neighbours = adjacent_vertices(i, G);
				//for (auto u : make_iterator_range(neighbours)) {
					//nextToNeighbor[(int)u] = true;
				//}
			}

		}
	}

	std::vector< genVoxels > generatorVoxels;
	priority_queue<genVoxels, std::vector<genVoxels>, CompareGen> genQ;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				int l = Label(labels[k][i + j * width]);
				if (l < n) {
					if (((int)nodes[l].type) == 2 || ((int)nodes[l].type) == 3) {
						if (!simplified[l]) {

							//if ( (rootMorphoSingle && nextToNeighbor[l]) || !rootMorphoSingle  ) {
							std::vector<Coordinate> voxels;
							findGenVoxels(labels, voxels, Label(labels[k][i + j * width]), width, height, numSlices, i, j, k);
							genVoxels gen; gen.n = nodes[Label(labels[k][i + j * width])]; gen.coords = voxels;
							genQ.push(gen);

							generatorVoxels.push_back(gen);

							simplified[l] = true;
							voxels.clear();
							//}

						}
					}
				}
			}
		}
	}
	int numNodes = nodes.size();
	float timeT = 0.0;
	float minigraphT = 0.0; float addETime = 0.0; float otherTime = 0.0;
	std::cout << "Generators to simplify " << genQ.size() << std::endl;

	std::vector<bool> visited(nodes.size(), false);
	int timer = 0; int overallIndex = 0;


	while (!genQ.empty()) {
		genVoxels gen = genQ.top();

		genQ.pop();
		if (nodes[gen.n.index].valid) {
			if (((int)nodes[gen.n.index].type) == 2 || ((int)nodes[gen.n.index].type) == 3) {
				std::vector<int> fNeighbors; std::vector<int> bNeighbors;
				auto neighbours = adjacent_vertices(gen.n.index, G);
				for (auto u : make_iterator_range(neighbours)) {
					if (nodes[u].valid) {
						if (((int)nodes[u].type) == 0 || ((int)nodes[u].type) == 2 || ((int)nodes[u].type) == 3) {
							if (edgeWt[{(int)gen.n.index, (int)u}] == 1) {
								fNeighbors.push_back((int)u);
							}
						}
						if (((int)nodes[u].type) == 1 || ((int)nodes[u].type) == 2 || ((int)nodes[u].type) == 3) {
							bNeighbors.push_back((int)u);
						}
					}
				}
				if (bNeighbors.size() == 0) {
					nodes[gen.n.index].valid = false;
					int nIndex = 0;
					for (int j = 0; j < fNeighbors.size(); j++) {
						if ((int)nodes[fNeighbors[j]].type == CORE) {
							nIndex = fNeighbors[j];
							edgeWt.erase({ gen.n.index, (int)fNeighbors[j] }); edgeWt.erase({ (int)fNeighbors[j] , gen.n.index });

						}
					}
					for (int p = 0; p < gen.coords.size(); p++) {
						labels[gen.coords[p].z][gen.coords[p].x + (width * (gen.coords[p].y))] = nIndex;
					}
					clear_vertex(gen.n.index, fgGWithFills); clear_vertex(gen.n.index, bgGWithCuts); clear_vertex(gen.n.index, G);

					continue;
				}
				if (fNeighbors.size() == 0) {
					nodes[gen.n.index].valid = false;
					int nIndex = 0;
					for (int j = 0; j < bNeighbors.size(); j++) {
						if ((int)nodes[bNeighbors[j]].type == N) {
							nIndex = bNeighbors[j];
							edgeWt.erase({ gen.n.index, (int)bNeighbors[j] }); edgeWt.erase({ (int)bNeighbors[j] , gen.n.index });

						}
					}
					for (int p = 0; p < gen.coords.size(); p++) {
						labels[gen.coords[p].z][gen.coords[p].x + (width * (gen.coords[p].y))] = nIndex;
					}
					clear_vertex(gen.n.index, fgGWithFills); clear_vertex(gen.n.index, bgGWithCuts); clear_vertex(gen.n.index, G);

					continue;
				}


				while (num_vertices(G) < nodes.size()) {
					add_vertex(G);
				}
				while (num_vertices(fgGWithFills) < nodes.size()) {
					add_vertex(fgGWithFills);
				}
				while (num_vertices(bgGWithCuts) < nodes.size()) {
					add_vertex(bgGWithCuts);
				}
				while (num_vertices(fgG) < nodes.size()) {
					add_vertex(fgG);
				}
				while (num_vertices(bgG) < nodes.size()) {
					add_vertex(bgG);
				}
				auto t1 = std::chrono::high_resolution_clock::now();
				if (((int)nodes[gen.n.index].type) == 2) {

					clear_vertex(gen.n.index, bgGWithCuts);
					for (int i = 0; i < nodes.size(); i++) {
						nodes[i].overallCompIndexBg = -1;
					}
					findOverallComponents(bgGWithCuts, nodes, 0, gen.n.index);


				}
				else {
					//if (nodes[gen.n.index].isArticulate) {
					clear_vertex(gen.n.index, fgGWithFills);
					for (int i = 0; i < nodes.size(); i++) {
						nodes[i].overallCompIndexFg = -1;
					}
					findOverallComponents(fgGWithFills, nodes, 1, gen.n.index);
					//}
				}
				auto t2 = std::chrono::high_resolution_clock::now();


				std::chrono::duration<double> elapsed = t2 - t1;
				timeT += elapsed.count();


				for (int i = 0; i < nodes.size(); i++) {
					nodes[i].isArticulate = false;
					nodes[i].isNew = false;
					nodes[i].tin = 0;
					nodes[i].low = 0;
					nodes[i].compIndex = -1;
				}
				if (((int)nodes[gen.n.index].type) == CUT) {
					fgLocalArtConnectivity = findTermComponents(coreG, nodes, CORE);//fgG
				}
				else {
					bgLocalArtConnectivity = findTermComponents(nG, nodes, N); //bgG
				}
				for (int i = 0; i < nodes.size(); i++) {
					if (nodes[i].type == CORE) {
						parentCompFg[nodes[i].compIndex] = nodes[i].compIndex;

					}
					if (nodes[i].type == N) {
						parentCompBg[nodes[i].compIndex] = nodes[i].compIndex;
					}
				}



				int numNodesInit = nodes.size(); bool generatorChanged = false;

				while (num_vertices(G) < nodes.size()) {
					add_vertex(G);
				}
				while (num_vertices(fgGWithFills) < nodes.size()) {
					add_vertex(fgGWithFills);
				}
				while (num_vertices(bgGWithCuts) < nodes.size()) {
					add_vertex(bgGWithCuts);
				}
				while (num_vertices(fgG) < nodes.size()) {
					add_vertex(fgG);
				}
				while (num_vertices(bgG) < nodes.size()) {
					add_vertex(bgG);
				}
				simplifyGenerator(G, labels, g_Image3D, fgLocalArtConnectivity, bgLocalArtConnectivity, numSlices, width, height, fgConn, nodes, gen.n.index, edgeWt, S, simpleDictionary3D, gen.coords, genMode,
					parentCompFg, parentCompBg, generatorChanged, fgGWithFills, bgGWithCuts, coreG, nG, minigraphT, addETime, otherTime, costType);



				map< std::vector<int>, int> newEdgeWt;
				if (nodes[gen.n.index].valid) {
					for (int i = 0; i < fNeighbors.size(); i++) {

						if (edgeWt.find({ gen.n.index, fNeighbors[i] }) != edgeWt.end()) {
							if (edgeWt[{ gen.n.index, fNeighbors[i] }] == 1) {
								add_edge(gen.n.index, fNeighbors[i], fgGWithFills);
							}
						}
					}
					for (int i = 0; i < bNeighbors.size(); i++) {
						if (edgeWt.find({ gen.n.index, bNeighbors[i] }) != edgeWt.end()) {
							add_edge(gen.n.index, bNeighbors[i], bgGWithCuts);
						}
					}
				}

				if (generatorChanged) {
					for (int i = 0; i < nodes.size(); i++) {
						nodes[i].isArticulateFg = false;
						nodes[i].isNew = false;
						nodes[i].tin = 0;
						nodes[i].low = 0;
						nodes[i].compIndex = -1;
						nodes[i].overallCompIndexFg = -1;
					}

					findArticulationPoints(fgGWithFills, nodes, 1, -1);

				}

				/**if (generatorChanged) {
					visited = std::vector<bool>(nodes.size(), false);
					for (int i = 0; i < nodes.size(); i++) {
						nodes[i].isArticulate = false;
					}
					timer = 0; overallIndex = 0;
					for (int i = 0; i < nodes.size(); i++) {
						if (((int)nodes[i].inFg) == 1 || ((int)nodes[i].type) == FILL) {
							if (!visited[i]) {
								int timer = 0; int timer1 = 0;
								dfs(fgGWithFills, i, visited, timer, -1, overallIndex, nodes);

								overallIndex += 1;
							}
						}
					}

					for (int i = 0; i < nodes.size(); i++) {
						nodes[i].overallCompIndexFg = -1;
					}
					findOverallComponents(fgGWithFills, nodes, 1, -1);
				}**/

			}
		}
	}
	cout << "minigraphTime " << minigraphT << " other time " << otherTime << " time outside " << timeT << endl;
}

int main(int argc, char** argv)
{
	string inFile, outFile;
	vector<float> shapes = { 0 };
	float vizEps;
	int hypernodeSize = INFINITY;
	int productThresh = INFINITY;
	int globalSteinerTime = 8;
	int localSteinerTime = 4;
	auto start = std::chrono::high_resolution_clock::now();
	params.heurbbtime = 3;
	int bbTime = 3;
	bool help = false;
	std::vector<float*> g_Image3D;
	int numSlices, width, height;
	bool shapeTopo = false;
	int beamSize = 1;
	int oneConflict = 1;
	float epsilon = INFINITY;
	bool propagate = false;
	int geomCost = 0;
	int distMode = 0;
	bool outputBB = false;
	string bbFile = "";
	bool findTotalTopo = false;
	int inFileType = 0;
	int outFileType = 0;
	bool createLevels = false;
	vector<int> shapeIndices;
	parseArguments(argc, inFile, outFile, argv, shapes, vizEps, hypernodeSize, productThresh,
		globalSteinerTime, localSteinerTime, bbTime, help, shapeTopo, beamSize,
		oneConflict, epsilon, propagate, geomCost, distMode, outputBB,
		bbFile, findTotalTopo, inFileType, outFileType, createLevels, shapeIndices);
	params.heurbbtime = bbTime;
	if (help) {
		return 1;
	}
	int maxVal = 0; int minVal = INFINITY; vector<int> maxCoord;
	loadImages(g_Image3D, numSlices, inFile, width, height, maxVal, minVal, maxCoord, inFileType);
	cout << "Slices loaded: " << numSlices << endl;
	if (numSlices == 0) {
		cout << "No images loaded. Check file directory " << endl;
		return 1;
	}
	std::ifstream simpleFile("simpleDictionaryFull.bin", ios::binary);
	std::vector<unsigned char> simpleDictionary3D(std::istreambuf_iterator<char>(simpleFile), {});
	if (!simpleFile) {
		std::cout << "Simplicity file 1 not read properly, " << simpleFile.gcount() << " bytes read " << endl;
		return 0;
	}
	if (createLevels) {
		vector<float*> layeredGImg;
		for (int s = 0; s < numSlices; s++) {
			float* layer = new float[width * height];

			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					layer[i + j * width] = 0;
					for (int t = 0; t < shapes.size(); t++) {
						if (g_Image3D[s][i + j * width] > shapes[t]) {
							layer[i + j * width] = shapes[t] + 1;
						}
					}
				}
			}
			layeredGImg.push_back(layer);
		}
		g_Image3D = layeredGImg;

	}
	if (findTotalTopo) {
		int numMin = 0;
		int numMax = 0;
		int numSaddle1 = 0;
		int numSaddle2 = 0;
		vector<vector<int>> xMask = { {-1,0,0}, {1,0,0} };
		vector<vector<int>> yMask = { {0,-1,0}, {0,1,0} };
		vector<vector<int>> zMask = { {0,0,-1}, {0,0,1} };
		int criticalPts = 0;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int s = 0; s < numSlices; s++) {
					//look at all structCube neighbors (subset is structCross)
					bool allNeighborGt = true;
					bool allNeighborLt = true;
					for (int k = 0; k < structCube.size(); k++) {
						vector<int> neighbor = { i + structCube[k][0], j + structCube[k][1], s + structCube[k][2] };
						if (neighbor[0] >= 0 && neighbor[0] < width && neighbor[1] >= 0 && neighbor[1] < height
							&& neighbor[2] >= 0 && neighbor[2] < numSlices) {
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] <= g_Image3D[s][i + j * width]) {
								allNeighborGt = false;
							}
							if (abs(neighbor[0] - i) + abs(neighbor[1] - j) + abs(neighbor[2] - s) == 1) {
								if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] >= g_Image3D[s][i + j * width]) {
									allNeighborLt = false;
								}
							}

						}
					}
					if (allNeighborGt) {
						numMin++;
						cout << "min " << g_Image3D[s][i + j * width] << " " << i << " " << j << " " << s << endl;
					}
					if (allNeighborLt) {
						numMax++;
						cout << "max " << g_Image3D[s][i + j * width] << " " << i << " " << j << " " << s << endl;
					}

					bool maxinX = true;
					bool mininX = true;
					bool maxinY = true;
					bool mininY = true;
					bool maxinZ = true;
					bool mininZ = true;
					for (int k = 0; k < xMask.size(); k++) {
						vector<int> neighbor = { i + xMask[k][0], j + xMask[k][1], s + xMask[k][2] };
						if (neighbor[0] >= 0 && neighbor[0] < width) {
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] <= g_Image3D[s][i + j * width]) {
								mininX = false;
							}
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] >= g_Image3D[s][i + j * width]) {
								maxinX = false;
							}
						}
						else {
							maxinX = false;
							mininX = false;
						}
						neighbor = { i + yMask[k][0], j + yMask[k][1], s + yMask[k][2] };
						if (neighbor[1] >= 0 && neighbor[1] < height) {
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] <= g_Image3D[s][i + j * width]) {
								mininY = false;
							}
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] >= g_Image3D[s][i + j * width]) {
								maxinY = false;
							}
						}
						else {
							maxinY = false;
							mininY = false;
						}

						neighbor = { i + zMask[k][0], j + zMask[k][1], s + zMask[k][2] };
						if (neighbor[2] >= 0 && neighbor[2] < numSlices) {
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] <= g_Image3D[s][i + j * width]) {
								mininZ = false;
							}
							if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] >= g_Image3D[s][i + j * width]) {
								maxinZ = false;
							}
						}
						else {
							maxinZ = false;
							mininZ = false;
						}
					}

					if ((mininX && maxinY && maxinZ) || (maxinX && mininY && maxinZ) || (maxinX && maxinY && mininZ)
						) {
						//cout << "s1 " << g_Image3D[s][i + j * width] << endl;
						numSaddle1++;
					}
					if ((maxinX && mininY && mininZ) || (mininX && mininY && maxinZ) || (mininX && maxinY && mininZ)
						) {

						numSaddle2++;
					}

					int minX = max(i - 1, 0); int minY = max(j - 1, 0); int minZ = max(s - 1, 0);
					int maxX = min(i + 2, width); int maxY = min(j + 2, height); int maxZ = min(s + 2, numSlices);


					std::vector< std::vector< std::vector<int>>> cubeN(3, std::vector<std::vector<int>>(3, std::vector<int>(3, 0)));
					for (int ni = minX; ni < maxX; ni++) {
						for (int nj = minY; nj < maxY; nj++) {
							for (int nk = minZ; nk < maxZ; nk++) {
								if (g_Image3D[nk][ni + nj * width] >= g_Image3D[s][i + j * width]) {
									cubeN[i - ni + 1][j - nj + 1][s - nk + 1] = 1;
								}
							}
						}
					}

					if (!(((int)simpleDictionary3D[neighborhoodToIndex(cubeN)]) == 49)) {
						//criticalPts++;
						//if (i == 118 && j == 146 && s == 20) {
						cout.precision(32);
						//cout << "currenpt " << i << " " << j << " " << s << " " << g_Image3D[s][i + j * width] << endl;
						int equal = false;
						for (int ni = minX; ni < maxX; ni++) {
							for (int nj = minY; nj < maxY; nj++) {
								for (int nk = minZ; nk < maxZ; nk++) {
									if (g_Image3D[nk][ni + nj * width] == g_Image3D[s][i + j * width]) {
										if (i != ni || j != nj || s != nk) {
											equal = true;
										}
										//cout.precision(32);
										//cout << ni << " " << nj << " " << nk << " " << i - ni + 1 << " " << j - nj + 1 << " " << s - nk + 1 << " " << 1 << " " << g_Image3D[nk][ni + nj * width] << endl;
									}

								}
							}
						}
						if (!equal) {
							criticalPts++;
							cout.precision(32);
							cout << "currenpt " << i << " " << j << " " << s << " " << g_Image3D[s][i + j * width] << endl;
							for (int ni = minX; ni < maxX; ni++) {
								for (int nj = minY; nj < maxY; nj++) {
									for (int nk = minZ; nk < maxZ; nk++) {
										cout.precision(32);
										cout << ni << " " << nj << " " << nk << " " << i - ni + 1 << " " << j - nj + 1 << " " << s - nk + 1 << " " << g_Image3D[nk][ni + nj * width] << endl;
									}
								}
							}
						}
						//}
					}
				}
			}
		}
		cout << "Maxima: " << numMax << " Minima: " << numMin << " 1-saddles: " << numSaddle1 << " 2-saddles: " << numSaddle2 << " critical points: " << criticalPts << endl;
		return 0;
	}

	vector<vector<float>> nonZeroCoords;
	queue<vector<int>> distQ;
	//vector<vector<bool>> inQ;
	vector<vector<float>> dist;
	for (int s = 0; s < numSlices; s++) {
		vector<float> layer(width * height, 0);
		//inQ.push_back(pqLayer);
		dist.push_back(layer);
	}
	int maxDist = 0;
	if (distMode == 1 || distMode == 2) {

		for (int s = 0; s < numSlices; s++) {
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					if (g_Image3D[s][(i)+(j) * (width)] > 0) {
						//nonZeroCoords.push_back({ (float)i,(float)j,(float)s });
						for (int k = 0; k < structCross3D.size(); k++) {
							int nx = i + structCross3D[k][0]; int ny = j + structCross3D[k][1]; int nz = s + structCross3D[k][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
								if (g_Image3D[nz][(nx)+(ny) * (width)] == 0) {
									if (dist[nz][nx + ny * width] == 0) {
										distQ.push({ nx,ny,nz });
										dist[nz][nx + ny * width] = dist[s][i + j * width] + 1;
									}
								}
							}
						}
					}
				}
			}
		}
		cout << "found initial queue " << endl;
		int addedToQ = distQ.size();
		int ctE = 1;
		while (addedToQ > 0) {
			cout << "count " << ctE << " " << distQ.size() << endl;
			ctE++;
			queue<vector<int>> nextQ;
			while (!distQ.empty()) {
				vector<int> c = distQ.front();
				distQ.pop();
				for (int k = 0; k < structCross3D.size(); k++) {
					int nx = c[0] + structCross3D[k][0]; int ny = c[1] + structCross3D[k][1]; int nz = c[2] + structCross3D[k][2];
					if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
						if (g_Image3D[nz][(nx)+(ny) * (width)] == 0) {
							if (dist[nz][nx + ny * width] == 0) {
								nextQ.push({ nx,ny,nz });
								dist[nz][nx + ny * width] = dist[c[2]][c[0] + c[1] * width] + 1;
							}
						}

					}
				}
			}
			distQ = nextQ;
			addedToQ = distQ.size();
		}
		maxDist = ctE;

		for (int s = 0; s < numSlices; s++) {
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					if (g_Image3D[s][(i)+(j) * (width)] > shapes[shapes.size() - 1]) {
						//nonZeroCoords.push_back({ (float)i,(float)j,(float)s });
						for (int k = 0; k < structCross3D.size(); k++) {
							int nx = i + structCross3D[k][0]; int ny = j + structCross3D[k][1]; int nz = s + structCross3D[k][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
								if (g_Image3D[nz][(nx)+(ny) * (width)] <= shapes[shapes.size() - 1]) {//shapes[shapes.size()-1]
									if (dist[s][i + j * width] == 0) {
										distQ.push({ i,j,s });
										dist[s][i + j * width] = -1;
									}
								}
							}
						}
					}
				}
			}
		}
		cout << "initialized inner q " << endl;
		addedToQ = distQ.size();
		ctE = 1;
		while (addedToQ > 0) {
			cout << "count " << ctE << " " << distQ.size() << endl;
			ctE++;
			queue<vector<int>> nextQ;
			while (!distQ.empty()) {
				vector<int> c = distQ.front();
				distQ.pop();
				for (int k = 0; k < structCross3D.size(); k++) {
					int nx = c[0] + structCross3D[k][0]; int ny = c[1] + structCross3D[k][1]; int nz = c[2] + structCross3D[k][2];
					if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
						if (g_Image3D[nz][(nx)+(ny) * (width)] > shapes[shapes.size() - 1]) {
							if (dist[nz][nx + ny * width] == 0) {
								nextQ.push({ nx,ny,nz });
								dist[nz][nx + ny * width] = dist[c[2]][c[0] + c[1] * width] - 1;
							}
						}

					}
				}
			}
			distQ = nextQ;
			addedToQ = distQ.size();
		}
		cout << "found initial queue " << endl;
	}
	std::vector<float*> g_Image3DOrig = g_Image3D;
	int origWidth = width;
	int origHeight = height;
	int numSlicesOrig = numSlices;
	width += 2;
	height += 2;
	numSlices += 2;
	g_Image3D.clear();
	float maxDim = 2 * max(max(width, height), numSlices);
	vector<float> maxC = { (float)maxCoord[0], (float)maxCoord[1], (float)maxCoord[2] };
	for (int s = 0; s < numSlices; s++) {
		float* layer = new float[width * height];
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				if (s <= 0 || s >= numSlices - 1 || i <= 0 || i >= width - 1 || j <= 0 || j >= height - 1) {
					if (distMode == 1 || distMode == 2) {
						layer[i + j * width] = -maxDist;
					}
					else {
						layer[i + j * width] = 0;
					}
					//cout << s << " " << i << " " << j << " " << numSlices << " " << width << " " << height << " " << origWidth * origHeight << endl;

				}
				else {
					//cout << s << " " << i << " " << j << " " << numSlices << " " << width << " " << height << " " << origWidth*origHeight << endl;
					if (g_Image3DOrig[s - 1][(i - 1) + (j - 1) * (origWidth)] == 0) {
						vector<float> coord = { (float)i, (float)j, (float)s };
						if (distMode == 1 || distMode == 2) {
							/**float minDist = 100000000;
							for (int n = 0; n < nonZeroCoords.size(); n++) {
								float dist = euclideanDistance(nonZeroCoords[n], coord);
								if (dist < minDist) {
									minDist = dist;
								}
							}**/
							//layer[i + j * width] = 1 + ((maxDim - minDist) / maxDim);
							layer[i + j * width] = -dist[s - 1][(i - 1) + (j - 1) * (origWidth)] + 1; // ((maxDim - dist[s - 1][(i - 1) + (j - 1) * (origWidth)]) / maxDim);
						}
						else {
							if (distMode == 0) {
								layer[i + j * width] = ((maxDim - euclideanDistance(maxC, coord)) / maxDim);
							}
						}
					}
					else {
						//if (g_Image3DOrig[sg_Image3D[maxCoord[2]][(maxCoord[0]) + (maxCoord[1]) * (width)]  - 1][(i - 1) + (j - 1) * (origWidth)] > shapes[shapes.size() - 1]) {
							//layer[i + j * width] = g_Image3DOrig[s - 1][(i - 1) + (j - 1) * (origWidth)]+1 -dist[s - 1][(i - 1) + (j - 1) * (origWidth)];
							//cout << "set to  " << layer[i + j * width] << endl;
						//}
						//else {
						layer[i + j * width] = g_Image3DOrig[s - 1][(i - 1) + (j - 1) * (origWidth)] + 1;
						//}
					}
				}
			}
		}
		g_Image3D.push_back(layer);
	}
	//g_Image3D[maxCoord[2] + 1][(maxCoord[0] + 1) + (maxCoord[1] + 1) * (width)] = maxVal+5;


	for (int i = 1; i < shapes.size(); i++) {
		shapes[i] = shapes[i] + 1;
	}

	int coreCt = 0;
	std::vector<vector<int> > coreLabels;
	for (int i = 0; i < g_Image3D.size(); i++) {
		coreLabels.push_back(vector<int>(width * height, 0));
	}
	int maxCoreSize = 0;
	int maxCoreLabel = -1;
	if (distMode == 2 || distMode == 1) {
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				for (int z = 0; z < numSlices; z++) {
					if (g_Image3D[z][x + y * width] > shapes[shapes.size() - 1]) {
						if (coreLabels[z][x + y * width] == 0) {
							coreCt++;
							queue<vector<int>> cQ;
							cQ.push({ x,y,z });
							coreLabels[z][x + y * width] = coreCt;
							int size = 1;
							while (!cQ.empty()) {
								vector<int> v = cQ.front();
								cQ.pop();
								for (int i = 0; i < structCross3D.size(); i++) {
									vector<int> neighbor = { v[0] + structCross3D[i][0], v[1] + structCross3D[i][1], v[2] + structCross3D[i][2] };
									if (neighbor[0] >= 0 && neighbor[0] < width && neighbor[1] >= 0 && neighbor[1] < height && neighbor[2] >= 0 && neighbor[2] < numSlices) {
										if (g_Image3D[neighbor[2]][neighbor[0] + neighbor[1] * width] > shapes[shapes.size() - 1]) {
											if (coreLabels[neighbor[2]][neighbor[0] + neighbor[1] * width] == 0) {
												coreLabels[neighbor[2]][neighbor[0] + neighbor[1] * width] = coreCt;
												cQ.push(neighbor);

												size++;
											}
										}
									}
								}
							}
							if (size > maxCoreSize) {
								maxCoreSize = size;
								maxCoreLabel = coreLabels[z][x + y * width];
								maxCoord = { x,y,z };

							}
						}
					}
				}
			}
		}

		//cout << "before pt q" << endl;
		queue<vector<int>> pointQ;
		pointQ.push(maxCoord);
		map<vector<int>, bool> visitedPt;
		visitedPt[maxCoord] = true;
		while (!pointQ.empty()) {
			vector<int> p = pointQ.front();
			pointQ.pop();
			//cout << dist[p[2] - 1][(p[0] - 1) + (p[1] - 1) * origWidth] << " " << dist[maxCoord[2] - 1][(maxCoord[0] - 1) + (maxCoord[1] - 1) * origWidth] << endl;
			if (dist[p[2] - 1][(p[0] - 1) + (p[1] - 1) * origWidth] < dist[maxCoord[2] - 1][(maxCoord[0] - 1) + (maxCoord[1] - 1) * origWidth]) {
				maxCoord = p;
			}
			for (int i = 0; i < structCross3D.size(); i++) {
				vector<int> neighbor = { p[0] + structCross3D[i][0], p[1] + structCross3D[i][1], p[2] + structCross3D[i][2] };
				if (neighbor[0] > 0 && neighbor[0] < width - 1 && neighbor[1] > 0 && neighbor[1] < height - 1 && neighbor[2] > 0 && neighbor[2] < numSlices - 1) {
					if (g_Image3DOrig[neighbor[2] - 1][(neighbor[0] - 1) + (neighbor[1] - 1) * (origWidth)] > 0) {
						if (visitedPt.find(neighbor) == visitedPt.end()) {
							visitedPt[neighbor] = true;
							pointQ.push(neighbor);
						}
					}
				}
			}
		}


		cout << "max core size " << maxCoreSize << " " << g_Image3D[maxCoord[2]][(maxCoord[0]) + (maxCoord[1]) * (width)] << " " << maxCoord[0] << " " << maxCoord[1] << " " << maxCoord[2] << " " << maxVal << endl;
	}
	//maxCoord = { 38,60,60 };
	/**if (distMode == 1) {
		maxCoord[0]++;
		maxCoord[1]++;
		maxCoord[2]++;
	}**/
	cout << "maxCoord " << maxCoord[0] << " " << maxCoord[1] << " " << maxCoord[2] <<  " " << g_Image3D[maxCoord[2]][(maxCoord[0]) + (maxCoord[1]) * (width)]  << endl;
	g_Image3D[maxCoord[2]][(maxCoord[0]) + (maxCoord[1]) * (width)] = maxVal + 5;
	
	if (outputBB) {
		g_Image3D[maxCoord[2]][(maxCoord[0]) + (maxCoord[1]) * (width)] = maxVal + 2;
		if (inFileType == 1 || outFileType == 1) {
			uint16 spp, bpp, photo;
			TIFF* out;
			int i, j;
			uint16 page;

			out = TIFFOpen(bbFile.c_str(), "w");
			if (!out)
			{
				fprintf(stderr, "Can't open %s for writing\n", argv[1]);
				return 1;
			}
			spp = 1;
			bpp = 32;
			photo = PHOTOMETRIC_MINISBLACK;

			for (int s = 0; s < numSlices; s++) {
				TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
				TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
				TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
				TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

				TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

				TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
				TIFFSetField(out, TIFFTAG_PAGENUMBER, s, numSlices);

				//float* dataV = g_Image3D[s];
				float* data = new float[width * height]();
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (i == 0 || j == 0 || s == 0 || i == width - 1 || j == height - 1 || s == numSlices - 1) {
							data[i + j * width] = 0;
						}
						else {
							if (g_Image3DOrig[s - 1][(i - 1) + (j - 1) * (origWidth)] == 0) {
								data[i + j * width] = 1;
							}
							else {
								data[i + j * width] = (int)g_Image3D[s][i + j * width];
							}
						}
						//data[i + j * width] = (float)data[i + j * width];
					}
				}

				for (j = 0; j < height; j++) {
					TIFFWriteScanline(out, &data[j * width], j, 0);
				}
				TIFFWriteDirectory(out);

			}
			TIFFClose(out);
		}
		else {
			for (int s = 0; s < numSlices; s++) {
				int digits = numDigits(s);
				string numStr = "";
				for (int n = 0; n < 4 - digits; n++) {
					numStr += "0";

				}
				numStr += std::to_string(s);
				string filename = bbFile + numStr + ".png";
				uint8_t* label8 = new uint8_t[width * height]();
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (i == 0 || j == 0 || s == 0 || i == width - 1 || j == height - 1 || s == numSlices - 1) {
							label8[i + j * width] = 0;
						}
						else {
							if (g_Image3DOrig[s - 1][(i - 1) + (j - 1) * (origWidth)] == 0) {
								label8[i + j * width] = 1;
							}
							else {
								label8[i + j * width] = (int)g_Image3D[s][i + j * width];
							}
						}
					}
				}
				int wrote = stbi_write_png(filename.c_str(), width, height, 1, label8, width);
			}
		}
		return 0;
	}
	cout << "maxVal " << maxVal << endl;
	shapes.push_back(maxVal + 2);

	if (distMode == 1 || distMode == 2) {
		shapes[0] = -maxDist;
	}

	for (int i = 0; i < shapes.size(); i++) {
		cout << "shape " << i << " is " << shapes[i] << endl;
	}

	//For each level, partition into kernel, cuts, fills, neighborhood
	std::vector< std::vector<uint32_t*> > levelLabels;
	for (int l = 0; l < shapes.size(); l++) {
		std::vector<uint32_t*> labels;
		for (int s = 0; s < numSlices; s++) {
			uint32_t* labelSlice = new uint32_t[width * height]();
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					labelSlice[i + j * width] = unvisited;
				}
			}
			labels.push_back(labelSlice);
		}
		levelLabels.push_back(labels);
	}

	int labelCt = -1;

	//Flood fill core and neighborhood
	std::cout << "Flood filling kernel and neighborhood components in one volume sweep " << endl;
	priority_queue<weightedCoord> corePQ, nPQ;

	float kernel = shapes[shapes.size() - 1];
	float shape = shapes[shapes.size() - 2];
	//flood fill core and neighborhood components in one sweep through volume
	std::vector< std::vector<node> > levelNodes;
	std::vector<node> nodes;
	labelCt += 1; //Found new component, increment new node index
	node n; n.type = 0; n.inFg = 1; n.index = labelCt;
	nodes.push_back(n);
	//((int)g_Image3D[maxCoord[2]][(maxCoord[0]) + (maxCoord[1]) * width]) - 1

	if (true) {
		vector<float*> distC;
		vector<bool*> visited;
		for (int s = 0; s < numSlices; s++) {
			float* layer = new float[width * height];
			bool* layerf = new bool[width * height];
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					layer[i + j * width] = 0;
					//visited[i + j * width] = false;
				}
			}
			distC.push_back(layer);
		}
		int step = 1;
		for (int s = 0; s < numSlices; s++) {
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) { //from shape of current level to kernel of last level
					if (g_Image3D[s][(i)+(j) * (width)] <= shape) {
						for (int k = 0; k < structCross3D.size(); k++) {
							int nx = i + structCross3D[k][0]; int ny = j + structCross3D[k][1]; int nz = s + structCross3D[k][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
								if (g_Image3D[nz][(nx)+(ny) * (width)] > shape) {
									if (distC[nz][nx + ny * width] == 0) {
										distQ.push({ nx,ny,nz });
										distC[nz][nx + ny * width] = step;
										step += 0.000001;
									}
								}
							}
						}
					}
				}
			}
		}
		int addedToQ = distQ.size();
		int ctE = 1;
		while (addedToQ > 0) {
			ctE++;
			queue<vector<int>> nextQ;
			while (!distQ.empty()) {
				vector<int> c = distQ.front();
				distQ.pop();
				for (int k = 0; k < structCross3D.size(); k++) {
					int nx = c[0] + structCross3D[k][0]; int ny = c[1] + structCross3D[k][1]; int nz = c[2] + structCross3D[k][2];
					if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
						//cout << "neighbor? " << endl;
						 //cout << g_Image3D[nz][(nx)+(ny) * (width)] << " " << distC[nz][(nx)+(ny) * (width)] << endl;
						if (g_Image3D[nz][(nx)+(ny) * (width)] > shape && distC[nz][(nx)+(ny) * (width)] == 0) {
							//cout << "label " << Label(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) << endl;
							if (Label(levelLabels[shapes.size() - 2][nz][nx + ny * width]) == unvisited) {
								nextQ.push({ nx,ny,nz });
								distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
							}
							else {
								if (((int)levelNodes[shapes.size() - 2][Label(levelLabels[shapes.size() - 2][nz][nx + ny * width])].type) != CORE) {
									nextQ.push({ nx,ny,nz });
									distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
								}
							}
						}

					}
				}
			}
			distQ = nextQ;
			addedToQ = distQ.size();
		}
		floodFillCore(g_Image3D, distC, distMode, levelLabels[levelLabels.size() - 1], numSlices, width, height, kernel, shape, maxCoord[0], maxCoord[1], maxCoord[2], labelCt, corePQ);
	}
	for (int i = 0; i < shapes.size(); i++) {
		levelNodes.push_back(nodes);
	}

	/**vector< uint8_t*> topoImg;
	for (int s = 0; s < numSlices; s++) {
		uint8_t* label8T = new uint8_t[width * height]();
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				uint32_t label = Label(levelLabels[levelLabels.size() - 1][s][i+j*width]);//levelNewToOldIndices[l][
				//if (((int)origLvlNodes[l][label].type) == 0 || ((int)origLvlNodes[l][label].type) == 2 || ((int)origLvlNodes[l][label].type) == 3) {
				if (label != unvisited) {
					if (nodes[label].type == CORE) {
						label8T[i + j * width] = 1;
					}
				}
				else {
					label8T[i + j * width] = 0;
				}
			}
		}
		topoImg.push_back(label8T);
	}
	std::vector<int> topoNums = getTopoFromBImg(topoImg, 1, width, height, numSlices);
	cout << "topo nums " << topoNums[0] << " " << topoNums[1] << " " << topoNums[2] << endl;**/
	/**for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				if (Label(levelLabels[levelLabels.size()-1][k][i + j * width]) == unvisited) {

					//Flood fill core component
					if ((float)g_Image3D[k][i + j * width] > kernel) {
						cout << "found core " << endl;
						labelCt += 1; //Found new component, increment new node index
						node n; n.type = 0; n.inFg = 1; n.index = labelCt;
						nodes.push_back(n);
						floodFillCore(g_Image3D, levelLabels[levelLabels.size()-1], numSlices, width, height, kernel, shape, i, j, k, labelCt, corePQ);
					}
				}
			}
		}
	}**/
	for (int s = 0; s < numSlices; s++) {
		float* layer = new float[width * height];
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				if (s > 0 && s < numSlices - 1 && i > 0 && i < width - 1 && j > 0 && j < height - 1) {
					if (g_Image3D[s][i + j * width] == -maxDist) {
						cout << "ne " << i << " " << j << " " << s << " " << g_Image3D[s][i + j * width] << endl;
					}
				}
			}
		}
	}
	cout << "done flood filling " << endl;
	for (int l = 0; l < shapes.size() - 1; l++) {
		float kernel = shapes[shapes.size() - l - 1];
		float shape = shapes[shapes.size() - l - 2];
		float nextShape = shapes[shapes.size() - l - 3];

		if (l == shapes.size() - 2) {
			nextShape = shape;
		}
		//levelLabels[shapes.size() - l - 2] = levelLabels[shapes.size() - l - 1];
		for (int s = 0; s < numSlices; s++) {
			uint32_t* labelSlice = new uint32_t[width * height]();
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					if (Label(levelLabels[shapes.size() - l - 1][s][i + j * width]) != unvisited) {
						int label = levelLabels[shapes.size() - l - 1][s][i + j * width];
						levelLabels[shapes.size() - l - 2][s][i + j * width] = label;
					}

				}
			}

		}

		priority_queue<weightedCoord> nextPQ;
		if (distMode == 2) {
			while (!distQ.empty()) {
				distQ.pop();
			}
			vector<float*> distC;
			vector<bool*> visited;
			for (int s = 0; s < numSlices; s++) {
				float* layer = new float[width * height];
				bool* layerf = new bool[width * height];
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						layer[i + j * width] = 0;
						//visited[i + j * width] = false;
					}
				}
				distC.push_back(layer);
			}
			int step = 1;
			for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) { //from shape of current level to kernel of last level
						if (g_Image3D[s][(i)+(j) * (width)] <= shape) {
							for (int k = 0; k < structCross3D.size(); k++) {
								int nx = i + structCross3D[k][0]; int ny = j + structCross3D[k][1]; int nz = s + structCross3D[k][2];
								if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
									if (g_Image3D[nz][(nx)+(ny) * (width)] > shape) {
										if (distC[nz][nx + ny * width] == 0) {
											distQ.push({ nx,ny,nz });
											distC[nz][nx + ny * width] = step;
											step += 0.000001;
										}
									}
								}
							}
						}
					}
				}
			}
			int addedToQ = distQ.size();
			int ctE = 1;
			while (addedToQ > 0) {
				ctE++;
				queue<vector<int>> nextQ;
				while (!distQ.empty()) {
					vector<int> c = distQ.front();
					distQ.pop();
					for (int k = 0; k < structCross3D.size(); k++) {
						int nx = c[0] + structCross3D[k][0]; int ny = c[1] + structCross3D[k][1]; int nz = c[2] + structCross3D[k][2];
						if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
							//cout << "neighbor? " << endl;
							 //cout << g_Image3D[nz][(nx)+(ny) * (width)] << " " << distC[nz][(nx)+(ny) * (width)] << endl;
							if (g_Image3D[nz][(nx)+(ny) * (width)] > shape && distC[nz][(nx)+(ny) * (width)] == 0) {
								//cout << "label " << Label(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) << endl;
								if (Label(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) == unvisited) {
									nextQ.push({ nx,ny,nz });
									distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
								}
								else {
									if (((int)levelNodes[shapes.size() - l - 2][Label(levelLabels[shapes.size() - l - 2][nz][nx + ny * width])].type) != CORE) {
										nextQ.push({ nx,ny,nz });
										distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
									}
								}
							}

						}
					}
				}
				distQ = nextQ;
				addedToQ = distQ.size();
			}
			/**
			for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (Visited(levelLabels[shapes.size() - l - 2][s][i + j * width])) {
							if (Label(levelLabels[shapes.size() - l - 2][s][i + j * width]) == unvisited) {
								setVisitedFlag(levelLabels[shapes.size() - l - 2][s][i + j * width], 0);
								cout << "visited set to not " << endl;
							}
						}
					}
				}
			}**/
			/**for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (Label(levelLabels[shapes.size() - l - 2][s][i + j * width]) != unvisited) {
							if (((int)levelNodes[shapes.size() - l - 2][Label(levelLabels[shapes.size() - l - 2][s][i + j * width])].type) == CORE) {
								//are neighbors not core and greater than shape?
								for (int k = 0; k < structCube.size(); k++) {
									int nx = i + structCube[k][0]; int ny = j + structCube[k][1]; int nz = s + structCube[k][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (g_Image3D[nz][nx + ny * width] > shape && Visited(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) == false) {//Visited(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) == false &&
											weightedCoord wc = { nx, ny, nz, ((float)(distC[nz][nx + width * ny])) };//((float)(distC[nz][nx + width * ny]))
											corePQ.push(wc);
											setVisitedFlag(levelLabels[shapes.size() - l - 2][nz][nx + ny * width], 1);
										}
									}
								}
							}
						}
					}
				}
			}**/
			cout << "corepq size " << corePQ.size() << endl;
			inflateCoreMulti(levelLabels[shapes.size() - l - 2], g_Image3D, g_Image3D, levelNodes[shapes.size() - l - 2], corePQ, nextPQ, numSlices, width, height, kernel, shape, nextShape, simpleDictionary3D, distMode);
			corePQ = nextPQ;
		}
		else {
			for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (Label(levelLabels[shapes.size() - l - 2][s][i + j * width]) != unvisited) {
							if (((int)levelNodes[shapes.size() - l - 2][Label(levelLabels[shapes.size() - l - 2][s][i + j * width])].type) != CORE) {
								setVisitedFlag(levelLabels[shapes.size() - l - 2][s][i + j * width], 0);
							}
						}
						else {
							setVisitedFlag(levelLabels[shapes.size() - l - 2][s][i + j * width], 0);
						}
					}
				}
			}
			for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (Label(levelLabels[shapes.size() - l - 2][s][i + j * width]) != unvisited) {
							if (((int)levelNodes[shapes.size() - l - 2][Label(levelLabels[shapes.size() - l - 2][s][i + j * width])].type) == CORE) {
								//are neighbors not core and greater than shape?
								for (int k = 0; k < structCross3D.size(); k++) {
									int nx = i + structCross3D[k][0]; int ny = j + structCross3D[k][1]; int nz = s + structCross3D[k][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (g_Image3D[nz][nx + ny * width] > shape && Visited(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) == false) {//Visited(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) == false && 
											weightedCoord wc = { nx, ny, nz, ((float)(g_Image3D[nz][nx + width * ny] - shape)) };
											corePQ.push(wc);
											setVisitedFlag(levelLabels[shapes.size() - l - 2][nz][nx + ny * width], 1);
										}
									}
								}
							}
						}
					}
				}
			}
			inflateCoreMulti(levelLabels[shapes.size() - l - 2], g_Image3D, g_Image3D, levelNodes[shapes.size() - l - 2], corePQ, nextPQ, numSlices, width, height, kernel, shape, nextShape, simpleDictionary3D, distMode);
			//corePQ = nextPQ;
		}
	}
	cout << "done inflating cores " << endl;

	if (true) {
		vector<float*> distC;
		for (int s = 0; s < numSlices; s++) {
			float* layer = new float[width * height];
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					layer[i + j * width] = 0;
				}
			}
			distC.push_back(layer);
		}

		for (int s = 0; s < numSlices; s++) {
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) { //from shape of current level to kernel of last level
					if (g_Image3D[s][(i)+(j) * (width)] > shape) {
						for (int k = 0; k < structCross3D.size(); k++) {
							int nx = i + structCross3D[k][0]; int ny = j + structCross3D[k][1]; int nz = s + structCross3D[k][2];
							if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
								if (g_Image3D[nz][(nx)+(ny) * (width)] <= shape) {
									if (distC[nz][nx + ny * width] == 0) {
										distQ.push({ nx,ny,nz });
										distC[nz][nx + ny * width] = 1;
									}
								}
							}
						}
					}
				}
			}
		}
		int addedToQ = distQ.size();
		int ctE = 1;
		while (addedToQ > 0) {
			ctE++;
			queue<vector<int>> nextQ;
			while (!distQ.empty()) {
				vector<int> c = distQ.front();
				distQ.pop();
				for (int k = 0; k < structCross3D.size(); k++) {
					int nx = c[0] + structCross3D[k][0]; int ny = c[1] + structCross3D[k][1]; int nz = c[2] + structCross3D[k][2];
					if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
						if (g_Image3D[nz][(nx)+(ny) * (width)] <= shape && distC[nz][(nx)+(ny) * (width)] == 0) {
							if (Label(levelLabels[1][nz][nx + ny * width]) == unvisited) {
								nextQ.push({ nx,ny,nz });
								distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
							}
							else {
								if (((int)levelNodes[1][Label(levelLabels[1][nz][nx + ny * width])].type) != N &&
									((int)levelNodes[1][Label(levelLabels[1][nz][nx + ny * width])].type) != CORE
									) {
									nextQ.push({ nx,ny,nz });
									distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
								}
							}
						}

					}
				}
			}
			distQ = nextQ;
			addedToQ = distQ.size();
		}

		cout << "begin flood filling neighborhood " << endl;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int k = 0; k < numSlices; k++) {
					if (Label(levelLabels[0][k][i + j * width]) == unvisited) {
						//if (k == 0) {

						//}
						//Flood fill neighborhood component
						if ((float)g_Image3D[k][i + j * width] <= shapes[0]) {
							labelCt += 1; //Found new component, increment new node index
							node n; n.type = 1; n.inFg = 0; n.index = labelCt;
							//cout << "new n shapes0 " << shapes[0] << " " << shapes[1] << endl;
							floodFillNeighborhood(g_Image3D, distC, distMode, levelLabels[0], numSlices, width, height, shapes[1], shapes[0], i, j, k, labelCt, nPQ, levelNodes[0]);
							levelNodes[0].push_back(n);

						}
					}
					if (Label(levelLabels[0][k][i + j * width]) != unvisited) {
						if (levelNodes[0][Label(levelLabels[0][k][i + j * width])].type == N) {
							if (g_Image3D[k][i + j * width] != -maxDist) {
								cout << "wrong n " << g_Image3D[k][i + j * width] << endl;
							}
						}
					}
					/**
					if (k <= 1 || k >= numSlices - 1 || i <= 1 || i >= width - 1 ||j <= 1 || j >= height - 1) {
						if (Label(levelLabels[0][k][i + j * width]) == unvisited) {
							cout << "unvisited n" << endl;
						}
					}
					else {
						if (Label(levelLabels[0][k][i + j * width]) != unvisited) {
							if (levelNodes[0][Label(levelLabels[0][k][i + j * width])].type == N) {
								cout << "N wrong place " << endl;
							}
						}
					}**/
				}
			}
		}
	}



	for (int i = 0; i < levelNodes.size(); i++) {
		levelNodes[i] = levelNodes[0];
	}
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			for (int z = 0; z < numSlices; z++) {
				if (Label(levelLabels[0][z][x + y * width]) != unvisited) {
					if (levelNodes[0][Label(levelLabels[0][z][x + y * width])].type == 1) {
						uint32_t label = levelLabels[0][z][x + y * width];
						levelLabels[0][z][x + y * width] = label;
					}
					else {
						setVisitedFlag(levelLabels[0][z][x + y * width], 0);
					}
				}
				else {
					setVisitedFlag(levelLabels[0][z][x + y * width], 0);
				}
			}
		}
	}
	cout << "finished flood filling neighborhood " << nPQ.size() << " " << levelNodes[0].size() << endl;
	cout << "begin deflating " << endl;
	for (int l = 0; l < shapes.size() - 1; l++) {
		float neighborhood = shapes[l];
		float shape = shapes[l + 1];
		float nextShape = shape;
		if (l != shapes.size() - 2) {
			nextShape = shapes[l + 2];
		}
		//levelLabels[l+1] = levelLabels[l];
		int neighborhoodAssigned = 0;
		int unvisitedCt = 0;
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				for (int z = 0; z < numSlices; z++) {

					if (Label(levelLabels[l][z][x + y * width]) != unvisited) {
						unvisitedCt++;
						if (levelNodes[l][Label(levelLabels[l][z][x + y * width])].type == 1) {
							uint32_t label = levelLabels[l][z][x + y * width];
							levelLabels[l + 1][z][x + y * width] = label;
							neighborhoodAssigned++;
						}
						else {
							setVisitedFlag(levelLabels[l + 1][z][x + y * width], 0);
						}
					}
					else {
						setVisitedFlag(levelLabels[l + 1][z][x + y * width], 0);
					}
				}
			}
		}

		priority_queue<weightedCoord> nextPQ;
		if (distMode == 2) {
			while (!distQ.empty()) {
				distQ.pop();
			}
			vector<float*> distC;
			for (int s = 0; s < numSlices; s++) {
				float* layer = new float[width * height];
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						layer[i + j * width] = 0;
					}
				}
				distC.push_back(layer);
			}

			for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) { //from shape of current level to kernel of last level
						if (g_Image3D[s][(i)+(j) * (width)] > shape) {
							for (int k = 0; k < structCube.size(); k++) {
								int nx = i + structCube[k][0]; int ny = j + structCube[k][1]; int nz = s + structCube[k][2];
								if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
									if (g_Image3D[nz][(nx)+(ny) * (width)] <= shape) {
										if (distC[nz][nx + ny * width] == 0) {
											distQ.push({ nx,ny,nz });
											distC[nz][nx + ny * width] = 1;
										}
									}
								}
							}
						}
					}
				}
			}
			int addedToQ = distQ.size();
			int ctE = 1;
			while (addedToQ > 0) {
				ctE++;
				queue<vector<int>> nextQ;
				while (!distQ.empty()) {
					vector<int> c = distQ.front();
					distQ.pop();
					for (int k = 0; k < structCube.size(); k++) {
						int nx = c[0] + structCube[k][0]; int ny = c[1] + structCube[k][1]; int nz = c[2] + structCube[k][2];
						if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
							//cout << "neighbor? " << endl;
							 //cout << g_Image3D[nz][(nx)+(ny) * (width)] << " " << distC[nz][(nx)+(ny) * (width)] << endl;
							if (g_Image3D[nz][(nx)+(ny) * (width)] <= shape && distC[nz][(nx)+(ny) * (width)] == 0) {
								//cout << "label " << Label(levelLabels[shapes.size() - l - 2][nz][nx + ny * width]) << endl;
								if (Label(levelLabels[l + 1][nz][nx + ny * width]) == unvisited) {
									nextQ.push({ nx,ny,nz });
									distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
								}
								else {
									if (((int)levelNodes[l + 1][Label(levelLabels[l + 1][nz][nx + ny * width])].type) != N &&
										((int)levelNodes[l + 1][Label(levelLabels[l + 1][nz][nx + ny * width])].type) != CORE
										) {
										nextQ.push({ nx,ny,nz });
										distC[nz][nx + ny * width] = distC[c[2]][c[0] + c[1] * width] + 1;
									}
								}
							}

						}
					}
				}
				distQ = nextQ;
				addedToQ = distQ.size();
				//cout << "addedToQ " << addedToQ << endl;

			}
			/**for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (Label(levelLabels[l + 1][s][i + j * width]) != unvisited) {
							if (((int)levelNodes[l + 1][Label(levelLabels[l + 1][s][i + j * width])].type) == N) {
								//are neighbors not N and less than shape?
								for (int k = 0; k < structCube.size(); k++) {
									int nx = i + structCube[k][0]; int ny = j + structCube[k][1]; int nz = s + structCube[k][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (Visited(levelLabels[l + 1][nz][nx + ny * width]) == false && g_Image3D[nz][nx + ny * width] <= shape) {
											weightedCoord wc = { nx, ny, nz, ((float)distC[nz][nx + width * ny]) };
											nPQ.push(wc);
											setVisitedFlag(levelLabels[l + 1][nz][nx + ny * width], 1);
										}
									}
								}
							}
						}
					}
				}
			}**/
			cout << neighborhood << " " << shape << " " << nextShape << " nPQ size " << nPQ.size() << " " << neighborhoodAssigned << " " << unvisitedCt << endl;

			deflateNeighborhoodMulti(levelLabels[l + 1], g_Image3D, g_Image3D, levelNodes[l + 1], nPQ, nextPQ, numSlices, width, height, neighborhood, shape, nextShape, simpleDictionary3D, distMode);
			//while (!nPQ.empty()) {
			///	nPQ.pop();
			//}
			nPQ = nextPQ;
		}
		else {
			for (int s = 0; s < numSlices; s++) {
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						if (Label(levelLabels[l + 1][s][i + j * width]) != unvisited) {
							if (((int)levelNodes[l + 1][Label(levelLabels[l + 1][s][i + j * width])].type) == N) {
								//are neighbors not N and less than shape?
								for (int k = 0; k < structCube.size(); k++) {
									int nx = i + structCube[k][0]; int ny = j + structCube[k][1]; int nz = s + structCube[k][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (Visited(levelLabels[l + 1][nz][nx + ny * width]) == false && g_Image3D[nz][nx + ny * width] <= shape) {
											weightedCoord wc = { nx, ny, nz, ((float)g_Image3D[nz][nx + width * ny] - shape) };
											nPQ.push(wc);
											setVisitedFlag(levelLabels[l + 1][nz][nx + ny * width], 1);
										}
									}
								}
							}
						}
					}
				}
			}
			cout << neighborhood << " " << shape << " " << nextShape << " nPQ size " << nPQ.size() << " " << neighborhoodAssigned << " " << unvisitedCt << endl;
			deflateNeighborhoodMulti(levelLabels[l + 1], g_Image3D, g_Image3D, levelNodes[l + 1], nPQ, nextPQ, numSlices, width, height, neighborhood, shape, nextShape, simpleDictionary3D, distMode);

			//nPQ = nextPQ;
		}
	}
	cout << "finished deflating " << endl;
	/**vector<vector<bool>> visited;
	for (int k = 0; k < numSlices; k++) {
		vector<bool> layer = vector<bool>(width * height, false);
		visited.push_back(layer);
	}
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int k = 0; k < numSlices; k++) {
				if (Label(levelLabels[0][k][i + j * width]) == unvisited) {


					if (g_Image3D[k][i + j * width] > shapes[0]) {
						vector< uint8_t*> topoImg;
						for (int z = 0; z < numSlices; z++) {
							uint8_t* label8T = new uint8_t[width * height]();
							for (int x = 0; x < width; x++) {
								for (int y = 0; y < height; y++) {
									uint32_t label = Label(levelLabels[0][z][x + y * width]);//levelNewToOldIndices[l][
									if (label != unvisited) {
										if (levelNodes[0][label].type == CORE) {
											label8T[x + y * width] = 1;
										}
									}
									else {
										label8T[x + y * width] = 0;
									}
								}
							}
							topoImg.push_back(label8T);
						}
						std::vector<int> topoNums = getTopoFromBImg(topoImg, 1, width, height, numSlices);

						std::queue<vector<int>> q;
						q.push({ i,j,k });
						visited[k][i + j * width] = true;
						topoImg[k][i + j * width] = 1;
						vector<vector<int>> changed = { {i,j,k} };
						while (!q.empty()) {
							vector<int> top = q.front();
							q.pop();
							for (int c = 0; c < structCube.size(); c++) {
								int nx = top[0] + structCube[c][0]; int ny = top[1] + structCube[c][0]; int nz = top[2] + structCube[c][2];
								if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
									uint32_t nlabel = Label(levelLabels[0][nz][nx + ny * width]);//levelNewToOldIndices[l][


									if (!visited[nz][nx + ny * width] && nlabel == unvisited && g_Image3D[nz][nx + ny * width] > shape) {
										visited[nz][nx + ny * width] = true;
										topoImg[nz][nx + ny * width] = 1;
										changed.push_back({ nx,ny,nz });
										q.push({ nx,ny,nz });
									}
								}
							}
						}
						std::vector<int> topoNums1 = getTopoFromBImg(topoImg, 1, width, height, numSlices);
						cout << topoNums[0] << " " << topoNums[1] << " " << topoNums[2] << " " << topoNums1[0] << " " << topoNums1[1] << " " << topoNums1[2] << endl;
					}
				}
			}
		}
	}**/
	Graph intersectingGraph;
	map< vector<int>, bool> edgeExists;
	int origLabelCt = labelCt;

	std::vector< map<vector<int>, int> > levelEdgeWts;
	std::vector<Graph> levelGraphs;



	for (int s = 0; s < shapes.size(); s++) {
		levelGraphs.push_back(Graph());
		map<vector<int>, int> edgeWts;
		int numK = 0; int numN = 0; int numC = 0; int numF = 0;
		int numUnvisited = 0;
		cout << "shape: " << shapes[s] << endl;
		labelCt = origLabelCt;
		//build edges across levels
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int k = 0; k < numSlices; k++) {
					if (Label(levelLabels[s][k][i + j * width]) == unvisited) {
						numUnvisited++;
						if (g_Image3D[k][i + j * width] > shapes[s]) {
							if (s == 0) {
								for (int c = 0; c < structCube.size(); c++) {
									int nx = i + structCube[c][0]; int ny = j + structCube[c][1]; int nz = s + structCube[c][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (Label(levelLabels[s][nz][nx + ny * width]) != unvisited) {
											if (levelNodes[s][Label(levelLabels[s][nz][nx + ny * width])].type == CORE) {
												uint32_t nlabel = levelLabels[s][nz][nx + ny * width];
												levelLabels[s][k][i + j * width] = nlabel;
											}
										}
									}
								}
							}
							else {
								labelCt += 1;
								node n; n.type = 2; n.inFg = 1; n.index = labelCt; n.level = s;
								identifyCutFromPixel3D(levelLabels, s, g_Image3D, levelGraphs[s], i, j, k, width,
									height, numSlices, labelCt, n, shapes[s], levelNodes, edgeWts, geomCost);
								//cout << (int)n.type << " added " << labelCt << endl;
								levelNodes[s].push_back(n);
							}
						}
						else {
							if (s == shapes.size() - 1) {
								for (int c = 0; c < structCube.size(); c++) {
									int nx = i + structCube[c][0]; int ny = j + structCube[c][1]; int nz = s + structCube[c][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (Label(levelLabels[s][nz][nx + ny * width]) != unvisited) {
											if (levelNodes[s][Label(levelLabels[s][nz][nx + ny * width])].type == N) {
												uint32_t nlabel = levelLabels[s][nz][nx + ny * width];
												levelLabels[s][k][i + j * width] = nlabel;
											}
										}
									}
								}
							}

							else {
								labelCt += 1;
								node n; n.type = 3; n.inFg = 0; n.index = labelCt; n.level = s;
								identifyFillFromPixel3D(levelLabels, s, g_Image3D, levelGraphs[s], i, j, k,
									width, height, numSlices, labelCt, n, shapes[s], levelNodes, edgeWts, geomCost);
								//cout << (int)n.type << " added " << labelCt << endl;
								levelNodes[s].push_back(n);
							}
						}
					}
					/**
					else {
						if (s == 0) {
							if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == CUT) {
								for (int c = 0; c < structCube.size(); c++) {
									int nx = i + structCube[c][0]; int ny = j + structCube[c][1]; int nz = s + structCube[c][2];
									if (nx >= 0 && nx < width && ny >= 0 && ny < height && nz >= 0 && nz < numSlices) {
										if (levelNodes[s][Label(levelLabels[s+1][nz][nx + ny * width])].type == CORE) {
											cout << i << " " << j << " " << k << endl;
										}
									}
								}
							}
						}
					}**/
					if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == 0) {
						numK++;
					}
					if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == 1) {
						numN++;
					}
					if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == 2) {
						numC++;
					}
					if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == 3) {
						numF++;
					}
				}
			}
		}
		cout << "Level " << s << " has " << levelNodes[s].size() << " nodes, " << numUnvisited << " voxels, numK " << numK << " numN " << numN << " numF " << numF << " numC " << numC << endl;
		levelEdgeWts.push_back(edgeWts);
	}

	/**for (int l = 0; l < shapes.size(); l++) {
		grapht fgG = grapht(); grapht bgG = grapht(); grapht coreG = grapht(); grapht  nG = grapht(); grapht fgGWithFills = grapht(); grapht bgGWithCuts = grapht();
		int numCuts = 0; int numF = 0;
		for (int i = 0; i < levelNodes[l].size(); i++) {
			add_vertex(fgG); add_vertex(fgGWithFills); add_vertex(bgG); add_vertex(bgGWithCuts); add_vertex(coreG); add_vertex(nG);

			if (((int)levelNodes[l][i].type) == CUT) {
				numCuts++;
			}
			if (((int)levelNodes[l][i].type) == FILL) {
				numF++;
			}
		}
		cout << "cuts before: " << numCuts << " num fills " << numF << endl;
		typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
		edge_iter ei, ei_end;
		for (tie(ei, ei_end) = edges(levelGraphs[l]); ei != ei_end; ++ei) {
			int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
			if ((((int)levelNodes[l][v1].type) == 0 || ((int)levelNodes[l][v1].type) == 2 || ((int)levelNodes[l][v1].type) == 3) && (((int)levelNodes[l][v2].type) == 0 || ((int)levelNodes[l][v2].type) == 2 || ((int)levelNodes[l][v2].type) == 3) && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) { //
				//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
				if (levelEdgeWts[l][{(int)ei->m_source, (int)ei->m_target}] == 1) { //If fg connectivity is cube or fg connectivity is structCross3D and edge is strong
					add_edge(v1, v2, 0, fgGWithFills);
				}
			}
			if ((((int)levelNodes[l][v1].type) == 0 || ((int)levelNodes[l][v1].type) == 2) && (((int)levelNodes[l][v2].type) == 0 || ((int)levelNodes[l][v2].type) == 2) && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) { //
				add_edge(v1, v2, 0, fgG);
			}
			if (((int)levelNodes[l][v1].type) == CORE && ((int)levelNodes[l][v2].type) == CORE && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) {
				//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
				if (levelEdgeWts[l][{(int)ei->m_source, (int)ei->m_target}] == 1) { //If fg connectivity is cube or fg connectivity is structCross3D and edge is strong
					add_edge(v1, v2, 0, coreG);
				}
			}

			if (((int)levelNodes[l][v1].inFg) == 0 && ((int)levelNodes[l][v2].inFg) == 0 && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) {
				//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
				add_edge(v1, v2, 0, bgG);
			}
			if ((((int)levelNodes[l][v1].type) == 1 || ((int)levelNodes[l][v1].type) == 2 || ((int)levelNodes[l][v1].type) == 3) && (((int)levelNodes[l][v2].type) == 1 || ((int)levelNodes[l][v2].type) == 2 || ((int)levelNodes[l][v2].type) == 3) && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) { //
				add_edge(v1, v2, 0, bgGWithCuts);
			}

			if (((int)levelNodes[l][v1].type) == N && ((int)levelNodes[l][v2].type) == N && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) {
				//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
				add_edge(v1, v2, 0, nG);
			}
		}

		cout << "Begin simplifying generators " << l << endl;
		simplifyGenerators(levelLabels[l], g_Image3D, numSlices, width, height, shapes[l], levelNodes[l], levelGraphs[l], fgG, fgGWithFills, bgGWithCuts, coreG, bgG, nG, levelEdgeWts[l], 1, simpleDictionary3D, 0, false, geomCost);
		numCuts = 0; numF = 0;
		for (int i = 0; i < levelNodes[l].size(); i++) {

			if (((int)levelNodes[l][i].type) == CUT) {
				numCuts++;
			}
			if (((int)levelNodes[l][i].type) == FILL) {
				numF++;
			}
		}
		cout << "after cuts: " << numCuts << " num fills " << numF << endl;
	}

	for (int s = 0; s < shapes.size(); s++) {
		int numC = 0;
		int numF = 0;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int k = 0; k < numSlices; k++) {

					if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == 2) {
						numC++;
					}
					if (levelNodes[s][Label(levelLabels[s][k][i + j * width])].type == 3) {
						numF++;
					}
				}
			}
		}
		cout << numC << " " << numF << " " << s << endl;
	}
	for (int l = 0; l < shapes.size(); l++) {
		for (int s = 0; s < numSlices; s++) {
			int digits = numDigits(s);
			string numStr = "";
			for (int n = 0; n < 4 - digits; n++) {
				numStr += "0";

			}
			numStr += std::to_string(s);
			string filenameC = outFile + "cuts/" + to_string(l) + "/" + numStr + ".png";
			string filenameF = outFile + "fills/" + to_string(l) + "/" + numStr + ".png";
			string filenameK = outFile + "kernels/" + to_string(l) + "/" + numStr + ".png";
			uint8_t* label8C = new uint8_t[width * height]();
			uint8_t* label8F = new uint8_t[width * height]();
			uint8_t* label8K = new uint8_t[width * height]();
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					label8C[i + j * width] = 0;
					label8F[i + j * width] = 0;
					label8K[i + j * width] = 0;
					uint32_t label = Label(levelLabels[l][s][i + j * width]);
					if (label != unvisited) {

						if (((int)levelNodes[l][label].type) == CUT) {
							label8C[i + j * width] = 1;
						}
						if (((int)levelNodes[l][label].type) == FILL) {
							label8F[i + j * width] = 1;
						}
						if (((int)levelNodes[l][label].type) == CORE) {
							//if (l == 1) {
							//	cout << "kernel " << endl;
							//}
							label8K[i + j * width] = 1;
						}
					}
				}
			}
			int wrote = stbi_write_png(filenameC.c_str(), width, height, 1, label8C, width);
			wrote = stbi_write_png(filenameF.c_str(), width, height, 1, label8F, width);
			wrote = stbi_write_png(filenameK.c_str(), width, height, 1, label8K, width);
		}
	}**/
	for (int l = 0; l < shapes.size(); l++) {
		for (int i = 0; i < levelNodes[l].size(); i++) {
			if ((int)levelNodes[l][i].type == CUT) {
				if (abs(levelNodes[l][i].floatCost) > 100000) {
					cout << "c " << levelNodes[l][i].floatCost << endl;
					//levelNodes[l][i].type = 0; //K
				}

			}
			
			if ((int)levelNodes[l][i].type == FILL) {
				if (abs(levelNodes[l][i].floatCost) > 100000) {
					cout << "f " << levelNodes[l][i].floatCost << endl;
					//levelNodes[l][i].type = 1; //Neighborhood
				}
				
			}
		}
	}

	if (epsilon < INFINITY) {
		for (int l = 0; l < shapes.size(); l++) {
			for (int i = 0; i < levelNodes[l].size(); i++) {
				if ((int)levelNodes[l][i].type == CUT) {
					if (abs(levelNodes[l][i].greatestDiff) >= epsilon) {
						levelNodes[l][i].type = CORE;
					}
				}
				if ((int)levelNodes[l][i].type == FILL) {
					//cout << levelNodes[l][i].greatestDiff << endl;
					if (abs(levelNodes[l][i].greatestDiff) >= epsilon) {
						levelNodes[l][i].type = 1; //Neighborhood
					}
				}
			}
		}
	}

	auto cutFillTime = std::chrono::high_resolution_clock::now();
	//Graph intersectionG;, vector< vector<int> > overallToLvlNodes
	std::chrono::duration<double> elapsedCf = cutFillTime - start;

	vector<int> newShapeIndices = { 0 };
	for (int i = 0; i < shapeIndices.size(); i++) {
		newShapeIndices.push_back(shapeIndices[i]);
	}
	newShapeIndices.push_back(shapes.size() - 1);
	shapeIndices = newShapeIndices;

	vector<vector<node>> lvlNodesNew;
	vector<Graph> graphsNew;
	std::vector< map<vector<int>, int> > levelEdgeWtsNew;
	std::vector< std::vector<uint32_t*> > levelLabelsNew;
	vector<float> shapesNew;
	for (int i = 0; i < shapeIndices.size(); i++) {
		//cout << "shape index " << shapeIndices[i] << " " << shapes.size() << endl;
		lvlNodesNew.push_back(levelNodes[shapeIndices[i]]);
		graphsNew.push_back(levelGraphs[shapeIndices[i]]);
		levelEdgeWtsNew.push_back(levelEdgeWts[shapeIndices[i]]);
		levelLabelsNew.push_back(levelLabels[shapeIndices[i]]);
		shapesNew.push_back(shapes[shapeIndices[i]]);
	}
	levelNodes = lvlNodesNew;
	levelGraphs = graphsNew;
	levelEdgeWts = levelEdgeWtsNew;
	levelLabels = levelLabelsNew;
	shapes = shapesNew;

	std::vector< std::vector<node>> origLvlNodes;
	State initialState;
	initialState.cost = 0; initialState.totalWtSum = 0;
	typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
	edge_iter ei, ei_end;
	vector< vector<vector<int>> > levelNewToOldComps;
	vector<Graph> levelGraphsInit = levelGraphs;
	vector<map<vector<int>, int> > levelEdgeWtsIn = levelEdgeWts;
	int h0S = 0; int h1S = 0; int h2S = 0;
	std::cout.precision(25);
	/**
	for (int l = 0; l < shapes.size(); l++) {
		for (int s = 0; s < numSlices; s++) {
			int digits = numDigits(s);
			string numStr = "";
			for (int n = 0; n < 4 - digits; n++) {
				numStr += "0";

			}
			numStr += std::to_string(s);
			string filenameC = outFile + "cuts/" + to_string(l) + "/" + numStr + ".png";
			string filenameF = outFile + "fills/" + to_string(l) + "/" + numStr + ".png";
			uint8_t* label8C = new uint8_t[width * height]();
			uint8_t* label8F = new uint8_t[width * height]();
			uint8_t* label8K = new uint8_t[width * height]();
			int sp = s + 1;
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					uint32_t label = Label(levelLabels[l][s][i + j * width]);
					if (label != unvisited) {

						if (((int)levelNodes[l][label].type) == CUT) {
							label8C[i + j * width] = 1;
						}
						if (((int)levelNodes[l][label].type) == FILL) {
							label8F[i + j * width] = 1;
						}
					}
				}
				int wrote = stbi_write_png(filenameC.c_str(), width, height, 1, label8C, width);
				wrote = stbi_write_png(filenameF.c_str(), width, height, 1, label8F, width);
			}
		}
	}**/

	if (shapeTopo) {
		for (int l = 0; l < shapes.size(); l++) {
			std::vector<int> eulerS = getTopoFromGImg(g_Image3D, shapes[l], 1, width, height, numSlices);
			cout << "Original Shape Topology for level " << shapes[l] << ": Components: " << eulerS[0] << " Cavities : " << eulerS[1] << " Cycles : " << eulerS[2] << endl;
			if (l > 0 && l < shapes.size() - 1) {
				h0S += eulerS[0];
				h1S += eulerS[2];
				h2S += eulerS[1];
			}
		}
	}
	for (int l = 0; l < shapes.size(); l++) { //shapes.size()
		cout << "solving level " << l << endl;
		//run sig. asia 2020 paper on each level
		//cout << "before get euler nums " << endl;
		std::vector<int> eulerNums = getEulerNumbers(levelNodes[l], levelLabels[l], width, height, numSlices);
		/**if (true) {
			grapht fgG = grapht(); grapht bgG = grapht(); grapht coreG = grapht(); grapht  nG = grapht(); grapht fgGWithFills = grapht(); grapht bgGWithCuts = grapht();
			int numCuts = 0; int numF = 0;
			for (int i = 0; i < levelNodes[l].size(); i++) {
				add_vertex(fgG); add_vertex(fgGWithFills); add_vertex(bgG); add_vertex(bgGWithCuts); add_vertex(coreG); add_vertex(nG);

				if (((int)levelNodes[l][i].type) == CUT) {
					numCuts++;
				}
				if (((int)levelNodes[l][i].type) == FILL) {
					numF++;
				}
			}
			cout << "cuts before: " << numCuts << " num fills " << numF << endl;
			typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
			edge_iter ei, ei_end;
			for (tie(ei, ei_end) = edges(levelGraphs[l]); ei != ei_end; ++ei) {
				int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
				if ((((int)levelNodes[l][v1].type) == 0 || ((int)levelNodes[l][v1].type) == 2 || ((int)levelNodes[l][v1].type) == 3) && (((int)levelNodes[l][v2].type) == 0 || ((int)levelNodes[l][v2].type) == 2 || ((int)levelNodes[l][v2].type) == 3) && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) { //
					//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
					if (levelEdgeWts[l][{(int)ei->m_source, (int)ei->m_target}] == 1) { //If fg connectivity is cube or fg connectivity is structCross3D and edge is strong
						add_edge(v1, v2, 0, fgGWithFills);
					}
				}
				if ((((int)levelNodes[l][v1].type) == 0 || ((int)levelNodes[l][v1].type) == 2) && (((int)levelNodes[l][v2].type) == 0 || ((int)levelNodes[l][v2].type) == 2) && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) { //
					add_edge(v1, v2, 0, fgG);
				}
				if (((int)levelNodes[l][v1].type) == CORE && ((int)levelNodes[l][v2].type) == CORE && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) {
					//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
					if (levelEdgeWts[l][{(int)ei->m_source, (int)ei->m_target}] == 1) { //If fg connectivity is cube or fg connectivity is structCross3D and edge is strong
						add_edge(v1, v2, 0, coreG);
					}
				}

				if (((int)levelNodes[l][v1].inFg) == 0 && ((int)levelNodes[l][v2].inFg) == 0 && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) {
					//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
					add_edge(v1, v2, 0, bgG);
				}
				if ((((int)levelNodes[l][v1].type) == 1 || ((int)levelNodes[l][v1].type) == 2 || ((int)levelNodes[l][v1].type) == 3) && (((int)levelNodes[l][v2].type) == 1 || ((int)levelNodes[l][v2].type) == 2 || ((int)levelNodes[l][v2].type) == 3) && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) { //
					add_edge(v1, v2, 0, bgGWithCuts);
				}

				if (((int)levelNodes[l][v1].type) == N && ((int)levelNodes[l][v2].type) == N && (levelNodes[l][v1].valid) && (levelNodes[l][v2].valid)) {
					//If src vtx or tgt vtx is core, cut, or fill, add to fg graph
					add_edge(v1, v2, 0, nG);
				}
			}

			cout << "Begin simplifying generators " << l << endl;
			simplifyGenerators(levelLabels[l], g_Image3D, numSlices, width, height, shapes[l], levelNodes[l], levelGraphs[l], fgG, fgGWithFills, bgGWithCuts, coreG, bgG, nG, levelEdgeWts[l], 1, simpleDictionary3D, 0, false, geomCost);
			numCuts = 0; numF = 0;
			for (int i = 0; i < levelNodes[l].size(); i++) {

				if (((int)levelNodes[l][i].type) == CUT) {
					numCuts++;
				}
				if (((int)levelNodes[l][i].type) == FILL) {
					numF++;
				}
			}
			cout << "after cuts: " << numCuts << " num fills " << numF << endl;
		}**/

		//cout << "after euler nums " << endl;
		initialState.levelEulerNums.push_back(eulerNums);
		//Check how simplified generators look like
		int costType = 0;
		if (l == 0 || l == shapes.size() - 1) {
			for (int i = 0; i < levelNodes[l].size(); i++) {
				if (levelNodes[l][i].type == 2) {
					levelNodes[l][i].type = 0;
				}
				if (levelNodes[l][i].type == 3) {
					levelNodes[l][i].type = 1;
				}
			}
		}
		//if (l == 1) {
		//	levelNodes[l][38].type = 0;
		//}
		//eulerNums = getEulerNumbers(levelNodes[l], levelLabels[l], width, height, numSlices);

		Graph fgG; Graph bgG;
		while (num_vertices(fgG) < levelNodes[l].size()) {
			add_vertex(fgG);
		}
		while (num_vertices(bgG) < levelNodes[l].size()) {
			add_vertex(bgG);
		}
		int64_t sumNodeCost = 0;
		int64_t cfCt = 0;
		/**
		Graph fgGt = Graph();
		Graph bgGt = Graph();
		for (int n = 0; n < levelNodes[l].size(); n++) {
			add_vertex(fgGt);
			add_vertex(bgGt);
		}
		typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
		edge_iter ei, ei_end;
		for (tie(ei, ei_end) = edges(levelGraphsInit[l]); ei != ei_end; ++ei) {
			int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
			if (!levelNodes[l][v1].valid || !levelNodes[l][v2].valid) { continue; }
			if ((int)levelNodes[l][v1].inFg == 1 && (int)levelNodes[l][v2].inFg == 1) {
				if (levelEdgeWtsIn[l][{v1, v2}] == 1) {
					add_edge(v1, v2, fgGt);
				}
			}

			if ((int)levelNodes[l][v1].inFg == 0 && (int)levelNodes[l][v2].inFg == 0) {
				add_edge(v1, v2, bgGt);
			}
		}

		vector<vector<int>> fgCompsIn = getComponents(fgGt, levelNodes[l], 1);
		vector<vector<int>> bgCompsIn = getComponents(bgGt, levelNodes[l], 0);
		int h0In = fgCompsIn.size();
		int h2In = bgCompsIn.size() - 1;
		int h1In = h0In + h2In - (eulerNums[0] - eulerNums[1] + eulerNums[2] - eulerNums[3]);
		**/
		for (int i = 0; i < levelNodes[l].size(); i++) {
			if (levelNodes[l][i].valid) {
				sumNodeCost += abs((int)levelNodes[l][i].floatCost);
				cfCt++;
				/**
				if ((int)levelNodes[l][i].type == FILL || ((int)levelNodes[l][i].type) == CUT) {
					Graph fgGn = Graph();
					Graph bgGn = Graph();
					int origLabel = levelNodes[l][i].type;
					char origFg = levelNodes[l][i].inFg;

					if ((int)levelNodes[l][i].inFg == 0) {
						levelNodes[l][i].inFg = 1;
						eulerNums[0] += levelNodes[l][i].v;
						eulerNums[1] += levelNodes[l][i].e;
						eulerNums[2] += levelNodes[l][i].f;
						eulerNums[3] += levelNodes[l][i].c;
					}
					else {
						levelNodes[l][i].inFg = 0;
						eulerNums[0] -= levelNodes[l][i].v;
						eulerNums[1] -= levelNodes[l][i].e;
						eulerNums[2] -= levelNodes[l][i].f;
						eulerNums[3] -= levelNodes[l][i].c;
					}
					for (int n = 0; n < levelNodes[l].size(); n++) {
						add_vertex(fgGn);
						add_vertex(bgGn);
					}
					for (tie(ei, ei_end) = edges(levelGraphsInit[l]); ei != ei_end; ++ei) {
						int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
						if (!levelNodes[l][v1].valid || !levelNodes[l][v2].valid) { continue; }
						if ((int)levelNodes[l][v1].inFg == 1 && (int)levelNodes[l][v2].inFg == 1) {
							if (levelEdgeWtsIn[l][{v1, v2}] == 1) {
								add_edge(v1, v2, fgGn);
							}
						}

						if ((int)levelNodes[l][v1].inFg == 0 && (int)levelNodes[l][v2].inFg == 0) {
							add_edge(v1, v2, bgGn);
						}
					}

					vector<vector<int>> fgCompsN = getComponents(fgGn, levelNodes[l], 1);
					vector<vector<int>> bgCompsN = getComponents(bgGn, levelNodes[l], 0);
					int h0N = fgCompsN.size();
					int h2N = bgCompsN.size() - 1;
					int h1N = h0N + h2N - (eulerNums[0] - eulerNums[1] + eulerNums[2] - eulerNums[3]);
					if (h0N + h2N + h1N >= h0In + h2In + h1In) {
						cout << "gt " << l << " " << origLabel << " " << h0N << " " << h2N << " " << h1N << " " << h0In << " " << h2In << " " << h1In << " " << levelNodes[l][i].floatCost << endl;

						if ((int)origFg == 1) {
							levelNodes[l][i].type = CORE;
							levelNodes[l][i].inFg = 1;
							eulerNums[0] += levelNodes[l][i].v;
							eulerNums[1] += levelNodes[l][i].e;
							eulerNums[2] += levelNodes[l][i].f;
							eulerNums[3] += levelNodes[l][i].c;
						}
						else {
							levelNodes[l][i].type = N;
							levelNodes[l][i].inFg = 0;
							eulerNums[0] -= levelNodes[l][i].v;
							eulerNums[1] -= levelNodes[l][i].e;
							eulerNums[2] -= levelNodes[l][i].f;
							eulerNums[3] -= levelNodes[l][i].c;
						}
					}
					else {
						levelNodes[l][i].inFg = origFg;
						if ((int)origFg == 1) {
							eulerNums[0] += levelNodes[l][i].v;
							eulerNums[1] += levelNodes[l][i].e;
							eulerNums[2] += levelNodes[l][i].f;
							eulerNums[3] += levelNodes[l][i].c;
						}
						else {
							eulerNums[0] -= levelNodes[l][i].v;
							eulerNums[1] -= levelNodes[l][i].e;
							eulerNums[2] -= levelNodes[l][i].f;
							eulerNums[3] -= levelNodes[l][i].c;
						}
						//levelNodes[l][i].labelCost = levelNodes[l][i].floatCost;
						sumNodeCost += abs((int)levelNodes[l][i].floatCost);
						cfCt++;
					}
				}**/
			}
		}

		//cout << "before add cost " << endl;
		for (int i = 0; i < levelNodes[l].size(); i++) {
			if (levelNodes[l][i].valid) {
				if (((int)levelNodes[l][i].type) == 2 || ((int)levelNodes[l][i].type) == 3) {
					int64_t eulerSum = levelNodes[l][i].v - levelNodes[l][i].e + levelNodes[l][i].f - levelNodes[l][i].c;
					if (((int)levelNodes[l][i].type) == 2 && eulerSum)
						levelNodes[l][i].intensity = levelNodes[l][i].floatCost;
					int64_t topoCost = eulerSum;
					levelNodes[l][i].labelCost = ((topoCost * cfCt * sumNodeCost) / 0.5) + levelNodes[l][i].floatCost;

					//if (abs(levelNodes[l][i].floatCost) > 150) {
						//cout << "level " << l << " type " << (int)levelNodes[l][i].type << " " << (int)levelNodes[l][i].inFg << " label cost " << levelNodes[l][i].labelCost << " geom cost " << levelNodes[l][i].floatCost << endl;
						//cout << "label cost " << levelNodes[l][i].labelCost << endl;
					//}
				}
			}
		}
		//cout << "after add cost " << endl;



		//Assign labeling for isolated cut-fill components and critical articulation points
		std::vector<bool> fillerB; std::vector<int> fillerI; map<int, int> fillerMap; //These variables are filler variables here; this function takes more arguments during the global steiner tree stage.
		std::vector<node> origNodes = levelNodes[l];

		origLvlNodes.push_back(origNodes);

		preprocessGraph(levelGraphsInit[l], levelNodes[l], levelEdgeWtsIn[l], fillerB, fillerI, fillerMap, false);
		//cout << "done preprocess" << endl;
		std::vector<node> nodesBeforeMerging = nodes; map<int, int> oldToNew;
		//Merge adjacent core terminals
		std::vector< std::vector<int> > newToOldComp = mergeAdjacentTerminals(levelGraphsInit[l], levelNodes[l], levelEdgeWtsIn[l], oldToNew);
		Graph origG = levelGraphs[l];
		//Find FG and BG graphs, to be used to find connected components of C and N in local graph
		tbb::concurrent_vector< hyperNode > globalHypernodes;
		std::vector<node> nodesGlobal; //Represents nodes after local stage: some have been assigned to core and neighborhood
		//Local steiner tree stage on clusters of cuts and fills
		//cout << "Original graph size: Nodes: " << nodes.size() << ", Edges: " << num_edges(G) << endl;
		int maxLocalGraphNodes = 0; int maxLocalGraphEdges = 0;
		int64_t wtSum = getWtSum(levelNodes[l]);

		//cout << "before local solve, wtSum: " << wtSum << endl;
		solveLocalGraphs(levelGraphsInit[l], levelNodes[l], levelNodes[l].size(), levelEdgeWtsIn[l], hypernodeSize,
			wtSum, productThresh, globalHypernodes, localSteinerTime, bbTime,
			maxLocalGraphNodes, maxLocalGraphEdges
		);
		//cout << "Max size of local hypergraph: Nodes: " << maxLocalGraphNodes << " Edges: " << maxLocalGraphEdges << 
		//	" wtSum " << wtSum <<
		//	endl;

		removeCAndNEdges(levelGraphsInit[l], levelNodes[l]); map<int, int> oldToNew2;

		std::vector< std::vector<int> > newToOldComp2 = mergeAdjacentTerminals(levelGraphsInit[l], levelNodes[l], levelEdgeWtsIn[l], oldToNew2);
		std::vector< std::vector<int> > newToOldTemp;
		for (int i = 0; i < newToOldComp2.size(); i++) {
			std::vector<int> combinedNewToOld;
			for (int j = 0; j < newToOldComp2[i].size(); j++) {
				int oldIndex = newToOldComp2[i][j];
				for (int k = 0; k < newToOldComp[oldIndex].size(); k++) {
					combinedNewToOld.push_back(newToOldComp[oldIndex][k]);
				}
			}
			newToOldTemp.push_back(combinedNewToOld);
		}

		newToOldComp = newToOldTemp;
		levelNewToOldComps.push_back(newToOldComp);
		for (int i = 0; i < globalHypernodes.size(); i++) {
			std::vector<int> subnodes;
			for (int j = 0; j < globalHypernodes[i].doubleSubnodes.size(); j++) {
				subnodes.push_back(oldToNew2[globalHypernodes[i].doubleSubnodes[j]]);
			}
			std::sort(subnodes.begin(), subnodes.end());
			subnodes.erase(unique(subnodes.begin(), subnodes.end()), subnodes.end());
			globalHypernodes[i] = hyperNode(subnodes, HYPERNODE, globalHypernodes[i].getSide());
		}

		//Global steiner tree stage
		int nodesToFix = 1; int numNodes = levelNodes[l].size();
		//cout << "before solve global " << endl;
		solveGlobalGraph(levelNodes[l], numNodes, levelGraphsInit[l], origG, globalHypernodes, wtSum, levelEdgeWtsIn[l], hypernodeSize, productThresh, globalSteinerTime, localSteinerTime,
			nodesToFix, bbTime);
		//cout << "after solve global " << endl;
		Level lvl;
		lvl.nodes = levelNodes[l];
		lvl.wtSum = wtSum;
		lvl.cfCt = cfCt;

		lvl.origNodes = origNodes;

		for (int j = 0; j < levelNodes[l].size(); j++) {
			for (int k = 0; k < newToOldComp[j].size(); k++) {
				int oldIndex = newToOldComp[j][k];
				lvl.origNodes[oldIndex].inFg = levelNodes[l][j].inFg;
			}
		}
		//cout << "found orig indices " << endl;
		//lvl.newToOldComp = newToOldComp;
		double geometryCost = 0;
		for (int j = 0; j < lvl.origNodes.size(); j++) {

			if (((int)origLvlNodes[l][j].type == FILL && (int)lvl.origNodes[j].inFg == 1)) {
				eulerNums[0] += lvl.origNodes[j].v;
				eulerNums[1] += lvl.origNodes[j].e;
				eulerNums[2] += lvl.origNodes[j].f;
				eulerNums[3] += lvl.origNodes[j].c;
				geometryCost += abs(origLvlNodes[l][j].floatCost);
			}

			if ((int)origLvlNodes[l][j].type == CUT && (int)lvl.origNodes[j].inFg == 0) {
				eulerNums[0] -= lvl.origNodes[j].v;
				eulerNums[1] -= lvl.origNodes[j].e;
				eulerNums[2] -= lvl.origNodes[j].f;
				eulerNums[3] -= lvl.origNodes[j].c;
				geometryCost += abs(origLvlNodes[l][j].floatCost);
			}
		}

		fgG = Graph();
		bgG = Graph();
		for (int i = 0; i < levelNodes[l].size(); i++) {
			add_vertex(fgG);
			add_vertex(bgG);
		}
		//typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
		//edge_iter ei, ei_end;
		for (tie(ei, ei_end) = edges(levelGraphsInit[l]); ei != ei_end; ++ei) {
			int v1 = (int)ei->m_source; int v2 = (int)ei->m_target;
			if (!levelNodes[l][v1].valid || !levelNodes[l][v2].valid) { continue; }
			if ((int)levelNodes[l][v1].inFg == 1 && (int)levelNodes[l][v2].inFg == 1) {
				if (levelEdgeWtsIn[l][{v1, v2}] == 1) {
					add_edge(v1, v2, fgG);
				}
			}

			if ((int)levelNodes[l][v1].inFg == 0 && (int)levelNodes[l][v2].inFg == 0) {
				add_edge(v1, v2, bgG);
			}


		}

		//swapLabelsGreedy(levelGraphsInit[l], fgG, bgG, levelNodes[l], levelEdgeWtsIn[l], 1, 0);

		//find lexicographical cost for this level
		vector<vector<int>> fgComps = getComponents(fgG, levelNodes[l], 1);
		vector<vector<int>> bgComps = getComponents(bgG, levelNodes[l], 0);
		lvl.h0 = fgComps.size();
		lvl.h2 = bgComps.size() - 1;
		int h1 = fgComps.size() + lvl.h2 - (eulerNums[0] - eulerNums[1] + eulerNums[2] - eulerNums[3]);
		lvl.h1 = h1;
		lvl.geomCost = geometryCost;
		initialState.totalWtSum += wtSum;
		initialState.levels.push_back(lvl);

	}

	for (int l = 1; l < initialState.levels.size() - 1; l++) {
		initialState.levels[l].cost = (initialState.levels[l].h0 + initialState.levels[l].h2 + initialState.levels[l].h1) * initialState.totalWtSum + initialState.levels[l].geomCost;
		if (l > 0 || l < shapes.size() - 1) {
			initialState.cost += initialState.levels[l].cost;
		}
	}


	if (shapeTopo) {
		cout << "Total shape topology: Components: " << h0S << " Cavities: " << h2S << " Cycles: " << h1S << endl;
	}
	cout << "done solving init state " << endl;
	auto initStateTime = std::chrono::high_resolution_clock::now();
	//Graph intersectionG;, vector< vector<int> > overallToLvlNodes
	std::chrono::duration<double> elapsedInit = initStateTime - cutFillTime;
	int totalNodeCt = 0;
	vector< vector<int> > overallToLvlNodes;
	for (int i = 0; i < origLvlNodes.size(); i++) {
		for (int j = 0; j < origLvlNodes[i].size(); j++) {
			origLvlNodes[i][j].totalNodeIndex = totalNodeCt;
			overallToLvlNodes.push_back({ i,j });
			totalNodeCt++;
		}
	}
	initialState.levelNewToOldComps = levelNewToOldComps;
	//map to original nodes
	for (int l = 0; l < initialState.levels.size(); l++) {
		std::vector<node> nodes = initialState.levels[l].nodes;
		initialState.levels[l].origNodes = origLvlNodes[l];
		for (int j = 0; j < nodes.size(); j++) {
			for (int k = 0; k < initialState.levelNewToOldComps[l][j].size(); k++) {
				int oldIndex = initialState.levelNewToOldComps[l][j][k];
				initialState.levels[l].origNodes[oldIndex].inFg = nodes[j].inFg;
			}
		}
	}
	//always want to check conflicts using origLvlNodes, update origLvlNodes to core or N

	map< vector<int>, bool> overallEdgeExists;
	Graph intersectionG;
	//find graphs across levels
	for (int i = 0; i < totalNodeCt; i++) {
		add_vertex(intersectionG);
	}

	for (int l = 1; l < shapes.size(); l++) {
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int k = 0; k < numSlices; k++) {

					int lowerIndex = Label(levelLabels[l - 1][k][i + j * width]);
					int upperIndex = Label(levelLabels[l][k][i + j * width]);
					if ((int)origLvlNodes[l - 1][lowerIndex].type == CUT || (int)origLvlNodes[l - 1][lowerIndex].type == FILL) {
						if ((int)origLvlNodes[l][upperIndex].type == CUT || (int)origLvlNodes[l][upperIndex].type == FILL) {
							int o1 = origLvlNodes[l - 1][lowerIndex].totalNodeIndex;
							int o2 = origLvlNodes[l][upperIndex].totalNodeIndex;
							if (overallEdgeExists.find({ o1, o2 }) == overallEdgeExists.end()) {
								overallEdgeExists[{o1, o2}] = true;
								overallEdgeExists[{o2, o1}] = true;
								add_edge(o1, o2, intersectionG);
							}
						}

					}
				}
			}
		}
	}

	double initialCost = 0;
	double initGeomCost1 = 0;
	for (int l = 1; l < origLvlNodes.size() - 1; l++) {
		double topoCost = initialState.levels[l].h0 + initialState.levels[l].h2 + initialState.levels[l].h1;
		initialCost += (topoCost * initialState.totalWtSum);

		float geomCost = 0;
		for (int i = 0; i < origLvlNodes[l].size(); i++) {
			if (origLvlNodes[l][i].type == CUT) {
				if (initialState.levels[l].origNodes[i].inFg == 0) {
					geomCost += abs(initialState.levels[l].origNodes[i].floatCost);
				}
			}
			if (origLvlNodes[l][i].type == FILL) {
				if (initialState.levels[l].origNodes[i].inFg == 1) {
					geomCost += abs(initialState.levels[l].origNodes[i].floatCost);
				}
			}
		}
		initGeomCost1 += geomCost;
		//cout << "init cost of level " << l << " is " << initialState.levels[l].h0 << " " << initialState.levels[l].h2 << " " << initialState.levels[l].h1 << " geom " << geomCost << " " << (topoCost * initialState.totalWtSum) <<  endl;
		initialCost += geomCost;
	}
	//cout << "initCost " << initialCost << endl;
	auto intersectionTime = std::chrono::high_resolution_clock::now();
	//cout << "find int graph " << endl;
	findContradictions(initialState, overallToLvlNodes, intersectionG, origLvlNodes);
	cout << "found contradictions " << initialState.contradictions.size() << endl;
	priority_queue<State, vector<State>, CompareState> stateQ;
	stateQ.push(initialState);
	State finalState = initialState;
	int ct = 0;
	map<std::vector<int>, bool> stateExists;
	stateExists[getStateEncoding(initialState)] = true;
	int statesPopped = 0;
	//Produce initial state contradictions
	map< vector<int>, int> lowerContradictions;
	map< vector<int>, int> upperContradictions;
	for (int i = 0; i < initialState.contradictions.size(); i++) {
		lowerContradictions[{initialState.contradictions[i].l1, initialState.contradictions[i].n1}] = 1;
		upperContradictions[{initialState.contradictions[i].l2, initialState.contradictions[i].n2}] = 1;
		//cout << "contra lvl " << top.contradictions[i].l1 << " " << top.contradictions[i].n1 << " " << top.contradictions[i].l2 << " " << top.contradictions[i].n2 << " " << endl;
	}
	int levelSolve = shapes.size() - 2;
	double bestCost = INFINITY;
	double geomCostInit = 0.0;
	for (int i = 0; i < initialState.levels.size(); i++) {
		geomCostInit += initialState.levels[i].geomCost;
	}

	if (propagate) {
		int bestStart = -1;
		for (int l = 1; l < shapes.size() - 1; l++) { //shapes.size() - 1
			cout << "fixing level " << l << endl;
			//if (l > 10) {
				//continue;
			//}
			//Fix all decisions at current level
			State start = initialState;
			for (int i = 1; i < shapes.size() - 1; i++) {
				for (int j = 0; j < origLvlNodes[i].size(); j++) {
					start.levels[i].origNodes[j].type = origLvlNodes[i][j].type;
				}
			}
			queue<int> ql;
			ql.push(l);
			//map<int, bool> lvlCovered;
			vector<int> lvlCovered(shapes.size(), false);
			int lower = l - 1;
			int upper = l + 1;
			while (lower >= 1 && upper <= shapes.size() - 2) {
				ql.push(lower);
				ql.push(upper);

				lower--;
				upper++;
			}
			while (lower >= 1) {
				ql.push(lower);
				lower--;
			}

			while (upper <= shapes.size() - 2) {
				ql.push(upper);
				upper++;
			}

			while (!ql.empty()) {
				int lvl = ql.front();
				start.contradictions.clear();
				findContradictions(start, overallToLvlNodes, intersectionG, origLvlNodes);
				map<vector<int>, bool> inContra;
				map<int, bool> inContraLvl;
				for (int i = 0; i < start.contradictions.size(); i++) {
					for (int j = 0; j < start.contradictions[i].lowerIndices.size(); j++) {
						int lowerLvl = overallToLvlNodes[start.contradictions[i].lowerIndices[j]][0];
						int lowerNode = overallToLvlNodes[start.contradictions[i].lowerIndices[j]][1];
						inContra[{lowerLvl, lowerNode}] = true;
					}
					for (int j = 0; j < start.contradictions[i].upperIndices.size(); j++) {
						int upperLvl = overallToLvlNodes[start.contradictions[i].upperIndices[j]][0];
						int upperNode = overallToLvlNodes[start.contradictions[i].upperIndices[j]][1];
						inContra[{upperLvl, upperNode}] = true;
					}
				}
				ql.pop();
				cout << "fixing sublevel " << lvl << endl;
				lvlCovered[lvl] = true;
				for (int i = 0; i < start.levels[lvl].origNodes.size(); i++) {
					if (start.levels[lvl].origNodes[i].inFg == 1) {
						if (start.levels[lvl].origNodes[i].type == CUT || start.levels[lvl].origNodes[i].type == FILL) {
							//propagate changes; everything below must be in FG

							vector<bool> visited(totalNodeCt + 1, false);
							queue<int> q;
							q.push(start.levels[lvl].origNodes[i].totalNodeIndex);
							visited[start.levels[lvl].origNodes[i].totalNodeIndex] = true;
							while (!q.empty()) { //search in lower levels
								int top = q.front();
								int cLvl = overallToLvlNodes[(int)top][0];
								int cIndex = overallToLvlNodes[(int)top][1];
								q.pop();
								if (inContra.find({ cLvl, cIndex }) != inContra.end()) {
									inContraLvl[cLvl] = true;
								}
								if (start.levels[cLvl].origNodes[cIndex].type == FILL) {
									auto neighboursp = adjacent_vertices(cIndex, levelGraphs[cLvl]);
									for (auto up : make_iterator_range(neighboursp)) {
										if (origLvlNodes[cLvl][(int)up].type == CUT) {
											start.levels[cLvl].origNodes[(int)up].inFg = 1;
											start.levels[cLvl].origNodes[(int)up].type = 0;

											auto neighboursLower = adjacent_vertices(start.levels[cLvl].origNodes[(int)up].totalNodeIndex, intersectionG);
											for (auto a : make_iterator_range(neighboursLower)) {
												int nLvlLower = overallToLvlNodes[(int)a][0];
												int nIndexLower = overallToLvlNodes[(int)a][1];
												if (nLvlLower < cLvl) {
													if ((start.levels[nLvlLower].origNodes[nIndexLower].type == CUT || start.levels[nLvlLower].origNodes[nIndexLower].type == FILL)
														&& !visited[(int)a]) {
														q.push((int)a);
														//cout << "lower pushed onto queue" << endl;
														visited[(int)a] = true;
													}
												}
											}
										}
									}
								}

								start.levels[cLvl].origNodes[cIndex].inFg = 1;
								start.levels[cLvl].origNodes[cIndex].type = 0;
								auto neighboursLower = adjacent_vertices(top, intersectionG);
								for (auto a : make_iterator_range(neighboursLower)) {
									int nLvlLower = overallToLvlNodes[(int)a][0];
									int nIndexLower = overallToLvlNodes[(int)a][1];
									if (nLvlLower < cLvl && !visited[(int)a]) {
										if ((start.levels[nLvlLower].origNodes[nIndexLower].type == CUT || start.levels[nLvlLower].origNodes[nIndexLower].type == FILL)
											&& !visited[(int)a]) {
											q.push((int)a);
											//cout << "lower pushed onto queue" << endl;
											visited[(int)a] = true;
										}
									}
								}

							}
						}


					}
					else {
						if (start.levels[lvl].origNodes[i].type == CUT || start.levels[lvl].origNodes[i].type == FILL) {
							//propagate changes; everything above must be in BG
							vector<bool> visited(totalNodeCt + 1, false);
							queue<int> q;
							q.push(start.levels[lvl].origNodes[i].totalNodeIndex);
							visited[start.levels[lvl].origNodes[i].totalNodeIndex] = true;
							while (!q.empty()) { //search in lower levels
								int top = q.front();
								int cLvl = overallToLvlNodes[(int)top][0];
								int cIndex = overallToLvlNodes[(int)top][1];
								q.pop();

								if (inContra.find({ cLvl, cIndex }) != inContra.end()) {
									inContraLvl[cLvl] = true;
								}
								if (start.levels[cLvl].origNodes[cIndex].type == CUT) {
									auto neighboursp = adjacent_vertices(cIndex, levelGraphs[cLvl]);
									for (auto up : make_iterator_range(neighboursp)) {
										if (origLvlNodes[cLvl][(int)up].type == FILL) {
											start.levels[cLvl].origNodes[(int)up].inFg = 0;
											start.levels[cLvl].origNodes[(int)up].type = 1;


											auto neighboursUpper = adjacent_vertices(start.levels[cLvl].origNodes[(int)up].totalNodeIndex, intersectionG);
											for (auto a : make_iterator_range(neighboursUpper)) {
												int nLvlUpper = overallToLvlNodes[(int)a][0];
												int nIndexUpper = overallToLvlNodes[(int)a][1];
												if (nLvlUpper > cLvl) {
													if ((start.levels[nLvlUpper].origNodes[nIndexUpper].type == CUT || start.levels[nLvlUpper].origNodes[nIndexUpper].type == FILL)
														&& !visited[(int)a]
														) {
														q.push((int)a);
														//cout << "upper pushed onto queue" << endl;
														visited[(int)a] = true;
													}
												}
											}
										}
									}
								}

								start.levels[cLvl].origNodes[cIndex].inFg = 0;
								start.levels[cLvl].origNodes[cIndex].type = 1;
								auto neighboursUpper = adjacent_vertices(top, intersectionG);
								for (auto a : make_iterator_range(neighboursUpper)) {
									int nLvlUpper = overallToLvlNodes[(int)a][0];
									int nIndexUpper = overallToLvlNodes[(int)a][1];
									if (nLvlUpper > cLvl) {
										if ((start.levels[nLvlUpper].origNodes[nIndexUpper].type == CUT || start.levels[nLvlUpper].origNodes[nIndexUpper].type == FILL
											) && !visited[(int)a]) {
											q.push((int)a);
											//cout << "upper pushed onto queue" << endl;
											visited[(int)a] = true;
										}
									}
								}

							}
						}


					}
				}
				for (int i = 1; i < shapes.size() - 1; i++) {
					if (!lvlCovered[i] && (inContraLvl.find(i) != inContraLvl.end())) {
						int cfCt = 0;
						for (int j = 0; j < start.levels[i].origNodes.size(); j++) {
							if (start.levels[i].origNodes[j].type == CUT || start.levels[i].origNodes[j].type == FILL) {
								cfCt++;
							}
						}
						if (cfCt > 0) {
							solveLevel(start, i, levelGraphs, levelEdgeWts, hypernodeSize, productThresh,
								localSteinerTime, bbTime, globalSteinerTime, origLvlNodes[i]
							);
							levelSolve++;
						}
						//cout << "solve level " << i << endl;
					}

				}

			}
			start.contradictions.clear();
			findContradictions(start, overallToLvlNodes, intersectionG, origLvlNodes);

			start.cost = 0;
			int h0Total = 0;
			int h1Total = 0;
			int h2Total = 0;
			double geomCost = 0;
			for (int i = 1; i < shapes.size() - 1; i++) {
				start.cost += start.levels[i].cost;
				h0Total += start.levels[i].h0;
				h1Total += start.levels[i].h1;
				h2Total += start.levels[i].h2;
				geomCost += abs(start.levels[i].geomCost);
			}
			cout << "has " << start.contradictions.size() << " contradictions, cost: " << start.cost <<
				" h0: " << h0Total << " h1: " << h1Total << " h2: " << h2Total << " " << geomCost << " initial cost: " << geomCostInit <<
				endl;
			if (start.cost < bestCost) {
				bestCost = start.cost;
				finalState = start;
				bestStart = l;
			}

		}
		cout << "best start level " << bestStart << endl;

	}
	else {

		while (!stateQ.empty()) {
			State top = stateQ.top();
			vector<int> topNums = getNums(top);
			//vector<int> topEnc = getStateEncoding(top);
			cout.precision(30);
			double geomCost = 0.0;
			for (int i = 0; i < top.levels.size(); i++) {
				geomCost += top.levels[i].geomCost;
			}
			cout << "top state cost " << top.cost << " geometric cost " << geomCost << " contradictions: " << top.contradictions.size() << " queue size " << stateQ.size() << " beam size " << beamSize << endl;
			//for (int i = 0; i < top.contradictions.size(); i++) {
				//cout << "contra lvl " << top.contradictions[i].l1 << " " << top.contradictions[i].n1 << " " << top.contradictions[i].l2 << " " << top.contradictions[i].n2 << " " << endl;
			//}
			stateQ.pop();
			statesPopped++;
			ct++;

			for (int l = 0; l < top.levels.size(); l++) {
				for (int i = 0; i < origLvlNodes[l].size(); i++) {
					if (origLvlNodes[l][i].type == FILL && top.levels[l].origNodes[i].inFg == 1) {
						auto neighbours = adjacent_vertices(i, levelGraphs[l]);
						for (auto u : make_iterator_range(neighbours)) {
							if (origLvlNodes[l][(int)u].type == CUT && top.levels[l].origNodes[(int)u].inFg == 0) {
								cout << "contradiction" << endl;
							}
						}
					}

				}

				//cout << "Level " << i << " topo " << top.levels[i].h0 << " " << top.levels[i].h2 << endl;
			}

			map<vector<int>, bool> inContra;

			int bestContra = 0;
			int bestContraScore = INFINITY;
			for (int i = 0; i < top.contradictions.size(); i++) {
				int diff = -abs((int)top.contradictions[i].lowerIndices.size() - (int)top.contradictions[i].upperIndices.size());
				if (diff < bestContraScore) {
					bestContra = i;
					bestContraScore = diff;
				}
				cout << top.contradictions[i].l1 << " " << top.contradictions[i].n1 << " " <<
					top.contradictions[i].l2 << " " << top.contradictions[i].n2 << endl;
				inContra[{top.contradictions[i].l1, top.contradictions[i].n1}] = true;
				inContra[{top.contradictions[i].l2, top.contradictions[i].n2}] = true;
			}
			map<vector<int>, bool> inBestContra;
			//cout << "in contra " << endl;
			if (top.contradictions.size() > 0) {
				Contradiction cTop = top.contradictions[bestContra];
				map<vector<int>, bool> inBestContra;
				inBestContra[{cTop.l1, cTop.n1}] = true;
				inBestContra[{cTop.l2, cTop.n2}] = true;
				cout << "best " << cTop.l1 << " " << cTop.n1 << endl;
			}
			cout << (inBestContra.find({ 2, 17 }) != inBestContra.end()) << endl;
			//cout << "best contra " << endl;
			for (int l = 0; l < shapes.size(); l++) {
				uint16 spp, bpp, photo;
				int i, j;
				uint16 page;
				string filenameState = outFile + to_string(statesPopped) + "_" + to_string(l) + "State.tif";
				TIFF* out = TIFFOpen(filenameState.c_str(), "w");
				if (!out)
				{
					fprintf(stderr, "Can't open %s for writing\n", argv[1]);
					return 1;
				}
				spp = 1;
				bpp = 32;
				photo = PHOTOMETRIC_MINISBLACK;
				float eps = 0.01;
				for (int s = 0; s < numSlices; s++) {
					TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
					TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
					TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
					TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
					TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
					TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
					TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

					TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

					TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
					TIFFSetField(out, TIFFTAG_PAGENUMBER, s, numSlices);

					int widthP = width;
					int heightP = height;
					int slicesP = numSlices;
					int sp = s;
					float* data = new float[width * height]();
					for (int i = 0; i < width; i++) {
						int ip = i;
						for (int j = 0; j < height; j++) {
							int jp = j;
							data[i + j * width] = 0;
							uint32_t label = Label(levelLabels[l][sp][ip + jp * widthP]);
							if (label != unvisited) {

								if (((int)origLvlNodes[l][label].type) == CORE) {
									//cout << 1 << endl;
									data[i + j * width] = 20;
									//cout << "1e" << endl;
								}
								else {
									if (((int)origLvlNodes[l][label].type) == CUT) {
										if (top.levels[l].origNodes[label].type == CORE) {
											//cout << "2" << endl;
											data[i + j * width] = 19;
											cout << "cut constrained to be 1 " << endl;
										}
										else {
											if (top.levels[l].origNodes[label].type == N) {
												//cout << "3" << endl;
												data[i + j * width] = 18;
												cout << "cut constrained to be 0 " << endl;
											}
											else {
												if (top.levels[l].origNodes[label].inFg == 1) {
													//cout << "4" << endl;
													data[i + j * width] = 17;
													//cout << "4e" << endl;
												}
												else {

													data[i + j * width] = 16;
												}
											}
										}
									}
									if (((int)origLvlNodes[l][label].type) == FILL) {
										//cout << "fill" << endl;
										if (top.levels[l].origNodes[label].type == CORE) {
											//cout << "5" << endl;
											data[i + j * width] = 15;
											cout << "fill constrained to be 1 " << endl;
										}
										else {
											//cout << "6" << endl;
											if (top.levels[l].origNodes[label].type == N) {

												data[i + j * width] = 14;
												cout << "fill constrained to be 0 " << endl;
											}
											else {
												if (top.levels[l].origNodes[label].inFg == 1) {
													data[i + j * width] = 13;
												}
												else {
													data[i + j * width] = 12;
												}
											}
											//cout << "7" << endl;
										}
									}
									//cout << "8" << endl;
									int label1 = (int)label;
									if (inContra.find({ l, label1 }) != inContra.end()) {
										cout << "9" << endl;
										if (top.levels[l].origNodes[label].inFg == 1) {
											if (((int)origLvlNodes[l][label].type) == CUT) {
												data[i + j * width] = 11;
											}
											else { //FILL
												data[i + j * width] = 10;
											}
										}
										else {
											if (((int)origLvlNodes[l][label].type) == FILL) {
												data[i + j * width] = 9;
											}
											else {
												data[i + j * width] = 8;
											}
										}
										//cout << "10" << endl;
									}
									//if (l == 2 && (int)label == 17) {
									//	cout << "in best contra? " << (inBestContra.find({ l, (int)label }) != inBestContra.end()) << endl;
									//}
									if (top.contradictions.size() > 0) {
										Contradiction cTop = top.contradictions[bestContra];

										if (l == 2 && (int)label == 17) {
											cout << "best contra " << cTop.l1 << " " << cTop.n1 << " " << (l == cTop.l1) << " " <<
												((int)label == cTop.n1) << 
												endl;
										}
										if(((l == cTop.l1) && ((int)label == cTop.n1)) || ((l == cTop.l2) && ((int)label == cTop.n2))) {
											//if (inBestContra.find({ l, (int)label }) != inBestContra.end()) {
												cout << "in best contra? " << l << " " << i << " " << j << endl;
												if (top.levels[l].origNodes[label].inFg == 1) {

													if (((int)origLvlNodes[l][label].type) == CUT) {
														data[i + j * width] = 7;
													}
													else { //FILL
														data[i + j * width] = 6;
													}
												}
												else {
													if (((int)origLvlNodes[l][label].type) == FILL) {
														data[i + j * width] = 5;
													}
													else {
														data[i + j * width] = 4;
													}
												}
											//}
											//cout << "12" << endl;
										}
									}
									//cout << "end " << endl;
								}
							}

						}
					}
					//cout << "write " << endl;
					for (j = 0; j < height; j++) {
						TIFFWriteScanline(out, &data[j * width], j, 0);
					}
					TIFFWriteDirectory(out);

				}
				TIFFClose(out);
			}

			if (top.contradictions.size() == 0 || ct == INFINITY) { //INFINITY
				finalState = top;
				break;
			}

			if (!oneConflict) {
				for (int i = 0; i < top.contradictions.size(); i++) {
					int lowerLvl = overallToLvlNodes[top.contradictions[i].lowerIndices[0]][0];
					int lowerNode = overallToLvlNodes[top.contradictions[i].lowerIndices[0]][1];
					if (top.levels[lowerLvl].origNodes[lowerNode].type == FILL || top.levels[lowerLvl].origNodes[lowerNode].type == CUT) {
						State lower = constrainLowerState(top, i, overallToLvlNodes, origLvlNodes,
							intersectionG, levelGraphsInit, levelEdgeWtsIn, hypernodeSize, productThresh, localSteinerTime,
							bbTime, globalSteinerTime, levelNewToOldComps, origLvlNodes, levelSolve);

						vector<int> lowerEnc = getStateEncoding(lower);
						if (stateExists.find(lowerEnc) == stateExists.end()) {
							stateQ.push(lower);
							stateExists[lowerEnc] = true;
						}
					}
					//vector<int> lowerNums = getNums(lower);
					int upperLvl = overallToLvlNodes[top.contradictions[i].upperIndices[0]][0];
					int upperNode = overallToLvlNodes[top.contradictions[i].upperIndices[0]][1];
					if (top.levels[upperLvl].origNodes[upperNode].type == FILL || top.levels[upperLvl].origNodes[upperNode].type == CUT) {
						State upper = constrainUpperState(top, i, overallToLvlNodes, origLvlNodes,
							intersectionG, levelGraphsInit, levelEdgeWtsIn, hypernodeSize, productThresh, localSteinerTime,
							bbTime, globalSteinerTime, levelNewToOldComps, origLvlNodes, levelSolve);
						vector<int> upperEnc = getStateEncoding(upper);
						if (stateExists.find(upperEnc) == stateExists.end()) {
							stateQ.push(upper);
							stateExists[upperEnc] = true;
						}
					}
				}
				if (beamSize < INFINITY) {
					priority_queue<State, vector<State>, CompareState> nextQ;
					while (!stateQ.empty() && nextQ.size() < beamSize) {
						nextQ.push(stateQ.top());
						stateQ.pop();
					}
					stateQ = nextQ;
				}
			}
			else { //one conflict approach
				//cout << "one conflict " << endl;
				
				//cout << "done " << endl;
	
				State lower = constrainLowerState(top, bestContra, overallToLvlNodes, origLvlNodes,
					intersectionG, levelGraphs, levelEdgeWts, hypernodeSize, productThresh, localSteinerTime,
					bbTime, globalSteinerTime, levelNewToOldComps, origLvlNodes, levelSolve);
				vector<int> lowerEnc = getStateEncoding(lower);
				vector<int> lowerNums = getNums(lower);

				//cout << "got lower" << endl;
				if (stateExists.find(lowerEnc) == stateExists.end()) {
					//cout << "lower exists? " << endl;
					stateQ.push(lower);
					stateExists[lowerEnc] = true;
				}
				State upper = constrainUpperState(top, bestContra, overallToLvlNodes, origLvlNodes,
					intersectionG, levelGraphs, levelEdgeWts, hypernodeSize, productThresh, localSteinerTime,
					bbTime, globalSteinerTime, levelNewToOldComps, origLvlNodes, levelSolve);
				vector<int> upperEnc = getStateEncoding(upper);
				if (stateExists.find(upperEnc) == stateExists.end()) {
					//cout << "upper exists? " << endl;
					stateQ.push(upper);
					stateExists[upperEnc] = true;
				}
				if (beamSize < INFINITY) {
					priority_queue<State, vector<State>, CompareState> nextQ;
					int pushCt = 0;
					while (!stateQ.empty() && nextQ.size() < beamSize) {
						nextQ.push(stateQ.top());
						stateQ.pop();
						pushCt++;
					}
					stateQ = nextQ;
				}
			}
		}
	}
	/**	double finalCost = 0;
		float finalGeomCost1 = 0;
		for (int l = 1; l < origLvlNodes.size()-1; l++) {
			float topoCost = finalState.levels[l].h0 + finalState.levels[l].h2 + finalState.levels[l].h1;
			finalCost += (topoCost * finalState.totalWtSum);

			float geomCost = 0;
			for (int i = 0; i < origLvlNodes[l].size(); i++) {
				if (origLvlNodes[l][i].type == CUT) {
					if (finalState.levels[l].origNodes[i].inFg == 0) {
						geomCost += abs(finalState.levels[l].origNodes[i].floatCost);
					}
				}
				if (origLvlNodes[l][i].type == FILL) {
					if (finalState.levels[l].origNodes[i].inFg == 1) {
						geomCost += abs(finalState.levels[l].origNodes[i].floatCost);
					}
				}
			}
			cout << "final cost of level " << l << " is " << finalState.levels[l].h0 << " " << finalState.levels[l].h2 << " " << finalState.levels[l].h1 << " geom " << geomCost << " " << (topoCost * initialState.totalWtSum) << endl;

			finalGeomCost1 += geomCost;
			finalCost += geomCost;
		}
		cout << "finalCost " << finalCost << endl;**/
	auto endSearchTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedEnd = endSearchTime - initStateTime;
	cout << "Initial state had " << initialState.contradictions.size() << " contradictions, states popped: " <<
		statesPopped <<
		endl;
	cout << "Time for generating cuts and fills: " << elapsedCf.count() << " seconds " << endl;
	cout << "Time for generating initial state: " << elapsedInit.count() << " seconds " << endl;
	//cout << "Time for other states: " << elapsedEnd.count() << " seconds " << endl;
	cout << "Time for search: " << elapsedInit.count() + elapsedEnd.count() << " seconds " << endl;
	cout << "Level solves: " << levelSolve << endl;
	cout << "Total time: " << elapsedCf.count() + elapsedInit.count() + elapsedEnd.count() << " seconds " << endl;
	//populate original indices
	vector<vector<node>> initNodes = origLvlNodes;
	if (propagate) {
		for (int l = 0; l < shapes.size(); l++) {
			for (int i = 0; i < finalState.levels[l].origNodes.size(); i++) {

				origLvlNodes[l][i].inFg = finalState.levels[l].origNodes[i].inFg;
			}
		}

	}
	else {
		int numNodes = 0;
		for (int l = 0; l < shapes.size(); l++) {
			for (int i = 0; i < finalState.levels[l].nodes.size(); i++) {
				if (finalState.levels[l].nodes[i].valid) {
					//cout << "in fg? " << (int)finalState.levels[l].nodes[i].type << " " << (int)finalState.levels[l].nodes[i].inFg << endl;
					vector<int> oldIndices = finalState.levelNewToOldComps[l][i];
					for (int j = 0; j < oldIndices.size(); j++) {
						origLvlNodes[l][oldIndices[j]].inFg = finalState.levels[l].nodes[i].inFg;
					}
				}
			}
			numNodes += origLvlNodes[l].size();
		}
		cout << "Number of original nodes: " << numNodes << endl;
	}

	for (int l = 0; l < shapes.size(); l++) {
		for (int i = 0; i < initialState.levels[l].nodes.size(); i++) {
			if (initialState.levels[l].nodes[i].valid) {
				//cout << "in fg? " << (int)finalState.levels[l].nodes[i].type << " " << (int)finalState.levels[l].nodes[i].inFg << endl;
				vector<int> oldIndices = initialState.levelNewToOldComps[l][i];
				for (int j = 0; j < oldIndices.size(); j++) {
					initNodes[l][oldIndices[j]].inFg = initialState.levels[l].nodes[i].inFg;
				}
			}
		}
	}
	numSlices = numSlicesOrig;
	width = origWidth;
	height = origHeight;
	int widthP = width + 2;
	int heightP = height + 2;
	int slicesP = numSlices + 2;
	for (int i = 1; i < shapes.size(); i++) {
		shapes[i] = shapes[i] - 1;
	}

	if (inFileType == 1 || outFileType == 1) {
		if (true) {
			uint16 spp, bpp, photo;
			TIFF* out;
			int i, j;
			uint16 page;
			string outFileTotal = outFile+".tif";
			out = TIFFOpen(outFileTotal.c_str(), "w");
			if (!out)
			{
				fprintf(stderr, "Can't open %s for writing\n", argv[1]);
				return 1;
			}
			spp = 1;
			bpp = 32;
			photo = PHOTOMETRIC_MINISBLACK;
			float eps = 0.01;
			for (int s = 0; s < numSlices; s++) {
				TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
				TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
				TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
				TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

				TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

				TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
				TIFFSetField(out, TIFFTAG_PAGENUMBER, s, numSlices);

				int sp = s + 1;
				float* data = new float[width * height]();
				for (int i = 0; i < width; i++) {
					int ip = i + 1;
					for (int j = 0; j < height; j++) {
						int jp = j + 1;
						for (int l = 1; l < shapes.size() - 1; l++) {
							uint32_t label = Label(levelLabels[l][sp][ip + jp * widthP]);
							if (label != unvisited) {
								if (((int)origLvlNodes[l][label].inFg) == 1) {
									data[i + j * width] = shapes[l] + eps;
								}
							}
						}
					}
				}

				for (j = 0; j < height; j++) {
					TIFFWriteScanline(out, &data[j * width], j, 0);
				}
				TIFFWriteDirectory(out);

			}
			TIFFClose(out);
		}
		for (int l = 0; l < shapes.size(); l++) {
			uint16 spp, bpp, photo;
			TIFF* out;
			int i, j;
			uint16 page;
			string filenameInit = outFile + to_string(l)+"init.tif";
			//string filenameF = outFile + to_string(l) + "F.tif";
			//string filenameC = outFile + to_string(l) + "C.tif";
			string filenameCF = outFile + to_string(l) + "CF.tif";
			out = TIFFOpen(filenameInit.c_str(), "w");
			//TIFF* outF = TIFFOpen(filenameF.c_str(), "w");
			//TIFF* outC = TIFFOpen(filenameC.c_str(), "w");
			TIFF* outCF = TIFFOpen(filenameCF.c_str(), "w");
			if (!out)
			{
				fprintf(stderr, "Can't open %s for writing\n", argv[1]);
				return 1;
			}
			spp = 1;
			bpp = 32;
			photo = PHOTOMETRIC_MINISBLACK;
			float eps = 0.01;
			for (int s = 0; s < numSlices; s++) {
				TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
				TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
				TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
				TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

				TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

				TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
				TIFFSetField(out, TIFFTAG_PAGENUMBER, s, numSlices);

				/**TIFFSetField(outF, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(outF, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(outF, TIFFTAG_BITSPERSAMPLE, bpp);
				TIFFSetField(outF, TIFFTAG_SAMPLESPERPIXEL, spp);
				TIFFSetField(outF, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(outF, TIFFTAG_PHOTOMETRIC, photo);
				TIFFSetField(outF, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

				TIFFSetField(outF, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

				TIFFSetField(outF, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
				TIFFSetField(outF, TIFFTAG_PAGENUMBER, s, numSlices);

				TIFFSetField(outC, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(outC, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(outC, TIFFTAG_BITSPERSAMPLE, bpp);
				TIFFSetField(outC, TIFFTAG_SAMPLESPERPIXEL, spp);
				TIFFSetField(outC, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(outC, TIFFTAG_PHOTOMETRIC, photo);
				TIFFSetField(outC, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

				TIFFSetField(outC, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

				TIFFSetField(outC, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
				TIFFSetField(outC, TIFFTAG_PAGENUMBER, s, numSlices);**/

				TIFFSetField(outCF, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(outCF, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(outCF, TIFFTAG_BITSPERSAMPLE, bpp);
				TIFFSetField(outCF, TIFFTAG_SAMPLESPERPIXEL, spp);
				TIFFSetField(outCF, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(outCF, TIFFTAG_PHOTOMETRIC, photo);
				TIFFSetField(outCF, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

				TIFFSetField(outCF, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

				TIFFSetField(outCF, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
				TIFFSetField(outCF, TIFFTAG_PAGENUMBER, s, numSlices);

				int sp = s + 1;
				float* data = new float[width * height]();
				float* dataF = new float[width * height]();
				float* dataC = new float[width * height]();
				float* dataCF = new float[width * height]();
				for (int i = 0; i < width; i++) {
					int ip = i + 1;
					for (int j = 0; j < height; j++) {
						int jp = j + 1;
						uint32_t label = Label(levelLabels[l][sp][ip + jp * widthP]);
						if (label != unvisited) {

							if (((int)initNodes[l][label].inFg) == 1) {
								data[i + j * width] = 1; //max(shapes[l] + 1, g_Image3DOrig[s][i + j * width]);
							}
							else {
								data[i + j * width] = 0;
							}

							if (((int)origLvlNodes[l][label].type) == CORE) {
								dataCF[i + j * width] = 4;
								
							}

							if(((int)origLvlNodes[l][label].type) == FILL) {
								dataF[i + j * width] = 1; //max(shapes[l] + 1, g_Image3DOrig[s][i + j * width]);
								dataCF[i + j * width] = 2;
							}

							if (((int)origLvlNodes[l][label].type) == CUT) {
								dataCF[i + j * width] = 3;
								dataC[i + j * width] = 1; //max(shapes[l] + 1, g_Image3DOrig[s][i + j * width]);
							}
						}
						
					}
				}

				for (j = 0; j < height; j++) {
					TIFFWriteScanline(out, &data[j * width], j, 0);
				}
				TIFFWriteDirectory(out);

				/**for (j = 0; j < height; j++) {
					TIFFWriteScanline(outF, &dataF[j * width], j, 0);
				}
				TIFFWriteDirectory(outF);

				for (j = 0; j < height; j++) {
					TIFFWriteScanline(outC, &dataC[j * width], j, 0);
				}
				TIFFWriteDirectory(outC);**/

				for (j = 0; j < height; j++) {
					TIFFWriteScanline(outCF, &dataCF[j * width], j, 0);
				}
				TIFFWriteDirectory(outCF);

			}
			TIFFClose(outCF);
			//TIFFClose(outF);
			//TIFFClose(outC);
		}
	}
	else {
		for (int l = 0; l < shapes.size(); l++) {

			for (int s = 0; s < numSlices; s++) {
				int digits = numDigits(s);
				string numStr = "";
				for (int n = 0; n < 4 - digits; n++) {
					numStr += "0";

				}
				numStr += std::to_string(s);
				string filenameInit = outFile + "init/" + to_string(l) + "/" + numStr + ".png";
				string filenameU = outFile + "upper_conflicts/" + to_string(l) + "/" + numStr + ".png";
				string filenameL = outFile + "lower_conflicts/" + to_string(l) + "/" + numStr + ".png";
				//cout << filenameL << endl;
				string filenameC = outFile + "cuts/" + to_string(l) + "/" + numStr + ".png";
				string filenameF = outFile + "fills/" + to_string(l) + "/" + numStr + ".png";
				string filenameK = outFile + "kernels/" + to_string(l) + "/" + numStr + ".png";
				uint8_t* label8Init = new uint8_t[width * height]();
				uint8_t* label8U = new uint8_t[width * height]();
				uint8_t* label8L = new uint8_t[width * height]();
				uint8_t* label8C = new uint8_t[width * height]();
				uint8_t* label8F = new uint8_t[width * height]();
				uint8_t* label8K = new uint8_t[width * height]();
				int sp = s + 1;
				for (int i = 0; i < width; i++) {
					int ip = i + 1;
					for (int j = 0; j < height; j++) {
						int jp = j + 1;
						uint32_t label = Label(levelLabels[l][sp][ip + jp * widthP]);
						if (label != unvisited) {
							if (((int)initNodes[l][label].inFg) == 1) {
								label8Init[i + j * width] = max(shapes[l] + 1, g_Image3DOrig[s][i + j * width]);
							}
							else {
								label8Init[i + j * width] = min(g_Image3DOrig[s][i + j * width], shapes[l]);
							}
							if (upperContradictions.find({ l, (int)label }) != upperContradictions.end()) {
								label8U[i + j * width] = shapes[l] + 1; // shapes[l] + 1;
							}
							if (lowerContradictions.find({ l, (int)label }) != lowerContradictions.end()) {
								label8L[i + j * width] = shapes[l] + 1; // shapes[l] + 1;
							}

							if (((int)origLvlNodes[l][label].type) == CUT) {
								if (l == 0) {
									cout << "cut 0" << endl;
								}
								label8C[i + j * width] = 1;
							}
							if (((int)origLvlNodes[l][label].type) == FILL) {
								label8F[i + j * width] = 1;
							}
							if (((int)origLvlNodes[l][label].type) == CORE) {
								label8K[i + j * width] = 1;
							}
						}
					}
				}
				int wrote = stbi_write_png(filenameInit.c_str(), width, height, 1, label8Init, width);
				wrote = stbi_write_png(filenameU.c_str(), width, height, 1, label8U, width);
				wrote = stbi_write_png(filenameL.c_str(), width, height, 1, label8L, width);
				wrote = stbi_write_png(filenameC.c_str(), width, height, 1, label8C, width);
				wrote = stbi_write_png(filenameF.c_str(), width, height, 1, label8F, width);
				wrote = stbi_write_png(filenameK.c_str(), width, height, 1, label8K, width);
			}
		}
		for (int s = 0; s < shapes.size(); s++) {
			cout << "shape " << s << " is " << shapes[s] << endl;
		}
		for (int s = 0; s < numSlices; s++) {

			int digits = numDigits(s);
			string numStr = "";
			for (int n = 0; n < 4 - digits; n++) {
				numStr += "0";

			}
			numStr += std::to_string(s);//+ to_string(l-1)
			//string filename = outFile + to_string(l - 1) + "/" + numStr + ".png";
			//string filename = outFile + "upper_conflicts/" + to_string(l - 1) + "/" + numStr + ".png";
			//string filename1 = outFile + "lower_conflicts/" + to_string(l - 1) + "/" + numStr + ".png";
			//string filenameCuts = outFile + "cuts/" + to_string(l - 1) + "/" + numStr + ".png";
			//string filenameFills = outFile + "fills/" + to_string(l - 1) + "/" + numStr + ".png";
			//string filenameInit = outFile + "init/" + to_string(l - 1) + "/" + numStr + ".png";

			string filename = outFile + "ours/" + numStr + ".png";  //outFile + numStr + ".png";

			uint8_t* label8 = new uint8_t[width * height]();
			uint8_t* label8T = new uint8_t[width * height]();
			uint8_t* label8Cuts = new uint8_t[width * height]();
			uint8_t* label8Fills = new uint8_t[width * height]();
			uint8_t* label8Init = new uint8_t[width * height]();
			int sp = s + 1;

			for (int i = 0; i < width; i++) {
				int ip = i + 1;
				for (int j = 0; j < height; j++) {
					//label8[i + j * width] = min(g_Image3DOrig[s][i + j * width], shapes[1]);
					label8Cuts[i + j * width] = 0;
					label8Fills[i + j * width] = 0;
					label8Init[i + j * width] = 0;
					int jp = j + 1;
					//int l = shapes.size() - 2;
					//int l = 9;

					for (int l = 1; l < shapes.size() - 1; l++) {
						uint32_t label = Label(levelLabels[l][sp][ip + jp * widthP]);
						if (label != unvisited) {
							//if (upperContradictions.find({l, (int)label}) != upperContradictions.end()) {
							//	label8[i + j * width] = shapes[l] + 1; // shapes[l] + 1;
							//}
							//if (lowerContradictions.find({ l, (int)label }) != lowerContradictions.end()) {
							//	label8T[i + j * width] = shapes[l] + 1; // shapes[l] + 1;
							//}
							//if ((int)origLvlNodes[l][label].type == 0) { // || ((int)origLvlNodes[l][label].type == FILL
								//label8[i + j * width] = 3;
							//}
							//else {
							///	if ((int)origLvlNodes[l][label].type == CUT) { // || ((int)origLvlNodes[l][label].type == FILL
								//	label8Cuts[i + j * width] = 1;
								//}
								//else {
								//	if ((int)origLvlNodes[l][label].type == FILL) { // || ((int)origLvlNodes[l][label].type == FILL
								//		label8Fills[i + j * width] = 1;
								//	}
									//else {
									//	label8[i + j * width] = 0;
									//}
								//}
							//}
							if (((int)origLvlNodes[l][label].inFg) == 1) {
								//if (l == shapes.size() - 2) {
									//label8[i + j * width] = max(g_Image3DOrig[s][i + j * width], shapes[l] + 1);
								//}
								//else {

								label8[i + j * width] = shapes[l] + 1;
								//}


								//label8Init[i + j * width] = shapes[l] + 1;
								/**if (l - 1 >= 0) {
									int label1 = Label(levelLabels[l - 1][sp][ip + jp * widthP]);
									if (((int)origLvlNodes[l - 1][label1].inFg) == 0) {
										//cout << "contra " << label << " " << label1 << " " << l << endl;
										//draw contradiction
										label8[i + j * width] = shapes[l] + 1;
									}
								}**/
							}
							//else {
								//if (l == 1) {
								//label8[i + j * width] = min(g_Image3DOrig[s][i + j * width], shapes[1]);
								//}
							//}

						//}
						}
						//if ((int)g_Image3DOrig[s][i + j * width] < shapes[1]) {
							//cout << (int)label8[i + j * width] << " " << g_Image3DOrig[s][i + j * width] << " " << endl;
						//}
					}
				}
			}
			int wrote = stbi_write_png(filename.c_str(), width, height, 1, label8, width);
			//int wroteL = stbi_write_png(filename1.c_str(), width, height, 1, label8T, width);
			//wrote = stbi_write_png(filenameCuts.c_str(), width, height, 1, label8Cuts, width);
			//wrote = stbi_write_png(filenameFills.c_str(), width, height, 1, label8Fills, width);
			//wrote = stbi_write_png(filenameInit.c_str(), width, height, 1, label8Init, width);
			//}
		}
		//}
		/**double finalGeomCost = 0;
		double initGeomCost = 0;
		for (int l = 0; l < origLvlNodes.size(); l++) {
			float levelICost = 0;
			float levelFCost = 0;
			for (int i = 0; i < origLvlNodes[l].size(); i++) {
				if ((int)origLvlNodes[l][i].type == FILL) {
					if (finalState.levels[l].origNodes[i].inFg == 1) {
						//finalCost += abs(origLvlNodes[l][i].intensity);
						finalGeomCost += abs(origLvlNodes[l][i].intensity);
						levelFCost += abs(origLvlNodes[l][i].intensity);
						//cout << "final fill with g " << origLvlNodes[l][i].intensity << endl;
					}
					if (initialState.levels[l].origNodes[i].inFg == 1) {
						initGeomCost += abs(origLvlNodes[l][i].intensity);
						//cout << "fill with g " << origLvlNodes[l][i].intensity << endl;
						levelICost += abs(origLvlNodes[l][i].intensity);
					}
				}
				if ((int)origLvlNodes[l][i].type == CUT) {
					if (finalState.levels[l].origNodes[i].inFg == 0) {
						//finalCost += abs(origLvlNodes[l][i].intensity);
						finalGeomCost += abs(origLvlNodes[l][i].intensity);
						levelFCost += abs(origLvlNodes[l][i].intensity);
						//cout << "final cut with g " << origLvlNodes[l][i].intensity << endl;
					}
					if (initialState.levels[l].origNodes[i].inFg == 0) {
						initGeomCost += abs(origLvlNodes[l][i].intensity);
						levelICost += abs(origLvlNodes[l][i].intensity);
					}
				}
			}
			//cout << "level " << l << " initial cost : " << levelICost << " final cost: " << levelFCost << endl;
		}**/
	}
	cout.precision(30);

	double igCost2 = 0;
	int initH0 = 0; int initH1 = 0; int initH2 = 0; double initGeomCost = 0;
	for (int l = 1; l < initialState.levels.size() - 1; l++) {
		initH0 += initialState.levels[l].h0;
		initH1 += initialState.levels[l].h1;
		initH2 += initialState.levels[l].h2;
		initGeomCost += initialState.levels[l].geomCost;
		for (int i = 0; i < initialState.levels[l].origNodes.size(); i++) {
			if ((origLvlNodes[l][i].type == CUT && initialState.levels[l].origNodes[i].inFg == 0) ||
				(origLvlNodes[l][i].type == FILL && initialState.levels[l].origNodes[i].inFg == 1)
				) {
				igCost2 += abs(origLvlNodes[l][i].floatCost);
			}
		}
		cout << "Initial topology of level " << l << ": Components: " << initialState.levels[l].h0 << " Cavities: " << initialState.levels[l].h2 << " Cycles: " << initialState.levels[l].h1 << endl; 
	}
	cout << "Initial state: Components: " << initH0 << " Cavities: " << initH2 << " Cycles: " << initH1 << " geom cost " << initGeomCost << " " << igCost2 << endl;
	int finalH0 = 0; int finalH1 = 0; int finalH2 = 0; double finalGeomCost = 0;
	double fgCost2 = 0;
	for (int l = 1; l < finalState.levels.size() - 1; l++) {
		finalH0 += finalState.levels[l].h0;
		finalH1 += finalState.levels[l].h1;
		finalH2 += finalState.levels[l].h2;
		finalGeomCost += finalState.levels[l].geomCost;
		for (int i = 0; i < finalState.levels[l].origNodes.size(); i++) {
			if ((origLvlNodes[l][i].type == CUT && finalState.levels[l].origNodes[i].inFg == 0) ||
				(origLvlNodes[l][i].type == FILL && finalState.levels[l].origNodes[i].inFg == 1)
				) {
				fgCost2 += abs(origLvlNodes[l][i].floatCost);
			}
		}
	}
	cout << "Final state: Components: " << finalH0 << " Cavities: " << finalH2 << " Cycles: " << finalH1 << " geom cost " << finalGeomCost << " " << fgCost2 << endl;
	int h0F = 0; int h1F = 0; int h2F = 0;
	//int l = 17;
	//finalCost = 0;
	for (int l = 1; l < shapes.size() - 1; l++) {
		//for (int l = 18; l < 19; l++) { //shapes.size();
			/**for (int j = 0; j < levelNodes[l].size(); j++) {
				if (levelNodes[l][j].valid) {
					cout << "Label " << l << " " << (int)levelNodes[l][j].type << " " << (int)levelNodes[l][j].inFg << endl;
				}
			}**/

		vector< uint8_t*> topoImg;
		for (int s = 0; s < numSlices; s++) {
			int sp = s + 1;
			uint8_t* label8T = new uint8_t[width * height]();
			for (int i = 0; i < width; i++) {
				int ip = i + 1;
				for (int j = 0; j < height; j++) {
					//if (((int)finalState.levels[l].nodes[Label(levelLabels[l][s][i + j * width])].type) == CUT) {
					//	cout << "cut in fg " << ((int)finalState.levels[l].nodes[Label(levelLabels[l][s][i + j * width])].inFg) << endl;
					//}
					int jp = j + 1;
					uint32_t label = Label(levelLabels[l][sp][ip + jp * widthP]);//levelNewToOldIndices[l][
					//if (((int)origLvlNodes[l][label].type) == 0 || ((int)origLvlNodes[l][label].type) == 2 || ((int)origLvlNodes[l][label].type) == 3) {
					if (((int)origLvlNodes[l][label].inFg) == 1) {
						label8T[i + j * width] = 1;
						if (l - 1 >= 0) {
							int label1 = Label(levelLabels[l - 1][sp][ip + jp * widthP]);
							if (((int)origLvlNodes[l - 1][label1].inFg) == 0) {
								//cout << "contra " << label << " " << label1 << " " << l << endl;
							}
						}
					}
					else {
						label8T[i + j * width] = 0;
					}
				}
			}
			topoImg.push_back(label8T);
		}
		std::vector<int> topoNums = getTopoFromBImg(topoImg, 1, width, height, numSlices);
		std::cout << "Final topology: Components: " << topoNums[0] << " Cavities: " << topoNums[1] << " Cycles: " << topoNums[2] << std::endl;
		//finalCost += ((topoNums[0] + topoNums[1] + topoNums[2]) * finalState.levels[l].wtSum);
		//finalCost += finalState.levels[l].geomCost;
		h0F += topoNums[0];
		h1F += topoNums[2];
		h2F += topoNums[1];
	}
	/**
		int ct1 = 0;
		std::vector< std::vector<int > > mask  = structCross3D;
		vector<vector<int> > comps;
		std::vector< std::vector< std::vector<bool>>> visited(width, std::vector<std::vector<bool>>(height, std::vector<bool>(numSlices, false)));
		unordered_map<int, int> compMap;
		for (int i = 0; i < width; i++) {
			int ip = i + 1;
			for (int j = 0; j < height; j++) {
				int jp = j + 1;
				for (int k = 0; k < numSlices; k++) {
					int sp = k + 1;
					//Unvisited foreground voxels
					uint32_t label = Label(levelLabels[18][sp][ip + jp * widthP]);///levelNewToOldIndices[l][
					if ((int)label == 1114) {
						for (int s = 0; s < structCube.size(); s++) {
							Coordinate np(ip + structCube[s][0], jp + structCube[s][1], sp + structCube[s][2]);
							uint32_t nlabel = Label(levelLabels[18][np.z][(np.x) + (np.y) * widthP]);
							cout << ip << " " << jp << " " << sp << " " << label << " neighbor " << np.x << " " << np.y << " " << np.z << " " <<
								nlabel << " in fg " << origLvlNodes[18][nlabel].inFg << " " << levelEdgeWts[18][{(int)label, (int)nlabel}] << " " << width << " " << height << " " << numSlices << endl;
						}
					}

					if (((int)origLvlNodes[18][label].inFg) == 1 && visited[i][j][k] == false) {

						ct1 += 1;
						cout << "component " << endl;
						std::queue<Coordinate> q;
						q.push(Coordinate(i, j, k));
						visited[i][j][k] = true;
						compMap[label] = ct1;
						vector<int> comp = { (int)label };
						while (!q.empty()) {
							Coordinate qp = q.front();
							q.pop();
							for (int s = 0; s < mask.size(); s++) {
								Coordinate np(qp.x + mask[s][0], qp.y + mask[s][1], qp.z + mask[s][2]);
								if (np.x >= 0 && np.x < width && np.y >= 0 && np.y < height && np.z >= 0 && np.z < numSlices) {
									uint32_t nlabel = Label(levelLabels[18][np.z+1][(np.x+1) + (np.y+1) * widthP]);
									if (((int)origLvlNodes[18][nlabel].inFg) == 1 && visited[np.x][np.y][np.z] == false) {
										visited[np.x][np.y][np.z] = true;
										q.push(np);
										//cout << "pushed into q " << endl;
										if (compMap.find(nlabel) != compMap.end()) {
											if (compMap[nlabel] != ct1) {
												cout << nlabel << " duplicated " << endl;
											}
										}
										else {
											comp.push_back(nlabel);
										}
										compMap[nlabel] = ct1;
									}
								}
							}
						}
						comps.push_back(comp);
					}
				}
			}
		}
		cout << "comps " << ct1 << endl;
		for (int i = 0; i < comps.size(); i++) {
			cout << "comp " << i << " has size " << comps[i].size() << endl;
		}
		auto neighbours = adjacent_vertices(comps[1][0], levelGraphs[18]);
		for (auto u : make_iterator_range(neighbours)) {
			cout << "between " << comps[1][0] << " of type " << (int)origLvlNodes[18][comps[1][0]].type << " in fg " << (int)origLvlNodes[18][comps[1][0]].inFg << " and " <<
				(int)u << " " << (int)origLvlNodes[18][(int)u].type << " in fg " << (int)origLvlNodes[18][u].inFg << " with edge wt " << levelEdgeWts[18][{comps[1][0], (int)u}] << endl;
		}**/

		//	std::cout << "Total final topology: Components: " << h0F << " Cavities: " << h2F << " Cycles: " << h1F << " geometry cost: " <<
			//	finalGeomCost << " initial geometry cost: " << initGeomCost << 
				//endl;

			/**
			cout << "begin exporting png " << endl;
			for (int l = 0; l < shapes.size(); l++) {
				for (int s = 0; s < numSlices; s++) {

					int digits = numDigits(s);
					string numStr = "";
					for (int n = 0; n < 4 - digits; n++) {
						numStr += "0";

					}
					numStr += std::to_string(s);
					string filename = outFile + numStr + "_"+to_string(l)+".png";

					uint8_t* label8 = new uint8_t[width * height]();
					uint8_t* label8T = new uint8_t[width * height]();

					for (int i = 0; i < width; i++) {
						for (int j = 0; j < height; j++) {
							label8[i + j * width] = 0;
							if (Label(levelLabels[l][s][i + j * width]) != unvisited) {
								if (((int)levelNodes[l][Label(levelLabels[l][s][i + j * width])].type) == 0) {
									label8[i + j * width] = 255;
								}
								else {
									if (((int)levelNodes[l][Label(levelLabels[l][s][i + j * width])].type) == 2) {
										label8[i + j * width] = 180;
									}
									else {
										if (((int)levelNodes[l][Label(levelLabels[l][s][i + j * width])].type) == 3) {
											label8[i + j * width] = 100;
										}
										else {
											if (((int)levelNodes[l][Label(levelLabels[l][s][i + j * width])].type) == 1) {
												label8[i + j * width] = 10;
											}
											else {
												//cout << levelNodes[l].size() << " " << Label(levelLabels[l][s][i + j * width]) << endl;

											}
										}

									}
								}
							}
						}
					}
					int wrote = stbi_write_png(filename.c_str(), width, height, 1, label8, width);
					}
					cout << "done " << endl;
				}**/

	return 0;
}