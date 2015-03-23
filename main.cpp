#include <iostream>
#include <fstream>	//files
#include <vector>
#include <algorithm>	//sort
//#include <stdlib.h>	//rand
#include <math.h>       // atan
#include <sstream>	//ostringstream
#include <iomanip>	//setfill, setw
#include <random>	//Random numbers with different distributions
#include <boost/concept_check.hpp>

#define PI 3.14159265

using namespace std;

typedef vector<float> floatVector;
typedef vector<floatVector> MatVector;

struct locationNode {
  int ID;
  float objectProbability;
  float priority;
  float area;
};


float DISTANCES_MEAN;
float DISTANCES_SIGMA;
int SELECTED_MAP;
string mapDataBase;


float optimalExpectedDistance;
vector<locationNode> optimalExpectedRoute;
float worstExpectedDistance;
vector<locationNode> worstExpectedRoute;
float optimalDistance;
vector<locationNode> optimalRoute;
float worstDistance;
vector<locationNode> worstRoute;
float optimalHeuristicDistance;
vector<locationNode> optimalHeuristicRoute;
float partialOptimalHeuristicDistance;
vector<locationNode> partialOptimalHeuristicRoute;


void generateRoomProbabilities(string distributionType, int nintervals, vector<double> &p){
  
  const bool showGraph=false;
  const int experiments=10000;             // number of experiments
  const int starsToShow=10*nintervals;     // maximum number of stars to distribute in the graph
  
  default_random_engine generator;
  
  uniform_real_distribution<double> distributionU(0.0,1.0);
  normal_distribution<double> distributionN(0.5,0.2);
  gamma_distribution<double> distributionG(2.0,0.2);
  exponential_distribution<double> distributionE(6.5);
  
  int notused=0;
  p.assign(nintervals,0);
  
  for (int i=0; i<experiments; ++i) {
    double number;
    if (distributionType=="Uniform") number = distributionU(generator);
    if (distributionType=="Normal") number = distributionN(generator);
    if (distributionType=="Gamma") number = distributionG(generator);
    if (distributionType=="Exponential") number = distributionE(generator);
    //cout << number << endl;
    if (number>=0 && number<1) {
      ++p[int(nintervals*number)];
    }
    else
      notused++;
  }
  
  cout << "Using distribution " << distributionType << endl;
  cout << "Percentage unused random numbers (out of range):" << float(100*notused)/experiments << endl << endl;
  
  cout << fixed; cout.precision(5);
  double totalP=0;
  for (int i=0; i<nintervals; ++i) {
    if (showGraph) {
      cout << float(i)/nintervals << "-" << float(i+1)/nintervals << ": ";
      cout << string(int(p[i])*starsToShow/experiments,'*');
      cout << "   " << p[i]/experiments << endl;
    }
    p[i]/=experiments;
    totalP+=p[i];
  }
  //Normalize to 1
  for (int i=0; i<nintervals; ++i) p[i]/totalP;
}


void loadMapData(int numNodes, int numMap, string mapType, MatVector &distancesTable) {
  string dataFileName;
  
  ostringstream strNodes;
  strNodes << setfill('0') << setw(2) << numNodes;
  ostringstream strNumMap;
  strNumMap << setfill('0') << setw(3) << numMap;
  dataFileName = mapDataBase + "map" + mapType + strNodes.str() + strNumMap.str() + ".txt";
  cout << "Loading data from " << dataFileName << endl;
  ifstream fileStream(dataFileName.c_str());
  
  string mapFileName;
  fileStream >> mapFileName;
  //cout << "Map image file: " << mapFileName << endl;
  
  int roomListSize;
  fileStream >> roomListSize;
  //cout << "Rooms: " << roomListSize << endl;
  
  for (int i=0; i<roomListSize; i++) {
    int x,y;
    fileStream >> x >> y;
  }
  for (int i=0; i<roomListSize; i++) {
    for (int j=0; j<roomListSize; j++) {
      fileStream >> distancesTable[i][j];
      //cout << distancesTable[i][j] << " ";
    }
    //cout << endl;
  }
  fileStream.close();
}

float maxDistance(MatVector distancesTable, vector<locationNode> locationList) {
  float maxDist=0;
  int numNodes=distancesTable.size();
  for (int i=0; i<numNodes; i++) {
    for (int j=0; j<numNodes; j++) {
      if (distancesTable[i][j]>maxDist) maxDist=distancesTable[i][j];
    }
  }
  float maxArea=0;
  numNodes=locationList.size();
  for (int i=0; i<numNodes; i++) {
    if (locationList[i].area>maxArea) maxArea=locationList[i].area;
  }
  return maxDist+maxArea;
}


string buildEnvironment(vector<locationNode> &locationList, MatVector &distancesTable, int numNodes, string mapType, string probabilityDistribution) {
  ostringstream stream2String;
  string msg;
  
  //Initialize random generators
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.0,1.0); 	//numbers in [0,1]
  
  //Pass the random numbers according the map number 
  //so, choosing different map generates different random numbers
  for (int i=0; i<SELECTED_MAP*numNodes*numNodes*numNodes; i++) { 
    float dump=distribution(generator);
  }
  
  vector<double> roomProbabilities;
  generateRoomProbabilities(probabilityDistribution, numNodes, roomProbabilities);
  for (int i=0; i<numNodes; i++) {
    locationNode location;
    location.ID=i;
    location.objectProbability=roomProbabilities[i];
    location.priority=0;
    location.area = ( DISTANCES_MEAN + DISTANCES_SIGMA * (distribution(generator)-0.5) )/2;  //Assume areas represent in general 1/2 of distances between rooms (Here area defined as robot movement inside rooms
    //Fix this for real maps where distances are not generated with random functions
    locationList.push_back(location);
  }
  
  floatVector rowMat;
  rowMat.assign(numNodes,0);
  for (int i=0; i<numNodes; i++) {
    distancesTable.push_back(rowMat);
  }
  
  if (mapType=="R") { 	//Type: House(H), Office (O), Random (R)
    //Generate distances randomly
    for (int i=1; i<numNodes; i++) {
      for (int j=0; j<i; j++) {
	distancesTable[i][j]= DISTANCES_MEAN + DISTANCES_SIGMA * (distribution(generator)-0.5);
	distancesTable[j][i]=distancesTable[i][j];
      }
    }
  }
  if (mapType=="H") { 	//Type: House(H), Office (O), Random (R) 
    //Load distances from files
    //int numMap = MAX_MAPS*(float(rand()) / RAND_MAX)+1; //Random house map
    int numMap = SELECTED_MAP;
    loadMapData(numNodes, numMap, mapType, distancesTable);
  }
  if (mapType=="O") { 	//Type: House(H), Office (O), Random (R) 
    //Not yet office maps
  }
  
  //Print node data
  cout << "Locations (ID,Priority,Probability,Area)" << endl;
  stream2String << "Locations (ID,Priority,Probability,Area)" << endl;
  for (int i=0; i<locationList.size(); i++) {
    cout << locationList[i].ID << ", ";
    cout << locationList[i].priority  << ", ";
    cout << locationList[i].objectProbability << ", ";
    cout << locationList[i].area << endl;
    
    stream2String << locationList[i].ID << ",";
    stream2String << locationList[i].priority  << ",";
    stream2String << locationList[i].objectProbability << ",";
    stream2String << locationList[i].area << endl;
  }
  
  //Print distances data
  cout << "Distances table" << endl;
  stream2String << "Distances table" << endl;
  for (int i=0; i<numNodes; i++) {
    for (int j=0; j<numNodes; j++) {
      cout << distancesTable[i][j] << ",";
      stream2String << distancesTable[i][j] << ",";
    }
    cout << endl;
    stream2String << endl;
  }
  msg = stream2String.str();
  return msg;
}

string msgLocationOrder(vector<locationNode> locationList, string sortingMethod, time_t initTime) {
  ostringstream stream2String;
  string msg;
  
  cout << "Locations ordered by " << sortingMethod << " in " << (time(NULL)-initTime) << " seconds" << endl;
  stream2String << "Locations ordered by " << sortingMethod << " in " << (time(NULL)-initTime) << " seconds" << endl;
  for (int i=0; i<locationList.size(); i++) {
    cout << locationList[i].ID << ", ";
    stream2String << locationList[i].ID << ",";
  }
  cout << endl;
  stream2String << endl;
  
  msg = stream2String.str();
  return msg;
}

bool wayToSortLocations(locationNode i, locationNode j) {
  return i.priority > j.priority;
}

void generateRouteMostProbable(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = locationList[i].objectProbability;
    locationList[startingNode].priority += locationList[i].priority;	//The highest priority
  }
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}

void generateRouteClosestFixed(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = 1/distancesTable[startingNode][i];
    locationList[startingNode].priority += locationList[i].priority;	//The highest priority
  }
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}

void generateRouteClosest(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  //The closest node from the previous. Not from the origin
  vector<locationNode> locationListSorted;
  locationListSorted.push_back(locationList[startingNode]);
  locationList.erase(locationList.begin()+startingNode);
  int nextLocation=startingNode;
  int totalLocations=locationList.size();
  for (int i=0; i<totalLocations; i++) {
    float nextLocationValue=0;
    for (int j=0; j<locationList.size(); j++) {
      float locationValue = 1/distancesTable[startingNode][locationList[j].ID];
      if (locationValue>nextLocationValue) {
	nextLocationValue=locationValue;
	nextLocation=j;
      }
    }
    startingNode=locationList[nextLocation].ID;
    locationListSorted.push_back(locationList[nextLocation]);
    locationList.erase(locationList.begin()+nextLocation);
  }
  locationList=locationListSorted;
}

void generateRouteRandom(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = float(rand()) / RAND_MAX;
    locationList[startingNode].priority += locationList[i].priority;	//The highest priority
  }
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}

void generateRouteExpectedDistanceFixed(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = locationList[i].objectProbability/distancesTable[startingNode][i];
    locationList[startingNode].priority += locationList[i].priority;	//The highest priority
  }
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}



void nonDominatedHull(int referenceNode, vector<locationNode> locationList, MatVector distancesTable, vector<locationNode> &locationListHull) {
  
  vector<locationNode> locationListHullTemp;
  float minDistance=distancesTable[referenceNode][locationList[0].ID];
  float maxProbabilityAtMinDistance=locationList[0].objectProbability;
  int minNode=0;
  
  //Generate first element of the partial convex hull
  for (int i=0; i<locationList.size(); i++) {
    if (distancesTable[referenceNode][locationList[i].ID]<minDistance) {
      minDistance=distancesTable[referenceNode][locationList[i].ID];
      minNode=i;
      maxProbabilityAtMinDistance=locationList[i].objectProbability;
    }
    else if (distancesTable[referenceNode][locationList[i].ID]==minDistance) {
      if (locationList[i].objectProbability > maxProbabilityAtMinDistance) {
	minNode=i;
	maxProbabilityAtMinDistance=locationList[i].objectProbability;
      }
    }
  }
  //Generate the rest of the elements of the partial convex hull
  
  float maxAngle=0.0; //any value for starting cicle of searching convex hull
  int maxAngleNode=0;
  while (maxAngle>=0) {
    locationListHullTemp.push_back(locationList[minNode]);
    locationList.erase(locationList.begin()+minNode);
    maxAngle=(-3)*PI; //Any lowest value for searching the maximun angle
    for (int i=0; i<locationList.size(); i++) {
      float angle=atan2( locationList[i].objectProbability-maxProbabilityAtMinDistance, distancesTable[referenceNode][locationList[i].ID]-minDistance);  //atan2(y,x)
      if (angle>maxAngle) {
	maxAngleNode=i;
	maxAngle=angle;
      }
      else if (angle==maxAngle) {
	if (distancesTable[referenceNode][locationList[i].ID] < distancesTable[referenceNode][locationList[maxAngleNode].ID]) {
	  maxAngleNode=i;
	  maxAngle=angle;
	}
      }
    }
    minDistance=distancesTable[referenceNode][locationList[maxAngleNode].ID];
    maxProbabilityAtMinDistance=locationList[maxAngleNode].objectProbability;
    minNode=maxAngleNode;
  }
  
  locationListHull=locationListHullTemp;
}


void notInList(vector<int> nodesInlist, vector<locationNode> locationList, vector<locationNode> &nodesNotInList) {
  nodesNotInList=locationList;
  for (int i=0; i<nodesInlist.size(); i++) {
    int j;
    for (j=0; j<nodesNotInList.size() && nodesNotInList[j].ID != nodesInlist[i]; j++) {
    };
    if (j<nodesNotInList.size()) {
      nodesNotInList.erase(nodesNotInList.begin()+j);
    }
  }
}


float optimalExpectedDistanceMurrieta;

void makeRouteMurrieta(MatVector distancesTable, vector< vector<int> > locationListPrevious, vector<locationNode> locationList,vector<int> &listSorted, int maxNumNodes, int numNodes) {
  vector< vector<int> > newLocationListPrevious;
  vector<int> optimalExpectedRouteMurrieta;
  //ver si ya la ruta es de longitud size of locationList
  if ( numNodes<maxNumNodes && locationListPrevious[0].size()<locationList.size() ) {
    
    for (int i=0; i<locationListPrevious.size(); i++) {
      //Make list with nodes not visited yet
      vector<locationNode> locationListTemp;
      notInList(locationListPrevious[i], locationList, locationListTemp);
      
      //Expandir el convex hull
      vector<locationNode> locationListHull;
      nonDominatedHull(locationListPrevious[i][locationListPrevious[i].size()-1], locationListTemp, distancesTable, locationListHull);
      for (int j=0; j<locationListHull.size(); j++) {
	vector<int> rowLocationListPrevious = locationListPrevious[i];
	rowLocationListPrevious.push_back(locationListHull[j].ID);
	newLocationListPrevious.push_back(rowLocationListPrevious);
	numNodes++; //Verificar aqui porque la condicion es de n nodos totales no n nodos en nivel final
	if (numNodes>maxNumNodes) { //Esto no parece logico que se corte aqui porque en este nivel no se han analizado todas las posibilidades
	  //evaluar todas las rutas en newLocationListPrevious
	}
      }
    }
  }
  else {
    newLocationListPrevious=locationListPrevious;
    float optimalExpectedDistanceLocal=optimalExpectedDistanceMurrieta;
    for (int i=0; i<newLocationListPrevious.size(); i++) {
      
      float accumulatedDistance=locationList[ newLocationListPrevious[i][0] ].area;
      float expectedDistance=accumulatedDistance * locationList[ newLocationListPrevious[i][0] ].objectProbability;
      for (int j=1; j<newLocationListPrevious[i].size(); j++) {
	accumulatedDistance+=distancesTable[ newLocationListPrevious[i][j-1] ][ newLocationListPrevious[i][j] ] + locationList[ newLocationListPrevious[i][j] ].area;
	expectedDistance+=accumulatedDistance * locationList[ newLocationListPrevious[i][j] ].objectProbability;
      }
      if (expectedDistance<optimalExpectedDistanceLocal) {
	optimalExpectedDistanceLocal=expectedDistance;
	optimalExpectedRouteMurrieta=newLocationListPrevious[i];
      }
    }
    listSorted=optimalExpectedRouteMurrieta;
    return;
  }
  makeRouteMurrieta(distancesTable, newLocationListPrevious, locationList, listSorted, maxNumNodes, numNodes);
}

void generateRouteExpectedDistanceMurrieta(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  vector<int> locationListSorted;
  
  locationListSorted.push_back(startingNode);
  int maxNodeSize=locationList.size();
  while ( locationListSorted.size()<locationList.size() ) {
    vector<locationNode> locationListTemp;
    notInList(locationListSorted, locationList, locationListTemp);
    int nextLocation=0;
    float nextLocationValue=0;
    for (int i=0; i<locationListTemp.size(); i++) {
      float locationValue = locationListTemp[i].objectProbability / distancesTable[locationListSorted[locationListSorted.size()-1] ] [ locationListTemp[i].ID ];
      
      if (locationValue>nextLocationValue) {
	nextLocationValue=locationValue;
	nextLocation=locationListTemp[i].ID;
      }
    }
    locationListSorted.push_back(nextLocation);
    
    vector< vector<int> > locationListInitial;
    locationListInitial.push_back(locationListSorted);
    makeRouteMurrieta(distancesTable, locationListInitial, locationList, locationListSorted, maxNodeSize, 1);
  }
  
  vector<locationNode> listSortedTemp;
  for (int i=0; i<locationListSorted.size(); i++) {
    //cout << locationListSorted[i] << " - " << locationList[locationListSorted[i]].ID << ",";
    listSortedTemp.push_back( locationList[locationListSorted[i]] );
  }
  //cout << endl;
  locationList=listSortedTemp;
}



//Tree manage for Monte Carlo Tree Search (MCTS) UCT
typedef struct nodeUCT nodeUCT;

struct nodeUCT {
  int ID;
  int N;
  float Q;
  nodeUCT * child; //Child
  nodeUCT * nextChild; //Brothers
  nodeUCT * parent;
};

float maxDistanceUCT;


nodeUCT* newChild(nodeUCT *parent, int ID, int N=0, float Q=0) {
  nodeUCT * chld = new nodeUCT;
  chld->ID = ID;
  chld->N = N;
  chld->Q = Q;
  chld->parent = parent;
  chld->child = NULL;
  chld->nextChild = NULL;
  if (parent->child == NULL) {
    parent->child = chld;
  }
  else {
    nodeUCT * tmp = parent->child;
    while(tmp->nextChild != NULL) {
      tmp = tmp->nextChild;
    }
    tmp->nextChild = chld;
  }
  return chld;
}

void showTree(nodeUCT* rootNode, int level=0) {
  if (rootNode!=NULL) {
    cout << " <" <<level<< "> " << rootNode->ID << " ";
    rootNode=rootNode->child;
    while(rootNode != NULL) {
      showTree(rootNode, level+1);
      rootNode = rootNode->nextChild;
    }
  }
}

void deleteTree(nodeUCT* rootNode) {
  if (rootNode!=NULL) {
    nodeUCT * tmp=rootNode;
    rootNode=rootNode->child;
    delete tmp;
    while(rootNode != NULL) {
      deleteTree(rootNode);
      rootNode = rootNode->nextChild;
    }
  }
}

void deleteNode(int ID,vector<locationNode> &locationList) {
  int i;
  for (i=0; i<locationList.size() && locationList[i].ID!=ID; i++) {
  }
  if (i<locationList.size()) {
    locationList.erase(locationList.begin()+i);
  }
}

bool fullyExpanded(nodeUCT* v, int branchSize) {
  //Check for deep level
  int deep=0;
  nodeUCT * parent = v;
  while (parent!=NULL) {
    deep++;
    parent = parent->parent;
  }
  //and now the number of children
  int children=0;
  nodeUCT * child=v->child;
  while (child!=NULL) {
    children++;
    child = child->nextChild;
  }
  return (branchSize-deep)==children;
}

nodeUCT* bestChild(nodeUCT *v, float c=1/sqrt(2)) {
  nodeUCT* child = v->child;
  float maxUCB1=-1;
  nodeUCT* maxUBC1Child;
  while (child!=NULL) {
    float UCB1=0;
    float N = (child->N!=0)?child->N:0.00001;  //Nodes never explored should be visited before but it's division by N=0
    UCB1= (child->Q)/N + c * sqrt( 2 * log( v->N)/N );
    if (UCB1>maxUCB1) {
      maxUCB1=UCB1;
      maxUBC1Child=child;
    }
    child = child->nextChild;
  }
  return maxUBC1Child;
}

nodeUCT* expand(nodeUCT *v,vector<locationNode> locationList) {
  vector<int> nodesInlist;
  nodeUCT* child = v->child;
  //Don't expand the same child
  while (child!=NULL) {
    nodesInlist.push_back(child->ID);
    child = child->nextChild;
  }
  //Neither its parents
  nodeUCT* parent = v;
  while (parent!=NULL) {
    nodesInlist.push_back(parent->ID);
    parent = parent->parent;
  }
  vector<locationNode> nodesNotInList;
  notInList(nodesInlist, locationList, nodesNotInList);
  nodeUCT* p = newChild(v, nodesNotInList[0].ID);
  return v;
}

void backup(nodeUCT *v, float delta) {
  while (v!=NULL) {
    v->N=v->N+1;
    //Best nodes with average min distance
    v->Q=v->Q+delta;
    v = v->parent;
  }
}

nodeUCT* treePolicy(nodeUCT *v, int deep, vector<locationNode> locationList) {
  while (deep>1) {
    if ( !fullyExpanded(v,locationList.size()) ) { 
      return expand(v, locationList);
    }
    else {
      if (v->child==NULL) {  //Tree shortest than max deep exploration
	deep=1;
      }
      else {
	v = bestChild(v);
      }
      deep--;
    }
  }
  return v;
}

float defaultPolicy(nodeUCT *v, vector<locationNode> locationList, MatVector distancesTable, vector<locationNode> locationListComplete) {
  vector<int> nodesInList;
  while (v!=NULL) {
    nodesInList.push_back(v->ID);
    v = v->parent;
  }
  reverse( nodesInList.begin(), nodesInList.end() );
  
  vector<locationNode> nodesNotInList;
  notInList(nodesInList, locationList, nodesNotInList);
  
  if (nodesNotInList.size()>0) {
    //Generate a random route until a leaf node
    for (int i=0; i<nodesNotInList.size(); i++) {
      nodesNotInList[i].priority = float(rand()) / RAND_MAX;
    }
    sort(nodesNotInList.begin(), nodesNotInList.end(), wayToSortLocations);
    
    for (int i=0; i<nodesNotInList.size(); i++) {
      nodesInList.push_back(nodesNotInList[i].ID);
    }
  }
  
  float accumulatedDistance=locationListComplete[ nodesInList[0] ].area;
  float expectedDistance=accumulatedDistance * (locationListComplete[ nodesInList[0] ].objectProbability);
  for (int i=1; i<nodesInList.size(); i++) {
    accumulatedDistance += distancesTable[ nodesInList[i-1] ][ nodesInList[i] ] +locationListComplete[ nodesInList[i] ].area;
    expectedDistance+=accumulatedDistance * (locationListComplete[ nodesInList[i] ].objectProbability);
  }
  return 1/(1+expectedDistance); //Shortest distance, biggest reward
}


void generateRouteUCT(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  //Parameters UCT
  const int maxIterations=100000;
  int maxDeep=5; //Min 2 deep beacause it must explore children
  if (locationList.size()>10) maxDeep=10; //If there are many rooms take more depth

  int maxRouteLength=locationList.size()-1; //-1?
  vector<locationNode> locationListBackup=locationList;
  vector<locationNode> locationListOrdered;
  vector<int> locationListOrderedIndex;
  locationListOrderedIndex.push_back(startingNode);
  int routeLength=0;
  while (routeLength<maxRouteLength) {

    nodeUCT* v0 = new nodeUCT;
    v0->ID = startingNode; v0->N = 0; v0->Q = 0; v0->child=NULL; v0->nextChild=NULL; v0->parent=NULL;
    int numIterations=0;
    while (numIterations<maxIterations) {
      nodeUCT* v1 = treePolicy(v0,maxDeep,locationList);
      float delta = defaultPolicy(v1,locationList, distancesTable, locationListBackup);
      
      backup(v1,delta);
      numIterations++;
    }
    //Save the complete route in the tree
    nodeUCT* v=v0;
    while (v->child!=NULL) {
      deleteNode(v->ID,locationList);
      v=bestChild(v,0);
      
      locationListOrderedIndex.push_back(v->ID);
      routeLength++;
    }
    startingNode=locationListOrderedIndex[locationListOrderedIndex.size()-1];
    
    deleteTree(v0);
  }
  for (int i=0; i<locationListOrderedIndex.size(); i++) {
    //cout << " Index " << locationListOrderedIndex[i] << "  ";
    locationListOrdered.push_back( locationListBackup[ locationListOrderedIndex[i] ] );
  }
  locationList=locationListOrdered;
}



void generateRouteExpectedDistanceSimple(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  vector<locationNode> locationListSorted;
  locationListSorted.push_back(locationList[startingNode]);
  locationList.erase(locationList.begin()+startingNode);
  int nextLocation=startingNode;
  int totalLocations=locationList.size();
  for (int i=0; i<totalLocations; i++) {
    float nextLocationValue=-1;
    for (int j=0; j<locationList.size(); j++) {
      //float locationValue = locationList[j].objectProbability/(distancesTable[startingNode][locationList[j].ID]); //Not in-room distance considered
      float locationValue = locationList[j].objectProbability/(distancesTable[startingNode][locationList[j].ID] * locationList[j].area);
      if (locationValue>nextLocationValue) {
	nextLocationValue=locationValue;
	nextLocation=j;
      }
    }
    startingNode=locationList[nextLocation].ID;
    locationListSorted.push_back(locationList[nextLocation]);
    locationList.erase(locationList.begin()+nextLocation);
  }
  locationList=locationListSorted;
}


//*******************************************************************
void makeRouteHeuristic(MatVector &distancesTable, int node, vector<locationNode> locationList, vector<locationNode> listSorted, int maxDepth) {
  if (maxDepth>1 && locationList.size()>1) {
    listSorted.push_back(locationList[node]);
    locationList.erase(locationList.begin()+node);
    for (int i=0; i<locationList.size(); i++) {
      makeRouteHeuristic(distancesTable, i, locationList, listSorted, maxDepth-1);
    }
  }
  else {
    listSorted.push_back(locationList[node]);
    float accumulatedDistance=listSorted[0].area;
    float heuristicDistance=0; // No heuristic for the first node because is the same for all
    for (int i=1; i<listSorted.size(); i++) {
      accumulatedDistance+=distancesTable[listSorted[i-1].ID][listSorted[i].ID] + listSorted[i].area;
      //heuristicDistance+=listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*listSorted[i].area );
      heuristicDistance+=listSorted[i].objectProbability/ (accumulatedDistance*listSorted[i].area );
    }
    if (heuristicDistance>partialOptimalHeuristicDistance) {
      partialOptimalHeuristicDistance=heuristicDistance;
      partialOptimalHeuristicRoute=listSorted;
    }
    /*
    cout << " Debug " << endl;
    for (int i=0; i<listSorted.size(); i++) {
      cout << listSorted[i].ID << ", ";
    }
    cout << endl;
    */
  }
}

void generateRoutePartialOptimalHeuristic(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  //Parameters Partial Optimal Heuritica
  int maxDepth=4;
  
  int maxRouteLength=locationList.size();
  
  vector<locationNode> locationListPartial=locationList; //Cambiar el nodo 0 por startingNode
  locationNode tmpNode = locationListPartial[startingNode];
  locationListPartial.erase(locationListPartial.begin()+startingNode);
  locationListPartial.insert(locationListPartial.begin(), tmpNode);
  
  vector<locationNode> listSorted;
  
  int k=0;
  int routeLength=0;
  while (routeLength<maxRouteLength) {
    vector<locationNode> listSortedTMP;
    partialOptimalHeuristicDistance=0;
    makeRouteHeuristic(distancesTable,0,locationListPartial,listSortedTMP, maxDepth);
    
    /*
    cout << "partialOptimalHeuristicRoute: ";
    for (int i=0; i<partialOptimalHeuristicRoute.size(); i++) cout << partialOptimalHeuristicRoute[i].ID << ",";
    cout << endl;
    */
    if (listSorted.size()>0) listSorted.erase(listSorted.end()-1);
    
    listSorted.insert(listSorted.end(), partialOptimalHeuristicRoute.begin(), partialOptimalHeuristicRoute.end());
    vector<int> listSortedIndex;
    for (int i=0; i<listSorted.size(); i++) listSortedIndex.push_back(listSorted[i].ID);
    notInList(listSortedIndex, locationList, locationListPartial);
    locationListPartial.insert(locationListPartial.begin(), partialOptimalHeuristicRoute[partialOptimalHeuristicRoute.size()-1]);
    
    /*
    cout << "locationListPartial: ";
    for (int i=0; i<locationListPartial.size(); i++) cout << locationListPartial[i].ID << ",";
    cout << endl;
    */
    //routeLength+=partialOptimalHeuristicRoute.size();
    routeLength=listSorted.size();
    
    /*
    cout << "ListSorted: ";
    for (int i=0; i<listSorted.size(); i++) cout << listSorted[i].ID << ",";
    cout << endl << routeLength << endl;
    */
    //k++;
    //if (k==7) terminate();
  }
  //terminate();
  locationList=listSorted;
  return;
}
//********************************************************************


void makeRoute(MatVector &distancesTable, int node, vector<locationNode> locationList,vector<locationNode> listSorted) {
  if (locationList.size()>1) {
    listSorted.push_back(locationList[node]);
    locationList.erase(locationList.begin()+node);
    for (int i=0; i<locationList.size(); i++) {
      makeRoute(distancesTable, i, locationList, listSorted);
    }
  }
  else {
    listSorted.push_back(locationList[node]);
    float accumulatedDistance=listSorted[0].area;
    float expectedDistance=accumulatedDistance*listSorted[0].objectProbability;
    float heuristicDistance=0; // No heuristic for the first node because is the same for all
    for (int i=1; i<listSorted.size(); i++) {
      accumulatedDistance+=distancesTable[listSorted[i-1].ID][listSorted[i].ID] + listSorted[i].area;
      expectedDistance+=accumulatedDistance*listSorted[i].objectProbability;
      //heuristicDistance+=listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*listSorted[i].area );
      heuristicDistance+=listSorted[i].objectProbability/ (accumulatedDistance*listSorted[i].area );
    }
    if (accumulatedDistance<optimalDistance) {
      optimalDistance=accumulatedDistance;
      optimalRoute=listSorted;
    }
    if (accumulatedDistance>worstDistance) {
      worstDistance=accumulatedDistance;
      worstRoute=listSorted;
    }
    if (expectedDistance<optimalExpectedDistance) {
      optimalExpectedDistance=expectedDistance;
      optimalExpectedRoute=listSorted;
    }
    if (expectedDistance>worstExpectedDistance) {
      worstExpectedDistance=expectedDistance;
      worstExpectedRoute=listSorted;
    }
    if (heuristicDistance>optimalHeuristicDistance) {
      optimalHeuristicDistance=heuristicDistance;
      optimalHeuristicRoute=listSorted;
    }
  }
}

void generateRouteOptimalAndWorst(vector<locationNode> &locationListOptimal, vector<locationNode> &locationListWorst, vector<locationNode> &locationListOptimalExpected, vector<locationNode> &locationListWorstExpected, vector<locationNode> &locationListOptimalHeuristic, MatVector distancesTable, int startingNode) {
  vector<locationNode> listSorted;
  
  makeRoute(distancesTable,startingNode,locationListOptimal,listSorted);
  locationListOptimal=optimalRoute;
  locationListWorst=worstRoute;
  locationListOptimalExpected=optimalExpectedRoute;
  locationListWorstExpected=worstExpectedRoute;
  locationListOptimalHeuristic=optimalHeuristicRoute;
  return;
}

void determineObjectPosition(vector<locationNode> locationList, int &objectPosition) {
  float weightedProbability = float(rand()) / RAND_MAX;
  float accumulatedProbability=0.0;
  int indexLocationObject=0;
  
  while (weightedProbability>accumulatedProbability+locationList[indexLocationObject].objectProbability) {
    accumulatedProbability=accumulatedProbability+locationList[indexLocationObject].objectProbability;
    indexLocationObject++;
  }
  objectPosition=locationList[indexLocationObject].ID;
}

void evaluateRoute(vector<locationNode> locationList, MatVector distancesTable, int objectPosition, int &rooms, float &distance) {
  rooms=1;
  distance=0;
  distance=locationList[0].area;  //If in-room distance is considered, add the starting room area
  for (int i=1; i<locationList.size() && locationList[i-1].ID!=objectPosition; i++) {
    //distance+=distancesTable[locationList[i-1].ID][locationList[i].ID]; //No in-room distance considered
    distance+=distancesTable[locationList[i-1].ID][locationList[i].ID]+locationList[i].area;
    rooms++;
  }
}

float evaluateRouteExpectedDistance(vector<locationNode> locationList, MatVector distancesTable) {
  float accumulatedDistance=locationList[0].area;
  float expectedDistance=accumulatedDistance*locationList[0].objectProbability;
  for (int i=1; i<locationList.size(); i++) {
    accumulatedDistance+=distancesTable[locationList[i-1].ID][locationList[i].ID] + locationList[i].area;
    expectedDistance+=accumulatedDistance*locationList[i].objectProbability;
  }
  return expectedDistance;
}


int main(int argc, char **argv) {
  
  mapDataBase="mapsDB/";
  DISTANCES_MEAN=1000.0;
  DISTANCES_SIGMA=1000.0;
  SELECTED_MAP=5; //Map to be evaluated
  string probabilityDistribution="Exponential"; // Uniform, Normal, Gamma, Exponential
  int minNodes=3;	//Minimum 3 Nodes
  int maxNodes=13;	//it has been tested a optimal policy to 14 nodes in 2 days of computing time. 13 nodes in 2 hours. 12 nodes 10 minutes
  int startingNode=0;
  int numRepetitions=10; //1000 is near the expected value
  string environmentType="H";	//Type: House(H), Office (O), Random (R)
  
  
  //Convert configuration in output filename
  string outputFileNameRunnings;
  string outputFileNameData;
  string outputFileNameExpected;
  ostringstream num2FormatString;
  num2FormatString << setfill('0') << setw(2) << minNodes << "-" << setfill('0') << setw(2) << maxNodes << "maps" << setfill('0') << setw(2) << SELECTED_MAP << "prob" << probabilityDistribution[0] << setfill('0') << setw(2);
  outputFileNameRunnings="experiment-" + environmentType + num2FormatString.str()+"-Running.txt";
  outputFileNameData="experiment-" + environmentType + num2FormatString.str()+"-Data.txt";
  outputFileNameExpected="experiment-" + environmentType + num2FormatString.str()+"-Expected.txt";
  cout << "Experiments robot runnings file " << outputFileNameRunnings << endl;
  cout << "Experiments maps and orders data file " << outputFileNameData << endl;
  cout << "Experiments expected distances file " << outputFileNameExpected << endl;
  
  ofstream fileStreamData(outputFileNameData.c_str());
  
  ofstream fileStream(outputFileNameRunnings.c_str());
  fileStream << "Nodes,Repetition,System Time,Object at Room,Starting At Node,Shortest Distance,Worst Distance,";
  fileStream << "Rooms Explored (Optimal),Traveled Distance (Optimal),";
  fileStream << "Rooms Explored (Worst),Traveled Distance (Worst),";
  fileStream << "Rooms Explored (OptimalExpected),Traveled Distance (OptimalExpected),";
  fileStream << "Rooms Explored (WorstExpected),Traveled Distance (WorstExpected),";
  fileStream << "Rooms Explored (OptimalExpectedHeuristic),Traveled Distance (OptimalExpectedHeuristic),";
  fileStream << "Rooms Explored (MostProbable),Traveled Distance (MostProbable),";
  fileStream << "Rooms Explored (Closest),Traveled Distance (Closest),";
  fileStream << "Rooms Explored (Random),Traveled Distance (Random),";
  fileStream << "Rooms Explored (MCTS-UCT),Traveled Distance (MCTS-UCT),";
  fileStream << "Rooms Explored (ExpectedDistanceFixed),Traveled Distance (ExpectedDistanceFixed),";
  fileStream << "Rooms Explored (ExpectedDistanceMurrieta),Traveled Distance (ExpectedDistanceMurrieta),";
  fileStream << "Rooms Explored (ExpectedDistanceSimple),Traveled Distance (ExpectedDistanceSimple),";
  fileStream << "Rooms Explored (PartialOptimalHeuristic),Traveled Distance (PartialOptimalHeuristic)";
  fileStream << endl;
  fileStream.close();
  
  ofstream fileStreamExpected(outputFileNameExpected.c_str());
  fileStreamExpected << "Map Type,Map Number,Nodes,Probability Distribution,Starting At Node,Shortest Distance,Worst Distance,";
  fileStreamExpected << "Optimal,";
  fileStreamExpected << "Worst,";
  fileStreamExpected << "OptimalExpected,";
  fileStreamExpected << "WorstExpected,";
  fileStreamExpected << "OptimalExpectedHeuristic,";
  fileStreamExpected << "MostProbable,";
  fileStreamExpected << "Closest,";
  fileStreamExpected << "Random,";
  fileStreamExpected << "MCTS-UCT,";
  fileStreamExpected << "ExpectedDistanceFixed,";
  fileStreamExpected << "ExpectedDistanceMurrieta,";
  fileStreamExpected << "ExpectedDistanceSimple,";
  fileStreamExpected << "PartialOptimalHeuristic";
  fileStreamExpected << endl;
  
  for (int numNodes=minNodes; numNodes<=maxNodes; numNodes++) {
    cout << "Environment with " << numNodes << endl;
    fileStreamData << "Environment with " << numNodes << endl;
    
    vector<locationNode> locationList;
    MatVector distancesTable;
    
    fileStreamData << buildEnvironment(locationList, distancesTable, numNodes, environmentType, probabilityDistribution);
    
    optimalDistance=maxDistance(distancesTable,locationList)*(numNodes+1);
    optimalExpectedDistance=optimalDistance;
    optimalExpectedDistanceMurrieta=optimalDistance; //Para el metodo del Dr. Murrieta
    maxDistanceUCT=optimalDistance; //Para standarizar en [0,1] las distancias
    worstDistance=0;
    worstExpectedDistance=0;
    optimalHeuristicDistance=0;
    partialOptimalHeuristicDistance=0;
    
    cout << endl << "Building strategies..." << endl;
    vector<locationNode> locationListOptimal, locationListWorst, locationListOptimalExpected, locationListWorstExpected, locationListOptimalHeuristic, locationListMostProbable, locationListClosest, locationListRandom, locationListUCT, locationListExpectedDistanceFixed, locationListExpectedDistanceMurrieta, locationListExpectedDistanceSimple, locationListPartialOptimalHeuristic;
    locationListOptimal=locationList;
    locationListWorst=locationList;
    locationListOptimalExpected=locationList;
    locationListWorstExpected=locationList;
    locationListOptimalHeuristic=locationList;
    locationListMostProbable=locationList;
    locationListClosest=locationList;
    locationListRandom=locationList;
    locationListUCT=locationList;
    locationListExpectedDistanceFixed=locationList;
    locationListExpectedDistanceMurrieta=locationList;
    locationListExpectedDistanceSimple=locationList;
    locationListPartialOptimalHeuristic=locationList;
    
    
    time_t initTime;
    initTime=time(NULL);
    generateRouteOptimalAndWorst(locationListOptimal, locationListWorst, locationListOptimalExpected, locationListWorstExpected, locationListOptimalHeuristic ,distancesTable, startingNode);
    
    //Save in file Optimal and Worst routes
    fileStreamExpected << environmentType << ",";
    fileStreamExpected << SELECTED_MAP << ",";
    fileStreamExpected << numNodes << ",";
    fileStreamExpected << probabilityDistribution << ",";
    fileStreamExpected << startingNode << ",";
    int tmpRooms; float tmpDistance;
    evaluateRoute(locationListOptimal, distancesTable, locationListOptimal[numNodes-1].ID, tmpRooms, tmpDistance);
    fileStreamExpected << optimalDistance << ",";
    evaluateRoute(locationListWorst, distancesTable, locationListWorst[numNodes-1].ID, tmpRooms, tmpDistance);
    fileStreamExpected << tmpDistance << ",";
    
    fileStreamExpected << evaluateRouteExpectedDistance(locationListOptimal, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListOptimal,"Optimal",initTime);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListWorst, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListWorst,"Worst",initTime);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListOptimalExpected, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListOptimalExpected,"Optimal Expected",initTime);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListWorstExpected, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListWorstExpected,"Worst Expected",initTime);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListOptimalHeuristic, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListOptimalHeuristic,"Optimal Heuristic",initTime);
    
    initTime=time(NULL);
    generateRouteMostProbable(locationListMostProbable, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListMostProbable, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListMostProbable,"Most Probable",initTime);
    initTime=time(NULL);
    generateRouteClosest(locationListClosest, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListClosest, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListClosest,"Closest",initTime);
    initTime=time(NULL);
    generateRouteRandom(locationListRandom, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListRandom, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListRandom,"Random",initTime);
    initTime=time(NULL);
    generateRouteUCT(locationListUCT, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListUCT, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListUCT,"Monte Carlo Tree Search - UCT",initTime);
    initTime=time(NULL);
    generateRouteExpectedDistanceFixed(locationListExpectedDistanceFixed, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListExpectedDistanceFixed, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListExpectedDistanceFixed,"Expected Distance Fixed",initTime);
    initTime=time(NULL);
    generateRouteExpectedDistanceMurrieta(locationListExpectedDistanceMurrieta, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListExpectedDistanceMurrieta, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListExpectedDistanceMurrieta,"Expected Distance Murrieta",initTime);
    initTime=time(NULL);
    generateRouteExpectedDistanceSimple(locationListExpectedDistanceSimple, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListExpectedDistanceSimple, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListExpectedDistanceSimple,"Expected Distance Simple",initTime);
    initTime=time(NULL);
    generateRoutePartialOptimalHeuristic(locationListPartialOptimalHeuristic, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListPartialOptimalHeuristic, distancesTable) << endl;
    fileStreamData << msgLocationOrder(locationListPartialOptimalHeuristic,"Partial Optimal Heuristic",initTime);
    
    cout << endl << "Begin of test!" << endl;
    for (int repetition=0; repetition<numRepetitions; repetition++) {
      
      cout << ".";
      //cout << "Repetition " << repetition << "  ";
      //startingNode=int((numNodes-1)*float(rand()) / RAND_MAX); //Change initial node
      
      int objectPosition=0;
      determineObjectPosition(locationList, objectPosition);
      
      //cout << "Object in room " << objectPosition << endl;
      
      fileStream.open(outputFileNameRunnings.c_str(), ios::app);
      fileStream << numNodes << ",";
      fileStream << repetition << ",";
      fileStream << time(NULL) << ",";
      fileStream << objectPosition << ",";
      fileStream << startingNode << ",";
      fileStream << optimalDistance << ",";
      fileStream << worstDistance << ",";
      
      int rooms;
      float distance;
      evaluateRoute(locationListOptimal, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListWorst, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListOptimalExpected, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListWorstExpected, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListMostProbable, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListClosest, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListRandom, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListUCT, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListExpectedDistanceFixed, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListExpectedDistanceMurrieta, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListExpectedDistanceSimple, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance << ",";
      evaluateRoute(locationListPartialOptimalHeuristic, distancesTable, objectPosition, rooms, distance);
      fileStream << rooms << "," << distance;
      
      fileStream << endl;
      fileStream.close();
    }
    
  }
  
  fileStreamExpected << endl;
  fileStreamExpected.close();
  
  fileStreamData << endl;
  fileStreamData.close();
  
  cout << endl << "End of test!" << endl;
  return 0;
}
