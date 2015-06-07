#include <vector>
#include <algorithm>	//sort

#define PI 3.14159265

using namespace std;

typedef vector<double> doubleVector;
typedef vector<doubleVector> MatVector;


struct locationNode {
  int ID;
  double objectProbability;
  double priority;
  double area;
};


double optimalExpectedDistance;
vector<locationNode> optimalExpectedRoute;
double worstExpectedDistance;
vector<locationNode> worstExpectedRoute;
double optimalDistance;
vector<locationNode> optimalRoute;
double worstDistance;
vector<locationNode> worstRoute;
double optimalHeuristicDistance;
vector<locationNode> optimalHeuristicRoute;
double partialOptimalHeuristicDistance;
vector<locationNode> partialOptimalHeuristicRoute;




//*****************************

bool wayToSortLocations(locationNode i, locationNode j) {
  return i.priority > j.priority;
}

double maxDistance(MatVector distancesTable, vector<locationNode> locationList) {
  double maxDist=0;
  int numNodes=distancesTable.size();
  for (int i=0; i<numNodes; i++) {
    for (int j=0; j<numNodes; j++) {
      if (distancesTable[i][j]>maxDist) maxDist=distancesTable[i][j];
    }
  }
  double maxArea=0;
  numNodes=locationList.size();
  for (int i=0; i<numNodes; i++) {
    if (locationList[i].area>maxArea) maxArea=locationList[i].area;
  }
  return maxDist+maxArea;
}


void generateRouteMostProbable(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  double maxPriority = 0;
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = locationList[i].objectProbability;
    maxPriority += locationList[i].priority;
  }
  locationList[startingNode].priority = maxPriority;	//The highest priority
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}

void generateRouteClosestFixed(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  double maxPriority = 0;
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = 1/distancesTable[startingNode][i];
    maxPriority += locationList[i].priority;
  }
  locationList[startingNode].priority = maxPriority;	//The highest priority
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
    double nextLocationValue=0;
    for (int j=0; j<locationList.size(); j++) {
      double locationValue = 1/distancesTable[startingNode][locationList[j].ID];
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
  double maxPriority = 0;
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = double(rand()) / RAND_MAX;
    maxPriority += locationList[i].priority;
  }
  locationList[startingNode].priority = maxPriority;	//The highest priority
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}

void generateRouteExpectedDistanceFixed(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  double maxPriority = 0;
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].priority = locationList[i].objectProbability/distancesTable[startingNode][i];
    maxPriority += locationList[i].priority;
  }
  locationList[startingNode].priority = maxPriority;	//The highest priority
  sort(locationList.begin(), locationList.end(), wayToSortLocations);
}



void nonDominatedHull(int referenceNode, vector<locationNode> locationList, MatVector distancesTable, vector<locationNode> &locationListHull) {
  
  vector<locationNode> locationListHullTemp;
  double minDistance=distancesTable[referenceNode][locationList[0].ID];
  double maxProbabilityAtMinDistance=locationList[0].objectProbability;
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
  
  double maxAngle=0.0; //any value for starting cicle of searching convex hull
  int maxAngleNode=0;
  while (maxAngle>=0) {
    locationListHullTemp.push_back(locationList[minNode]);
    locationList.erase(locationList.begin()+minNode);
    maxAngle=(-3)*PI; //Any lowest value for searching the maximun angle
    for (int i=0; i<locationList.size(); i++) {
      double angle=atan2( locationList[i].objectProbability-maxProbabilityAtMinDistance, distancesTable[referenceNode][locationList[i].ID]-minDistance);  //atan2(y,x)
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


double optimalExpectedDistanceMurrieta;

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
    double optimalExpectedDistanceLocal=optimalExpectedDistanceMurrieta;
    for (int i=0; i<newLocationListPrevious.size(); i++) {
      
      double accumulatedDistance=sqrt(locationList[ newLocationListPrevious[i][0] ].area);
      double expectedDistance=accumulatedDistance * locationList[ newLocationListPrevious[i][0] ].objectProbability;
      for (int j=1; j<newLocationListPrevious[i].size(); j++) {
	accumulatedDistance+=distancesTable[ newLocationListPrevious[i][j-1] ][ newLocationListPrevious[i][j] ] + sqrt(locationList[ newLocationListPrevious[i][j] ].area);
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
    double nextLocationValue=0;
    for (int i=0; i<locationListTemp.size(); i++) {
      double locationValue = locationListTemp[i].objectProbability / distancesTable[locationListSorted[locationListSorted.size()-1] ] [ locationListTemp[i].ID ];
      
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
  double Q;
  nodeUCT * child; //Child
  nodeUCT * nextChild; //Brothers
  nodeUCT * parent;
};

double maxDistanceUCT;


nodeUCT* newChild(nodeUCT *parent, int ID, int N=0, double Q=0) {
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

nodeUCT* bestChild(nodeUCT *v, double c=1/sqrt(2)) {
  nodeUCT* child = v->child;
  double maxUCB1=-1;
  nodeUCT* maxUBC1Child;
  while (child!=NULL) {
    double UCB1=0;
    double N = (child->N!=0)?child->N:0.00001;  //Nodes never explored should be visited before but it's division by N=0
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

void backup(nodeUCT *v, double delta) {
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

double defaultPolicy(nodeUCT *v, vector<locationNode> locationList, MatVector distancesTable, vector<locationNode> locationListComplete) {
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
      nodesNotInList[i].priority = double(rand()) / RAND_MAX;
    }
    sort(nodesNotInList.begin(), nodesNotInList.end(), wayToSortLocations);
    
    for (int i=0; i<nodesNotInList.size(); i++) {
      nodesInList.push_back(nodesNotInList[i].ID);
    }
  }
  
  double accumulatedDistance=sqrt(locationListComplete[ nodesInList[0] ].area);
  double expectedDistance=accumulatedDistance * (locationListComplete[ nodesInList[0] ].objectProbability);
  for (int i=1; i<nodesInList.size(); i++) {
    accumulatedDistance += distancesTable[ nodesInList[i-1] ][ nodesInList[i] ] + sqrt(locationListComplete[ nodesInList[i] ].area);
    expectedDistance+=accumulatedDistance * (locationListComplete[ nodesInList[i] ].objectProbability);
  }
  return 1/(1+expectedDistance); //Shortest distance, biggest reward
}


void generateRouteUCT(vector<locationNode> &locationList, MatVector distancesTable, int startingNode) {
  //Parameters UCT
  const int maxIterations=100; //10000 $$$$###################################$$$$$$$$$%%%%%%%%%%%%%%%%&&&&&&&&&&&&&
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
      double delta = defaultPolicy(v1,locationList, distancesTable, locationListBackup);
      
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
    double nextLocationValue=-1;
    for (int j=0; j<locationList.size(); j++) {
      double locationValue = locationList[j].objectProbability/(distancesTable[startingNode][locationList[j].ID]); //Not in-room distance considered
      //double locationValue = locationList[j].objectProbability/(distancesTable[startingNode][locationList[j].ID] * (locationList[j].area) ); //Not root squared room area
      //double locationValue = locationList[j].objectProbability/(distancesTable[startingNode][locationList[j].ID] * sqrt(locationList[j].area) );
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
    double accumulatedDistance=sqrt(listSorted[0].area);
    double accumulatedProbability=listSorted[0].objectProbability;  //To test the new heuristic
    double heuristicDistance=1; // No heuristic for the first node because is the same for all
    for (int i=1; i<listSorted.size(); i++) {
      accumulatedDistance+=distancesTable[listSorted[i-1].ID][listSorted[i].ID] + sqrt(listSorted[i].area);
      heuristicDistance+=listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*sqrt(listSorted[i].area) );
      
      //Old heuristic
      //heuristicDistance+=listSorted[i].objectProbability/ (accumulatedDistance*sqrt(listSorted[i].area) );
      
      //To test the new heuristic
      //accumulatedProbability+=listSorted[i].objectProbability;
      //heuristicDistance+=accumulatedProbability/ (accumulatedDistance*sqrt(listSorted[i].area) );
      
      //heuristicDistance*=listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*sqrt(listSorted[i].area) );
      //heuristicDistance*=listSorted[i].objectProbability/ (accumulatedDistance*sqrt(listSorted[i].area) );
      
    }
    if (heuristicDistance>partialOptimalHeuristicDistance) {
      partialOptimalHeuristicDistance=heuristicDistance;
      partialOptimalHeuristicRoute=listSorted;
    }
    /*
     *   cout << " Debug " << endl;
     *   for (int i=0; i<listSorted.size(); i++) {
     *     cout << listSorted[i].ID << ", ";
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
     *   cout << "partialOptimalHeuristicRoute: ";
     *   for (int i=0; i<partialOptimalHeuristicRoute.size(); i++) cout << partialOptimalHeuristicRoute[i].ID << ",";
     *   cout << endl;
     */
    if (listSorted.size()>0) listSorted.erase(listSorted.end()-1);
    
    listSorted.insert(listSorted.end(), partialOptimalHeuristicRoute.begin(), partialOptimalHeuristicRoute.end());
    vector<int> listSortedIndex;
    for (int i=0; i<listSorted.size(); i++) listSortedIndex.push_back(listSorted[i].ID);
    notInList(listSortedIndex, locationList, locationListPartial);
    locationListPartial.insert(locationListPartial.begin(), partialOptimalHeuristicRoute[partialOptimalHeuristicRoute.size()-1]);
    
    /*
     *   cout << "locationListPartial: ";
     *   for (int i=0; i<locationListPartial.size(); i++) cout << locationListPartial[i].ID << ",";
     *   cout << endl;
     */
    //routeLength+=partialOptimalHeuristicRoute.size();
    routeLength=listSorted.size();
    
    /*
     *   cout << "ListSorted: ";
     *   for (int i=0; i<listSorted.size(); i++) cout << listSorted[i].ID << ",";
     *   cout << endl << routeLength << endl;
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
    double accumulatedDistance=sqrt(listSorted[0].area);  //Aproximatly the perimeter, not the area
    double expectedDistance=accumulatedDistance*listSorted[0].objectProbability;
    double heuristicDistance=0; // No heuristic for the first node because is the same for all
    //double factor=10000000000000000000.0;
    
    for (int i=1; i<listSorted.size(); i++) {
      accumulatedDistance+=distancesTable[listSorted[i-1].ID][listSorted[i].ID] + sqrt(listSorted[i].area);
      expectedDistance+=accumulatedDistance*listSorted[i].objectProbability;
      heuristicDistance+=listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*sqrt(listSorted[i].area) );
      //heuristicDistance+=listSorted[i].objectProbability/ (accumulatedDistance*listSorted[i].area );
      //heuristicDistance*= listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*listSorted[i].area );
      //heuristicDistance+= listSorted[i].objectProbability/ (accumulatedDistance*listSorted[i].area );
      //heuristicDistance*= factor * listSorted[i].objectProbability/ (accumulatedDistance*listSorted[i].area);
      //heuristicDistance*= factor * listSorted[i].objectProbability/ (distancesTable[listSorted[i-1].ID][listSorted[i].ID]*listSorted[i].area);
      
      
      //factor/=10;
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
      //cout << endl << "Heuristic " << heuristicDistance << "Optimal " << optimalHeuristicDistance << "   Size " << optimalHeuristicRoute.size() << endl;
    }
    //cout << endl << "Heuristic " << heuristicDistance << "   Optimal " << optimalHeuristicDistance << "   Size " << optimalHeuristicRoute.size() << endl;
  }
}

void generateRouteOptimalAndWorst(vector<locationNode> &locationListOptimal, vector<locationNode> &locationListWorst, vector<locationNode> &locationListOptimalExpected, vector<locationNode> &locationListWorstExpected, vector<locationNode> &locationListOptimalHeuristic, MatVector distancesTable, int startingNode) {
  vector<locationNode> listSorted;
  
  //CHANGE THIS FOR MAKING FAST EXPERIMENTS. OPTIMAL SOLUTIONS CAN BE DISABLED
  //return; //Comented ENABLE else DISABLED
  
  makeRoute(distancesTable,startingNode,locationListOptimal,listSorted);
  locationListOptimal=optimalRoute;
  locationListWorst=worstRoute;
  locationListOptimalExpected=optimalExpectedRoute;
  locationListWorstExpected=worstExpectedRoute;
  locationListOptimalHeuristic=optimalHeuristicRoute;
  return;
}