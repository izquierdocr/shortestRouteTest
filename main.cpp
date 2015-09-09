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
#include "explorationMethods.cpp"

using namespace std;



double DISTANCES_MEAN;
double DISTANCES_SIGMA;
int selectedMap;
string mapDataBase;

vector<double> areas;


vector <double> extractProbability(string str, int maxRoom) {
  //string str must follow format objectNumber-orderRoom1,orderRoom2,orderRoom3,etc
  vector <double> A;
  size_t separatorPos = str.find(',',0);
  while (separatorPos!=string::npos && separatorPos<str.length()) {
    size_t separatorPos2 = str.find(',',separatorPos+1);
    if (separatorPos2==string::npos) separatorPos2=str.length();
    string numberSTR=str.substr(separatorPos+1,separatorPos2-(separatorPos+1));
    double number = atof(numberSTR.c_str());
    A.push_back( number );
    separatorPos=separatorPos2;
  }
  return A;
}

void loadUserProbabilities(string object, int numRooms, vector<double> &p){
  int numObjects=72;  // ######################################## Beacuse those are in the file
  const char *objectvector[] = {"apple","shoe","coffee","clothe","laptop","bread","pen","book","bed sheet","cellphone","spoon","fork","glass of water","handbag","food","towel","chair","medicine","tv remote","coke","broom","mop","key","scissor","comb","bible","ac remote","cell charger","tablet","plate","backpack","baby bottle","headphone","nail clipper","jacket","hand cream","inhaler","cosmetic bag","fly swatter","pillow","blanket","milk","shirt","sock","cup","glasses","knife","soap","coat","pumpkin","orange","paddle","ball","dinosaur","bottle","toy car","frying pan","cd","dvd","videogame","toy","potato chips","cracker","cookie","extinguisher","phone","printer","potty","bookshelf","trash","fridge","softener"};
  vector<string> objects(objectvector, objectvector+numObjects);
  
  vector<string>::iterator it;
  it = find (objects.begin(), objects.end(), object);
  if (it == objects.end()) {cout << "Object not found in file" << endl; exit (EXIT_FAILURE);}
  size_t tmpPos = it - objects.begin();
  int objectPosition=tmpPos+1;
  
  ifstream inputFile("/home/izquierdocr/projects/internetQuery/data/bordaUserPb.txt");
  string lineA;
  for (int i=0; i<objectPosition; i++) {
    getline(inputFile, lineA); //This do not check if the file is incomplete or unstructured
  }
  p=extractProbability(lineA, numRooms);
  reverse(p.begin(),p.end());
}


//Load from file
void loadRoomProbabilities(string object, int numRooms, vector<double> &p, string fileName){
  int numObjects=72;  // ######################################## Beacuse those are in the file
  const char *objectvector[] = {"apple","shoe","coffee","clothe","laptop","bread","pen","book","bed sheet","cellphone","spoon","fork","glass of water","handbag","food","towel","chair","medicine","tv remote","coke","broom","mop","key","scissor","comb","bible","ac remote","cell charger","tablet","plate","backpack","baby bottle","headphone","nail clipper","jacket","hand cream","inhaler","cosmetic bag","fly swatter","pillow","blanket","milk","shirt","sock","cup","glasses","knife","soap","coat","pumpkin","orange","paddle","ball","dinosaur","bottle","toy car","frying pan","cd","dvd","videogame","toy","potato chips","cracker","cookie","extinguisher","phone","printer","potty","bookshelf","trash","fridge","softener"};
  vector<string> objects(objectvector, objectvector+numObjects);
  
  vector<string>::iterator it;
  it = find (objects.begin(), objects.end(), object);
  if (it == objects.end()) {cout << "Object not found in file" << endl; exit (EXIT_FAILURE);}
  size_t tmpPos = it - objects.begin();
  int objectPosition=tmpPos+1;
  
  ifstream inputFile(fileName);
    string lineA;
  for (int i=0; i<objectPosition; i++) {
    getline(inputFile, lineA); //This do not check if the file is incomplete or unstructured
  }
  p=extractProbability(lineA, numRooms);
}

//Generate syntethic (random) probabilities
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
  cout << "Percentage unused random numbers (out of range):" << double(100*notused)/experiments << endl << endl;
  
  cout << fixed; cout.precision(5);
  double totalP=0;
  for (int i=0; i<nintervals; ++i) {
    if (showGraph) {
      cout << double(i)/nintervals << "-" << double(i+1)/nintervals << ": ";
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
  
  areas.clear();
  for (int i=0; i<roomListSize; i++) {
    int x,y;
    fileStream >> x >> y;
    double area;
    fileStream >> area;
    areas.push_back(area);
    //Load the number of border vertex for each room
    int borderCount=0;
    fileStream >> borderCount;
    //Load the vertexs
    for (int j=0; j<borderCount; j++) {
      int bx, by;
      fileStream >> bx >> by;
    }
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



string buildEnvironment(vector<locationNode> &locationList, MatVector &distancesTable, int numNodes, string mapType, string probabilityDistribution, string object) {
  ostringstream stream2String;
  string msg;
  
  //Initialize random generators
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.0,1.0); 	//numbers in [0,1]
  
  //Pass the random numbers according the map number 
  //so, choosing different map generates different random numbers
  for (int i=0; i<selectedMap*numNodes*numNodes*numNodes; i++) { 
    double dump=distribution(generator);
  }
  
  vector<double> roomProbabilities;
  if (probabilityDistribution!="File")
    generateRoomProbabilities(probabilityDistribution, numNodes, roomProbabilities); //Generated with random
  else {
    loadRoomProbabilities(object, numNodes, roomProbabilities, "/home/izquierdocr/projects/internetQuery/data/experimentQueryWeightedProbabilityWebPb.txt"); //Loaded from files
    //loadRoomProbabilities(object, numNodes, roomProbabilities, "/home/izquierdocr/projects/internetQuery/data/bordaUserPb.txt"); //Loaded from users
  }
  for (int i=0; i<numNodes; i++) {
    locationNode location;
    location.ID=i;
    location.objectProbability=roomProbabilities[numNodes-i-1]; //Only for evaluating cross probabilities. This can be changed to [i]. Choosen only for having other order than 0,1,2...
    location.priority=0;
    locationList.push_back(location);
  }
  
  doubleVector rowMat;
  rowMat.assign(numNodes,0);
  for (int i=0; i<numNodes; i++) {
    distancesTable.push_back(rowMat);
  }
  
  if (mapType=="R") { 	//Type: House(H), Office (O), Random (R)
    for (int i=0; i<numNodes; i++) {
      locationList[i].area = ( DISTANCES_MEAN + DISTANCES_SIGMA * (distribution(generator)-0.5) )/2;  //Assume areas represent in general 1/2 of distances between rooms (Here area defined as robot movement inside rooms
    }
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
    //int numMap = MAX_MAPS*(double(rand()) / RAND_MAX)+1; //Random house map
    int numMap = selectedMap;
    loadMapData(numNodes, numMap, mapType, distancesTable);
    for (int i=0; i<numNodes; i++) {
      locationList[i].area = areas[i];
      //locationList[i].area = 1; //Testing with equal areas
    }
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

//#########################

void determineObjectPosition(vector<locationNode> locationList, int &objectPosition) {
  double weightedProbability = double(rand()) / RAND_MAX;
  double accumulatedProbability=0.0;
  int indexLocationObject=0;
  
  while (weightedProbability>accumulatedProbability+locationList[indexLocationObject].objectProbability) {
    accumulatedProbability=accumulatedProbability+locationList[indexLocationObject].objectProbability;
    indexLocationObject++;
  }
  objectPosition=locationList[indexLocationObject].ID;
}

void evaluateRoute(vector<locationNode> locationList, MatVector distancesTable, int objectPosition, int &rooms, double &distance) {
  rooms=1;
  distance=0;
  distance=sqrt(locationList[0].area);  //If in-room distance is considered, add the starting room area
  for (int i=1; i<locationList.size() && locationList[i-1].ID!=objectPosition; i++) {
    //distance+=distancesTable[locationList[i-1].ID][locationList[i].ID]; //No in-room distance considered
    distance+=distancesTable[locationList[i-1].ID][locationList[i].ID]+sqrt(locationList[i].area);
    rooms++;
  }
}


//evaluateRouteExpectedDistance
double evaluateRouteExpectedDistance(vector<locationNode> locationList, MatVector distancesTable) {
  /*
  //Change probabilities to those generated by users
  //Only for evaluating cross probabilities - Comment in normal use
  vector<double> p;
  loadUserProbabilities("cup", 10, p);
  for (int i=0; i<locationList.size(); i++) {
    locationList[i].objectProbability=p[locationList[i].ID];
  }
  */
  
  
  double accumulatedDistance=sqrt(locationList[0].area);
  double expectedDistance=accumulatedDistance*locationList[0].objectProbability;
  for (int i=1; i<locationList.size(); i++) {
    accumulatedDistance+=distancesTable[locationList[i-1].ID][locationList[i].ID] + sqrt(locationList[i].area);
    expectedDistance+=accumulatedDistance*locationList[i].objectProbability;
  }
  cout << "###Expected Distance: " << expectedDistance << "  ";
  return expectedDistance;
}

/*
//evaluateRouteHeuristic
double evaluateRouteExpectedDistance(vector<locationNode> locationList, MatVector distancesTable) {
  double accumulatedDistance=0;
  //double expectedDistance=accumulatedDistance*locationList[0].objectProbability;
  for (int i=1; i<locationList.size(); i++) {
    accumulatedDistance+=locationList[i].objectProbability/ ( distancesTable[locationList[i-1].ID][locationList[i].ID] * sqrt(locationList[i].area) );
    //expectedDistance+=accumulatedDistance*locationList[i].objectProbability;
  }
  cout << "###Expected Distance: " << accumulatedDistance << "  ";
  return accumulatedDistance;
}
*/

string msgLocationOrder(vector<locationNode> locationList, string sortingMethod, clock_t initTime) {
  ostringstream stream2String;
  string msg;
  
  cout << "Locations ordered by " << sortingMethod << " in " << ((float)(clock()-initTime))/CLOCKS_PER_SEC << " seconds" << endl;
  stream2String << "Locations ordered by " << sortingMethod << " in " << ((float)(clock()-initTime))/CLOCKS_PER_SEC << " seconds" << endl;
  for (int i=0; i<locationList.size(); i++) {
    cout << locationList[i].ID << ", ";
    stream2String << locationList[i].ID << ",";
  }
  cout << endl;
  stream2String << endl;
  
  msg = stream2String.str();
  return msg;
}

int main(int argc, char **argv) {
  
  mapDataBase="mapsDB/";
  DISTANCES_MEAN=1000.0;
  DISTANCES_SIGMA=1000.0;
  int minMap=1;
  int maxMap=5;
  selectedMap=5; //Map to be evaluated
  string probabilityDistribution="Exponential"; // Uniform, Normal, Gamma, Exponential, File
  string object="coke"; // Only useful with File
  int minNodes=10;	//Minimum 3 Nodes
  int maxNodes=10;	//it has been tested a optimal policy to 14 nodes in 2 days of computing time. 13 nodes in 2 hours. 12 nodes 10 minutes
  int minStartingNode=0;
  int maxStartingNode=0;
  int startingNode=0;
  int numRepetitions=10; //1000 is near the expected value
  string environmentType="H";	//Type: House(H), Office (O), Random (R)
  
  //Convert configuration in output filename
  string outputFileNameRunnings;
  string outputFileNameData;
  string outputFileNameExpected;
  ostringstream num2FormatString;
  //"coke","knife","milk","tv remote","cup"
  num2FormatString << "cokeWEB-WEB" << setfill('0') << setw(2) << minNodes << "-" << setfill('0') << setw(2) << maxNodes << "maps" << setfill('0') << setw(2) << minMap << "-" << setfill('0') << setw(2) << maxMap  << "strt" << setfill('0') << setw(2) << minStartingNode << "-" << setfill('0') << setw(2) << maxStartingNode << "prob" << probabilityDistribution[0];
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
  for (selectedMap=minMap; selectedMap<=maxMap; selectedMap++) {
  for (startingNode=minStartingNode; startingNode<=maxStartingNode && startingNode<numNodes; startingNode++) {
    cout << "Environment with " << numNodes << endl;
    fileStreamData << "Environment with " << numNodes << " rooms" << endl;
    
    vector<locationNode> locationList;
    MatVector distancesTable;
    
    fileStreamData << buildEnvironment(locationList, distancesTable, numNodes, environmentType, probabilityDistribution, object);
    
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
    
    
    clock_t initTime;
    initTime=clock();
    generateRouteOptimalAndWorst(locationListOptimal, locationListWorst, locationListOptimalExpected, locationListWorstExpected, locationListOptimalHeuristic ,distancesTable, startingNode);
    
    //Save in file Optimal and Worst routes
    fileStreamExpected << environmentType << ",";
    fileStreamExpected << selectedMap << ",";
    fileStreamExpected << numNodes << ",";
    fileStreamExpected << probabilityDistribution << ",";
    fileStreamExpected << startingNode << ",";
    int tmpRooms; double tmpDistance;
    evaluateRoute(locationListOptimal, distancesTable, locationListOptimal[numNodes-1].ID, tmpRooms, tmpDistance);
    fileStreamExpected << tmpDistance << ",";
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
    
    initTime=clock();
    generateRouteMostProbable(locationListMostProbable, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListMostProbable, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListMostProbable,"Most Probable",initTime);
    initTime=clock();
    generateRouteClosest(locationListClosest, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListClosest, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListClosest,"Closest",initTime);
    initTime=clock();
    generateRouteRandom(locationListRandom, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListRandom, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListRandom,"Random",initTime);
    initTime=clock();
    generateRouteUCT(locationListUCT, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListUCT, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListUCT,"Monte Carlo Tree Search - UCT",initTime);
    initTime=clock();
    generateRouteExpectedDistanceFixed(locationListExpectedDistanceFixed, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListExpectedDistanceFixed, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListExpectedDistanceFixed,"Expected Distance Fixed",initTime);
    initTime=clock();
    generateRouteExpectedDistanceMurrieta(locationListExpectedDistanceMurrieta, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListExpectedDistanceMurrieta, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListExpectedDistanceMurrieta,"Expected Distance Murrieta",initTime);
    initTime=clock();
    generateRouteExpectedDistanceSimple(locationListExpectedDistanceSimple, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListExpectedDistanceSimple, distancesTable) << ",";
    fileStreamData << msgLocationOrder(locationListExpectedDistanceSimple,"Expected Distance Simple",initTime);
    initTime=clock();
    generateRoutePartialOptimalHeuristic(locationListPartialOptimalHeuristic, distancesTable, startingNode);
    fileStreamExpected << evaluateRouteExpectedDistance(locationListPartialOptimalHeuristic, distancesTable) << endl;
    fileStreamData << msgLocationOrder(locationListPartialOptimalHeuristic,"Partial Optimal Heuristic",initTime);
    
    cout << endl << "Begin of test!" << endl;
    for (int repetition=0; repetition<numRepetitions; repetition++) {
      
      cout << ".";
      //cout << "Repetition " << repetition << "  ";
      //startingNode=int((numNodes-1)*double(rand()) / RAND_MAX); //Change initial node
      
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
      double distance;
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
  }
  }
  
  fileStreamExpected << endl;
  fileStreamExpected.close();
  
  fileStreamData << endl;
  fileStreamData.close();
  
  cout << endl << "End of test!" << endl;
  return 0;
}

