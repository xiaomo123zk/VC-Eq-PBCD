#include <iostream>
// #include <string>
#include <vector>
// #include <math.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdio.h>
#include <unordered_set>
#define INF 1.0e13
#define ZERO 1.0e-13
#define INVALID -1      /* Represents an invalid value. */
#define WAS_IN_QUEUE -7 /* Shows that the node was in the queue before. (7 is for luck.) */
using namespace std;
// Formatting has been reorganized.
// OpenMP-enabled version.
class CNode // Node
{
public:
    int ID;
    string Name; // STL variable-length string.
    int Origin_ID = -2;
    // int haveOrigin = -1; // Indicates whether this node has been used as an origin.
    vector<int> IncomingLink;
    vector<int> OutgoingLink;
    // CNode() : Origin_ID(-1) {} // Constructor.
};
class CLink // Link
{
public:
    int ID;
    CNode *pInNode;
    CNode *pOutNode;
    double FreeFlowTravelTime;
    double Capacity;
    double Alpha = 0.15;
    double Power = 4.0;
};
class CPath
{
public:
    int *plinkArray;
    int n = 0;
};
class COrigin
{
public:
    int ID;
    CNode *pOriginNode;
    vector<int> DestinationNode;
    vector<double> ODDemand;
    vector<int> Block; // 20240703
    vector<CPath> *PathSet;
    vector<double> *PathFlow;
};
class CODPair
{
public:
    int ori, des, desID, oriID;
};
vector<CODPair> m_CODPair;
vector<CNode> m_Node;
vector<CLink> m_Link;
vector<COrigin> m_Origin;
double *LinkFlow;
double *LinkCost;
double *LinkCostDiff;
double *LinkFreeTravelTime;
double *LinkCostCoef;
double *LinkCostDiffCoef;
// Parallel configuration
int PG_num_threads = omp_get_max_threads(); // Use all available threads.
int PG_num_threads_else = omp_get_max_threads();
// int elseparalevel = 500; // 16; 24; 32; 128; 256;
// int paralevel = 128;     // 16; 24; 32; 128; 256;
// OD and link buffers
int *ODIncludedLinkIndex;
int ODIncludedLinkNum;
int *ViolationODPairIndex1D = NULL;
int m_nViolation1D = 0;
int SumViolation;
double **ODPairGap;
double UEGap;
int largeBlock = 0;            // 202407
vector<vector<int>> groups;    // 202407
vector<vector<int>> groupFlag; // 202407
// vector<vector<int>> bad_groups;  // 202407
// Shortest-path settings
int FIRST_THRU_NODE = 1;
// UE algorithm parameters
int InnerIterNum = 1000; // 1000
int GapCheckIterNum = 100;
double step = 0.3; // 0.09
double MaxUEGap = 1.0e-12;
double MaxODGap;
// PBCD similarity settings
double threshold = 0.4;
bool bDrop;
double **FloydMatrix;
double Similarity;
// Timing statistics
double CPUTime = 0, CPUtimeE = 0, GapCheckperiod = 0;
// Function declarations
void ReadData(string fileName, string simName);
void NetworkInitialization();
void GetShortestPath(int orinode, double *ShortestPathCost, int *ShortestPathParent);
void UpdateLinkTravelCostAndDiff();
double GetUEGap();
double ObjectiveValue();
void GPUpdatePathSet(bool bFlowUpdate);
void GPMethodInitilization();
void GPFlowAdjustforOneODPair(COrigin *pOrigin, int des, bool bWithoutCheck);

void GPFlowAdjustforOneODPairWOCost(COrigin *pOrigin, int des, double *LinkflowAdj);
void GreedyFlowAdjustforOneODPairWOCost(COrigin *pOrigin, int des, double *LinkflowAdj);

void UE_PBCD();
void BCUpdateLinkCostDiff();
void BCUpdateLinkCostDiff2(double **LinkFlowAdj_thr);
// Utility declarations
void shuffleArray(int array[], int length);
double random(double start, double end); // Generate a random number within the given interval.

// Utility declarations
#pragma region "Basic"
// Data loader for the VC-PBCD experiments. The OD file contains the vertex-coloring block ID
// generated from OD path-overlap information, and the assignment algorithm consumes those blocks.
void ReadData(string fileName, string simName)
{
    string lineStr;
    string inlineStr;
    stringstream stream;
    int intData, intData1;
    double doubleData;
    int count;

    // read fileName_node.tntp
    ifstream nodeFile(fileName + "_node.tntp");
    if (nodeFile.is_open())
    {
        count = 0;
        while (getline(nodeFile, lineStr)) // Continue until the end of the file is reached.
        {
            // pNode = new CNode();
            CNode Node;
            stringstream ss(lineStr);
            getline(ss, inlineStr, '\t');
            Node.Name = inlineStr;
            Node.ID = count;
            count++;
            m_Node.push_back(Node);
        }
        nodeFile.close();
        printf("Read node data finish!\n");
    }
    else
    {
        printf("Read node data wrong!\n");
    }
    // read fileName_link.csv
    ifstream linkFile(fileName + "_link.csv");
    count = 0;
    if (linkFile.is_open())
    {
        while (getline(linkFile, lineStr)) // Continue until the end of the file is reached.
        {
            CLink Link;
            Link.ID = count;
            stringstream ss(lineStr);

            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData;
            Link.pInNode = &m_Node[intData - 1];
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData;
            Link.pOutNode = &m_Node[intData - 1];
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> doubleData;
            Link.Capacity = doubleData;
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> doubleData;
            // Link.FreeFlowTravelTime = doubleData;
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> doubleData;
            Link.FreeFlowTravelTime = doubleData;
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> doubleData;
            Link.Alpha = doubleData;
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> doubleData;
            Link.Power = doubleData;
            stream.clear();
            m_Link.push_back(Link);
            count++;
        }
        linkFile.close();
        printf("Read link data finish!\n");
    }
    else
    {
        printf("Read link data wrong!\n");
    }
    // read index-OD pairs
    int desCounter = 0;
    int oriCounter = -1;
    // int ori = 1;          Birm_od_VertexColor_equ_sim09
    ifstream odFile1(fileName + "_od_VertexColor_equ_sim" + simName + ".csv");
    CNode *pNode;
    if (odFile1.is_open())
    {
        while (getline(odFile1, lineStr)) // Continue until the end of the file is reached.
        {
            CODPair pODPair;
            stringstream ss(lineStr);

            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData;
            pNode = &m_Node[intData - 1];
            stream.clear();
            // if (intData == ori)
            if (pNode->Origin_ID == -2)
            {
                pNode->Origin_ID = -1;
                desCounter = 0;
                oriCounter++;
            }
            else
            {
                // ori = intData;
                // desCounter = 0;
                // oriCounter++;
            }
            pODPair.ori = intData;
            pODPair.oriID = oriCounter; // oriID records the index in the origin set, which is not always ID - 1.
            // stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData;
            pODPair.des = intData;
            pODPair.desID = desCounter;
            desCounter++;
            stream.clear();

            m_CODPair.push_back(pODPair);
        }
        odFile1.close();
        printf("Read index-od pair data finish!\n");
    }
    else
    {
        printf("Read od data wrong!\n");
    }
    // read fileName_od.csv
    ifstream odFile(fileName + "_od_VertexColor_equ_sim" + simName + ".csv");
    // CNode *pNode;
    COrigin *pOrigin;
    // int line = -1;
    if (odFile.is_open())
    {
        while (getline(odFile, lineStr)) // Continue until the end of the file is reached.
        {
            // line++;
            stringstream ss(lineStr);
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData;
            pNode = &m_Node[intData - 1];
            stream.clear();

            if (pNode->Origin_ID < 0)
            {
                // pOrigin = new COrigin();
                COrigin Origin;
                Origin.ID = m_Origin.size();
                Origin.pOriginNode = pNode;
                pNode->Origin_ID = Origin.ID;
                m_Origin.push_back(Origin);
                pOrigin = &m_Origin[pNode->Origin_ID];
            }
            else
            {
                pOrigin = &m_Origin[pNode->Origin_ID];
            }
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData;
            pNode = &m_Node[intData - 1];
            stream.clear();
            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> doubleData;
            stream.clear();

            getline(ss, inlineStr, ',');
            stream << inlineStr;
            stream >> intData1;
            stream.clear();
            // getline(ss, inlineStr, ',');
            // stream << inlineStr;
            // stream >> doubleData;
            // stream.clear();
            // cout << "OD demand:" << doubleData << endl;
            if (doubleData > 0)
            {
                pOrigin->DestinationNode.push_back(pNode->ID);
                pOrigin->ODDemand.push_back(doubleData);
                pOrigin->Block.push_back(intData1);
                // pOrigin->ODPairIndex.push_back(line);
                if (largeBlock < intData1)
                    largeBlock = intData1;
            }
            stream.clear();
        }
        odFile.close();
        printf("Read origin data finish!\n");
    }
    else
    {
        printf("Read origin data wrong!\n");
    }
}
// Initialize network adjacency, BPR coefficients, and reusable link/OD buffers.
void NetworkInitialization()
{
    CLink *pLink;
    double temp;
    //********************************************
    LinkCostDiff = new double[m_Link.size()];
    LinkFlow = new double[m_Link.size()];
    LinkCost = new double[m_Link.size()];
    LinkFreeTravelTime = new double[m_Link.size()];
    LinkCostCoef = new double[m_Link.size()];
    LinkCostDiffCoef = new double[m_Link.size()];
    //********************************************
    ODIncludedLinkIndex = new int[m_Link.size()];
    for (int link = 0; link < m_Link.size(); link++)
    {
        pLink = &m_Link.at(link);
        pLink->pInNode->OutgoingLink.push_back(link); // Add this link index to the corresponding node.
        pLink->pOutNode->IncomingLink.push_back(link);
        LinkFreeTravelTime[link] = pLink->FreeFlowTravelTime; // Extract the free-flow travel time from the link record.
        temp = pow(pLink->Capacity, pLink->Power);
        LinkCostCoef[link] = LinkFreeTravelTime[link] * pLink->Alpha / temp;
        LinkCostDiffCoef[link] = LinkFreeTravelTime[link] * pLink->Alpha * pLink->Power / temp;
    }
    //***************************************************
}

// BPR link-cost update used by the UE objective and path marginal-cost calculations.
void UpdateLinkTravelCostAndDiff()
{
    double temp, temp1;
    for (int link = 0; link < m_Link.size(); link++)
    {
        temp1 = LinkFlow[link];
        temp = temp1 * temp1 * temp1;
        LinkCost[link] = LinkFreeTravelTime[link] + LinkCostCoef[link] * temp * temp1;
        LinkCostDiff[link] = LinkCostDiffCoef[link] * temp;
    }
}
// Merge thread-local link-flow adjustments after each block pass and update only changed links.
void BCUpdateLinkCostDiff2(double **LinkFlowAdj_thr)
{
    int link_size = m_Link.size();
#pragma omp parallel for num_threads(PG_num_threads)
    for (int i = 0; i < link_size; i++)
    {
        bool adj = false;
        for (int j = 0; j < PG_num_threads; j++)
        {
            if (LinkFlowAdj_thr[j][i] != 0)
            {
                LinkFlow[i] += LinkFlowAdj_thr[j][i];
                adj = true;
                LinkFlowAdj_thr[j][i] = 0; // Reset here so the buffer is ready for the next use.
            }
        }
        if (adj)
        {
            // if (LinkFlow[i] < zero)
            //     LinkFlow[i] = 0.0;
            double temp1 = LinkFlow[i];
            double temp = temp1 * temp1 * temp1;
            if (LinkFreeTravelTime[i] == 0)
            {
                LinkCost[i] = 0;
                LinkCostDiff[i] = 0;
            }
            else
            {
                LinkCost[i] = LinkFreeTravelTime[i] + LinkCostCoef[i] * temp * temp1;
                LinkCostDiff[i] = LinkCostDiffCoef[i] * temp;
            }
        }
    }
}

// Relative gap metric: compares current link-based system cost with shortest-path OD cost.
double GetUEGap() // Parallelized relative-gap calculation.
{
    double link_systemcost = 0.0, od_systemcost = 0.0;
    int ori, link, des;
    for (link = 0; link < m_Link.size(); link++)
    {
        link_systemcost += LinkFlow[link] * LinkCost[link];
    }
#pragma omp parallel num_threads(PG_num_threads)
    {
#pragma omp for reduction(+ : od_systemcost)
        for (int ori = 0; ori < m_Origin.size(); ori++)
        {
            double *ShortestPathCost = new double[m_Node.size()];
            int *ShortestPathParent = new int[m_Node.size()];
            COrigin *pOrigin = &m_Origin.at(ori);
            GetShortestPath(pOrigin->pOriginNode->ID, ShortestPathCost, ShortestPathParent);
            for (int des = 0; des < pOrigin->DestinationNode.size(); des++)
                od_systemcost += pOrigin->ODDemand.at(des) * ShortestPathCost[pOrigin->DestinationNode.at(des)];
            delete[] ShortestPathParent;
            delete[] ShortestPathCost;
        }
    }
    // printf("  od_systemcost:%f, link_systemcost:%f \n", od_systemcost, link_systemcost);
    return 1 - od_systemcost / link_systemcost;
}
double ObjectiveValue()
{
    double obj = 0.0;
    double link_integral = 0.0;
    CLink *pLink;
    for (int link = 0; link < m_Link.size(); link++)
    {
        pLink = &m_Link.at(link);
        link_integral = LinkFreeTravelTime[link] * LinkFlow[link] + LinkCostCoef[link] * pow(LinkFlow[link], 5) / 5;
        obj += link_integral;
    }
    return obj;
}
void GetShortestPath(int orinode, double *ShortestPathCost, int *ShortestPathParent)
{
    CNode *pNode;
    CLink *pLink;
    int m_nNode = m_Node.size();
    int now, NewNode;
    double NewCost;
    int QueueFirst, QueueLast;
    int *QueueNext =
        new int[m_nNode](); // Queue of nodes to check; each value is the next node ID or a marker that the node has already entered the queue.
    for (int node = 0; node < m_nNode; node++)
    {
        QueueNext[node] = INVALID;
        ShortestPathCost[node] = INF;
        ShortestPathParent[node] = INVALID;
    }
    now = orinode;
    QueueNext[now] = WAS_IN_QUEUE;
    ShortestPathParent[now] = INVALID;
    ShortestPathCost[now] = 0.0;
    QueueFirst = QueueLast = INVALID;
    while ((now != INVALID) && (now != WAS_IN_QUEUE))
    {
        pNode = &m_Node[now]; // Corresponding node in the network.
        if (now >= FIRST_THRU_NODE - 1 || now == orinode)
        {
            for (int index = 0; index < pNode->OutgoingLink.size(); index++)
            {
                /* For every link that terminate at "now": */
                pLink = &m_Link[pNode->OutgoingLink[index]]; // Find the link between the corresponding node and index.
                NewNode = pLink->pOutNode->ID;
                NewCost = ShortestPathCost[now] + LinkCost[pLink->ID];
                if (ShortestPathCost[NewNode] > NewCost)
                {
                    /* If the new lable is better than the old one, correct it, and make sure that the new node to the
                     * queue. */
                    ShortestPathCost[NewNode] = NewCost;
                    ShortestPathParent[NewNode] = pLink->ID;
                    /* If the new node was in the queue before, add it as the first in the queue. */
                    if (QueueNext[NewNode] == WAS_IN_QUEUE)
                    {
                        QueueNext[NewNode] = QueueFirst;
                        QueueFirst = NewNode;
                        if (QueueLast == INVALID)
                            QueueLast = NewNode;
                    }
                    /* If the new node is not in the queue, and wasn't there before, add it at the end of the queue. */
                    else if (QueueNext[NewNode] == INVALID && NewNode != QueueLast)
                    {
                        if (QueueLast != INVALID)
                        {
                            /*Usually*/
                            QueueNext[QueueLast] = NewNode;
                            QueueLast = NewNode;
                        }
                        else
                        {
                            /* If the queue is empty, initialize it. */
                            QueueFirst = QueueLast = NewNode;
                            QueueNext[QueueLast] = INVALID;
                        }
                    }
                    /* If the new node is in the queue, just leave it there. (Do nothing) */
                }
            }
        }
        /* Get the first node out of the queue, and use it as the current node. */
        now = QueueFirst;
        if ((now == INVALID) || (now == WAS_IN_QUEUE))
            break;
        QueueFirst = QueueNext[now];
        QueueNext[now] = WAS_IN_QUEUE;
        if (QueueLast == now)
            QueueLast = INVALID;
    }
    delete[] QueueNext;
}

#pragma endregion
#pragma region "Specific-GP"
// Builds the initial feasible solution, performs all-or-nothing assignment, and maps OD pairs
// into vertex-coloring blocks for the PBCD iterations.
void GPMethodInitilization() // Algorithm initialization.
{
    double start = 0, finish = 0;
    COrigin *pOrigin;
    start = omp_get_wtime(); // Start timing.
    for (int link = 0; link < m_Link.size(); link++)
        LinkCost[link] = LinkFreeTravelTime[link];
    ODPairGap = new double *[m_Origin.size()];
    ViolationODPairIndex1D = new int[m_CODPair.size()];
    for (int ori = 0; ori < m_Origin.size(); ori++)
    {
        pOrigin = &m_Origin.at(ori);
        ODPairGap[ori] = new double[pOrigin->DestinationNode.size()];
        pOrigin->PathSet = new vector<CPath>[pOrigin->DestinationNode.size()];
        pOrigin->PathFlow = new vector<double>[pOrigin->DestinationNode.size()];
    }
    // All-or-nothing assignment.
    GPUpdatePathSet(false);
    CPath Path;
    for (int link = 0; link < m_Link.size(); link++)
        LinkFlow[link] = 0.0;
    groups.resize(largeBlock); // Set the first dimension of groups to the block size.
    groupFlag.resize(largeBlock);
    // bad_groups.resize(largeBlock); //20241013

    int ODID = 0;
    for (int ori = 0; ori < m_Origin.size(); ori++)
    {
        pOrigin = &m_Origin.at(ori);
        for (int des = 0; des < pOrigin->DestinationNode.size(); des++)
        {
            pOrigin->PathFlow[des].at(0) = pOrigin->ODDemand.at(des);
            for (int pathindex = 0; pathindex < pOrigin->PathSet[des].size(); pathindex++)
            {
                Path = pOrigin->PathSet[des].at(pathindex);
                for (int pathlink = 0; pathlink < Path.n; pathlink++)
                    LinkFlow[Path.plinkArray[pathlink]] += pOrigin->PathFlow[des].at(pathindex);
            }
            for (int color = 1; color < largeBlock + 1; color++)
            {
                if (color == pOrigin->Block[des])
                {
                    groups[color - 1].push_back(ODID);
                    groupFlag[color - 1].push_back(0);
                    // bad_groups[color - 1].push_back(ODID);
                    break;
                }
            }
            ODID++;
        }
    }

    UpdateLinkTravelCostAndDiff();
    finish = omp_get_wtime(); // End timing.
    UEGap = GetUEGap();
    double obj = ObjectiveValue();
    cout << setprecision(20) << obj << ",0," << setprecision(10) << UEGap << ", " << (finish - start) << endl;
}
// Serial GP restricted OD-subproblem update. This version immediately applies link-flow changes
// and refreshes affected link costs for sequential Gauss-Seidel-style updates.
void GPFlowAdjustforOneODPair(COrigin *pOrigin, int des, bool bWithoutCheck)
{
    int ori = pOrigin->ID;
    if (pOrigin->PathSet[des].size() > 1)
    {
        // ODcomputCountElse += 1;
        ODIncludedLinkNum = 0;
        int spathindex = 0;
        int link = 0;
        CPath Path;
        //*********************************************
        // Update path information
        bool *bLink = new bool[m_Link.size()](); // Whether this link has been included.
        int m_nPath = pOrigin->PathSet[des].size();
        double *ODPathCost = new double[m_nPath]();
        double *ODPathCostDiff = new double[m_nPath]();
        double *ODPathFlow = new double[m_nPath]();
        int *bRemoved = new int[m_nPath]();
        double ODTotalCost = 0.0, minpathcost = INF;

        // update path cost
        for (int path = 0; path < m_nPath; path++)
        {
            Path = pOrigin->PathSet[des].at(path);
            for (int pathlink = 0; pathlink < Path.n; pathlink++)
            {
                link = Path.plinkArray[pathlink];
                ODPathCost[path] += LinkCost[link];
                ODPathCostDiff[path] += LinkCostDiff[link];
                if (!bLink[link])
                {
                    ODIncludedLinkIndex[ODIncludedLinkNum] = link; // Record links used by the OD pair for later link-flow and link-cost updates.
                    bLink[link] = true;
                    ODIncludedLinkNum++;
                }
            }
            if (minpathcost > ODPathCost[path])
            {
                minpathcost = ODPathCost[path]; // Find the minimum-cost path ID.
                spathindex = path;
            }
            ODTotalCost += ODPathCost[path] * pOrigin->PathFlow[des].at(path);
            if (ODPathCostDiff[path] < ZERO)
                ODPathCostDiff[path] = ZERO;
        }
        //**************************************************
        // // Find the shortest path among all paths.

        ODPairGap[ori][des] = 1 - minpathcost * pOrigin->ODDemand.at(des) / ODTotalCost;
        double spTempFlow = pOrigin->ODDemand.at(des);
        double tempPF = 0;
        if (ODPairGap[ori][des] > MaxODGap || bWithoutCheck) // Always execute in the outer loop; check conditionally in the inner loop.
        {
            // update path flow
            for (int path = 0; path < m_nPath; path++)
            {
                if (path != spathindex)
                {
                    tempPF = pOrigin->PathFlow[des].at(path) - step * (ODPathCost[path] - ODPathCost[spathindex]) / ODPathCostDiff[path];
                    if (tempPF <= 0.0) // ZERO
                    {
                        bRemoved[path] = 1;
                        ODPathFlow[path] = 0.0;
                    }
                    else
                    {
                        spTempFlow -= tempPF;
                        ODPathFlow[path] = tempPF;
                    }
                }
            }
            ODPathFlow[spathindex] = spTempFlow;
            //***************************************************
            // update link flow, cost, and cost diffential
            for (int path = m_nPath - 1; path >= 0; path--)
            {
                Path = pOrigin->PathSet[des].at(path);
                for (int pathlink = 0; pathlink < Path.n; pathlink++)
                {
                    LinkFlow[Path.plinkArray[pathlink]] += (ODPathFlow[path] - pOrigin->PathFlow[des].at(path)); // GaussSeidel
                    if (LinkFlow[Path.plinkArray[pathlink]] < 1e-12)                                             // abs(LinkFlow[Path.plinkArray[pathlink]])<1e-12
                    {
                        LinkFlow[Path.plinkArray[pathlink]] = 0.0;
                    }
                    if (LinkFlow[Path.plinkArray[pathlink]] < 0.0)
                    {
                        cout << " XX " << endl;
                    }
                }
                if (bRemoved[path])
                {
                    delete[] pOrigin->PathSet[des].at(path).plinkArray;
                    pOrigin->PathFlow[des].erase(pOrigin->PathFlow[des].begin() + path);
                    pOrigin->PathSet[des].erase(pOrigin->PathSet[des].begin() + path);
                }
                else
                    pOrigin->PathFlow[des][path] = ODPathFlow[path];
            }
            double temp = 0, temp1 = 0;
            for (int linkindex = 0; linkindex < ODIncludedLinkNum; linkindex++)
            {
                link = ODIncludedLinkIndex[linkindex];
                temp1 = LinkFlow[link];
                temp = temp1 * temp1 * temp1;
                LinkCost[link] = LinkFreeTravelTime[link] + LinkCostCoef[link] * temp * temp1;
                LinkCostDiff[link] = LinkCostDiffCoef[link] * temp;
            }
            SumViolation++;
        }
        delete[] bLink;
        delete[] ODPathCost;
        delete[] ODPathCostDiff;
        delete[] ODPathFlow;
        delete[] bRemoved;
    }
    else
        ODPairGap[ori][des] = 0.0;
}

// Greedy restricted OD-subproblem update for larger path sets. Link-flow deltas are accumulated
// in the thread-local buffer and merged after the block update.
void GreedyFlowAdjustforOneODPairWOCost(COrigin *pOrigin, int des, double *LinkflowAdj)
{
  int ori = pOrigin->ID;
  // if (pOrigin->PathSet[des].size() > 1)
  // {
  int temppathindex, pathindex, tempindex;
  double temp_B = 0.0, temp_C = 0.0, temp_wrs = INF;
  double tempvalue;
  double odflow = pOrigin->ODDemand[des];
  double tempflow;
  int path, pathlink;
  CPath Path;
  int link, linkindex;
  // bool *bLink = new bool[m_Link.size()]();
  // ODIncludedLinkNum = 0;
  //*********************************************
  // Update path information
  int m_nPath = pOrigin->PathSet[des].size();
  double *ODPathCost = new double[m_nPath]();
  double *ODPathCostDiff = new double[m_nPath]();
  double *ODPathCostAuxiliary = new double[m_nPath]();
  double *ODPathFlow = new double[m_nPath]();
  int *bRemoved = new int[m_nPath]();
  int *ODPathOrder = new int[m_nPath]();
  double ODTotalCost = 0.0, minpathcost = INF;
  // update path cost
  for (path = 0; path < m_nPath; path++)
  {
    Path = pOrigin->PathSet[des][path];
    ODPathOrder[path] = path;
    for (pathlink = 0; pathlink < Path.n; pathlink++)
    {
      link = Path.plinkArray[pathlink];
      ODPathCost[path] += LinkCost[link];
      ODPathCostDiff[path] += LinkCostDiff[link];
      // if (!bLink[link])
      // {
      //   ODIncludedLinkIndex[ODIncludedLinkNum] = link; // Record links used by the OD pair for later link-flow and link-cost updates.
      //   bLink[link] = true;
      //   ODIncludedLinkNum++;
      // }
    }
    if (minpathcost > ODPathCost[path])
      minpathcost = ODPathCost[path]; // Find the minimum cost.
    ODTotalCost += ODPathCost[path] * pOrigin->PathFlow[des][path];
    if (ODPathCostDiff[path] < ZERO)
      ODPathCostDiff[path] = ZERO;
    ODPathCostAuxiliary[path] = ODPathCost[path] - ODPathCostDiff[path] * pOrigin->PathFlow[des][path];
  }
  //**************************************************
  ODPairGap[ori][des] = 1 - minpathcost * pOrigin->ODDemand[des] / ODTotalCost;
  if (ODPairGap[ori][des] > MaxODGap) // Always execute in the outer loop; check conditionally in the inner loop.
  {
    // update path flow
    temp_B = 0.0;
    temp_C = 0.0;
    temp_wrs = INF;
    // Sort in ascending order.
    for (pathindex = 0; pathindex < m_nPath; pathindex++)
    {
      for (temppathindex = pathindex + 1; temppathindex < m_nPath; temppathindex++)
      {
        if (ODPathCostAuxiliary[ODPathOrder[temppathindex]] < ODPathCostAuxiliary[ODPathOrder[pathindex]])
        {
          tempindex = ODPathOrder[temppathindex]; // Find the smaller path index.
          ODPathOrder[temppathindex] = ODPathOrder[pathindex];
          ODPathOrder[pathindex] = tempindex; // Swap positions.
        }
      }
      // Find the minimum ODPathCostAuxiliary.
      path = ODPathOrder[pathindex];
      if (ODPathCostAuxiliary[path] >= temp_wrs)
        break; // Later paths must have no flow.
      tempvalue = 1.0 / (ODPathCostDiff[path] * odflow);
      temp_C += ODPathCostAuxiliary[path] * tempvalue;
      temp_B += tempvalue;
      temp_wrs = (1.0 + temp_C) / temp_B;
    }

    tempflow = odflow;
    for (temppathindex = 0; temppathindex < pathindex - 1; temppathindex++)
    {
      path = ODPathOrder[temppathindex];
      ODPathFlow[path] = (temp_wrs - ODPathCostAuxiliary[path]) / ODPathCostDiff[path];
      tempflow -= ODPathFlow[path];
    }
    ODPathFlow[ODPathOrder[temppathindex]] = tempflow;
    for (temppathindex = pathindex; temppathindex < m_nPath; temppathindex++)
    {
      ODPathFlow[ODPathOrder[temppathindex]] = 0.0;
      bRemoved[ODPathOrder[temppathindex]] = true;
    }
    //***************************************************
    // update link flow, cost, and cost diffential
    for (path = m_nPath - 1; path >= 0; path--)
    {
      Path = pOrigin->PathSet[des][path];
      for (pathlink = 0; pathlink < Path.n; pathlink++)
      {
        LinkflowAdj[Path.plinkArray[pathlink]] += (ODPathFlow[path] - pOrigin->PathFlow[des].at(path)); // GaussSeidel
        // LinkFlow[Path.plinkArray[pathlink]] += (ODPathFlow[path] - pOrigin->PathFlow[des][path]); // GaussSeidel
      }
      if (bRemoved[path])
      {
        delete[] pOrigin->PathSet[des].at(path).plinkArray;
        pOrigin->PathFlow[des].erase(pOrigin->PathFlow[des].begin() + path);
        pOrigin->PathSet[des].erase(pOrigin->PathSet[des].begin() + path);
      }
      else
        pOrigin->PathFlow[des][path] = ODPathFlow[path];
    }
  }
  // delete[] bLink;
  delete[] ODPathCost;
  delete[] ODPathCostDiff;
  delete[] ODPathCostAuxiliary;
  delete[] ODPathFlow;
  delete[] bRemoved;
  delete[] ODPathOrder;
  // }
  // else
  //   ODPairGap[ori][des] = 0.0;
}

// Adjust OD path flows. pOrigin is the origin object address; des is the destination index in PathSet, ODDemand, PathFlow, and related origin-owned containers; LinkflowAdj temporarily stores network adjustment information. Link data is not updated. No return value.
// Parallel GP update for smaller OD path sets. The caller chooses this branch by path-set size.
void GPFlowAdjustforOneODPairWOCost(COrigin *pOrigin, int des, double *LinkflowAdj)
{
  int ori = pOrigin->ID;
  // if (pOrigin->PathSet[des].size() > 1)
  // {
  int spathindex = 0;
  int link = 0;
  CPath Path;
  //*********************************************
  // Update path information
  int m_nPath = pOrigin->PathSet[des].size();
  double *ODPathCost = new double[m_nPath]();
  double *ODPathCostDiff = new double[m_nPath]();
  double *ODPathFlow = new double[m_nPath]();
  int *bRemoved = new int[m_nPath]();
  double ODTotalCost = 0.0, minpathcost = INF;

  // update path cost
  for (int path = 0; path < m_nPath; path++)
  {
    Path = pOrigin->PathSet[des].at(path);
    for (int pathlink = 0; pathlink < Path.n; pathlink++)
    {
      link = Path.plinkArray[pathlink];
      ODPathCost[path] += LinkCost[link];
      ODPathCostDiff[path] += LinkCostDiff[link];
    }
    if (minpathcost > ODPathCost[path])
    {
      minpathcost = ODPathCost[path]; // Find the minimum-cost path ID.
      spathindex = path;
    }
    ODTotalCost += ODPathCost[path] * pOrigin->PathFlow[des].at(path);
    if (ODPathCostDiff[path] < ZERO)
      ODPathCostDiff[path] = ZERO;
  }
  //**************************************************
  // // Find the shortest path among all paths.

  ODPairGap[ori][des] = 1 - minpathcost * pOrigin->ODDemand.at(des) / ODTotalCost;
  double spTempFlow = pOrigin->ODDemand.at(des);
  if (ODPairGap[ori][des] > MaxODGap) // Always execute in the outer loop; check conditionally in the inner loop.
  {
    // update path flow
    for (int path = 0; path < m_nPath; path++)
    {
      if (path != spathindex)
      {
        double tempPF = pOrigin->PathFlow[des].at(path) - step * (ODPathCost[path] - ODPathCost[spathindex]) / ODPathCostDiff[path];
        if (tempPF <= 0.0) // ZERO
        {
          bRemoved[path] = 1;
          ODPathFlow[path] = 0.0;
        }
        else
        {
          spTempFlow -= tempPF;
          ODPathFlow[path] = tempPF;
        }
      }
    }
    ODPathFlow[spathindex] = spTempFlow;
    // update link flow, cost, and cost diffential
    for (int path = m_nPath - 1; path >= 0; path--)
    {
      Path = pOrigin->PathSet[des].at(path);
      for (int pathlink = 0; pathlink < Path.n; pathlink++)
      {
        LinkflowAdj[Path.plinkArray[pathlink]] += (ODPathFlow[path] - pOrigin->PathFlow[des].at(path)); // GaussSeidel
        // if (LinkflowAdj[Path.plinkArray[pathlink]] < 1e-14)
        //   LinkflowAdj[Path.plinkArray[pathlink]] = 0.0;
      }
      if (bRemoved[path])
      {
        delete[] pOrigin->PathSet[des].at(path).plinkArray;
        pOrigin->PathFlow[des].erase(pOrigin->PathFlow[des].begin() + path);
        pOrigin->PathSet[des].erase(pOrigin->PathSet[des].begin() + path);
      }
      else
        pOrigin->PathFlow[des][path] = ODPathFlow[path];
    }
    // No need to update cost or cost differential.
  }
  delete[] ODPathCost;
  delete[] ODPathCostDiff;
  delete[] ODPathFlow;
  delete[] bRemoved;
  // }
  // else
  //   ODPairGap[ori][des] = 0.0;
}

void BCUpdateLinkCostDiff()
{
    int link = 0;
    double temp = 0, temp1 = 0;
    for (int linkindex = 0; linkindex < ODIncludedLinkNum; linkindex++)
    {
        link = ODIncludedLinkIndex[linkindex];
        temp1 = LinkFlow[link];
        temp = temp1 * temp1 * temp1;
        LinkCost[link] = LinkFreeTravelTime[link] + LinkCostCoef[link] * temp * temp1;
        LinkCostDiff[link] = LinkCostDiffCoef[link] * temp;
    }
}

// Column-generation step: recompute shortest paths under current link costs and add any missing
// shortest paths to each OD path set.
void GPUpdatePathSet(bool bFlowUpdate) // Update the path set.
{
#pragma omp parallel num_threads(PG_num_threads)
    {
#pragma omp for
        for (int ori = 0; ori < m_Origin.size(); ori++)
        {
            int *pReversePath = new int[m_Node.size()]; // Path buffer that needs to be reversed.
            int pathlinknum, templink, tempnode;
            CLink *pTempLink;
            bool bPathExist, bPathEqual;
            CPath *pPath; // Path represented by links.
            CPath TempPath;
            COrigin *pOrigin = &m_Origin.at(ori);
            double *ShortestPathCost = new double[m_Node.size()];
            int *ShortestPathParent = new int[m_Node.size()];
            GetShortestPath(pOrigin->pOriginNode->ID, ShortestPathCost, ShortestPathParent);
            for (int des = 0; des < pOrigin->DestinationNode.size(); des++)
            {
                //****************************
                pathlinknum = 0;
                tempnode = pOrigin->DestinationNode[des];
                templink = ShortestPathParent[tempnode]; // templink is the link ID.
                while (templink > -1)
                {
                    pReversePath[pathlinknum] = templink;
                    pathlinknum++;
                    pTempLink = &m_Link[templink];
                    tempnode = pTempLink->pInNode->ID;
                    templink = ShortestPathParent[tempnode];
                }
                pPath = new CPath();
                pPath->n = pathlinknum;
                pPath->plinkArray = new int[pathlinknum];
                for (int pathlink = 0; pathlink < pathlinknum; pathlink++)
                    pPath->plinkArray[pathlink] = pReversePath[pathlinknum - pathlink - 1]; // Build the new shortest path for the current OD pair.
                //*****************************************
                // Check the existence of the route
                bPathExist = false;
                for (int path = 0; path < pOrigin->PathSet[des].size(); path++)
                {
                    TempPath = pOrigin->PathSet[des][path];
                    bPathEqual = true;
                    if (TempPath.n != pathlinknum) // Check path length.
                        bPathEqual = false;
                    else
                    {
                        for (int pathlink = 0; pathlink < pathlinknum; pathlink++)
                        {
                            if (pPath->plinkArray[pathlink] != TempPath.plinkArray[pathlink])
                            {
                                bPathEqual = false;
                                break;
                            }
                        }
                    }
                    if (bPathEqual)
                    {
                        bPathExist = true;
                        break;
                    }
                }
                if (!bPathExist)
                {
                    pOrigin->PathSet[des].push_back(*pPath); // Add the new path.
                    pOrigin->PathFlow[des].push_back(
                        0.0); // Add the flow for the new path.
                              //  if (bFlowUpdate)
                              //  {
                              //      GPFlowAdjustforOneODPair(pOrigin, des, true); // Function variant 3: true.
                              //  }
                }
            }
            delete[] ShortestPathParent;
            delete[] ShortestPathCost;
            delete[] pReversePath;
        }
    }
}

double random(double start, double end) // Generate a random number within the given interval.
{
    return start + (end - start) * rand() / (RAND_MAX + 1.0);
}

// Main VC-PBCD algorithm: alternate shortest-path generation with block coordinate flow updates.
// Vertex-coloring blocks reduce path overlap among OD pairs solved in parallel.
void UE_PBCD() // Solve UE traffic assignment using the GP path-based method.
{
    COrigin *pOrigin;
    int nOutLoop = 0;
    double obj = 0;
    int des;
    MaxODGap = 1.0e-3;
    double GapRate = 4.0;
    double start, finish, CPUTime, PG_period = 0, PFA_period = 0, RG_period = 0; // Define timing variables.
    NetworkInitialization();
    start = omp_get_wtime(); // Start timing.
    GPMethodInitilization();
    // Floyd();

    double If_Time = 0, If_PFA_Time = 0, If_SumAdj_Time = 0;
    double Else_Time = 0, Else_PFA_Time = 0, Else_SumAdj_Time = 0, Else_Rem_Time = 0;

    //////////// Thread initialization //////////
    double **LinkFlowAdj_thr = new double *[PG_num_threads]();
#pragma omp parallel for num_threads(PG_num_threads)
    for (int j = 0; j < PG_num_threads; j++)
    {
        LinkFlowAdj_thr[j] = new double[m_Link.size()];
        // ViolationODIndex_thr[j] = new int[m_CODPair.size() + 1000];
        double link_size = m_Link.size();
        for (int m = 0; m < link_size; m++)
            LinkFlowAdj_thr[j][m] = 0;
    }
    while (UEGap > MaxUEGap && nOutLoop <= 30 && CPUTime < 1000) // 500
    {
        nOutLoop++;
        double PG_beg = omp_get_wtime();
        GPUpdatePathSet(true);
        double PG_end = omp_get_wtime();
        PG_period += (PG_end - PG_beg);
        double PFA_beg = omp_get_wtime();
        // printf("nOutLoop:%d \n", nOutLoop);

        // Copy directly to the temporary bad_groups container.
        vector<vector<int>> bad_groups;
        // vector<vector<int>> bad_groupFlag;

        for (int iter = 0; iter < InnerIterNum; iter++)
        {
            // Full block sweep: evaluate each vertex-coloring block and keep only OD pairs above the gap threshold.
            if (iter % GapCheckIterNum == 0)
            {
                bad_groups = groups; // Store the current full block grouping.
                // bad_groupFlag = groupFlag;
                m_nViolation1D = 0;
                for (int blockID = 1; blockID < largeBlock + 1; blockID++)
                { // Read OD pairs block by block.
                    vector<bool> toDelete(groups[blockID - 1].size(), false);
                    ODIncludedLinkNum = 0; // bLink should also be defined outside and passed in by pointer.

#pragma omp parallel for num_threads(PG_num_threads)
                    for (int i = 0; i < groups[blockID - 1].size(); ++i) // Keep a backup group containing only violating OD pairs.
                    {
                        int current_ODID = groups[blockID - 1][i];
                        int ori = m_CODPair[current_ODID].oriID;
                        int des = m_CODPair[current_ODID].desID;
                        int destrueid = m_CODPair[current_ODID].des - 1;

                        int tid = omp_get_thread_num();
                        COrigin *pOrigin = &m_Origin.at(ori);
                        if (pOrigin->PathSet[des].size() == 1)
                        {
                            ODPairGap[ori][des] = 0.0;
                        }
                        if (pOrigin->PathSet[des].size() == 2 || pOrigin->PathSet[des].size() == 3)
                        {
                            GPFlowAdjustforOneODPairWOCost(pOrigin, des, LinkFlowAdj_thr[tid]); //   LinkFlowAdj_thr[tid][link]
                        }
                        else if (pOrigin->PathSet[des].size() > 3)
                        {
                            GreedyFlowAdjustforOneODPairWOCost(pOrigin, des, LinkFlowAdj_thr[tid]);
                        }
                        if (ODPairGap[ori][des] < MaxODGap) // Remove OD pairs that already satisfy the gap threshold.
                        {
                            toDelete[i] = true; // Mark for deletion in the parallel loop.
                        }
                    }

                    // Serially remove marked elements.
                    auto &bad_group = bad_groups[blockID - 1];
                    auto it = bad_group.begin();
                    for (int i = 0; i < toDelete.size(); ++i)
                    {
                        if (toDelete[i])
                        {
                            it = bad_group.erase(it);
                        }
                        else
                        {
                            ++it;
                        }
                    }
                    BCUpdateLinkCostDiff2(LinkFlowAdj_thr);
                }
                // printf("3 \n");
                // if (m_nViolation1D == 0)
                //     break;
            }
            else
            {
                // Restricted sweep: revisit only the OD pairs retained in bad_groups.
                // printf("4 \n");
                for (int blockID = 1; blockID < largeBlock + 1; blockID++)
                {
                    vector<bool> ElsetoDelete(bad_groups[blockID - 1].size(), false);

                    if (bad_groups[blockID - 1].size() == 0)
                    {
                    }
                    else
                    {

#pragma omp parallel for num_threads(PG_num_threads_else)
                        for (int i = 0; i < bad_groups[blockID - 1].size(); ++i)
                        {
                            // if (bad_groupFlag[blockID - 1][i] == 1)
                            // {
                            int current_ODID = bad_groups[blockID - 1][i];
                            // int next_ODID_value = groups[blockID - 1][i + 1];
                            int ori = m_CODPair[current_ODID].oriID;
                            int des = m_CODPair[current_ODID].desID;
                            int destrueid = m_CODPair[current_ODID].des - 1;

                            int tid = omp_get_thread_num();
                            COrigin *pOrigin = &m_Origin.at(ori);
                            if (pOrigin->PathSet[des].size() == 1)
                            {
                                ODPairGap[ori][des] = 0.0;
                            }
                            if (pOrigin->PathSet[des].size() == 2)
                            {
                                GPFlowAdjustforOneODPairWOCost(pOrigin, des, LinkFlowAdj_thr[tid]); //   LinkFlowAdj_thr[tid][link]
                            }
                            else if (pOrigin->PathSet[des].size() > 2)
                            {
                                GreedyFlowAdjustforOneODPairWOCost(pOrigin, des, LinkFlowAdj_thr[tid]);
                            }
                            if (ODPairGap[ori][des] < MaxODGap) // Remove OD pairs that already satisfy the gap threshold.
                            {
                                // bad_groups[blockID - 1].erase(bad_groups[blockID - 1].begin() + i);
                                ElsetoDelete[i] = true;
                            }
                            // printf("4.2, bad_groups_bs[%d]:%d  \n", blockID, bad_groups[blockID - 1].size());
                        }
                        // Serially remove marked elements.
                        auto &bad_group = bad_groups[blockID - 1];
                        auto it = bad_group.begin();
                        for (int i = 0; i < ElsetoDelete.size(); ++i)
                        {
                            if (ElsetoDelete[i])
                            {
                                it = bad_group.erase(it);
                            }
                            else
                            {
                                ++it;
                            }
                        }

                        BCUpdateLinkCostDiff2(LinkFlowAdj_thr);
                    }
                }
                // if (m_nViolation1D == 0)
                //     break;
            }
        }

        double PFA_end = omp_get_wtime();
        PFA_period += PFA_end - PFA_beg;
        double RG_beg = omp_get_wtime();
        UEGap = GetUEGap();
        double RG_end = omp_get_wtime();
        RG_period += RG_end - RG_beg;
        MaxODGap = UEGap / GapRate;
        finish = omp_get_wtime(); // End timing.
        CPUTime = finish - start; // Compute runtime.
        obj = ObjectiveValue();
        cout << setprecision(20) << obj << ", " << nOutLoop << "," << setprecision(10) << UEGap << ", " << CPUTime
             << ", " << PG_period << ", " << PFA_period << ", " << RG_period << endl;
    }

    delete[] LinkFlowAdj_thr;
    //********************************
    printf("%s\t", "wtime:");
    printf("%f\n", CPUTime);
    printf("%s\t", "nOutLoop:");
    printf("%d\n", nOutLoop);
    printf("%s\t", "UEGap:");
    printf("%e\n", UEGap);
    printf("%s\t", "LinkFlow[0]:");
    printf("%f\n", LinkFlow[0]);
}

#pragma endregion

int main(int argc, char *argv[])
{
    // At least two arguments, argv[1] and argv[2], are required, so argc should be checked against 3.
    if (argc < 3)
    {
        std::cout << "Please provide an integer argument." << std::endl;
        return 1;
    }

    int number = std::stoi(argv[1]);
    int number2 = std::stoi(argv[2]);
    // Phi 1526; Birm 899; CS 1; Anaheim:39 SF:1
    if (number == 0) // SiouxFall
    {
        int first_node = 1;
        string fileName("F:\\VC-Eq-PBCD\\data\\SF");
        string simName("07");
        ReadData(fileName, simName);
        // elseparalevel = 64;
        UE_PBCD();
    }
    else if (number == 1) //  Anaheim
    {
        int first_node = 39;
        PG_num_threads = number2; // omp_get_max_threads()
        PG_num_threads_else = PG_num_threads;
        string fileName("F:\\VC-Eq-PBCD\\data\\Ana");
        string simName("09");
        ReadData(fileName, simName);

        UE_PBCD();
    }
    else if (number == 2) //  CS
    {
        int first_node = 1;
        PG_num_threads = number2; // omp_get_max_threads()
        PG_num_threads_else = PG_num_threads;
        string fileName("F:\\VC-Eq-PBCD\\data\\CS");
        string simName("07");
        ReadData(fileName, simName);

        UE_PBCD();
    }
    else if (number == 3) //  Birmingham
    {
        int first_node = 899;
        PG_num_threads = number2; // omp_get_max_threads()
        PG_num_threads_else = PG_num_threads;
        string fileName("F:\\VC-Eq-PBCD\\data\\Birm");
        string simName("08");
        ReadData(fileName, simName);

        UE_PBCD();
    }
    else if (number == 4) //  Philadelphia
    {
        int first_node = 1526;
        PG_num_threads = number2; // omp_get_max_threads()
        PG_num_threads_else = PG_num_threads;
        string fileName("F:\\VC-Eq-PBCD\\data\\Phi");


        string simName("09");
        ReadData(fileName, simName);

        UE_PBCD();
    }

    return 0;
}
