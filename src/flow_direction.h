#ifndef FLOW_DIRECTION_H
#define FLOW_DIRECTION_H

#include <algorithm>
#include <vector>

#include "commonLib.h"
#include "linearpart.h"
#include "sparsepartition.h"

#include "arp.h"

using std::vector;

// Checks if cells cross
bool cellsCross(int k, int i, int j, linearpart<short>& flowDir)
{
    int n1, c1, n2, c2;

    switch(k) {
    case 2:
        n1=1;
        c1=4;
        n2=3;
        c2=8;
        break;
    case 4:
        n1=3;
        c1=6;
        n2=5;
        c2=2;
        break;
    case 6:
        n1=7;
        c1=4;
        n2=5;
        c2=8;
        break;
    case 8:
        n1=1;
        c1=6;
        n2=7;
        c2=2;
        break;
    default:
        return 0;
    }

    int in1=i+d1[n1];
    int jn1=j+d2[n1];
    int in2=i+d1[n2];
    int jn2=j+d2[n2];

    if (flowDir.getData(in1,jn1) == c1 || flowDir.getData(in2,jn2) == c2)
    {
        return true;
    }

    return false;
}

// FIXME: D-inf, should we care if cells cross?
bool cellsCross(int k, int i, int j, linearpart<float>& flowDir)
{
    return false;
}

template<typename Algo, typename E>
int markPits(E& elev, linearpart<typename Algo::FlowType>& flowDir, vector<vector<node>>& islands, SparsePartition<int>& inc);

template<typename Algo>
size_t propagateIncrements(linearpart<typename Algo::FlowType>& flowDir, SparsePartition<int>& inc, vector<node>& queue) {
    size_t numInc = 0;
    int st = 1;
    
    vector<node> newFlats;
    while (!queue.empty()) {
        for(node flat : queue) {
            // Duplicate. already set
            if (inc.getData(flat.x, flat.y) > 0)
                continue;

            for (int k = 1; k <= 8; k++) {
                if (cellsCross(k, flat.x, flat.y, flowDir))
                    continue;

                int in = flat.x + d1[k];
                int jn = flat.y + d2[k];

                if (!flowDir.isInPartition(in, jn)) 
                    continue;

                typename Algo::FlowType flow = flowDir.getData(in,jn);

                if (!Algo::HasFlow(flow) && inc.getData(in, jn) == 0) {
                    newFlats.push_back(node(in, jn));
                    inc.setData(in, jn, -1);
                }
            }

            numInc++;
            inc.setData(flat.x, flat.y, st);
        }

        queue.clear();
        queue.swap(newFlats);
        st++;
    }

    if (st < 0) {
        printf("WARNING: increment overflow during propagation (st = %d)\n", st);
    }

    return numInc;
}

template<typename Algo>
size_t propagateBorderIncrements_async(linearpart<typename Algo::FlowType>& flowDir, SparsePartition<int>& inc, vector<node>& changes, bool border_changed[4])
{
    int nx = flowDir.getnx();
    int ny = flowDir.getny();

    vector<node> queue, newFlats;

    // Find the starting nodes at the edge of the raster
    int ignoredGhostCells = 0;

    auto checkFlat = [&](int in, int jn, int st) {
        if (!flowDir.isInPartition(in, jn))
            return;

        bool noFlow = !Algo::HasFlow(flowDir.getData(in, jn));
        auto neighInc = inc.getData(in, jn);

        if (noFlow && (neighInc == 0 || std::abs(neighInc) > st + 1)) {
            // If neighbor increment is positive, we are overriding a larger increment
            // and it is not yet in the queue
            if (neighInc >= 0) {
                queue.emplace_back(in, jn);
            }

            // Here we set a negative increment if it's still not set
            //
            // Another flat might be neighboring the same cell with a lower increment,
            // which has to override the higher increment (that hasn't been set yet but is in the queue).
            inc.setData(in, jn, -(st + 1));
        }
    };

    for (auto flat : changes) {
        int st = inc.getData(flat.x, flat.y);

        if (st == 0)
            continue;

        if (st == INT_MAX) {
            ignoredGhostCells++;
            continue;
        }

        bool y_border = flat.x == -1 || flat.x == nx;
        bool x_border = flat.y == -1 || flat.y == ny;
        auto in = flat.x == -1 ? 0 : nx - 1;
        auto jn = flat.y == -1 ? 0 : ny - 1;

        if (x_border && y_border) {
            // single pixel neighbor at a corner
            checkFlat(in, jn, st);
        } else if (x_border) {
            for (auto in : {flat.x - 1, flat.x, flat.x + 1}) { 
                checkFlat(in, jn, st);
            }
        } else if (y_border) {
            for (auto jn : {flat.y - 1, flat.y, flat.y + 1}) { 
                checkFlat(in, jn, st);
            }
        } else {
            printf("error: received cell change not at border (%d %d)\n", flat.x, flat.y);
        }
    }

    if (ignoredGhostCells > 0) {
       printf("warning: ignored %d ghost cells which were at upper limit (%d)\n", ignoredGhostCells, INT_MAX);
    }

    size_t numChanged = 0;
    size_t abandonedCells = 0;

    while (!queue.empty()) {
        for(node flat : queue) {
            // Increments are stored as negative for the cells that have been added to the queue
            // (to signify that they need to be explored)
            auto st = -inc.getData(flat.x, flat.y);

            // I don't think this is possible anymore, but just in case.
            if (st <= 0) {
                printf("warning: unexpected non-negative increment @ (%d, %d) - %d\n", flat.x, flat.y, -st);
                continue;
            }

            inc.setData(flat.x, flat.y, st);
            numChanged++;

            border_changed[0] |= flat.y == 0;
            border_changed[1] |= flat.x == 0;
            border_changed[2] |= flat.x == nx - 1;
            border_changed[3] |= flat.y == ny - 1;

            if (st == INT_MAX) {
                abandonedCells++;
                continue;
            }

            for (int k = 1; k <= 8; k++) {
                if (cellsCross(k, flat.x, flat.y, flowDir))
                    continue;

                int in = flat.x + d1[k];
                int jn = flat.y + d2[k];

                if (!flowDir.isInPartition(in, jn))
                    continue;

                bool noFlow = !Algo::HasFlow(flowDir.getData(in, jn));
                auto neighInc = inc.getData(in, jn);

                if (noFlow && (neighInc == 0 || std::abs(neighInc) > st + 1)) {
                    // If neighbor increment is positive, we are overriding a larger increment
                    // and it is not yet in the queue
                    if (neighInc >= 0) {
                        newFlats.emplace_back(in, jn);
                    }

                    inc.setData(in, jn, -(st + 1));
                }
            }
        }

        queue.clear();
        queue.swap(newFlats);
    }

    if (abandonedCells > 0) {
        printf("warning: gave up propagating %zu cells because they were at upper limit (%d)\n", abandonedCells, INT_MAX);
    }

    return numChanged;
}

template<typename Algo, typename E>
void flowTowardsLower(E& elev, linearpart<typename Algo::FlowType>& flowDir, vector<vector<node>>& islands, SparsePartition<int>& inc)
{
    long nx = flowDir.getnx();
    long ny = flowDir.getny();

    vector<node> lowBoundaries;

    // Find low boundaries. 
    for(auto& island : islands) {
        for(node flat : island) {
            float flatElev = elev.getData(flat.x, flat.y);

            for (int k = 1; k <= 8; k++) {
                if (cellsCross(k, flat.x, flat.y, flowDir))
                    continue;

                int in = flat.x + d1[k];
                int jn = flat.y + d2[k];

                if (!flowDir.hasAccess(in, jn)) {
                    continue;
                }

                auto elevDiff = flatElev - elev.getData(in,jn);
                typename Algo::FlowType flow = flowDir.getData(in,jn);

                // Adjacent cell drains and is equal or lower in elevation so this is a low boundary
                if (elevDiff >= 0 && Algo::HasFlow(flow)) {
                    lowBoundaries.push_back(flat);
                    inc.setData(flat.x, flat.y, -1);

                    // No need to check the other neighbors
                    break;
                } 
            }
        }
    }

    size_t numInc = propagateIncrements<Algo>(flowDir, inc, lowBoundaries);
}

template<typename Algo, typename E>
void flowFromHigher(E& elev, linearpart<typename Algo::FlowType>& flowDir, vector<vector<node>>& islands, SparsePartition<int>& inc)
{
    long nx = flowDir.getnx();
    long ny = flowDir.getny();

    vector<node> highBoundaries;

    // Find high boundaries
    for (auto& island : islands) {
        for (node flat : island) {
            float flatElev = elev.getData(flat.x, flat.y);
            bool highBoundary = false;

            for (int k = 1; k <= 8; k++) {
                if (cellsCross(k, flat.x, flat.y, flowDir))
                    continue;

                int in = flat.x + d1[k];
                int jn = flat.y + d2[k];

                if (!flowDir.hasAccess(in, jn))
                    continue;

                auto elevDiff = flatElev - elev.getData(in, jn);

                if (elevDiff < 0) {
                    // Adjacent cell has higher elevation so this is a high boundary
                    highBoundary = true;
                    break;
                }
            }

            if (highBoundary) {
                inc.setData(flat.x, flat.y, -1);
                highBoundaries.push_back(flat);
            }
        }
    }

    propagateIncrements<Algo>(flowDir, inc, highBoundaries);
}

template<typename Algo, typename E>
int markPits(E& elev, linearpart<typename Algo::FlowType>& flowDir, vector<vector<node>>& islands, SparsePartition<int>& inc)
{
    int nx = flowDir.getnx();
    int ny = flowDir.getny();

    int numPits = 0;

    //There are pits remaining - set direction to no data
    for (auto& island : islands) {
        for (node flat : island) {
            bool skip = false;

            for (int k=1; k<=8; k++) {
                if (cellsCross(k, flat.x, flat.y, flowDir))
                    continue;

                int in = flat.x + d1[k];
                int jn = flat.y + d2[k];

                if (!flowDir.hasAccess(in, jn)) 
                    continue;

                auto elevDiff = elev.getData(flat.x, flat.y) - elev.getData(in, jn);
                typename Algo::FlowType flow = flowDir.getData(in,jn);

                // Adjacent cell drains and is equal or lower in elevation so this is a low boundary
                if (elevDiff >= 0 && Algo::HasFlow(flow)) {
                    skip = true;
                    break;
                } else if (!Algo::HasFlow(flow)) {
                    // If neighbor is in flat

                    if (inc.getData(in,jn) >= 0){ 
                        skip = true;
                        break;
                    }
                }
            }
            
            // mark pit
            if (!skip) {
                numPits++;
                flowDir.setToNodata(flat.x, flat.y);
            }  
        }
    }

    return numPits;
}

template<typename Algo>
void findIslands(linearpart<typename Algo::FlowType>& flowDir, vector<vector<node>>& localIslands, vector<vector<node>>& sharedIslands)
{
    int globalWidth = flowDir.gettotalx();
    int globalHeight = flowDir.gettotaly();
    int width = flowDir.getnx();
    int height = flowDir.getny();

    int numIslands = 0;
    SparsePartition<int> islandLabel(globalWidth, globalHeight, 0);
    vector<node> q, tempVector;

    for (int j=0; j<height; j++) {
        for(int i=0; i<width; i++) {
            if (Algo::HasFlow(flowDir.getData(i, j))) {
                continue;
            }

            node flat(i, j);

            if (islandLabel.getData(flat.x, flat.y) != 0) {
                continue;
            }

            bool shared = false;
            int label = ++numIslands;
            
            q.push_back(flat);

            while(!q.empty()) {
                node flat = q.back();
                q.pop_back();

                if (islandLabel.getData(flat.x, flat.y) != 0) {
                    continue;
                }

                islandLabel.setData(flat.x, flat.y, label);
                tempVector.push_back(flat);

                for (int k=1; k<=8; k++) {
                    int in = flat.x + d1[k];
                    int jn = flat.y + d2[k];

                    bool ghostCell = (in == -1) || (in == width) || (jn == -1) || (jn == height);

                    if (ghostCell && flowDir.hasAccess(in, jn)) {
                        if (!Algo::HasFlow(flowDir.getData(in, jn)))
                        {
                            shared = true;
                        }
                    }

                    if (!flowDir.isInPartition(in, jn))
                        continue;

                    if (!Algo::HasFlow(flowDir.getData(in, jn)))
                        q.push_back(node(in, jn));
                }
            }

            if (!shared) {
                localIslands.push_back(vector<node>(tempVector));
            } else {
                sharedIslands.push_back(vector<node>(tempVector));
            }

            tempVector.clear();
        }
    }
}

template<typename Algo, typename E>
long resolveFlats(E& elev, SparsePartition<int>& inc, linearpart<typename Algo::FlowType>& flowDir, std::vector<std::vector<node>>& islands)
{
    int globalWidth = flowDir.gettotalx();
    int globalHeight = flowDir.gettotaly();
    
    int rank;
    MPI_Comm_rank(MCW, &rank);
    
    if (rank==0) {
        fprintf(stderr,"Resolving flats\n");
        fflush(stderr);
    }

    SparsePartition<int> s(globalWidth, globalHeight, 0);
    
    flowTowardsLower<Algo>(elev, flowDir, islands, inc);

    // Not all grid cells were resolved - pits remain
    // Remaining grid cells are unresolvable pits
    markPits<Algo>(elev, flowDir, islands, inc);

    // Drain flats away from higher adjacent terrain
    flowFromHigher<Algo>(elev, flowDir, islands, s);

    // High flow must be inverted before it is combined
    //
    // higherFlowMax has to be greater than all of the increments
    // higherFlowMax can be maximum value of the data type but it will cause overflow problems if more than one iteration is needed
    int higherFlowMax = 0;

    for (auto& island : islands) {
        for (node flat : island) {    
            int val = s.getData(flat.x, flat.y);

            if (val > higherFlowMax)
                higherFlowMax = val;
        }
    }

    for (auto& island : islands) {
        for (auto flat : island) {
            inc.addToData(flat.x, flat.y, higherFlowMax - s.getData(flat.x, flat.y));
        }
    }

    if (rank==0) {
        fprintf(stderr,"Setting directions\n");
        fflush(stderr);
    }

    long flatsRemaining = 0;
    for (auto& island : islands) {
        for (node flat : island) {
            Algo::SetFlow(flat.x, flat.y, flowDir, elev, inc);

            if (!Algo::HasFlow(flowDir.getData(flat.x, flat.y))) {
                flatsRemaining++;
            }
        }
    }

    auto hasFlowDirection = [&](const node& n) { return Algo::HasFlow(flowDir.getData(n.x, n.y)); };
    auto isEmpty = [&](const std::vector<node>& i) { return i.empty(); };
    
    // Remove flats which have flow direction set
    for (auto& island : islands) {
        island.erase(std::remove_if(island.begin(), island.end(), hasFlowDirection), island.end());
    }

    // Remove empty islands
    islands.erase(std::remove_if(islands.begin(), islands.end(), isEmpty), islands.end());

    return flatsRemaining;
}

template<typename Algo, typename E>
long resolveFlats_parallel_async(E& elev, SparsePartition<int>& inc, linearpart<typename Algo::FlowType>& flowDir, vector<vector<node>>& islands)
{
    int globalWidth = flowDir.gettotalx();
    int globalHeight = flowDir.gettotaly();
    
    
    int rank;
    MPI_Comm_rank(MCW, &rank);

    SparsePartition<int> higherGradient(globalWidth, globalHeight, 0);

    flowTowardsLower<Algo>(elev, flowDir, islands, inc);
    flowFromHigher<Algo>(elev, flowDir, islands, higherGradient);

    {
        AsyncRasterProcessor arp;

        arp.add(&inc, [&flowDir, &inc](std::vector<node>& border_diff, bool borders_updated[4]) {
            //printf("rank %d: got low gradient update - %d cells\n", rank, border_diff.size());

            propagateBorderIncrements_async<Algo>(flowDir, inc, border_diff, borders_updated);
        });

        arp.add(&higherGradient, [&flowDir, &higherGradient](std::vector<node>& border_diff, bool borders_updated[4]) {
            //printf("rank %d: got high gradient update - %d cells\n", rank, border_diff.size());

            propagateBorderIncrements_async<Algo>(flowDir, higherGradient, border_diff, borders_updated);
        });

        arp.run();
     }

    // Not all grid cells were resolved - pits remain
    // Remaining grid cells are unresolvable pits
    markPits<Algo>(elev, flowDir, islands, inc);

    // High flow must be inverted before it is combined
    //
    // higherFlowMax has to be greater than all of the increments
    // higherFlowMax can be maximum value of the data type (e.g. 65535) but it will cause overflow problems if more than one iteration is needed
    int higherFlowMax = 0;

    for (auto& island : islands) {
        for (auto& flat : island) {
            int val = higherGradient.getData(flat.x, flat.y);
        
            if (val > higherFlowMax)
                higherFlowMax = val;
        }
    }

    // FIXME: Is this needed? would it affect directions at the border?
    // It is local to a flat area, but can that be reduced further to minimize comm?
    int globalHigherFlowmax = 0;
    MPI_Allreduce(&higherFlowMax, &globalHigherFlowmax, 1, MPI_INT, MPI_MAX, MCW);

    size_t badCells = 0;

    for (auto& island : islands) {
        for (auto flat : island) {
            auto val = inc.getData(flat.x, flat.y);
            auto highFlow = higherGradient.getData(flat.x, flat.y);

            inc.setData(flat.x, flat.y, val + (globalHigherFlowmax - highFlow));

            if (val < 0 || val == INT_MAX || highFlow < 0 || highFlow == INT_MAX) {
                badCells++;
            }
        }
    }

    if (badCells > 0) {
        printf("warning rank %d: %d increment values either incorrect or overflown\n", rank, badCells);
    }

    inc.share();

    if (rank==0) {
        fprintf(stderr,"\nPRL: Setting directions\n");
        fflush(stderr);
    }

    uint64_t localFlatsRemaining = 0, globalFlatsRemaining = 0;

    for (auto& island : islands) {
        for (node flat : island) {
            Algo::SetFlow(flat.x, flat.y, flowDir, elev, inc);
    
            if (!Algo::HasFlow(flowDir.getData(flat.x, flat.y))) {
                localFlatsRemaining++;
            }
        }
    }

    flowDir.share();
    MPI_Allreduce(&localFlatsRemaining, &globalFlatsRemaining, 1, MPI_UINT64_T, MPI_SUM, MCW);

    auto hasFlowDirection = [&](const node& n) { return Algo::HasFlow(flowDir.getData(n.x, n.y)); };
    auto isEmpty = [&](const std::vector<node>& i) { return i.empty(); };
    
    // Remove flats which have flow direction set
    for (auto& island : islands) {
        island.erase(std::remove_if(island.begin(), island.end(), hasFlowDirection), island.end());
    }

    // Remove empty islands
    islands.erase(std::remove_if(islands.begin(), islands.end(), isEmpty), islands.end());

    return globalFlatsRemaining;
}

#endif
