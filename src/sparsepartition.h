#ifndef SPARSEPARTITION_H
#define SPARSEPARTITION_H

// Sparse partition that internally stores the raster as tiles. 
// 
// It greatly reduces the memory footprint for temporary
// rasters that are used sparsely during computation.

#include <cstring>

#include <algorithm>
#include <vector>
#include <memory>

#include <mpi.h>

#include "const.h"
#include "commonLib.h"
#include "linearpart.h"

using std::unique_ptr;
using std::vector;

const int BLOCK_SIZE_BITS = 8;

const int BLOCK_SIZE = 1 << BLOCK_SIZE_BITS;
const int BLOCK_MASK = ~(BLOCK_SIZE - 1);

class AsyncPartition { 
    public:
        virtual size_t asyncGetBufferSize() = 0;
        virtual size_t asyncBytesToSend(int count) = 0;
        virtual int asyncGetChanges(int neighbor, void* buffer) = 0;
        virtual void asyncStoreChanges(int neighbor, std::vector<node>& changes) = 0;

        virtual int neighborRank(int neighbor) = 0;
        virtual int queueRecv(int neighbor, int tag, MPI_Request* req) = 0;
};

constexpr int getBlockCount(int pixels) {
    return ((pixels + BLOCK_SIZE - 1) & BLOCK_MASK) >> BLOCK_SIZE_BITS;
}

template<typename T>
struct Cell {
    int position;
    T value;
};

// FIXME: get rid of this and use MPI message sizes
template<typename T>
struct BorderUpdate {
    int count;
    Cell<T> cells[];
};

template<typename T>
class SparsePartition : public AsyncPartition {
    public:
        SparsePartition(T noData) : noDataValue(noData) {}

        SparsePartition(int globalWidth, int globalHeight, T noData)
            : globalWidth(globalWidth), globalHeight(globalHeight), noDataValue(noData)
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            
            findNeighbors();

            widthBlocks = getBlockCount(width + 2*BLOCK_SIZE);
            int blocksNeeded = widthBlocks * getBlockCount(height + 2*BLOCK_SIZE);
            blocks.resize(blocksNeeded);

            allocBuffers();
        }

        SparsePartition(const SparsePartition& bp) = delete;
        SparsePartition(SparsePartition&& bp) = default;

        SparsePartition& operator=(const SparsePartition&) = delete;
        SparsePartition& operator=(SparsePartition&& bp) = default;
        
        T& getPixel(int gx, int gy) {
            T* blockData = getBlock(gx, gy, true);

            int x = gx & ~BLOCK_MASK;
            int y = gy & ~BLOCK_MASK;

            return blockData[x + y*BLOCK_SIZE];
        }

        T getData(int gx, int gy) {
            if (getBlock(gx, gy) == nullptr)
                return noDataValue;

            return getPixel(gx, gy);
        }

        void setData(int gx, int gy, T val) {
            getPixel(gx, gy) = val;
        }

        void addToData(int gx, int gy, T val) {
            getPixel(gx, gy) += val;
        }

        void share() {
            if (size == 1) return;

            const int dims[9] = { 1, width, 1, height, 0, height, 1, width, 1 };

            MPI_Request requests[9 * 2];
            std::fill_n(requests, 9*2, MPI_REQUEST_NULL);

            for (int i = 0; i < 9; i++) {
                if (neighbors[i] == -1)
                    continue;

                MPI_Irecv(receiveBuf[i].get(), dims[i] * sizeof(T), MPI_BYTE,
                    neighbors[i], 0, MPI_COMM_WORLD, &requests[i]);
            }

            for (int i = 0; i < 9; i++) {
                if (neighbors[i] == -1)
                    continue;

                // Load the data
                loadBorderBuf(diffBuf[i].get(), i); 

                MPI_Isend(diffBuf[i].get(), dims[i] * sizeof(T), MPI_BYTE,
                    neighbors[i], 0, MPI_COMM_WORLD, &requests[i + 9]);
            }

            MPI_Waitall(2 * 9, requests, MPI_STATUSES_IGNORE);

            for (int i = 0; i < 9; i++) {
                if (neighbors[i] == -1)
                    continue;
                
                storeBorderBuf((T*) receiveBuf[i].get(), i);
            }
        }

        // Async stuff
        size_t asyncGetBufferSize() {
            // int count + width * (int pos + T value)
            return sizeof(int) + std::max(height, width) * (sizeof(int) + sizeof(T));
        }

        size_t asyncBytesToSend(int count) {
            // int + count * (int + T)
            return sizeof(int) + count * (sizeof(int) + sizeof(T));
        }
        
        // rename asyncSerializeChanges
        int asyncGetChanges(int neighbor, void* buffer) {
            if (neighbors[neighbor] == -1)
                return 0;

            BorderUpdate<T>* update = (BorderUpdate<T>*) buffer;
           
            update->count = 0;

            switch(neighbor) {
                case 1:
                    diffHorizontal(update, 1, 0);
                    break;

                case 3:
                    diffVertical(update, 3, 0);
                    break;

                case 5:
                    diffVertical(update, 5, width - 1);
                    break;

                case 7:
                    diffHorizontal(update, 7, height - 1);
                    break;

                case 0:
                case 2:
                case 6:
                case 8:
                    diffCorner(update, neighbor);
                    break;
            }

            return update->count;
        }

        void diffHorizontal(BorderUpdate<T>* update, int neighbor, int y) {
            for (int x = 0; x < width; x += BLOCK_SIZE) {
                T* block = getStride(x, y, false);

                if (block == nullptr) {
                    continue;
                }
                    
                int sz = std::min(width - x, BLOCK_SIZE);
                for(int i = 0; i < sz; i++) {
                    if (block[i] == diffBuf[neighbor][x + i]) {
                        continue;
                    }

                    auto& cell = update->cells[update->count++];
                    cell.position = x + i;
                    cell.value = block[i];

                    diffBuf[neighbor][x + i] = block[i];
                }
            }
        }
        
        void diffVertical(BorderUpdate<T>* update, int neighbor, int x) {
            int localX = x & ~BLOCK_MASK;

            for (int y = 0; y < height; y += BLOCK_SIZE) {
                T* block = getBlock(x, y, false);

                if (block == nullptr) {
                    continue;
                }
                    
                int sz = std::min(height - y, BLOCK_SIZE);
                for(int i = 0; i < sz; i++) {
                    auto val = block[i * BLOCK_SIZE + localX];

                    if (val == diffBuf[neighbor][y + i]) {
                        continue;
                    }

                    auto& cell = update->cells[update->count++];
                    cell.position = y + i;
                    cell.value = val;

                    diffBuf[neighbor][y + i] = val;
                }
            }
        }
        
        void diffCorner(BorderUpdate<T>* update, int neighbor) {
            int x, y;

            switch(neighbor) {
                case 0:
                    x = 0;
                    y = 0;
                    break;

                case 2:
                    x = width - 1;
                    y = 0;
                    break;

                case 6:
                    x = 0;
                    y = height - 1;
                    break;

                case 8:
                    x = width - 1;
                    y = height - 1;
                    break;

                default: printf("rank %d: unknown corner neighbor %d\n", rank, neighbor);
            }

            T currentValue = getPixel(x, y);

            if (diffBuf[neighbor][0] != currentValue) {
                update->cells[update->count++] = { 0, currentValue };
                diffBuf[neighbor][0] = currentValue;
            }
        }

        void asyncStoreChanges(int neighbor, std::vector<node>& changes) {
            if (neighbors[neighbor] == -1) {
                printf("rank %d: storing changes for unknown neighbor %d\n", rank, neighbor);
                return;
            }

            BorderUpdate<T>* update = (BorderUpdate<T>*) receiveBuf[neighbor].get();

            using wtf = void (*)(int, int, int, int&, int&);

            const wtf get_xy[] = {
                [](int pos, int width, int height, int& x, int& y) { x = -1; y = -1; }, 
                [](int pos, int width, int height, int& x, int& y) { x = pos; y = -1; }, 
                [](int pos, int width, int height, int& x, int& y) { x = width; y = -1; }, 
                
                [](int pos, int width, int height, int& x, int& y) { x = -1; y = pos; },
                [](int pos, int width, int height, int& x, int& y) {  },
                [](int pos, int width, int height, int& x, int& y) { x = width; y = pos; },

                [](int pos, int width, int height, int& x, int& y) { x = -1; y = height; }, 
                [](int pos, int width, int height, int& x, int& y) { x = pos; y = height; }, 
                [](int pos, int width, int height, int& x, int& y) { x = width; y = height; }, 
            };

            for (int n = 0; n < update->count; n++) {
                auto& c = update->cells[n];

                int x, y;
                get_xy[neighbor](c.position, width, height, x, y);

                setData(x, y, c.value);
                changes.emplace_back(x, y);
            }
        }

        int neighborRank(int neighbor) {
            return neighbors[neighbor];
        }

        int queueRecv(int neighbor, int tag, MPI_Request* req) {
            if (neighbors[neighbor] == -1)
                return -1;

            return MPI_Irecv(receiveBuf[neighbor].get(), receiveBufSizes[neighbor],
                    MPI_BYTE, neighbors[neighbor], tag, MPI_COMM_WORLD, req); 
        }
    
    private:
        void findNeighbors() {
            // Find number of partitions in both dimensions
            int sizeX, sizeY;
            if (tdpartition::decompType == DECOMP_BLOCK) {
                int firstFactor, secondFactor;
                // find the closest pair of factors of processor size
                findClosestFactors(size, firstFactor, secondFactor);

                // largest factor goes to largest dimension
                if (globalHeight > globalWidth) {
                    sizeX = std::min(firstFactor, secondFactor);
                    sizeY = std::max(firstFactor, secondFactor);
                } else {
                    sizeX = std::max(firstFactor, secondFactor);
                    sizeY = std::min(firstFactor, secondFactor);
                }
            } else if (tdpartition::decompType == DECOMP_ROW) {
                sizeX = 1;
                sizeY = size;
            } else {
                sizeX = size;
                sizeY = 1;
            }

            int coordX = rank % sizeX;
            int coordY = rank / sizeX;

            int cx = (globalWidth / sizeX) * coordX;
            int cy = (globalHeight / sizeY) * coordY;
            
            // Set widths
            width = coordX != (sizeX - 1) ? globalWidth / sizeX : globalWidth - cx;
            height = coordY != (sizeY - 1) ? globalHeight / sizeY : globalHeight - cy;

            //printf("rank %d is at %d,%d - w x h = %d %d\n", rank, coordX, coordY, width, height);
            
            int startX = std::max(0, coordX - 1);
            int startY = std::max(0, coordY - 1);
            int endX = std::min(sizeX, coordX + 2);
            int endY = std::min(sizeY, coordY + 2);

            for(int x = startX; x < endX; x++) {
                for(int y = startY; y < endY; y++) {
                    // Skip current rank
                    if (x == coordX && y == coordY)
                        continue;

                    int n = (y-coordY+1) * 3 + (x-coordX+1);
                    neighbors[n] = y * sizeX + x;
                   
                    //printf("rank %d: neighbor %d %d (%d) is at %d\n", rank, x, y, n, neighbors[n]);
                }
            }
        }

        void allocBuffers() {
            const int dims[9] = { 1, width, 1, height, 0, height, 1, width, 1 };
            int totalSz = 0;

            for(int i = 0; i < 9; i++) {
                if (neighbors[i] == -1)
                    continue;

                int sz = sizeof(int) + dims[i] * (sizeof(int) + sizeof(T));
                totalSz += sz + dims[i] * sizeof(T);
    
                diffBuf[i] = unique_ptr<T[]>(new T[dims[i]]);
                std::fill(diffBuf[i].get(), diffBuf[i].get() + dims[i], noDataValue);

                receiveBuf[i] = unique_ptr<uint8_t[]>(new uint8_t[sz]);
                receiveBufSizes[i] = sz;
            }

            //printf("Allocated %s bytes for SparsePartition\n", humanReadableSize(totalSz).c_str());
        }

        void loadBorderBuf(T* buf, int neighbor) {
            switch(neighbor) {
                case 0:
                    *buf = getPixel(0, 0);
                    break;

                case 1:
                    loadHorizontalBuf(buf, 0);
                    break;

                case 2:
                    *buf = getPixel(width - 1, 0);
                    break;

                case 3:
                    loadVerticalBuf(buf, 0);
                    break;

                case 5:
                    loadVerticalBuf(buf, width - 1);
                    break;

                case 6:
                    *buf = getPixel(0, height - 1);
                    break;

                case 7:
                    loadHorizontalBuf(buf, height - 1);
                    break;

                case 8:
                    *buf = getPixel(width - 1, height - 1);
                    break;
            
                default:
                    printf("WARNING: loadBorderBuf unknown neighbor %d\n", neighbor);
                    break;
            }
        }

        void loadHorizontalBuf(T* buf, int y) {
            for (int x = 0; x < width; x += BLOCK_SIZE) {
                T* block = getStride(x, y, false);

                if (block == nullptr) {
                    for(int i = 0; i < BLOCK_SIZE; i++) {
                        if (x + i >= width)
                            break;

                        buf[x + i] = noDataValue;
                    }
                } else {
                    int sz = std::min(width - x, BLOCK_SIZE) * sizeof(T);

                    memcpy(&buf[x], block, sz);
                }
            }    
        }
        
        void loadVerticalBuf(T* buf, int x) {
            int localX = x & ~BLOCK_MASK;

            for (int y = 0; y < height; y += BLOCK_SIZE) {
                T* block = getBlock(x, y, false);

                int sz = std::min(height - y, BLOCK_SIZE);

                if (block == nullptr) {
                    for(int i = 0; i < sz; i++) {
                        buf[y + i] = noDataValue;
                    }
                } else {
                    for(int i = 0; i < sz; i++) {
                        buf[y + i] = block[i * BLOCK_SIZE + localX];
                    }
                }
            }    
        }

        void storeBorderBuf(T* buf, int neighbor) {
            switch(neighbor) {
                case 0:
                    getPixel(-1, -1) = *buf;
                    break;

                case 1:
                    storeHorizontalBuf(buf, -1);
                    break;

                case 2:
                    getPixel(width, -1) = *buf;
                    break;

                case 3:
                    storeVerticalBuf(buf, -1);
                    break;

                case 5:
                    storeVerticalBuf(buf, width);
                    break;

                case 6:
                    getPixel(-1, height) = *buf;
                    break;

                case 7:
                    storeHorizontalBuf(buf, height);
                    break;

                case 8:
                    getPixel(width, height) = *buf;
                    break;

                default:
                    printf("WARNING: storeBorderBuf unknown neighbor %d\n", neighbor);
                    break;
            }
        }

        void storeHorizontalBuf(T* buf, int y) {
            for (int x = 0; x < width; x += BLOCK_SIZE) {
                T* block = getStride(x, y, false);

                if (block == nullptr) {
                    // Check if all values are no data to save memory
                    for(int i = 0; i < BLOCK_SIZE; i++) {
                        if (buf[x + i] != noDataValue)
                        {
                            block = getStride(x, y, true);
                            break;
                        }
                    }
                }

                if (block != nullptr) {
                    memcpy(block, &buf[x], sizeof(T)*BLOCK_SIZE);
                }
            }
        }
        
        void storeVerticalBuf(T* buf, int x) {
            int localX = x & ~BLOCK_MASK;
            for (int y = 0; y < height; y += BLOCK_SIZE) {
                T* block = getBlock(x, y, false);

                if (block == nullptr) {
                    // Check if all values are no data to save memory
                    for(int i = 0; i < BLOCK_SIZE; i++) {
                        if (buf[localX + i*BLOCK_SIZE] != noDataValue)
                        {
                            block = getBlock(x, y, true);
                            break;
                        }
                    }
                }

                if (block != nullptr) {
                    int sz = std::min(height - y, BLOCK_SIZE);

                    for(int i = 0; i < sz; i++) {
                        block[i*BLOCK_SIZE + localX] = buf[y + i];
                    }
                }
            }
        }

        T* getBlock(int gx, int gy, bool create=false) { 
            if (blocks.empty())
                return nullptr;

            gx += BLOCK_SIZE;
            gy += BLOCK_SIZE;

            int bx = (gx & BLOCK_MASK) >> BLOCK_SIZE_BITS;
            int by = (gy & BLOCK_MASK) >> BLOCK_SIZE_BITS;

            int id = bx + by * widthBlocks;

            T* blockData = blocks[id].get();

            if (blockData == nullptr && create) {
                // Create new block
                blocks[id] = unique_ptr<T[]>(new T[BLOCK_SIZE * BLOCK_SIZE]);
                blockData = blocks[id].get();

                std::fill(blockData, blockData + BLOCK_SIZE*BLOCK_SIZE, noDataValue);
            }

            return blockData;
        }

        T* getStride(int gx, int gy, bool create=false) {
            T* block = getBlock(gx, gy, create);

            if (block == nullptr)
                return nullptr;

            int y = gy & ~BLOCK_MASK;

            return &block[y*BLOCK_SIZE];
        }

        int rank = -1;
        int size = 0;

        MPI_Comm cartComm = MPI_COMM_NULL;
        int neighbors[9] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };

        int globalWidth = 0;
        int globalHeight = 0;
        T noDataValue;

        int width = 0;
        int widthBlocks = 0;
        int height = 0;

        std::vector<unique_ptr<T[]>> blocks;
        unique_ptr<uint8_t[]> receiveBuf[9];
        int receiveBufSizes[9] = {};
        unique_ptr<T[]> diffBuf[9];
};

#endif //SPARSEPARTITION_H
