#ifndef ARP_H
#define ARP_H

#include <algorithm>
#include <memory>
#include <vector>

#include <mpi.h>

#include "commonLib.h"
#include "sparsepartition.h"

using std::vector;
using std::unique_ptr;

class AsyncRasterProcessor
{
    public:
        using update_fn = std::function<void(std::vector<node>&, bool[4])>;

        AsyncRasterProcessor() {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
       
            if (rank == 0) {
                outstanding_updates = std::unique_ptr<int[]>(new int[size]);
            }
        }

        void add(AsyncPartition* raster, update_fn fn) {
            rasters.push_back(raster);
            raster_update_fns.push_back(fn);
        }

        void run() {
            int num_borders = rasters.size() * 2;

            size_t max_border_size = 0;
            for (auto* r : rasters) {
                max_border_size = std::max(max_border_size, r->asyncGetBufferSize());
            }

            // FIXME: pick reasonable size
            // FIXME: actually, get rid of buffered sends
            size_t buffer_size = size * num_borders * (max_border_size + MPI_BSEND_OVERHEAD);
            mpi_buffer = unique_ptr<uint8_t[]>(new uint8_t[buffer_size]);
            MPI_Buffer_attach(mpi_buffer.get(), buffer_size);

            if (rank == 0) {
                printf("ASYNC: Allocating %s for MPI buffer\n", humanReadableSize(buffer_size).c_str());
                std::fill(outstanding_updates.get(), outstanding_updates.get() + size, rasters.size());
            }

            receive_requests.resize(rasters.size() * 9, MPI_REQUEST_NULL); 

            // Start receiving updates
            for (int x = 0; x < rasters.size(); x++) {
                auto* r = rasters[x];
                int req_idx = x * 9;

                for(int i = 0; i < 9; i++) {
                    int tag = ((12 + x) << 4) | i;

                    r->queueRecv(i, tag, &receive_requests[req_idx + i]);
                }
            }

            // for control messages
            if (rank == 0) {
                update_receive_buffer = unique_ptr<int[]>(new int[9 * size]);

                // Make sure there is no reallocation that invalidates MPI_Request pointers passed to Irecv
                receive_requests.reserve(receive_requests.size() + size);

                for (int i = 0; i < size; i++) {
                    receive_requests.push_back(MPI_REQUEST_NULL);
                    MPI_Irecv(&update_receive_buffer[i * 9], 9, MPI_INT, i, 100, MPI_COMM_WORLD, &receive_requests.back()); 
                }
            } else {
                receive_requests.push_back(MPI_REQUEST_NULL);
                MPI_Irecv(completion_buffer, 9, MPI_INT, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &receive_requests.back());
            }

            // Tags:
            //  100 - control message
            //  192+ - border change messages

            // Send initial border data if it is modified (or for all cases?)
            for (int x = 0; x < rasters.size(); x++) {
                auto* r = rasters[x];

                int updates[9];
                std::fill_n(updates, 9, -1);
                updates[4] = 1;

                for (int i = 0; i < 9; i++) {
                    send_border_changes(r, x, i, updates);
                }

                MPI_Bsend(updates, 9, MPI_INT, 0, 100, MCW);
            }

            vector<vector<node>> border_changes(rasters.size());
            vector<int> pending_updates(rasters.size());

            vector<int> indices(receive_requests.size());
            vector<MPI_Status> statuses(receive_requests.size());

            while (true) {
                int count = 0;
                int err = MPI_Waitsome(receive_requests.size(), &receive_requests[0], &count, &indices[0], &statuses[0]);

                if (err != MPI_SUCCESS) {
                    printf("MPI_Waitany failed on rank %d - %d\n", rank, err);
                }

                if (count > 1) {
                    //printf("rank %d: got %d messages\n", rank, count);
                }

                for (int i = 0; i < count; i++) {
                    MPI_Status& status = statuses[i];
                    int index = indices[i];

                    //printf("rank %d: got update for idx %d - %d %d\n", rank, index, status.MPI_SOURCE, status.MPI_TAG);
                    if (status.MPI_TAG == 100) {
                        int* updates = &update_receive_buffer[status.MPI_SOURCE * 9];
                        bool done = process_status(status, updates);

                        if (done) {
                            goto all_done;
                        } else {
                            MPI_Irecv(updates, 9, MPI_INT, status.MPI_SOURCE, 100, MPI_COMM_WORLD, &receive_requests[index]);
                        }
                    } else if (status.MPI_TAG >= 192) {
                        // Border update
                        total_comms++;

                        int bytes_recv;
                        MPI_Get_count(&status, MPI_BYTE, &bytes_recv);
                        comm_bytes_used += bytes_recv;

                        int raster_n = (status.MPI_TAG >> 4) - 12;
                        int neighbor = status.MPI_TAG & 15;

                        //printf("%d: Got border upd to %d %d\n", rank, raster_n, neighbor);
                        pending_updates[raster_n]++;

                        auto r = rasters[raster_n];
                        r->asyncStoreChanges(neighbor, border_changes[raster_n]);
                        r->queueRecv(neighbor, status.MPI_TAG, &receive_requests[index]);
                    } else {
                        printf("%d: unrecognized tag %d from %d\n", rank, status.MPI_TAG, status.MPI_SOURCE);
                    }
                }

                for (int raster_n = 0; raster_n < rasters.size(); raster_n++) {
                    auto r = rasters[raster_n];
                    auto& changes = border_changes[raster_n];

                    if (pending_updates[raster_n] == 0)
                        continue;

                    bool borders_updated[4] = {};
                    raster_update_fns[raster_n](changes, borders_updated);
                    changes.clear();

                    int updates[9]; 
                    std::fill_n(updates, 9, -1); 
                    updates[4] = pending_updates[raster_n];

                    if (borders_updated[0] && borders_updated[1]) {
                        send_border_changes(r, raster_n, 0, updates);
                    }

                    if (borders_updated[0]) {
                        send_border_changes(r, raster_n, 1, updates);
                    } 

                    if (borders_updated[1] && borders_updated[2]) {
                        send_border_changes(r, raster_n, 2, updates);
                    }

                    if (borders_updated[1]) {
                        send_border_changes(r, raster_n, 3, updates);
                    }

                    if (borders_updated[2]) {
                        send_border_changes(r, raster_n, 5, updates);
                    }

                    if (borders_updated[1] && borders_updated[3]) {
                        send_border_changes(r, raster_n, 6, updates);
                    }

                    if (borders_updated[3]) {
                        send_border_changes(r, raster_n, 7, updates);
                    }

                    if (borders_updated[2] && borders_updated[3]) {
                        send_border_changes(r, raster_n, 8, updates);
                    }

                    MPI_Bsend(updates, 9, MPI_INT, 0, 100, MCW);
                    pending_updates[raster_n] = 0;
                }
            }

all_done:
            for(MPI_Request& req : receive_requests) {
                if (req != MPI_REQUEST_NULL) {
                    MPI_Cancel(&req);
                }
            }

            // Detach our MPI buffer
            void* buf_addr;
            int buf_size;
            MPI_Buffer_detach(&buf_addr, &buf_size);

            int global_num_comms = 0;
            int global_bytes_used = 0;

            MPI_Reduce(&total_comms, &global_num_comms, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&comm_bytes_used, &global_bytes_used, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                printf("ASYNC: took %d border transfers (%s)\n", global_num_comms, humanReadableSize(global_bytes_used).c_str());
            }
        }

        void send_border_changes(AsyncPartition* r, int raster_n, int neighbor, int updates[]) {
            int dest_rank = r->neighborRank(neighbor);

            if (dest_rank == -1)
                return;

            auto buffer = unique_ptr<uint8_t[]>(new uint8_t[r->asyncGetBufferSize()]);
            int num_changes = r->asyncGetChanges(neighbor, buffer.get());

            if (num_changes == 0)
                return;

            //printf("rank %d: sent %d cells to %d\n", rank, num_changes, neighbor);

            updates[neighbor] = dest_rank;

            size_t nb = r->asyncBytesToSend(num_changes);
            int receivingNeighborIndex = 8 - neighbor;
            int tag = ((12 + raster_n) << 4) | receivingNeighborIndex;

            MPI_Bsend(buffer.get(), nb, MPI_BYTE, dest_rank, tag, MPI_COMM_WORLD);
        }

        bool process_status(MPI_Status& status, int updates[9]) {
            if (rank == 0) {
                // Process status notifications from other ranks
                outstanding_updates[status.MPI_SOURCE] -= updates[4];
 
                for (int i = 0; i < 9; i++) {
                    if (i == 4) continue;

                    if (updates[i] != -1) { outstanding_updates[updates[i]]++; }
                }

                bool done = true;
                for (int x = 0; x < size; x++) {
                    if (outstanding_updates[x] != 0) {
                        done = false;
                        break;
                    }
                }
                
                if (done) {
                    int unused[9] = {};

                    for (int x = 1; x < size; x++) {
                        MPI_Bsend(unused, 9, MPI_INT, x, 100, MPI_COMM_WORLD);
                    }

                    return true;
                }
            } else {
                // Root signaled global completion - we are done.
                return true;
            }

            return false;
        }

    private:
        int rank, size;
        int total_comms = 0;
        int comm_bytes_used = 0;

        unique_ptr<uint8_t[]> mpi_buffer;
        unique_ptr<int[]> outstanding_updates;
        unique_ptr<int[]> update_receive_buffer;

        vector<AsyncPartition*> rasters;
        vector<update_fn> raster_update_fns;

        vector<MPI_Request> receive_requests;

        // FIXME: move to std::array
        int completion_buffer[9];
};

#endif
