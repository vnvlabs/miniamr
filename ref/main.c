// ************************************************************************
//
// miniAMR: stencil computations with boundary exchange and AMR.
//
// Copyright (2014) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
// Questions? Contact Courtenay T. Vaughan (ctvaugh@sandia.gov)
//                    Richard F. Barrett (rfbarre@sandia.gov)
//
// ************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define MA_MAIN
#include "block.h"
#include "comm.h"
#include "timer.h"
#include "proto.h"

#include "VnV.h"

#include "param.h"

INJECTION_EXECUTABLE(MINIAMR)


int main(int argc, char** argv)
{
   int i, ierr, object_num;
   double *objs, t1;
   

   ierr = MPI_Init(&argc, &argv);
   ierr = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_pe);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_pes);
   
   INJECTION_INITIALIZE(MINIAMR, &argc, &argv, "./vnv-input.json");
   
   counter_malloc = 0;
   size_malloc = 0.0;
   num_objects = object_num = 0;

   /* set initial values */
   if (!my_pe) {
      for (i = 1; i < argc; i++)
         
         
         else if (!strcmp(argv[i], "--num_objects")) {
            num_objects = atoi(argv[++i]);
            objects = (object *) ma_malloc(num_objects*sizeof(object),
                                           __FILE__, __LINE__);
            object_num = 0;
         } else if (!strcmp(argv[i], "--object")) {
            
            if (object_num >= num_objects) {
               printf("object number greater than num_objects\n");
               MPI_Abort(MPI_COMM_WORLD, -1);
            }
            objects[object_num].type = atoi(argv[++i]);
            objects[object_num].bounce = atoi(argv[++i]);
            objects[object_num].cen[0] = atof(argv[++i]);
            objects[object_num].cen[1] = atof(argv[++i]);
            objects[object_num].cen[2] = atof(argv[++i]);
            objects[object_num].move[0] = atof(argv[++i]);
            objects[object_num].move[1] = atof(argv[++i]);
            objects[object_num].move[2] = atof(argv[++i]);
            objects[object_num].size[0] = atof(argv[++i]);
            objects[object_num].size[1] = atof(argv[++i]);
            objects[object_num].size[2] = atof(argv[++i]);
            objects[object_num].inc[0] = atof(argv[++i]);
            objects[object_num].inc[1] = atof(argv[++i]);
            objects[object_num].inc[2] = atof(argv[++i]);
            object_num++;
         

   if (object_num != num_objects) {
      printf("Error - number of objects less than specified");
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   if (reorder == -1)
     if (!lb_method)
         reorder = 1;
     else
         reorder = 0;

   if (check_input())
      MPI_Abort(MPI_COMM_WORLD, -1);

   if (!block_change)
      block_change = num_refine;

   for (object_num = 0; object_num < num_objects; object_num++) {
     for (i = 0; i < 3; i++) {
         objects[object_num].orig_cen[i] = objects[object_num].cen[i];
         objects[object_num].orig_move[i] = objects[object_num].move[i];
         objects[object_num].orig_size[i] = objects[object_num].size[i];
     }
   }

   allocate();

   if (lb_method >= 2) {
      init_x = npx*init_block_x;
      init_y = npy*init_block_y;
      init_z = npz*init_block_z;
   } else {
      // above will be done later for this case
      init_x = npx;
      init_y = npy;
      init_z = npz;
   }

   timer_main = timer() - t1;
   driver();

   INJECTION_POINT("MINIAMR", VSELF, "MY_FIRST_INJECTION_POINT", argc);


   profile();

   deallocate();

   fflush(NULL);

   MPI_Barrier(MPI_COMM_WORLD);

   INJECTION_FINALIZE("MINIAMR");
   MPI_Finalize();

   exit(0);
}

// =================================== print_help_message ====================

void print_help_message(void)
{
   printf("(Optional) command line input is of the form: \n\n");

   printf("--nx - block size x (even && > 0)\n");
   printf("--ny - block size y (even && > 0)\n");
   printf("--nz - block size z (even && > 0)\n");
   printf("--init_x - initial blocks in x (> 0)\n");
   printf("--init_y - initial blocks in y (> 0)\n");
   printf("--init_z - initial blocks in z (> 0)\n");
   printf("--reorder - ordering of blocks if initial number > 1\n");
   printf("--npx - (0 < npx <= num_pes)\n");
   printf("--npy - (0 < npy <= num_pes)\n");
   printf("--npz - (0 < npz <= num_pes)\n");
   printf("--max_blocks - maximun number of blocks per core\n");
   printf("--num_refine - (>= 0) number of levels of refinement\n");
   printf("--block_change - (>= 0) number of levels a block can change in a timestep\n");
   printf("--uniform_refine - if 1, then grid is uniformly refined\n");
   printf("--refine_freq - frequency (in timesteps) of checking for refinement\n");
   printf("--inbalance - percentage inbalance to trigger inbalance\n");
   printf("--lb_opt - load balancing - 0 = none, 1 = each refine, 2 = each refine phase\n");
   printf("--num_vars - number of variables (> 0)\n");
   printf("--comm_vars - number of vars to communicate together\n");
   printf("--num_tsteps - number of timesteps (> 0)\n");
   printf("--time - time to run problem with delta by object speed (> 0.0)\n");
   printf("--stages_per_ts - number of comm/calc stages per timestep\n");
   printf("--checksum_freq - number of stages between checksums\n");
   printf("--stencil - 0 (variable work) or 7 or 27 point (27 will not work with refinement (except uniform))\n");
   printf("--error_tol - (e^{-error_tol} ; >= 0) \n");
   printf("--report_diffusion - report check sums for each variable\n");
   printf("--report_perf - 0, 1, 2\n");
   printf("--plot_freq - frequency (timesteps) of plotting (0 for none)\n");
   printf("--code - closely minic communication of different codes\n");
   printf("         0 minimal sends, 1 send ghosts, 2 send ghosts and process on send\n");
   printf("--permute - altenates directions in communication\n");
   printf("--blocking_send - use blocking sends instead of nonblocking\n");
   printf("--refine_ghost - use full extent of block (including ghosts) to determine if block is refined\n");
   printf("--change_dir - allow the RCB algorithm to change the directions of the cuts each load balance step\n");
   printf("--group_blocks - change the RCB algorithm so that a group of blocks with the same center all get put onto the same side of a cut\n");
   printf("--limit_move - limit the number of blocks that can be moved during load balance (number that is a percentage of the total number of blocks)\n");
   printf("--send_faces - send each face individually instead of packing all faces going to a rank together\n");
   printf("--rcb - use RCB algorithm for load balancing (default)\n");
   printf("--morton - use Morton Space Filling Curve algorithm for load balancing\n");
   printf("--hilbert - use Hilbert like Space Filling Curve algorithm for load balancing\n");
   printf("--num_objects - (>= 0) number of objects to cause refinement\n");
   printf("--object - type, position, movement, size, size rate of change\n");

   printf("All associated settings are integers except for objects\n");
}

// =================================== allocate ==============================

void allocate(void)
{
   int i, j, k, m, n;

   num_blocks = (num_sz *) ma_malloc((num_refine+1)*sizeof(num_sz),
                                  __FILE__, __LINE__);
   num_blocks[0] = num_pes*init_block_x*init_block_y*init_block_z;
   local_num_blocks = (num_sz *) ma_malloc((num_refine+1)*sizeof(num_sz),
                                        __FILE__, __LINE__);
   local_num_blocks[0] = init_block_x*init_block_y*init_block_z;

   blocks = (block *) ma_malloc(max_num_blocks*sizeof(block),
                                __FILE__, __LINE__);

   for (n = 0; n < max_num_blocks; n++) {
      blocks[n].number = -1;
      blocks[n].array = (double ****) ma_malloc(num_vars*sizeof(double ***),
                                                __FILE__, __LINE__);
      for (m = 0; m < num_vars; m++) {
         blocks[n].array[m] = (double ***)
                              ma_malloc((x_block_size+2)*sizeof(double **),
                                        __FILE__, __LINE__);
         for (i = 0; i < x_block_size+2; i++) {
            blocks[n].array[m][i] = (double **)
                                   ma_malloc((y_block_size+2)*sizeof(double *),
                                             __FILE__, __LINE__);
            for (j = 0; j < y_block_size+2; j++)
               blocks[n].array[m][i][j] = (double *)
                                     ma_malloc((z_block_size+2)*sizeof(double),
                                               __FILE__, __LINE__);
         }
      }
   }

   sorted_list = (sorted_block *)ma_malloc(max_num_blocks*sizeof(sorted_block),
                                           __FILE__, __LINE__);
   sorted_index = (int *) ma_malloc((num_refine+2)*sizeof(int),
                                    __FILE__, __LINE__);

   max_num_parents = max_num_blocks;  // Guess at number needed
   parents = (parent *) ma_malloc(max_num_parents*sizeof(parent),
                                  __FILE__, __LINE__);
   for (n = 0; n < max_num_parents; n++)
      parents[n].number = -1;

   max_num_dots = 3*max_num_blocks;     // Guess at number needed
   if (!lb_method) {
      dots = (dot *) ma_malloc(max_num_dots*sizeof(dot), __FILE__, __LINE__);
      for (n = 0; n < max_num_dots; n++)
         dots[n].number = -1;
   } else {
      spots = (spot *) ma_malloc(max_num_dots*sizeof(spot), __FILE__, __LINE__);
      for (n = 0; n < max_num_dots; n++)
         spots[n].number = -1;
   }

   grid_sum = (double *)ma_malloc(num_vars*sizeof(double), __FILE__, __LINE__);

   p8 = (num_sz *) ma_malloc((num_refine+2)*sizeof(num_sz), __FILE__, __LINE__);
   p2 = (int *) ma_malloc((num_refine+2)*sizeof(int), __FILE__, __LINE__);
   block_start = (num_sz *) ma_malloc((num_refine+1)*sizeof(num_sz),
                                      __FILE__, __LINE__);

   from = (int *) ma_malloc(num_pes*sizeof(int), __FILE__, __LINE__);
   to   = (int *) ma_malloc(num_pes*sizeof(int), __FILE__, __LINE__);

   // first try at allocating comm arrays
   for (i = 0; i < 3; i++) {
      if (num_refine)
         max_comm_part[i] = 20;
      else
         max_comm_part[i] = 2;
      comm_partner[i] = (int *) ma_malloc(max_comm_part[i]*sizeof(int),
                                          __FILE__, __LINE__);
      send_size[i] = (int *) ma_malloc(max_comm_part[i]*sizeof(int),
                                       __FILE__, __LINE__);
      recv_size[i] = (int *) ma_malloc(max_comm_part[i]*sizeof(int),
                                       __FILE__, __LINE__);
      comm_index[i] = (int *) ma_malloc(max_comm_part[i]*sizeof(int),
                                        __FILE__, __LINE__);
      comm_num[i] = (int *) ma_malloc(max_comm_part[i]*sizeof(int),
                                      __FILE__, __LINE__);
      if (num_refine)
         max_num_cases[i] = 100;
      else if (i == 0)
         max_num_cases[i] = 2*init_block_y*init_block_z;
      else if (i == 1)
         max_num_cases[i] = 2*init_block_x*init_block_z;
      else
         max_num_cases[i] = 2*init_block_x*init_block_y;
      comm_block[i] = (int *) ma_malloc(max_num_cases[i]*sizeof(int),
                                        __FILE__, __LINE__);
      comm_face_case[i] = (int *) ma_malloc(max_num_cases[i]*sizeof(int),
                                            __FILE__, __LINE__);
      comm_pos[i] = (int *) ma_malloc(max_num_cases[i]*sizeof(int),
                                      __FILE__, __LINE__);
      comm_pos1[i] = (int *)ma_malloc(max_num_cases[i]*sizeof(int),
                                      __FILE__, __LINE__);
      comm_send_off[i] = (int *) ma_malloc(max_num_cases[i]*sizeof(int),
                                           __FILE__, __LINE__);
      comm_recv_off[i] = (int *) ma_malloc(max_num_cases[i]*sizeof(int),
                                           __FILE__, __LINE__);
   }

   if (num_refine) {
      par_b.max_part = 10;
      par_b.max_cases = 100;
      par_p.max_part = 10;
      par_p.max_cases = 100;
      par_p1.max_part = 10;
      par_p1.max_cases = 100;
   } else {
      par_b.max_part = 1;
      par_b.max_cases = 1;
      par_p.max_part = 1;
      par_p.max_cases = 1;
      par_p1.max_part = 1;
      par_p1.max_cases = 1;
   }
   par_b.comm_part = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                       __FILE__, __LINE__);
   par_b.comm_num = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                      __FILE__, __LINE__);
   par_b.index = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                   __FILE__, __LINE__);
   par_b.comm_b = (num_sz *) ma_malloc(par_b.max_cases*sizeof(num_sz),
                                    __FILE__, __LINE__);
   par_b.comm_p = (num_sz *) ma_malloc(par_b.max_cases*sizeof(num_sz),
                                    __FILE__, __LINE__);
   par_b.comm_c = (int *) ma_malloc(par_b.max_cases*sizeof(int),
                                    __FILE__, __LINE__);

   par_p.comm_part = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                       __FILE__, __LINE__);
   par_p.comm_num = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                      __FILE__, __LINE__);
   par_p.index = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                   __FILE__, __LINE__);
   par_p.comm_b = (num_sz *) ma_malloc(par_b.max_cases*sizeof(num_sz),
                                    __FILE__, __LINE__);
   par_p.comm_p = (num_sz *) ma_malloc(par_b.max_cases*sizeof(num_sz),
                                    __FILE__, __LINE__);
   par_p.comm_c = (int *) ma_malloc(par_b.max_cases*sizeof(int),
                                    __FILE__, __LINE__);

   par_p1.comm_part = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                       __FILE__, __LINE__);
   par_p1.comm_num = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                      __FILE__, __LINE__);
   par_p1.index = (int *) ma_malloc(par_b.max_part*sizeof(int),
                                   __FILE__, __LINE__);
   par_p1.comm_b = (num_sz *) ma_malloc(par_b.max_cases*sizeof(num_sz),
                                    __FILE__, __LINE__);
   par_p1.comm_p = (num_sz *) ma_malloc(par_b.max_cases*sizeof(num_sz),
                                    __FILE__, __LINE__);
   par_p1.comm_c = (int *) ma_malloc(par_b.max_cases*sizeof(int),
                                    __FILE__, __LINE__);

   if (num_refine) {
      s_buf_size = (int) (0.10*((double)max_num_blocks))*comm_vars*
                   (x_block_size+2)*(y_block_size+2)*(z_block_size+2);
      if (s_buf_size < (num_vars*x_block_size*y_block_size*z_block_size + 49))
         s_buf_size = num_vars*x_block_size*y_block_size*z_block_size + 49;
      r_buf_size = 5*s_buf_size;
   } else {
      i = init_block_x*(x_block_size+2);
      j = init_block_y*(y_block_size+2);
      k = init_block_z*(z_block_size+2);
      if (i > j)         // do not need ordering just two largest
         if (j > k)      // i > j > k
            s_buf_size = i*j;
         else            // i > j && k > j
            s_buf_size = i*k;
      else if (i > k)    // j > i > k
            s_buf_size = i*j;
         else            // j > i && k > i
            s_buf_size = j*k;
      r_buf_size = 2*s_buf_size;
   }
   send_buff = (double *) ma_malloc(s_buf_size*sizeof(double),
                                    __FILE__, __LINE__);
   recv_buff = (double *) ma_malloc(r_buf_size*sizeof(double),
                                    __FILE__, __LINE__);

   if (!stencil)
      a0 = (double *) ma_malloc((num_vars/4)*sizeof(double),
                                __FILE__, __LINE__);
}

// =================================== deallocate ============================

void deallocate(void)
{
   int i, j, m, n;

   for (n = 0; n < max_num_blocks; n++) {
      for (m = 0; m < num_vars; m++) {
         for (i = 0; i < x_block_size+2; i++) {
            for (j = 0; j < y_block_size+2; j++)
               free(blocks[n].array[m][i][j]);
            free(blocks[n].array[m][i]);
         }
         free(blocks[n].array[m]);
      }
      free(blocks[n].array);
   }
   free(blocks);

   free(sorted_list);
   free(sorted_index);

   free(objects);

   free(grid_sum);

   free(p8);
   free(p2);

   free(from);
   free(to);

   for (i = 0; i < 3; i++) {
      free(comm_partner[i]);
      free(send_size[i]);
      free(recv_size[i]);
      free(comm_index[i]);
      free(comm_num[i]);
      free(comm_block[i]);
      free(comm_face_case[i]);
      free(comm_pos[i]);
      free(comm_pos1[i]);
      free(comm_send_off[i]);
      free(comm_recv_off[i]);
   }

   free(send_buff);
   free(recv_buff);
}

int check_input(void)
{
   int error = 0;

   if (init_block_x < 1 || init_block_y < 1 || init_block_z < 1) {
      printf("initial blocks on processor must be positive\n");
      error = 1;
   }
   if (max_num_blocks < init_block_x*init_block_y*init_block_z) {
      printf("max_num_blocks not large enough\n");
      error = 1;
   }
   if (x_block_size < 1 || y_block_size < 1 || z_block_size < 1) {
      printf("block size must be positive\n");
      error = 1;
   }
   if (((x_block_size/2)*2) != x_block_size) {
      printf("block size in x direction must be even\n");
      error = 1;
   }
   if (((y_block_size/2)*2) != y_block_size) {
      printf("block size in y direction must be even\n");
      error = 1;
   }
   if (((z_block_size/2)*2) != z_block_size) {
      printf("block size in z direction must be even\n");
      error = 1;
   }
   if (num_refine < 0) {
      printf("number of refinement levels must be non-negative\n");
      error = 1;
   }
   if (block_change < 0) {
      printf("number of refinement levels must be non-negative\n");
      error = 1;
   }
   if (num_vars < 1) {
      printf("number of variables must be positive\n");
      error = 1;
   }
   if (num_pes != npx*npy*npz) {
      printf("number of processors used does not match number allocated\n");
      error = 1;
   }
   if (stencil != 0 && stencil != 7 && stencil != 27) {
      printf("illegal value for stencil\n");
      error = 1;
   }
   if (stencil == 0 && num_vars < 8) {
      printf("if stencil is 0, num_vars must be more than 8\n");
      error = 1;
   }
   if (stencil == 27 && num_refine && !uniform_refine)
      printf("WARNING: 27 point stencil with non-uniform refinement: answers may diverge\n");
   if (comm_vars == 0 || comm_vars > num_vars)
      comm_vars = num_vars;
   if (code < 0 || code > 2) {
      printf("code must be 0, 1, or 2\n");
      error = 1;
   }
   if (lb_opt < 0 || lb_opt > 2) {
      printf("lb_opt must be 0, 1, or 2\n");
      error = 1;
   }

   return (error);
}
