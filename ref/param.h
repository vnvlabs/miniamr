// ************************************************************************
//
// miniAMR: stencil computations with boundary exchange and AMR.
//
// Copyright (2014) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
//
// This library is free software) you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation) either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY) without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library) if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
// Questions? Contact Courtenay T. Vaughan (ctvaugh@sandia.gov)
//                    Richard F. Barrett (rfbarre@sandia.gov)
//
// ************************************************************************

#define PARAMSLIST \
X(max_num_blocks , integer,  500, , "description goes here" )\
X(num_refine , integer ,   5, "description goes here" )\
X(uniform_refine , integer ,   0, "description goes here" )\
X(x_block_size , integer ,   10, "description goes here" )\
X(y_block_size , integer ,   10, "description goes here" )\
X(z_block_size , integer ,   10, "description goes here" )\
X(num_vars , integer ,   40, "description goes here" )\
X(comm_vars , integer ,   0, "description goes here" )\
X(init_block_x , integer ,   1, "description goes here" )\
X(init_block_y , integer ,   1, "description goes here" )\
X(init_block_z , integer ,   1, "description goes here" )\
X(reorder , integer ,   -1, "description goes here" )\
X(npx , integer ,   1, "description goes here" )\
X(npy , integer ,   1, "description goes here" )\
X(npz , integer ,   1, "description goes here" )\
X(inbalance , integer ,   0, "description goes here" )\
X(refine_freq , integer ,   5, "description goes here" )\
X(report_diffusion , integer ,   0, "description goes here" )\
X(error_tol , integer ,   8, "description goes here" )\
X(use_tsteps , integer ,   0, "description goes here" )\
X(num_tsteps , integer ,   20, "description goes here" )\
X(use_time , integer ,   0, "description goes here" )\
X(end_time , integer ,   0.0, "description goes here" )\
X(stages_per_ts , integer ,   20, "description goes here" )\
X(checksum_freq , integer ,   5, "description goes here" )\
X(stencil , integer ,   7, "description goes here" )\
X(report_perf , integer ,   12, "description goes here" )\
X(plot_freq , integer ,   0, "description goes here" )\
X(num_objects , integer ,   0, "description goes here" )\
X(lb_opt , integer ,   1, "description goes here" )\
X(block_change , integer ,   0, "description goes here" )\
X(code , integer ,   0, "description goes here" )\
X(permute , integer ,   0, "description goes here" )\
X(nonblocking , integer ,   1, "description goes here" )\
X(refine_ghost , integer ,   0, "description goes here" )\
X(change_dir , integer ,   0, "description goes here" )\
X(group_blocks , integer ,   0, "description goes here" )\
X(limit_move , integer ,   0, "description goes here" )\
X(send_faces , integer ,   0, "description goes here" )\
X(lb_method , integer ,   0, "description goes here" )

#define X(name,type,default, desc) type name = default;
struct _Params {
    PARAMSLIST; 
}
typedef _Params Params;
#undef X

#define X(name,type,default,desc) "\"" #name "\" : {\"type\" : \"" #type "\" , \"description\" : \""#desc" \" , \"default\" :" #default "}," }

const char* param_schema() {
    return "{ \"type\" : \"object\", \"properties\" : { " PARAMSLIST "\"objects\" : { \"type\" : \"object\" } } }";
}
#undef X


#define X(type,name, def) \
  if(cjson_object_get(json,#name,&result) && cjson_##type(result,&##type##_)) { \
      params->##name = type##_; \
  }


/**
 * The options for MINI AMR
 * ------------------------
 */
INJECTION_OPTIONS(MINIAMR, ""){
   Params* p = (Params*) malloc(sizeof(Params));
   int integer_;
   cjson result;
   PARAMSLIST;


   if (cjson_object_contains(json,"num-tsteps")) p->use_tsteps = 1;
   if (cjson_object_contains(json,"time")) p->use_time = 1;
   
   if (cjson_object_contains(json,"group_blocks")) p->group_blocks = 1;
   else if (cjson_object_contains(json,"break_ties")) p->group_blocks = 2;
   else {p->group_blocks = 0;}
   
   if (cjson_object_contains(json,"rcb")) p->lb_method = 0;
   else if (cjson_object_contains(json,"morton")) p->lb_method = 1;
   else if (cjson_object_contains(json,"hilbert")) p->lb_method = 2;
   else if (cjson_object_contains(json,"trunc_hilbert")) p->lb_method = 3;
   
   return (void*) p;

  


}
