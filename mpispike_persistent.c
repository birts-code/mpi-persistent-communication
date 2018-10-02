#include <../../nrnconf.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* do not want the redef in the dynamic load case */
#include <nrnmpiuse.h>   

#if NRNMPI_DYNAMICLOAD
#include <nrnmpi_dynam.h>
#endif

#include <nrnmpi.h>
#include <hocdec.h>

#if NRNMPI
#include "nrnmpidec.h"
#include "nrnmpi_impl.h"
#include "mpispike.h"
#include <mpi.h>



typedef struct {
	int nspike;
	int gid[nrn_spikebuf_size];
	double spiketime[nrn_spikebuf_size];
} NRNMPI_Spikebuf;





extern void nrnbbs_context_wait();
extern int srplog;
extern int srpcount;

//sonra sil
extern double srpprintintov;
double srpprintbuf;
//
extern MPI_Request *srprequests;
extern MPI_Request *srprequests2;
MPI_Request srprequestarray[40];
extern int srpnout_;
extern int *srpnin_;
int srpnin_arr[2000];
extern double srpexectime;
extern double srpexectimedbl;


static int np;
static int* displs;
static int* byteovfl; /* for the compressed transfer method */
static MPI_Datatype spike_type;

static void pgvts_op(double* in, double* inout, int* len, MPI_Datatype* dptr);
static MPI_Op mpi_pgvts_op;

static void make_spike_type() {
	NRNMPI_Spike s;
	int block_lengths[2];
	MPI_Aint displacements[2];
	MPI_Aint addresses[3];
	MPI_Datatype typelist[2];

	typelist[0] = MPI_INT;
	typelist[1] = MPI_DOUBLE;

	block_lengths[0] = block_lengths[1] = 1;

	MPI_Address(&s, &addresses[0]);
	MPI_Address(&(s.gid), &addresses[1]);
	MPI_Address(&(s.spiketime), &addresses[2]);

	displacements[0] = addresses[1] - addresses[0];
	displacements[1] = addresses[2] - addresses[0];

	MPI_Type_struct(2, block_lengths, displacements, typelist, &spike_type);
	MPI_Type_commit(&spike_type);

	MPI_Op_create((MPI_User_function*)pgvts_op, 1, &mpi_pgvts_op);
}

void nrnmpi_spike_initialize() {
	make_spike_type();
}

#if nrn_spikebuf_size > 0

static MPI_Datatype spikebuf_type;

static void make_spikebuf_type() {
	NRNMPI_Spikebuf s;
	int block_lengths[3];
	MPI_Aint displacements[3];
	MPI_Aint addresses[4];
	MPI_Datatype typelist[3];

	typelist[0] = MPI_INT;
	typelist[1] = MPI_INT;
	typelist[2] = MPI_DOUBLE;

	block_lengths[0] = 1;
	block_lengths[1] = nrn_spikebuf_size;
	block_lengths[2] = nrn_spikebuf_size;

	MPI_Address(&s, &addresses[0]);
	MPI_Address(&(s.nspike), &addresses[1]);
	MPI_Address(&(s.gid[0]), &addresses[2]);
	MPI_Address(&(s.spiketime[0]), &addresses[3]);

	displacements[0] = addresses[1] - addresses[0];
	displacements[1] = addresses[2] - addresses[0];
	displacements[2] = addresses[3] - addresses[0];

	MPI_Type_struct(3, block_lengths, displacements, typelist, &spikebuf_type);
	MPI_Type_commit(&spikebuf_type);
}
#endif



int request_creator(int *srpsendbuf, int srpsendcount, MPI_Datatype srpsendtype,int srprecvbuf[],int srprecvcount, MPI_Datatype srprecvtype, MPI_Comm srpcomm,MPI_Request reqarr[]){
	
	 int dst_tree_root,my_tree_root,srpsend_offset,srprecv_offset,perststeps,srpcommsize,srpindex,srpmyrank,dst,srprecv_bytes,srpsend_bytes;
	   
       	 MPI_Comm_size( srpcomm, &srpcommsize);
	 MPI_Comm_rank( srpcomm, &srpmyrank );
	 
	 srplog = 0;
	 int srpsize2 = srpcommsize;
	 while(srpsize2){
		srpsize2 >>= 1;
		srplog++;
	 }
	 srplog --;
	 
	 int mask = 0x1;
	 int srpiter = 0;
	 
	 perststeps = 0;
	 srpindex = 0;	
     
         int srpcurr_cnt = srprecvcount;
	 
	 //MPI_Type_size(srprecvtype,&srprecv_bytes);
	 //MPI_Type_size(srpsendtype,&srpsend_bytes);
         //int *srpptr = srprecvbuf;
	 //int *srpptr2 = srpsendbuf;
	 srprecvbuf[srpmyrank] = srpsendbuf[0];
	
	 
	 while(perststeps < srplog){
		dst= srpmyrank ^ mask;

		dst_tree_root = dst >> srpiter;
                dst_tree_root <<= srpiter;
                
           	my_tree_root = srpmyrank >> srpiter;
           	my_tree_root <<= srpiter;

		srpsend_offset = my_tree_root * srprecvcount;
          	srprecv_offset = dst_tree_root * srprecvcount;

			//printf("srpmyrank %d srpiter %d sendoff %d recvoff %d\n",srpmyrank,srpiter,srpsend_offset,srprecv_offset);
			
		if(srpiter == 0){
			MPI_Send_init(&srprecvbuf[srpsend_offset],srpsendcount,srpsendtype, dst, 1,srpcomm, &reqarr[srpindex]);
			MPI_Recv_init(&srprecvbuf[srprecv_offset],srprecvcount,srprecvtype,dst,1,srpcomm,&reqarr[srpindex+1]);
			//printf("process %d step %d sends %d to process %d\n",srpmyrank,perststeps,srprecvbuf[srpsend_offset],dst);
			//printf("process %d srpiter %d receives to nin_[%d] from process %d\n", srpmyrank,srpiter,srprecv_offset,dst);
		}
		else{
			MPI_Send_init(&srprecvbuf[srpsend_offset],srpcurr_cnt,srprecvtype, dst, 1,srpcomm, &reqarr[srpindex]);
			//printf("process %d srpiter %d sends %d to process %d\n",srpmyrank,srpiter,srprecvbuf[srpsend_offset],dst);
    			MPI_Recv_init(&srprecvbuf[srprecv_offset],srpcurr_cnt,srprecvtype,dst,1,srpcomm,&reqarr[srpindex+1]);
			//printf("process %d srpiter %d receives to rbuf1[%d] from process %d\n", srpmyrank,srpiter,srprecv_offset,dst);
		}
			
				//int j;
				//for(j=0;j<MESSAGESIZE*srpcommsize;j++)
               // printf("process %d srpiter %d rbuf1[%d] %lf \n",srpmyrank,srpiter,j,srprecvbuff[j]);
        
			perststeps++;
			mask <<= 1;
			srpcurr_cnt += srpcurr_cnt;
			srpiter++;
			srpindex += 2;
	}


	return 0;	
}


MPI_Request *request_creator2(const void *srpsendbuf, int srpsendcount, MPI_Datatype srpsendtype,void *srprecvbuf,int srprecvcount, MPI_Datatype srprecvtype, MPI_Comm srpcomm){
	
	 int dst_tree_root,my_tree_root,send_offset,recv_offset,perststeps,srpcommsize,srpindex,srpmyrank,dst,recv_bytes,send_bytes;
	 
	 MPI_Comm_size( srpcomm, &srpcommsize);
	 MPI_Comm_rank( srpcomm, &srpmyrank );
	 
	 srplog = 0;
	 int srpsize2 = srpcommsize;
	 while(srpsize2){
		srpsize2 >>= 1;
		srplog++;
	 }
	 srplog --;
	 
	 int srpmask = 0x1;
	 int srpiter = 0;
	 
	 perststeps = 0;
	 srpindex = 0;	
     
         int srpcurr_cnt = srprecvcount;
	 MPI_Request *srpreqs2 = (MPI_Request*)malloc(sizeof(MPI_Request)*srplog*2);

	 
	 while(perststeps < srplog){
			dst= srpmyrank ^ srpmask;

			dst_tree_root = dst >> srpiter;
            dst_tree_root <<= srpiter;
                
           	my_tree_root = srpmyrank >> srpiter;
           	my_tree_root <<= srpiter;

		send_offset = my_tree_root * srprecvcount;
          	recv_offset = dst_tree_root * srprecvcount;
			
			if(srpiter == 0){
				((NRNMPI_Spikebuf*)srprecvbuf)[srpmyrank*srprecvcount].nspike = ((NRNMPI_Spikebuf*)srpsendbuf)->nspike;
				int t;
					
				for(t = 0; t < nrn_spikebuf_size; t++){
					((NRNMPI_Spikebuf*)srprecvbuf)[srpmyrank*srprecvcount].gid[t] =((NRNMPI_Spikebuf*)srpsendbuf)->gid[t];
					((NRNMPI_Spikebuf*)srprecvbuf)[srpmyrank*srprecvcount].spiketime[t] = ((NRNMPI_Spikebuf*)srpsendbuf)->spiketime[t];
				}

		
                               MPI_Send_init(&((NRNMPI_Spikebuf*)srprecvbuf)[send_offset],srpsendcount,srpsendtype,dst, 1,srpcomm, &srpreqs2[srpindex]);
				MPI_Recv_init(&((NRNMPI_Spikebuf*)srprecvbuf)[recv_offset],srprecvcount,srprecvtype,dst,1,srpcomm,&srpreqs2[srpindex+1]);	

				
			}
			else{
				MPI_Send_init(&((NRNMPI_Spikebuf*)srprecvbuf)[send_offset],srpcurr_cnt,srprecvtype, dst, 1,srpcomm, &srpreqs2[srpindex]);
    			MPI_Recv_init(&((NRNMPI_Spikebuf*)srprecvbuf)[recv_offset],srpcurr_cnt,srprecvtype,dst,1,srpcomm,&srpreqs2[srpindex+1]);
			}
        
			perststeps++;
			srpmask <<= 1;
			srpcurr_cnt += srpcurr_cnt;
			srpiter++;
			srpindex += 2;
	}


	return srpreqs2;	
}





int request_initializer(MPI_Request *srpreqs, MPI_Comm srpcomm){
     int perststeps,srpindex;
     /*
	 MPI_Comm_size( srpcomm, &srpcommsize);
	 
	 srplog = 0;
	 srpsize2 = srpcommsize;
	 while(srpsize2){
		srpsize2 >>= 1;
		srplog++;
	 }
	 srplog--;
     */

     // optimization issues first send
     MPI_Start(&srpreqs[0]);

     // issues MPI_recv 
     perststeps = 0;
     srpindex = 0;	
     while(perststeps<srplog){
       MPI_Start(&srpreqs[srpindex+1]);
       perststeps++;
       srpindex += 2;
     }


	
	 perststeps = 0;
	 srpindex = 0;	
	 while(perststeps<srplog){
	   if(srpindex != 0) MPI_Start(&srpreqs[srpindex]);
		//MPI_Start(&srpreqs[srpindex+1]);
				
		MPI_Wait(&srpreqs[srpindex+1],MPI_STATUS_IGNORE);

#if 0
		int i, flag=0;		
		for(i=0;i<600;i++){
		  MPI_Test(&srpreqs[srpindex+1], &flag, MPI_STATUS_IGNORE);
		  if(flag) break;
		}		
		if(!flag) MPI_Wait(&srpreqs[srpindex+1],MPI_STATUS_IGNORE);
#endif

		perststeps++;
		srpindex += 2;
				
	}

	/*perststeps = 0;
	srpindex = 0;	
	while(perststeps<srplog){
	 	MPI_Request_free(&srpreqs[srpindex]);
		MPI_Request_free(&srpreqs[srpindex+1]);	
		srpindex += 2;	
		perststeps++;
	}*/
	return 0;
}

#if 0
int request_initializer(MPI_Request *srpreqs, MPI_Comm srpcomm){
     int srpcommsize,srpsize2,perststeps,srpindex;
 

	 perststeps = 0;
	 srpindex = 0;	
	 while(perststeps<srplog){
		MPI_Start(&srpreqs[srpindex+1]);
		//	MPI_Start(&srpreqs[srpindex]);
		if(perststeps == 0 ) MPI_Start(&srpreqs[0]);
				
		//MPI_Wait(&srpreqs[srpindex+1],MPI_STATUS_IGNORE);
		perststeps++;
		srpindex += 2;
				
	}
         int i, flag=0;
	 perststeps = 0;
	 srpindex = 0;	
	 while(perststeps<srplog){
	   //MPI_Start(&srpreqs[srpindex+1]);
	   if(perststeps != 0 ) MPI_Start(&srpreqs[srpindex]);

		
	   //	MPI_Wait(&srpreqs[srpindex+1],MPI_STATUS_IGNORE);
	   //#if 0

		flag=0;		
		for(i=0;i<60000000;i++){
		  MPI_Test(&srpreqs[srpindex+1], &flag, MPI_STATUS_IGNORE);
		  if(flag) break;
		}		
		if(!flag) MPI_Wait(&srpreqs[srpindex+1],MPI_STATUS_IGNORE);
		//#endif
		perststeps++;
		srpindex += 2;
				
	}

	return 0;
}
#endif

int nrnmpi_spike_exchange() {
	int i, n, novfl, n1;
	if (!displs) {
		np = nrnmpi_numprocs;
		displs = (int*)hoc_Emalloc(np*sizeof(int)); hoc_malchk();
		displs[0] = 0;
#if nrn_spikebuf_size > 0		
		make_spikebuf_type();
#endif
	}
	nrnbbs_context_wait();
#if nrn_spikebuf_size == 0
	int srpsize,srprank;
	
	MPI_Comm_size( nrnmpi_comm, &srpsize );
   	MPI_Comm_rank( nrnmpi_comm, &srprank );

	////MPI_Allgather(&nout_, 1, MPI_INT, nin_, 1, MPI_INT, nrnmpi_comm);
	
	//srpnout_ = nout_;
	if(srpprintintov < 2){
			fprintf(stderr,"First allgather\n");
			srpprintintov++;
	}
	if(srpcount == 0){
		//srpnin_ = (int*)malloc(sizeof(int)*nrnmpi_numprocs);
		//request_creator(&nout_,1,MPI_INT,srpnin_arr,1,MPI_INT,nrnmpi_comm,srprequestarray);
	  request_creator(&nout_,1,MPI_INT,nin_,1,MPI_INT,nrnmpi_comm,srprequestarray);
	  // srpcount++;
	  srpexectime = 0;
	}
	srpnin_arr[srprank] = nout_;
	double srptime = MPI_Wtime();
	request_initializer(srprequestarray, nrnmpi_comm);
	srptime = MPI_Wtime() - srptime;
	 srpcount++;
	srpexectime += srptime;
#if 0
        int srpindex;
	for(srpindex = 0; srpindex < nrnmpi_numprocs;srpindex++){
	  nin_[srpindex] = srpnin_arr[srpindex];	
	}
#endif	
// printf("proc %d nout was %d iter  %d \n",srprank,nout_,srpcount);

	/*for(srpindex = 0; srpindex < nrnmpi_numprocs;srpindex++){
		printf("proc %d iter %d nin_[%d] %d \n",srprank,srpcount,srpindex,nin_[srpindex]);
	}*/

	n = nin_[0];
	for (i=1; i < np; ++i) {
		displs[i] = n;
		n += nin_[i];
	}
	if (n) {
		if (icapacity_ < n) {
			icapacity_ = n + 10;
			free(spikein_);
			spikein_ = (NRNMPI_Spike*)hoc_Emalloc(icapacity_ * sizeof(NRNMPI_Spike)); hoc_malchk();
		}
		if(srpprintintov < 2){
			fprintf(stderr,"First allgather overflow\n");
			srpprintintov++;
		}
		MPI_Allgatherv(spikeout_, nout_, spike_type, spikein_, nin_, displs, spike_type, nrnmpi_comm);
	}
	//printf("came to the else section \n");
#else

	int srpsize,srprank,srpi,srpj;
	
	MPI_Comm_size( nrnmpi_comm, &srpsize );
   	MPI_Comm_rank( nrnmpi_comm, &srprank );

	///////MPI_Allgather(spbufout_, 1, spikebuf_type, spbufin_, 1, spikebuf_type, nrnmpi_comm); idi
	if(srpcount2 == 0){
		srpspbufin_ = (NRNMPI_Spikebuf*)malloc(sizeof(NRNMPI_Spikebuf)*nrnmpi_numprocs);
		srprequests2 = request_creator2(&spbufout_[0],1,spikebuf_type,&srpspbufin_[0],1,spikebuf_type,nrnmpi_comm);
		srpcount2++;
		srpexectimedbl = 0.0;
	}
	double srptimedbl = MPI_Wtime();
	request_initializer(srprequests2, nrnmpi_comm);
	srptimedbl = MPI_Wtime() - srptimedbl;
	srpexectimedbl += srptimedbl;
	for(srpi=0;srpi<nrnmpi_numprocs;srpi++){
		spbufin_[srpi].nspike = srpspbufin_[srpi].nspike;
		for(srpj=0;srpj<nrn_spikebuf_size;srpj++){
			spbufin_[srpi].gid[srpj] = srpspbufin_[srpi].gid[srpj];
			spbufin_[srpi].spiketime[srpj] = srpspbufin_[srpi].spiketime[srpj];
		}
	}
	printf("inside the second allgather nrn_spikebuf_size %d\n",nrn_spikebuf_size);
	novfl = 0;
	n = spbufin_[0].nspike;
	if (n > nrn_spikebuf_size) {
		nin_[0] = n - nrn_spikebuf_size;
		novfl += nin_[0];
	}else{
		nin_[0] = 0;
	}
	for (i=1; i < np; ++i) {
		displs[i] = novfl;
		n1 = spbufin_[i].nspike;
		n += n1;
		if (n1 > nrn_spikebuf_size) {
			nin_[i] = n1 - nrn_spikebuf_size;
			novfl += nin_[i];
		}else{
			nin_[i] = 0;
		}
	}
	if (novfl) {
		if (icapacity_ < novfl) {
			icapacity_ = novfl + 10;
			free(spikein_);
			spikein_ = (NRNMPI_Spike*)hoc_Emalloc(icapacity_ * sizeof(NRNMPI_Spike)); hoc_malchk();
		}
		n1 = (nout_ > nrn_spikebuf_size) ? nout_ - nrn_spikebuf_size : 0;
		fprintf(stderr,"inside the second allgather nrn_spikebuf_size overflow %d\n",nrn_spikebuf_size);
		MPI_Allgatherv(spikeout_, n1, spike_type, spikein_, nin_, displs, spike_type, nrnmpi_comm);
	}
	ovfl_ = novfl;
#endif
	return n;
}




/*
The compressed spike format is restricted to the fixed step method and is
a sequence of unsigned char. 
nspike = buf[0]*256 + buf[1]
a sequence of spiketime, localgid pairs. There are nspike of them.
	spiketime is relative to the last transfer time in units of dt.
	note that this requires a mindelay < 256*dt.
	localgid is an unsigned int, unsigned short,
	or unsigned char in size depending on the range and thus takes
	4, 2, or 1 byte respectively. To be machine independent we do our
	own byte coding. When the localgid range is smaller than the true
	gid range, the gid->PreSyn are remapped into
	hostid specific	maps. If there are not many holes, i.e just about every
	spike from a source machine is delivered to some cell on a
	target machine, then instead of	a hash map, a vector is used.
The allgather sends the first part of the buf and the allgatherv buffer
sends any overflow.
*/
int nrnmpi_spike_exchange_compressed() {
	fprintf(stderr,"in compressed exchange\n");
	int i, novfl, n, ntot, idx, bs, bstot; /* n is #spikes, bs is #byte overflow */
	if (!displs) {
		np = nrnmpi_numprocs;
		displs = (int*)hoc_Emalloc(np*sizeof(int)); hoc_malchk();
		displs[0] = 0;
		byteovfl = (int*)hoc_Emalloc(np*sizeof(int)); hoc_malchk();
	}
	nrnbbs_context_wait();

	MPI_Allgather(spfixout_, ag_send_size_, MPI_BYTE, spfixin_, ag_send_size_, MPI_BYTE, nrnmpi_comm);
	novfl = 0;
	ntot = 0;
	bstot = 0;
	for (i=0; i < np; ++i) {
		displs[i] = bstot;
		idx = i*ag_send_size_;
		n = spfixin_[idx++]*256;
		n += spfixin_[idx++];
		ntot += n;
		nin_[i] = n;
		if (n > ag_send_nspike_) {
			bs = 2 + n*(1 + localgid_size_) - ag_send_size_;
			byteovfl[i] = bs;
			bstot += bs;
			novfl += n - ag_send_nspike_;
		}else{
			byteovfl[i] = 0;
		}
	}
	if (novfl) {
		if (ovfl_capacity_ < novfl) {
			ovfl_capacity_ = novfl + 10;
			free(spfixin_ovfl_);
			spfixin_ovfl_ = (unsigned char*)hoc_Emalloc(ovfl_capacity_ * (1 + localgid_size_)*sizeof(unsigned char)); hoc_malchk();
		}
		bs = byteovfl[nrnmpi_myid];
		/*
		note that the spfixout_ buffer is one since the overflow
		is contiguous to the first part. But the spfixin_ovfl_ is
		completely separate from the spfixin_ since the latter
		dynamically changes its size during a run.
		*/
		MPI_Allgatherv(spfixout_ + ag_send_size_, bs, MPI_BYTE, spfixin_ovfl_, byteovfl, displs, MPI_BYTE, nrnmpi_comm);
	}
	ovfl_ = novfl;
	return ntot;
}

double nrnmpi_mindelay(double m) {
	double result;
	if (!nrnmpi_use) { return m; }
	nrnbbs_context_wait();
	MPI_Allreduce(&m, &result, 1, MPI_DOUBLE, MPI_MIN, nrnmpi_comm);
	return result;
}

int nrnmpi_int_allmax(int x) {
	int result;
	if (nrnmpi_numprocs < 2) { return x; }
	nrnbbs_context_wait();
	MPI_Allreduce(&x, &result, 1, MPI_INT, MPI_MAX, nrnmpi_comm);
	return result;
}

extern void nrnmpi_int_gather(int* s, int* r, int cnt, int root) {
	fprintf(stderr,"in nrnmpi_int_gather\n");
	MPI_Gather(s, cnt, MPI_INT, r, cnt, MPI_INT, root, nrnmpi_comm);
}

extern void nrnmpi_int_gatherv(int* s, int scnt,
    int* r, int* rcnt, int* rdispl, int root) {
	fprintf(stderr,"in nrnmpi_int_gatherv\n");
	MPI_Gatherv(s, scnt, MPI_INT,
		r, rcnt, rdispl, MPI_INT, root, nrnmpi_comm);
}

extern void nrnmpi_int_alltoallv(int* s, int* scnt, int* sdispl,
    int* r, int* rcnt, int* rdispl) {
	fprintf(stderr,"in nrnmpi_int_alltoallv\n");
	MPI_Alltoallv(s, scnt, sdispl, MPI_INT,
		r, rcnt, rdispl, MPI_INT, nrnmpi_comm);
}

extern void nrnmpi_dbl_alltoallv(double* s, int* scnt, int* sdispl,
    double* r, int* rcnt, int* rdispl) {
	fprintf(stderr,"in nrnmpi_dbl_alltoallv\n");
	MPI_Alltoallv(s, scnt, sdispl, MPI_DOUBLE,
		r, rcnt, rdispl, MPI_DOUBLE, nrnmpi_comm);
}

extern void nrnmpi_char_alltoallv(char* s, int* scnt, int* sdispl,
    char* r, int* rcnt, int* rdispl) {
	fprintf(stderr,"in nrnmpi_char_alltoallv\n");
	MPI_Alltoallv(s, scnt, sdispl, MPI_CHAR,
		r, rcnt, rdispl, MPI_CHAR, nrnmpi_comm);
}

/* following are for the partrans */

void nrnmpi_int_allgather(int* s, int* r, int n) {
	fprintf(stderr,"in nrnmpi_int_allgather\n");
	MPI_Allgather(s, n,  MPI_INT, r, n, MPI_INT, nrnmpi_comm);
}

void nrnmpi_int_allgatherv(int* s, int* r, int* n, int* dspl) {
	fprintf(stderr,"in nrnmpi_int_allgatherv\n");
	MPI_Allgatherv(s, n[nrnmpi_myid],  MPI_INT,
		r, n, dspl, MPI_INT, nrnmpi_comm);
}

void nrnmpi_dbl_allgatherv(double* s, double* r, int* n, int* dspl) {
	fprintf(stderr,"in nrnmpi_dbl_allgatherv\n");
	MPI_Allgatherv(s, n[nrnmpi_myid],  MPI_DOUBLE,
		r, n, dspl, MPI_DOUBLE, nrnmpi_comm);
}

void nrnmpi_dbl_broadcast(double* buf, int cnt, int root) {
	fprintf(stderr,"in nrnmpi_dbl_broadcast\n");
	MPI_Bcast(buf, cnt,  MPI_DOUBLE, root, nrnmpi_comm);
}

void nrnmpi_int_broadcast(int* buf, int cnt, int root) {
	fprintf(stderr,"in nrnmpi_int_broadcast\n");
	MPI_Bcast(buf, cnt,  MPI_INT, root, nrnmpi_comm);
}

void nrnmpi_char_broadcast(char* buf, int cnt, int root) {
	fprintf(stderr,"in nrnmpi_char_broadcast\n");
	MPI_Bcast(buf, cnt,  MPI_CHAR, root, nrnmpi_comm);
}

int nrnmpi_int_sum_reduce(int in) {
	int result;
	MPI_Allreduce(&in, &result, 1, MPI_INT, MPI_SUM, nrnmpi_comm);
	return result;
}

void nrnmpi_assert_opstep(int opstep, double t) {
	/* all machines in comm should have same opstep and same t. */
	double buf[2];
	if (nrnmpi_numprocs < 2) { return; }
	buf[0] = (double)opstep;
	buf[1] = t;
	MPI_Bcast(buf, 2, MPI_DOUBLE, 0, nrnmpi_comm);
	if (opstep != (int)buf[0]  || t != buf[1]) {
		printf("%d opstep=%d %d  t=%g t-troot=%g\n", nrnmpi_myid, opstep,
			(int)buf[0], t, t-buf[1]);
		hoc_execerror("nrnmpi_assert_opstep failed", (char*)0);		
	}
}

double nrnmpi_dbl_allmin(double x) {
	double result;
	if (nrnmpi_numprocs < 2) { return x; }
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_MIN, nrnmpi_comm);
	return result;
}

static void pgvts_op(double* in, double* inout, int* len, MPI_Datatype* dptr){
	int i, r=0;
	assert(*dptr == MPI_DOUBLE);
	assert(*len == 4);
	if (in[0] < inout[0]) {
 		/* least time has highest priority */
 		r = 1;
	}else if (in[0] == inout[0]) {
		/* when times are equal then */
		if (in[1] < inout[1]) {
			/* NetParEvent done last */
			r = 1;
		}else if (in[1] == inout[1]) {
			/* when times and ops are equal then */
			if (in[2] < inout[2]) {
				/* init done next to last.*/
				r = 1;
			}else if (in[2] == inout[2]) {
				/* when times, ops, and inits are equal then */
				if (in[3] < inout[3]) {
					/* choose lowest rank */
					r = 1;
				}
			}
		}
	}
	if (r) {
		for (i=0; i < 4; ++i) { inout[i] = in[i]; }	
	}
}

int nrnmpi_pgvts_least(double* t, int* op, int* init) {
	int i;
	double ibuf[4], obuf[4];
	ibuf[0] = *t;
	ibuf[1] = (double)(*op);
	ibuf[2] = (double)(*init);
	ibuf[3] = (double)nrnmpi_myid;
	for (i=0; i < 4; ++i) {
		obuf[i] = ibuf[i];
	}
	MPI_Allreduce(ibuf, obuf, 4, MPI_DOUBLE, mpi_pgvts_op, nrnmpi_comm);
	assert(obuf[0] <= *t);
	if (obuf[0] == *t) {
	  assert((int)obuf[1] <= *op);
	  if ((int)obuf[1] == *op) {
	    assert((int)obuf[2] <= *init);
	    if ((int)obuf[2] == *init) {
	      assert((int)obuf[3] <= nrnmpi_myid);
	    }
	  }
	}
	*t = obuf[0];
	*op = (int)obuf[1];
	*init = (int)obuf[2];
	if (nrnmpi_myid == (int)obuf[3]) {
		return 1;
	}
	return 0;
}

/* following for splitcell.cpp transfer */
void nrnmpi_send_doubles(double* pd, int cnt, int dest, int tag) {
	fprintf(stderr,"splitcell nrnmpi_send_doubles\n");
	MPI_Send(pd, cnt, MPI_DOUBLE, dest, tag, nrnmpi_comm);
}

void nrnmpi_recv_doubles(double* pd, int cnt, int src, int tag) {
	MPI_Status status;
	fprintf(stderr,"splitcell nrnmpi_recv_doubles\n");
	MPI_Recv(pd, cnt, MPI_DOUBLE, src, tag, nrnmpi_comm, &status);
}

void nrnmpi_postrecv_doubles(double* pd, int cnt, int src, int tag, void** request) {
	fprintf(stderr,"splitcell nrnmpi_postrecv_doubles\n");
	MPI_Irecv(pd, cnt, MPI_DOUBLE, src, tag, nrnmpi_comm, (MPI_Request*)request);
}

void nrnmpi_wait(void** request) {
	MPI_Status status;
	MPI_Wait((MPI_Request*)request, &status);
}

void nrnmpi_barrier() {
	if (nrnmpi_numprocs < 2) { return; }
	MPI_Barrier(nrnmpi_comm);
}

double nrnmpi_dbl_allreduce(double x, int type) {
	double result;
	MPI_Op t;
	if (nrnmpi_numprocs < 2) { return x; }
	if (type == 1) {
		t = MPI_SUM;
	}else if (type == 2) {
		t = MPI_MAX;
	}else{
		t = MPI_MIN;
	}
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, t, nrnmpi_comm);
	return result;
}

void nrnmpi_dbl_allreduce_vec(double* src, double* dest, int cnt, int type) {
	int i;
	MPI_Op t;
	assert(src != dest);
	if (nrnmpi_numprocs < 2) {
		for (i = 0; i < cnt; ++i) {
			dest[i] = src[i];
		}
		return;
	}
	if (type == 1) {
		t = MPI_SUM;
	}else if (type == 2) {
		t = MPI_MAX;
	}else{
		t = MPI_MIN;
	}
	MPI_Allreduce(src, dest, cnt, MPI_DOUBLE, t, nrnmpi_comm);
	return;
}

void nrnmpi_longdbl_allreduce_vec(longdbl* src, longdbl* dest, int cnt, int type) {
	int i;
	MPI_Op t;
	assert(src != dest);
	if (nrnmpi_numprocs < 2) {
		for (i = 0; i < cnt; ++i) {
			dest[i] = src[i];
		}
		return;
	}
	if (type == 1) {
		t = MPI_SUM;
	}else if (type == 2) {
		t = MPI_MAX;
	}else{
		t = MPI_MIN;
	}
	MPI_Allreduce(src, dest, cnt, MPI_LONG_DOUBLE, t, nrnmpi_comm);
	return;
}

void nrnmpi_long_allreduce_vec(long* src, long* dest, int cnt, int type) {
	int i;
	MPI_Op t;
	assert(src != dest);
	if (nrnmpi_numprocs < 2) {
		for (i = 0; i < cnt; ++i) {
			dest[i] = src[i];
		}
		return;
	}
	if (type == 1) {
		t = MPI_SUM;
	}else if (type == 2) {
		t = MPI_MAX;
	}else{
		t = MPI_MIN;
	}
	MPI_Allreduce(src, dest, cnt, MPI_LONG, t, nrnmpi_comm);
	return;
}

void nrnmpi_dbl_allgather(double* s, double* r, int n) {
	fprintf(stderr,"in nrnmpi_dbl_allgather\n");
	MPI_Allgather(s, n,  MPI_DOUBLE, r, n, MPI_DOUBLE, nrnmpi_comm);
}

#if BGPDMA

static MPI_Comm bgp_comm;

void nrnmpi_bgp_comm() {
	if (!bgp_comm) {
		MPI_Comm_dup(MPI_COMM_WORLD, &bgp_comm);
	}
}

void nrnmpi_bgp_multisend(NRNMPI_Spike* spk, int n, int* hosts) {
	int i;
	MPI_Request r;
	MPI_Status status;
	fprintf(stderr,"in nrnmpi_bgp_multisend\n");
	for (i=0; i < n; ++i) {
		MPI_Isend(spk, 1, spike_type, hosts[i], 1, bgp_comm, &r);
		MPI_Request_free(&r);
	}
}

int nrnmpi_bgp_single_advance(NRNMPI_Spike* spk) {
	int flag = 0;
	MPI_Status status;
	MPI_Iprobe(MPI_ANY_SOURCE, 1, bgp_comm, &flag, &status);
	if (flag) {
		MPI_Recv(spk, 1, spike_type, MPI_ANY_SOURCE, 1, bgp_comm, &status);
	}
	return flag;
}

static int iii;
int nrnmpi_bgp_conserve(int nsend, int nrecv) {
	int tcnts[2];
	tcnts[0] = nsend - nrecv;
	MPI_Allreduce(tcnts, tcnts+1, 1, MPI_INT, MPI_SUM, bgp_comm);
	return tcnts[1];
}

#endif /*BGPDMA*/

#endif /*NRNMPI*/
