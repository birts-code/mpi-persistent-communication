#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h> 
#define NITER 200000

main(int argc, char *argv[]) {
	int commsize,myrank,i;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    	MPI_Comm adjdist;
	
	
	int indegree = commsize;
	int outdegree = commsize;
	int sources[indegree],destinations[outdegree],index,reorder,sendcount,recvcount,sourceweights[indegree],destweights[outdegree];
	reorder = index = 0;
	for(i=0;i<indegree;i++){

			sources[index] = i;
			destinations[index] = i;
			sourceweights[i] =0;
			destweights[i] =0;
			index++;
		
	}


	MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,indegree,sources,sourceweights,outdegree,destinations,destweights,MPI_INFO_NULL,reorder,&adjdist);

	//MPI_Dist_graph_neighbors(dist,indegree,sources,MPI_UNWEIGHTED,outdegree,destinations,MPI_UNWEIGHTED);

	int sendbuf = myrank;
	int sendbuf2 = myrank;
	sendcount = recvcount = 1;
	double neightime,normaltime;
	int *recvbuf,*recvbuf2;

	recvbuf = malloc(commsize * sizeof(int));
	recvbuf2 = malloc(commsize * sizeof(int));

	for(i=0;i<commsize*recvcount;i++){
		recvbuf[i] = recvbuf2[i] = -11;
	}
	
	neightime = MPI_Wtime();
	for(i=0;i<NITER;i++){
		MPI_Neighbor_allgather(&sendbuf, sendcount, MPI_INT, &recvbuf[0], recvcount, MPI_INT,adjdist);
	}
	neightime = MPI_Wtime() - neightime;

	normaltime = MPI_Wtime();
	for(i=0;i<NITER;i++){
		MPI_Allgather(&sendbuf2, sendcount, MPI_INT, &recvbuf2[0], recvcount, MPI_INT,MPI_COMM_WORLD);
	}
	normaltime = MPI_Wtime() - normaltime;

	/*for(i=0;i<commsize;i++)
		printf("rank %d %d %d\n",myrank,recvbuf[i],recvbuf2[i]);*/

	int isequal = 1;
	for(i=0;i<commsize*recvcount;i++){
		if(recvbuf[i] != recvbuf2[i]){
			isequal = 0;
			break;			
		}
		
	}
	if(isequal)
	printf("process %d neightime %lf normaltime %lf\n",myrank,neightime,normaltime);
	else{	
		printf("process %d wrong allgather\n",myrank);
		for(i=0;i<commsize*recvcount;i++){
			printf("process %d recvbuf[i] = %d recvbuf2[i] = %d\n",myrank,recvbuf[i],recvbuf2[i]);
		}
	}

	MPI_Finalize();
	return;
}
