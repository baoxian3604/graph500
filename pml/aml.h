/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/
#include <pthread.h>
#include <mpi.h>
#define NBUFS 4
#define NBUFR 4
#define BUFSIZE 1024*32
typedef struct pml_thread{
	int size;
	int rank;
	void *arg;
	pthread_t *tid;
	char *buf;
	MPI_Request req[NBUFR];
	pthread_barrier_t *barrierp;
}pml_thread;
typedef struct pml_comm{
	int parse;//0,creating thread;1, all thread created
	int size;
	pthread_t *tids;
	pml_thread *pths;
	char *glob_rbuf;
	pthread_barrier_t barrier;
	pthread_mutex_t mutex;
}pml_comm;
	//MPI-like init,finalize calls
	extern int  aml_init(int *,char***);
	extern int	pml_comm_create(int size,void*(*entry)(void* p),void **arg,pml_comm *comm);
	extern int	pml_init();
	extern void	pml_wait(pml_comm *comm);
	extern void aml_finalize(void);
	//barrier which ensures that all AM sent before the barrier are completed everywhere after the barrier
	extern int aml_barrier( void );
	extern int pml_barrier( void );
	//register active message function(collective call)
	extern void aml_register_handler(void(*f)(int,void*,int),int n);
	//send AM to another(myself is ok) node
	//execution of AM might be delayed till next aml_barrier() call
	extern int aml_send(void *srcaddr, int type,int length, int node );
	extern int pml_send(void *srcaddr, int type,int length, int node );
	
	// rank and size
	extern int aml_my_pe( void );
	extern int aml_n_pes( void );
	extern int pml_long_allsum(int n);
#define my_pe aml_my_pe
#define num_pes aml_n_pes

#define aml_time() MPI_Wtime()

#define aml_long_allsum(p) MPI_Allreduce(MPI_IN_PLACE,p,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD)
#define aml_long_allmin(p) MPI_Allreduce(MPI_IN_PLACE,p,1,MPI_LONG_LONG,MPI_MIN,MPI_COMM_WORLD)
#define aml_long_allmax(p) MPI_Allreduce(MPI_IN_PLACE,p,1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD)
