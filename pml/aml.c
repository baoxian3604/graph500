/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// AML: active messages library v 1.0
// MPI-3 passive transport
// transparent message aggregation greatly increases message-rate for loosy interconnects
// shared memory optimization used
// Implementation basic v1.0
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <unistd.h>
#include <mpi.h>
#include <pthread.h>

#define MAXGROUPS 65536		//number of nodes (core processes form a group on a same node)
#define AGGR (1024*32) //aggregation buffer size per dest in bytes : internode
#define AGGR_intra (1024*32) //aggregation buffer size per dest in bytes : intranode
#define NRECV 4 // number of preposted recvs internode
#define NRECV_intra 4 // number of preposted recvs intranode
#define NSEND 4 // number of available sends internode
#define NSEND_intra 4 // number of send intranode
#define SOATTR __attribute__((visibility("default")))

#define SENDSOURCE(node) ( sendbuf+(AGGR*nbuf[node]))
#define SENDSOURCE_intra(node) ( sendbuf_intra+(AGGR_intra*nbuf_intra[node]) )

#define ushort unsigned short
static int myproc,num_procs;
static int mygroup,num_groups;
static int mylocal,group_size;
#ifndef PROCS_PER_NODE_NOT_POWER_OF_TWO
static int loggroup,groupmask;
#define PROC_FROM_GROUPLOCAL(g,l) ((l)+((g)<<loggroup))
#define GROUP_FROM_PROC(p) ((p) >> loggroup)
#define LOCAL_FROM_PROC(p) ((p) & groupmask)
#else
#define PROC_FROM_GROUPLOCAL(g,l) ((g)*group_size+(l))
#define GROUP_FROM_PROC(p) ((p)/group_size)
#define LOCAL_FROM_PROC(p) ((p)%group_size)
#endif
volatile static int ack=0;

volatile static int inbarrier=0;

static void (*aml_handlers[256]) (int,void *,int); //pointers to user-provided AM handlers

//internode comm (proc number X from each group)
//intranode comm (all cores of one nodegroup)
MPI_Comm comm, comm_intra;

// MPI stuff for sends
static char *sendbuf; //coalescing buffers, most of memory is allocated is here
static int *sendsize; //buffer occupacy in bytes
static ushort *acks; //aggregated acks
static ushort *nbuf; //actual buffer for each group/localcore
static ushort activebuf[NSEND];// N_buffer used in transfer(0..NSEND{_intra}-1)
static MPI_Request rqsend[NSEND];
// MPI stuff for recv
static char recvbuf[AGGR*NRECV];
static MPI_Request rqrecv[NRECV];

unsigned long long nbytes_sent,nbytes_rcvd;

static char *sendbuf_intra;
static int *sendsize_intra;
static ushort *acks_intra;
static ushort *nbuf_intra;
static ushort activebuf_intra[NSEND_intra];
static MPI_Request rqsend_intra[NSEND_intra];
static char recvbuf_intra[AGGR_intra*NRECV_intra];
static MPI_Request rqrecv_intra[NRECV_intra];
volatile static int ack_intra=0;
inline void aml_send_intra(void *srcaddr, int type, int length, int local ,int from);

void aml_finalize(void);
void aml_barrier(void);

SOATTR void aml_register_handler(void(*f)(int,void*,int),int n) { aml_barrier(); aml_handlers[n]=f; aml_barrier(); }

struct __attribute__((__packed__)) hdr { //header of internode message
	ushort sz;
	char hndl;
	char routing;
};
//process internode messages
static void process(int fromgroup,int length ,char* message) {
	int i = 0;
	int from = PROC_FROM_GROUPLOCAL(fromgroup,mylocal);
	while ( i < length ) {
		void* m = message+i;
		struct hdr *h = (struct hdr *)m;
		int hsz=h->sz;
		int hndl=h->hndl;
		int destlocal = LOCAL_FROM_PROC(h->routing);
		if(destlocal == mylocal)
			aml_handlers[hndl](from,m+sizeof(struct hdr),hsz);
		else
			aml_send_intra(m+sizeof(struct hdr),hndl,hsz,destlocal,from);
		i += hsz + sizeof(struct hdr);
	}
}
struct __attribute__((__packed__)) hdri { //header of internode message
	ushort routing;
	ushort sz;
	char hndl;
};

//process intranode messages
static void process_intra(int fromlocal,int length ,char* message) {
	int i=0;
	while ( i < length ) {
		void*m = message+i;
		struct hdri *h = (struct hdri *)m;
		int hsz=h->sz;
		int hndl=h->hndl;
		aml_handlers[hndl](PROC_FROM_GROUPLOCAL((int)(h->routing),fromlocal),m+sizeof(struct hdri),hsz);
		i += sizeof(struct hdri) + hsz;
	}
}

// poll intranode message
const int one=1;
inline void aml_poll_intra(void) {
	int flag, from, length,index;
	MPI_Status status;
	MPI_Testany( NRECV_intra,rqrecv_intra, &index, &flag, &status );
	if ( flag ) {
		MPI_Get_count( &status, MPI_CHAR, &length );
		if(length>0){
			ack_intra -= *(int *)(recvbuf_intra +AGGR_intra*index);
			length-=sizeof(int);
		}
		if(length>0) { //no confirmation & processing for ack only messages
			from = status.MPI_SOURCE;
			if(inbarrier){
				MPI_Send(&one, 4, MPI_CHAR,from, 1, comm_intra); //ack now
			}
			else
			{
				acks_intra[from]++; //normally we have delayed ack
			}
			process_intra( from, length,recvbuf_intra +AGGR_intra*index+sizeof(int));
		}
		MPI_Start( rqrecv_intra+index);
	}
}
// poll internode message
static void aml_poll(void) {
	int flag, from, length,index;
	MPI_Status status;

	aml_poll_intra();

	MPI_Testany( NRECV,rqrecv,&index, &flag, &status );
	if ( flag ) {
		MPI_Get_count( &status, MPI_CHAR, &length );
		nbytes_rcvd+=length;
		if(length>0){
			ack -= *(int *)(recvbuf +AGGR_intra*index);
			length-=sizeof(int);
		}
		if(length>0) { //no confirmation & processing for ack only messages
			from = status.MPI_SOURCE;
			if(inbarrier)
				MPI_Send(&one, 1, MPI_INT,from, 1, comm); //ack now
			else
				acks[from]++; //normally we have delayed ack
			process( from, length,recvbuf+AGGR*index+sizeof(int) );
		}
		MPI_Start( rqrecv+index );
	}
}

//flush internode buffer to destination node
inline void flush_buffer( int node ) {
	MPI_Status stsend;
	int flag=0,index,tmp;
	if (sendsize[node] == 0 && acks[node]==0 ) return;
	while (!flag) {
		aml_poll();
		MPI_Testany(NSEND,rqsend,&index,&flag,&stsend);
	}
	((int *)(SENDSOURCE(node)))[0]=acks[node];
	MPI_Isend(SENDSOURCE(node), sendsize[node], MPI_CHAR,node, node, comm, rqsend+index );
	nbytes_sent+=sendsize[node];
	if (sendsize_intra[node] > sizeof(int)) ack++;
	sendsize[node] = sizeof(int);
	acks[node] = 0;
	tmp=activebuf[index]; activebuf[index]=nbuf[node]; nbuf[node]=tmp; //swap bufs

}
//flush intranode buffer, NB:node is local number of pe in group
inline void flush_buffer_intra( int node ) {
	MPI_Status stsend;
	int flag=0,index,tmp;
	if (sendsize_intra[node] == 0 && acks_intra[node]==0 ) return;
	while (!flag) {
		aml_poll_intra();
		MPI_Testany(NSEND_intra,rqsend_intra,&index,&flag,&stsend);
	}
	((int *)(SENDSOURCE_intra(node)))[0]=acks_intra[node];
	MPI_Isend( SENDSOURCE_intra(node), sendsize_intra[node], MPI_CHAR,
			node, node, comm_intra, rqsend_intra+index );
	if (sendsize_intra[node] > sizeof(int)) ack_intra++;
	sendsize_intra[node] = sizeof(int);
	acks_intra[node] = 0;
	tmp=activebuf_intra[index]; activebuf_intra[index]=nbuf_intra[node]; nbuf_intra[node]=tmp; //swap bufs

}

inline void aml_send_intra(void *src, int type, int length, int local, int from) {
	//send to _another_ process from same group
	int nmax = AGGR_intra - sendsize_intra[local] - sizeof(struct hdri);
	if ( nmax < length ) {
		flush_buffer_intra(local);
	}
	char* dst = (SENDSOURCE_intra(local)+sendsize_intra[local]);
	struct hdri *h=(struct hdri *)dst;
	h->routing = GROUP_FROM_PROC(from);
	h->sz=length;
	h->hndl = type;
	sendsize_intra[local] += length+sizeof(struct hdri);

	memcpy(dst+sizeof(struct hdri),src,length);
}
SOATTR void aml_send(void *src, int type,int length, int node ) {
	if ( node == myproc )
		return aml_handlers[type](myproc,src,length);

	int group = GROUP_FROM_PROC(node);
	int local = LOCAL_FROM_PROC(node);

	//send to another node in my group
	if ( group == mygroup )
		return aml_send_intra(src,type,length,local,myproc);

	//send to another group
	int nmax = AGGR - sendsize[group]-sizeof(struct hdr);
	if ( nmax < length ) {
		flush_buffer(group);
	}
	char* dst = (SENDSOURCE(group)+sendsize[group]);
	struct hdr *h=(struct hdr*)dst;
	h->routing = local;
	h->hndl = type;
	h->sz=length;
	sendsize[group] += length+sizeof(struct hdr);
	memcpy(dst+sizeof(struct hdr),src,length);
}


int stringCmp( const void *a, const void *b)
{ return strcmp((const char*)a,(const char*)b);  }

// Should be called by user instead of MPI_Init()
SOATTR int aml_init( int *argc, char ***argv ) {
	int r, i, j,tmpmax;

	int provided=-1;
	r=MPI_Init_thread(argc,argv,MPI_THREAD_MULTIPLE,&provided);
	if(provided!=MPI_THREAD_MULTIPLE){
		printf("pml init failed\n");
	}
	if ( r != MPI_SUCCESS ) return r;

	MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myproc );

	//split communicator
	char host_name[MPI_MAX_PROCESSOR_NAME];
	char (*host_names)[MPI_MAX_PROCESSOR_NAME];
	int namelen,bytes,n,color;
	MPI_Get_processor_name(host_name,&namelen);

	bytes = num_procs * sizeof(char[MPI_MAX_PROCESSOR_NAME]);
	host_names = (char (*)[MPI_MAX_PROCESSOR_NAME]) malloc(bytes);
	strcpy(host_names[myproc], host_name);
	for (n=0; n<num_procs; n++)
		MPI_Bcast(&(host_names[n]),MPI_MAX_PROCESSOR_NAME, MPI_CHAR, n, MPI_COMM_WORLD);
	qsort(host_names, num_procs, sizeof(char[MPI_MAX_PROCESSOR_NAME]), stringCmp);
	color = 0;
	for (n=0; n<num_procs; n++)  {
		if(n>0 && strcmp(host_names[n-1], host_names[n])) color++;
		if(strcmp(host_name, host_names[n]) == 0) break;
	}
	
	free(host_names);
	MPI_Comm_split(MPI_COMM_WORLD, color, myproc, &comm_intra);

	//find intranode numbers and make internode communicator
	MPI_Comm_size( comm_intra, &group_size );
	MPI_Comm_rank( comm_intra, &mylocal );

	MPI_Comm_split(MPI_COMM_WORLD, mylocal, myproc, &comm);

	MPI_Comm_size( comm, &num_groups );
	MPI_Comm_rank( comm, &mygroup );

	//first nonblocking barriers are blocking,so we call them now
	MPI_Request hndl;
	MPI_Ibarrier(comm,&hndl);
	MPI_Wait(&hndl,MPI_STATUS_IGNORE);
	MPI_Ibarrier(comm_intra,&hndl);
	MPI_Wait(&hndl,MPI_STATUS_IGNORE);

#ifndef PROCS_PER_NODE_NOT_POWER_OF_TWO
	groupmask=group_size-1;
	if((group_size&groupmask)) { printf("AML: Fatal: non power2 groupsize unsupported. Define macro PROCS_PER_NODE_NOT_POWER_OF_TWO to override\n");return -1;}
	for (loggroup = 0; loggroup < group_size; loggroup++)
		if ((1 << loggroup) == group_size) break;
#endif
	if(myproc!=PROC_FROM_GROUPLOCAL(mygroup,mylocal)) {printf("AML: Fatal: Strange group rank assignment scheme.\n");return -1;}
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);

	CPU_SET(mylocal,&cpuset); //FIXME ? would it work good enough on all architectures?
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#ifdef DEBUGSTATS
	if(myproc==0) printf ("AML: multicore, num_groups %d group_size %d\n",num_groups,group_size);
#ifdef PROCS_PER_NODE_NOT_POWER_OF_TWO
	if(myproc==0) printf ("AML: multicore, PROCS_PER_NODE_NOT_POWER_OF_TWO defined\n");
#else
	if(myproc==0) printf ("AML: multicore, loggroup=%d groupmask=%d\n",loggroup,groupmask);
#endif
	if(myproc==0) printf ("NRECV=%d NRECVi=%d NSEND=%d  NSENDi=%d AGGR=%dK AGGRi=%dK\n",NRECV,NRECV_intra,NSEND,NSEND_intra,AGGR>>10,AGGR_intra>>10);
#endif
	if(num_groups>MAXGROUPS) { if(myproc==0) printf("AML:v1.0 reference:unsupported num_groups > MAXGROUPS=%d\n",MAXGROUPS); exit(-1); }
	fflush(NULL);
	//init preposted recvs: NRECV internode
	for(i=0;i<NRECV;i++)  {
		r = MPI_Recv_init( recvbuf+AGGR*i, AGGR, MPI_CHAR,MPI_ANY_SOURCE, MPI_ANY_TAG, comm,rqrecv+i );
		if ( r != MPI_SUCCESS ) return r;
	}
	sendbuf = (char *)malloc( AGGR*(num_groups+NSEND));
	if ( !sendbuf ) return -1;
	memset(sendbuf,0,AGGR*(num_groups+NSEND));
	sendsize = (int *)malloc( num_groups*sizeof(*sendsize) );
	if (!sendsize) return -1;
	acks = (short unsigned int*)malloc( num_groups*sizeof(*acks) );
	if (!acks) return -1;
	nbuf = (short unsigned int*)malloc( num_groups*sizeof(*nbuf) );
	if (!nbuf) return -1;


	for(i=0;i<NRECV_intra;i++)  {
		r = MPI_Recv_init( recvbuf_intra+AGGR_intra*i, AGGR_intra, MPI_CHAR,MPI_ANY_SOURCE, MPI_ANY_TAG, comm_intra,rqrecv_intra+i );
		if ( r != MPI_SUCCESS ) return r;
	}
	sendbuf_intra = (char *)malloc( AGGR_intra*(group_size+NSEND_intra));
	if ( !sendbuf_intra ) return -1;
	memset(sendbuf_intra,0,AGGR_intra*(group_size+NSEND_intra));
	sendsize_intra = (int *)malloc( group_size*sizeof(*sendsize_intra) );
	if (!sendsize_intra) return -1;
	acks_intra = (short unsigned int*)malloc( group_size*sizeof(*acks_intra) );
	if (!acks_intra) return -1;
	nbuf_intra = (short unsigned int*)malloc( group_size*sizeof(*nbuf_intra) );
	if (!nbuf_intra) return -1;
	for ( j = 0; j < group_size; j++ ) {
		sendsize_intra[j] = sizeof(int); nbuf_intra[j] = j; acks_intra[j]=0;
	}
	for(i=0;i<NRECV_intra;i++)
		MPI_Start(rqrecv_intra+i);

	for ( j = 0; j < NSEND_intra; j++ ) {
		MPI_Isend( NULL, 0, MPI_CHAR, MPI_PROC_NULL, 0, comm_intra, rqsend_intra+j );
		activebuf_intra[j]=group_size+j;
	}

	for ( j = 0; j < num_groups; j++ ) {
		sendsize[j] = sizeof(int); nbuf[j] = j;  acks[j]=0;
	}
	for(i=0;i<NRECV;i++)
		MPI_Start( rqrecv+i );
	for ( j = 0; j < NSEND; j++ ) {
		MPI_Isend( NULL, 0, MPI_CHAR, MPI_PROC_NULL, 0, comm, rqsend+j );
		activebuf[j]=num_groups+j;
	}
	return 0;
}
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
#define BUF_OF(buf,i) (((char *)(buf))+BUFSIZE*i)
typedef struct pml_comm{
	int parse;//0,creating thread;1, all thread created
	int size;
	pthread_t *tids;
	pml_thread *pths;
	char *glob_rbuf;
	pthread_barrier_t barrier;
	pthread_mutex_t mutex;
}pml_comm;
pml_comm glob_comm;
MPI_Comm pml_comm_intra;
SOATTR int pml_comm_create(int size,void*(*entry)(void* p),void **arg,pml_comm *comm){
	glob_comm.size=size;
	glob_comm.parse=0;
	pthread_mutex_init(&glob_comm.mutex, NULL); 
	glob_comm.tids=(pthread_t *)malloc(sizeof(pthread_t)*size);
	glob_comm.pths=(pml_thread *)malloc(sizeof(pml_thread)*size);
	glob_comm.glob_rbuf=(char *)malloc(BUFSIZE*size*NBUFR);
	pml_thread *pths=glob_comm.pths;
	MPI_Comm_split(comm_intra, 0, myproc, &pml_comm_intra);
	pthread_attr_t attr;
	
	pthread_barrier_init(&glob_comm.barrier,NULL,size);
	int r;
	cpu_set_t cpuset;
	for(int i=0;i<size;i++){
		pths[i].size=size;
		pths[i].rank=size*myproc+i;
		if(arg)pths[i].arg=arg[i];
		else pths[i].arg=NULL;
		pths[i].tid=&(glob_comm.tids[i]);
		pths[i].buf=BUF_OF(glob_comm.glob_rbuf,i*NBUFR);
		pths[i].barrierp=&(glob_comm.barrier);
		for(int j=0;j<NBUFR;j++){
			r = MPI_Recv_init( BUF_OF(pths[i].buf,j), BUFSIZE, MPI_CHAR,MPI_ANY_SOURCE, pths[i].rank, pml_comm_intra,pths[i].req+j );
			MPI_Start(pths[i].req+j);
			if ( r != MPI_SUCCESS ) return r;
		}
		pthread_attr_init(&attr);
		CPU_ZERO(&cpuset);
		CPU_SET(pths[i].rank,&cpuset);
		pthread_attr_setaffinity_np(&attr, sizeof(cpuset), &cpuset);
		pthread_create(&glob_comm.tids[i],&attr,entry,&pths[i]);
		
	}
	
	glob_comm.parse=1;
	if(comm)*comm=glob_comm;
	return 0;
}

static __thread int thread_rank;
static __thread int self;
SOATTR int pml_init(pml_thread *pth){
	self=pth->rank%pth->size;
	thread_rank=pth->rank;
	return thread_rank;
}
SOATTR void pml_wait(pml_comm *comm){
	glob_comm.parse=2;
	for(int i=0;i<comm->size;i++){
		pthread_join(comm->tids[i],NULL);
	}
	glob_comm.parse=0;
}
static void pml_process(int from,int length ,char* message) {
	int i=0;
	while ( i < length ) {
		void*m = message+i;
		struct hdri *h = (struct hdri *)m;
		int hsz=h->sz;
		int hndl=h->hndl;
		aml_handlers[hndl](h->routing,m+sizeof(struct hdri),hsz);
		i += sizeof(struct hdri) + hsz;
	}
}
//int poll_time=0;
#ifndef p2p
inline void pml_poll(void) {//recev any tag
	//poll_time++;
	int flag, from, length,index;
	static pthread_mutex_t mutex2=PTHREAD_MUTEX_INITIALIZER;
	MPI_Status status;	
	//printf("???? %d %d\n",ack_intra,myproc);
	//printf("!!!!\n");
	if(pthread_mutex_trylock(&mutex2))return;
	MPI_Testany( NRECV_intra,rqrecv_intra, &index, &flag, &status );
	if ( flag ) {
		char *bufhead=recvbuf_intra +AGGR_intra*index;
		MPI_Get_count( &status, MPI_CHAR, &length );
		if(length>0){
			__sync_fetch_and_sub(&ack_intra,*(int *)bufhead);
			length-=sizeof(int);
		}
		if(length>0) { 
			from = status.MPI_SOURCE;
			if(inbarrier){
				MPI_Send(&one, 4, MPI_CHAR,from, ((struct hdri*)(bufhead+sizeof(int)))->routing, comm_intra); //ack now
			}
			else
			{
				__sync_fetch_and_add(&acks_intra[from],1); //normally we have delayed ack
			}
			pml_process( 0, length,bufhead+sizeof(int));
		}
		MPI_Start( rqrecv_intra+index);
	}
		pthread_mutex_unlock(&mutex2);
}
#else
inline void pml_poll(void) {
	int flag, from, length,index;
	MPI_Status status;
	MPI_Testany( NBUFR,glob_comm.pths[self].req, &index, &flag, &status );
	if ( flag ) {
		char *bufhead=BUF_OF(glob_comm.pths[self].buf,index);
		MPI_Get_count( &status, MPI_CHAR, &length );
		if(length>0){
			__sync_fetch_and_sub(&ack_intra,*(int *)bufhead);
			length-=sizeof(int);
			//if(*(int *)bufhead)printf("recv ack=%d from %d\n",*(int *)BUF_OF(glob_comm.pths[self].buf,index),status.MPI_SOURCE);
		}
		if(length>0) { //no confirmation & processing for ack only messages
		////	printf("recv msg in %d from %d inbarrier?%d\n",thread_rank,((struct hdri*)(bufhead+sizeof(int)))->routing,inbarrier);
			from = status.MPI_SOURCE;
			if(inbarrier){
				MPI_Send(&one, 4, MPI_CHAR,from, ((struct hdri*)(bufhead+sizeof(int)))->routing, pml_comm_intra); //ack now
			}
			else
			{
				__sync_fetch_and_add(&acks_intra[from],1); //normally we have delayed ack
			}
			pml_process( 0, length,bufhead+sizeof(int));
		}
		MPI_Start( glob_comm.pths[self].req+index);
	}
}
#endif
inline void pml_flush_buffer( int group ) {
	MPI_Status stsend;
	int flag=0,index,tmp;
	if (sendsize_intra[group] == 0 && acks_intra[group]==0 ) return;
	while (!flag) {
		pml_poll();
		MPI_Testany(NSEND_intra,rqsend_intra,&index,&flag,&stsend);
	}
	((int *)(SENDSOURCE_intra(group)))[0]=__sync_lock_test_and_set(&acks_intra[group],0);
	
#ifdef p2p
	MPI_Isend( SENDSOURCE_intra(group), sendsize_intra[group], MPI_CHAR,group, group*glob_comm.size, pml_comm_intra, rqsend_intra+index );
#else
	MPI_Isend( SENDSOURCE_intra(group), sendsize_intra[group], MPI_CHAR,group, group*glob_comm.size, comm_intra, rqsend_intra+index );
	//if (sendsize_intra[group] > sizeof(int))printf("myproc=%d %d\n",myproc,group);
#endif
	if (sendsize_intra[group] > sizeof(int)) __sync_fetch_and_add(&ack_intra,1);
	sendsize_intra[group] = sizeof(int);
	tmp=activebuf_intra[index]; activebuf_intra[index]=nbuf_intra[group]; nbuf_intra[group]=tmp; 

}
int writing=0;
//int send_time=0;
SOATTR void pml_send(void *src, int type,int length, int node ) {
	//send_time++;
#ifdef p2p
	if ( node == thread_rank )
		return aml_handlers[type](thread_rank,src,length);
#else 
	if ( node == thread_rank/glob_comm.size )
		return aml_handlers[type](thread_rank,src,length);
#endif
	int group = node;
	while(pthread_mutex_trylock(&glob_comm.mutex)){sleep(0);}
	int nmax = AGGR_intra - sendsize_intra[group] - sizeof(struct hdri);
	if ( nmax < length ) {
		while(writing){sleep(0);}
		pml_flush_buffer(group);
	}
	char* dst = (SENDSOURCE_intra(group)+sendsize_intra[group]);
	__sync_fetch_and_add(&writing,1);
	__sync_fetch_and_add(&sendsize_intra[group],length+sizeof(struct hdri));
	pthread_mutex_unlock(&glob_comm.mutex);
	struct hdri *h=(struct hdri *)dst;
	h->routing = thread_rank;
	//printf("thread_rank=%d\n",thread_rank);
	h->sz=length;
	h->hndl = type;
	memcpy(dst+sizeof(struct hdri),src,length);
	__sync_fetch_and_sub(&writing,1);
	
}
volatile int sync1=0,sync2=0;
SOATTR int pml_barrier( void ) {
	//printf("!!!ack_intra=%d %d %d\n",ack_intra,myproc,num_procs);
	int __pthread_num=glob_comm.size;
	while(sync2>=__pthread_num){
		//pml_poll();
		sleep(0);
		};
	int r;
	if(__pthread_num==__sync_add_and_fetch(&sync1,1)){
		sync1=0;
			
		int i,flag;
		MPI_Request hndl;
		inbarrier++;
		for ( i = 0; i < num_procs; i++ ) {
			pml_flush_buffer(i);
		}
		while(ack_intra!=0){ 
		pml_poll();}
		MPI_Ibarrier(comm_intra,&hndl);
		flag=0;
		while(flag==0) {
			MPI_Test(&hndl,&flag,MPI_STATUS_IGNORE); pml_poll(); }

		inbarrier--;
		MPI_Barrier(MPI_COMM_WORLD);
		if(__pthread_num==1)return 0;
		r=__sync_fetch_and_add(&sync2,__pthread_num);
	}else {
		__sync_fetch_and_add(&sync2,1);
		while(sync2<__pthread_num){
			sleep(0);
			//pml_poll();
		};
		r=__sync_fetch_and_sub(&sync2,1);
		if((r-1)%__pthread_num==0){
			__sync_sub_and_fetch(&sync2,__pthread_num);
		}
	}
	return r;
}
volatile int __sync3=0;
int __sum=0;
SOATTR int pml_long_allsum(int n){
	__sync_fetch_and_add(&__sum,n);
	pthread_barrier_wait(&glob_comm.barrier);
	if(pthread_mutex_trylock(&glob_comm.mutex)==0){
		MPI_Allreduce(MPI_IN_PLACE,&__sum,1,MPI_INT,MPI_SUM,pml_comm_intra);
		pthread_barrier_wait(&glob_comm.barrier);
		pthread_mutex_unlock(&glob_comm.mutex);
	}else 
	pthread_barrier_wait(&glob_comm.barrier);
	int r=__sum;
	if(__sync_add_and_fetch(&__sync3,1)==glob_comm.size){
		__sync_sub_and_fetch(&__sum,r);
		__sync3=0;
	}
	return r;
}
SOATTR void aml_barrier( void ) {
	int i,flag;
	MPI_Request hndl;
	inbarrier++;
	//1. flush internode buffers
	for ( i = 1; i < num_groups; i++ ) {
		int group=(mygroup+i)%num_groups;
		flush_buffer(group);
	}
	//2. wait for all internode being acknowledged
	while(ack!=0) aml_poll();
	//3. notify everybody that all my internode messages were received
	MPI_Ibarrier(comm,&hndl);
	//4. receive internode until barrier done
	flag=0;
	while(flag==0) {
		MPI_Test(&hndl,&flag,MPI_STATUS_IGNORE); aml_poll(); }
	// NB: All internode received here. I can receive some more intranode.

	//5. Flush all intranode buffers
	for ( i = 1; i < group_size; i++ ) {
		int localproc=LOCAL_FROM_PROC(mylocal+i);
		flush_buffer_intra(localproc);
	}
	//inbarrier=2;
	//6. wait for all intranode being acknowledged
	while(ack_intra!=0) aml_poll_intra();
	//7. notify everybody that all my intranode messages were received
	MPI_Ibarrier(comm_intra,&hndl);
	//8. receive internode until barrier done
	flag=0;
	while(flag==0) {
		MPI_Test(&hndl,&flag,MPI_STATUS_IGNORE); aml_poll_intra(); }
	inbarrier--;
	MPI_Barrier(MPI_COMM_WORLD);
}

SOATTR void aml_finalize( void ) {
	int i;
	aml_barrier();
	for(i=0;i<NRECV;i++)
		MPI_Cancel(rqrecv+i);
#ifndef NOINTRA
	for(i=0;i<NRECV_intra;i++)
		MPI_Cancel(rqrecv_intra+i);
#endif
	MPI_Finalize();
}

SOATTR int aml_my_pe(void) { return myproc; }
SOATTR int aml_n_pes(void) { return num_procs; }
