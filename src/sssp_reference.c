/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 3 SSSP
// Simple parallel delta-stepping with relaxations as Active Messages

#include "aml.h"
#include "common.h"
#include "csr_reference.h"
#include "bitmap_reference.h"
#include <string.h>

#ifdef DEBUGSTATS
extern int64_t nbytes_sent,nbytes_rcvd;
#endif
// variables shared from bfs_reference
extern oned_csr_graph g;
extern int qc,q2c;
extern int* q1,*q2;
extern int* rowstarts;
extern int64_t* column,*pred_glob,visited_size;
extern unsigned long * visited;
#ifdef SSSP
//global variables as those accesed by active message handler
float *glob_dist;
float glob_maxdelta, glob_mindelta; //range for current bucket
float *weights;
volatile int lightphase;
#define TH_N 2
//Relaxation data type 
typedef struct  __attribute__((__packed__)) relaxmsg {
	float w; //weight of an edge
	int dest_vloc; //local index of destination vertex
	int src_vloc; //local index of source vertex
} relaxmsg;

extern int __thread poll_time;
extern int __thread send_time;
// Active message handler for relaxation
void relaxhndl(int from, void* dat, int sz) {
	relaxmsg* m = (relaxmsg*) dat;
	send_time++;
	//printf("recv %d->%d\n",VERTEX_TO_GLOBAL(from/2,m->src_vloc),VERTEX_TO_GLOBAL(rank,m->dest_vloc));
	int vloc = m->dest_vloc;
	float w = m->w;
	float *dest_dist = &glob_dist[vloc];
	//check if relaxation is needed: either new path is shorter or vertex not reached earlier
	if (*dest_dist < 0 || *dest_dist > w) {
		*dest_dist = w; //update distance
		//printf("from/2=%d\n",from/2);
		pred_glob[vloc]=VERTEX_TO_GLOBAL(from/TH_N,m->src_vloc); //update path

		if(lightphase && !TEST_VISITEDLOC(vloc)) //Bitmap used to track if was already relaxed with light edge
		{
			if(w < glob_maxdelta) { //if falls into current bucket needs further reprocessing
				q2[__sync_fetch_and_add(&q2c,1)] = vloc;
				SET_VISITEDLOC(vloc);
			}
		}
	}
}

//Sending relaxation active message
void send_relax(int64_t glob, float weight,int fromloc) {
	
	//printf("send %d->%d\n",VERTEX_TO_GLOBAL(rank,fromloc),glob);
	relaxmsg m = {weight,VERTEX_LOCAL(glob),fromloc};
	pml_send(&m,1,sizeof(relaxmsg),VERTEX_OWNER(glob));
}
int rt;
void *entry(void *arg){
	poll_time=0;
	send_time=0;
	//printf("1 %f\n",aml_time());
	pml_init(arg);
	unsigned int i,j;
	long sum=0;
	int localrank=((struct pml_thread*)arg)->rank%((struct pml_thread*)arg)->size;
	float delta = 0.1;
	float *dist=glob_dist;
	int64_t* pred=pred_glob;
	int root=rt;
	if(localrank==0){
		glob_mindelta=0.0;
		glob_maxdelta=delta;
		weights=g.weights;
		qc=0;q2c=0;
		if (VERTEX_OWNER(root) == my_pe()) {
			q1[0]=VERTEX_LOCAL(root);
			qc=1;
			dist[VERTEX_LOCAL(root)]=0.0;
			pred[VERTEX_LOCAL(root)]=root;
		}	
	}
	pml_barrier();
	sum=1;
	int64_t lastvisited=1;
	while(sum!=0) {
	//printf("2 %f\n",aml_time());
#ifdef DEBUGSTATS
		double t0 = aml_time();
		nbytes_sent=0;
#endif
		//1. iterate over light edges
		while(sum!=0) {
			CLEAN_VISITED();
			lightphase=1;
			pml_barrier();
			
			for(i=0;i<qc;i++)
				for(j=rowstarts[q1[i]];j<rowstarts[q1[i]+1];j++)
					if(weights[j]<delta)
						send_relax(COLUMN(j),dist[q1[i]]+weights[j],q1[i]);
			pml_barrier();
			if(localrank==0){
				qc=q2c;q2c=0;int *tmp=q1;q1=q2;q2=tmp;
				sum=qc;	
			}else sum=0;
			sum=pml_long_allsum(sum);
		}
		lightphase=0;
		pml_barrier();

	//printf("4 %f\n",aml_time());
		//2. iterate over S and heavy edges
		for(i=0;i<g.nlocalverts;i++)
			if(dist[i]>=glob_mindelta && dist[i] < glob_maxdelta) {
				for(j=rowstarts[i];j<rowstarts[i+1];j++)
					if(weights[j]>=delta)
						send_relax(COLUMN(j),dist[i]+weights[j],i);
			}
		pml_barrier();
	//printf("5 %f\n",aml_time());
	if(localrank==0){
		glob_mindelta=glob_maxdelta;
		glob_maxdelta+=delta;
		qc=0;
	}
		sum=0;
		pml_barrier();
		//3. Bucket processing and checking termination condition
		int64_t lvlvisited=0;
	if(localrank==0){
		for(i=0;i<g.nlocalverts;i++)
			if(dist[i]>=glob_mindelta) {
				sum++; //how many are still to be processed
				if (dist[i] < glob_maxdelta)
					q1[__sync_fetch_and_add(&qc,1)]=i; //this is lowest bucket
			} else if(dist[i]!=-1.0) lvlvisited++;
	}
	
		sum=pml_long_allsum(sum);
#ifdef DEBUGSTATS
		t0-=aml_time();
		lvlvisited=pml_long_allsum(lvlvisited);
		nbytes_sent=pml_long_allsum(nbytes_sent);
		if(!my_pe()) printf("--lvl[%1.2f..%1.2f] visited %lld (total %llu) in %5.2fs, network aggr %5.2fGb/s\n",glob_mindelta,glob_maxdelta,lvlvisited-lastvisited,lvlvisited,-t0,-(double)nbytes_sent*8.0/(1.e9*t0));
		lastvisited = lvlvisited;
#endif
	}
	printf("polltime=%d send_time=%d\n",poll_time,send_time);
	return NULL;
}
int __pthread_num=TH_N;
volatile extern int sync1,sync2;
void *test(void *arg){
	
	//printf("1 %f\n",aml_time());
	pml_init(arg);
	int tosend=1;
	if(rank==0)tosend=1;
	else tosend=0;
	for(int i=0;i<1138364/((struct pml_thread*)arg)->size;i++){send_relax(tosend,0.2,8);if(i%1000==0)q2c=0;}
	for(int i=0;i<88;i++)pml_barrier();
	/*{while(sync2>=__pthread_num){sleep(0);};
	int r;
	if(__pthread_num==__sync_add_and_fetch(&sync1,1)){
		sync1=0;
		MPI_Barrier(MPI_COMM_WORLD);
		r=__sync_fetch_and_add(&sync2,__pthread_num);
	}else {
		__sync_fetch_and_add(&sync2,1);
		while(sync2<__pthread_num){sleep(0);};
		r=__sync_fetch_and_sub(&sync2,1);
		if((r-1)%__pthread_num==0){
			__sync_sub_and_fetch(&sync2,__pthread_num);
		}
	}}*/
	return NULL;
}
void run_sssp(int64_t root,int64_t* pred,float *dist) {
	rt=root;
	glob_dist=dist;
	pred_glob=pred;
	aml_register_handler(relaxhndl,1);
	pml_comm comm;
	pml_comm_create(TH_N,entry,NULL,&comm);
	pml_wait(&comm);

}

void clean_shortest(float* dist) {
	int i;
	for(i=0;i<g.nlocalverts;i++) dist[i]=-1.0;
}
#endif
