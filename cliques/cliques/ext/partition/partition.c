/* Emulates partition_prune_fast.py */
/* Only works on a machine with sizeof(unsigned long) == 8 */
/* for graphs with at most 64 vertices */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define MAX_LINE_LENGTH 1000
#define SET_MAX_SIZE 63

int num_partitions = 0;

typedef struct graph{
  int n;
  int m;
  int **links;
  int *degrees;
  int *capacities;
} graph;

/* size is the max number of sets */
/* current is the index of the first free set */
typedef struct partition{
  int size;
  int current;
  unsigned long *sets;
} partition;

/* Global variables */
graph *g;
unsigned long nodes_to_place;
partition *part;
unsigned long part_union;

void report_error(char *s){
  fprintf(stderr,"%s\n",s);
  exit(-1);
}


void free_graph(graph *g){
  if (g!=NULL) {
    if (g->links!=NULL) {
      if (g->links[0]!=NULL)
        free(g->links[0]);
      free(g->links);
    }
    if (g->capacities!=NULL)
      free(g->capacities);
    if (g->degrees!=NULL)
      free(g->degrees);
    free(g);
  }
}

graph *graph_from_file(FILE *f){
  char line[MAX_LINE_LENGTH];
  int i, u, v;
  graph *g;

  assert( (g=(graph *)malloc(sizeof(graph))) != NULL );
  
  /* read n */
  if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
    report_error("graph_from_file: read error (fgets) 1");
  if( sscanf(line, "%d\n", &(g->n)) != 1 )
    report_error("graph_from_file: read error (sscanf) 2");
  
  /* read the degree sequence */
  if( (g->capacities=(int *)malloc(g->n*sizeof(int))) == NULL )
    report_error("graph_from_file: malloc() error 2");
  if( (g->degrees=(int *)calloc(g->n,sizeof(int))) == NULL )
    report_error("graph_from_file: calloc() error");
  for(i=0;i<g->n;i++){
    if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
      report_error("graph_from_file; read error (fgets) 2");
    if( sscanf(line, "%d %d\n", &v, &(g->capacities[i])) != 2 )
      report_error("graph_from_file; read error (sscanf) 2");
    if( v != i ){
      fprintf(stderr,"Line just read : %s\n i = %d; v = %d\n",line,i,v);
      report_error("graph_from_file: error while reading degrees");
    }
  }
  
  /* compute the number of links */
  g->m=0;
  for(i=0;i<g->n;i++)
    g->m += g->capacities[i];
  g->m /= 2;
  
  /* create contiguous space for links */
  if (g->n==0){
    g->links = NULL; g->degrees = NULL; g->capacities = NULL;
  }
  else {
    if( (g->links=(int **)malloc(g->n*sizeof(int*))) == NULL )
      report_error("graph_from_file: malloc() error 3");
    if( (g->links[0]=(int *)malloc(2*g->m*sizeof(int))) == NULL )
      report_error("graph_from_file: malloc() error 4");
    for(i=1;i<g->n;i++)
      g->links[i] = g->links[i-1] + g->capacities[i-1];
  }

  /* read the links */
  for(i=0;i<g->m;i++) {
    if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
      report_error("graph_from_file; read error (fgets) 3");
    if( sscanf(line, "%d %d\n", &u, &v) != 2 ){
      fprintf(stderr,"Attempt to scan link #%d failed. Line read:%s\n", i, line);
      report_error("graph_from_file; read error (sscanf) 3");
    }
    if ( (u>=g->n) || (v>=g->n) || (u<0) || (v<0) ) {
      fprintf(stderr,"Line just read: %s",line);
      report_error("graph_from_file: bad node number");
    }
    if ( (g->degrees[u]>=g->capacities[u]) ||
         (g->degrees[v]>=g->capacities[v]) ){
      fprintf(stderr, "reading link %s\n", line);
      report_error("graph_from_file: too many links for a node");
    }
    g->links[u][g->degrees[u]] = v;
    g->degrees[u]++;
    g->links[v][g->degrees[v]] = u;
    g->degrees[v]++;
  }
  for(i=0;i<g->n;i++)
    if (g->degrees[i]!=g->capacities[i])
      report_error("graph_from_file: capacities <> degrees");
  /*  printf("%s\n", line);*/
  /* horrible hack */
/*   if( fgets(line,MAX_LINE_LENGTH,f) != NULL ) */
/*     printf("%s\n", line); */
/*     report_error("graph_from_file; too many lines"); */
  
  return(g);
}

/* set functions */
unsigned long empty_set(){
  return(0);
}

/* fails if e is already in set */
unsigned long add_to_set(int e, unsigned long s){
  unsigned long tmp=1;

  tmp=tmp<<e;
  assert( (s & tmp) == 0 );
  s += tmp;
  return(s);
}

unsigned long add_to_set_no_fail(int e, unsigned long s){
  unsigned long tmp=1;

  tmp=tmp<<e;
  if( (s & tmp) == 0 )
    return(s+tmp);
  return(s);
}

/* checks whether e is in the set */
/* fails if it is not */
unsigned long remove_from_set(int e, unsigned long s){
  unsigned long tmp=1;

  tmp=tmp<<e;
  assert( (s & tmp) != 0 );
  s -= tmp;
  return(s);
}

unsigned long set_union(unsigned long s1, unsigned long s2){
  return(s1 | s2);
}

/* returns a set with elements of s1 which are not in s2 */
unsigned long set_difference(unsigned long s1, unsigned long s2){
  return(s1 - (s1 & s2));
}

/* n is the largest possible element in the set */
void print_setold(unsigned long s, int n){
  int i;
  unsigned long test;
  
  test=1;
  i=0;
  while( i < n ){
    if( (test & s) != 0 )
      printf("%d ",i);
    i+=1;
    test=test<<1;
/*      printf("i : %d, test : %ld\n", i, test);  */
  }
  if( (test & s) != 0 )
    printf("%d", i);
  printf("\n");
}

/* returns smallest element from non-empty set */
int set_get_first(unsigned long s){
  int i;
  unsigned long test;

  assert( s != 0 );
  test=1;
  i=0;
  while(1){
    if( (test & s) != 0 )
      return(i);
    i+=1;
    test=test<<1;
  }
}

/* set has length one is corresponding number is a power of two */
/* http://www.cprogramming.com/tutorial/powtwosol.html */
int has_length_one(unsigned long x)
{
  return !((x-1) & x);
}

/* returns the number of elements in the set */
/* must be smaller than SET_MAX_SIZE */
/* int cardinal(unsigned long s){ */
/*   int i; */
/*   unsigned long test; */

/*   assert( s != 0 ); */
/*   test=1; */
/*   i=0; */
/*   while(1){ */
/*     if( (test & s) != 0 ) */
/*       return(i); */
/*     i+=1; */
/*     test=test<<1; */
/*   } */
/* } */

void print_set_no_newline(unsigned long s, int n){
  int i;
  unsigned long test;
  
  test=1;
  i=0;
  while( i < n ){
    if( (test & s) != 0 )
      printf("%d ",i);
    i+=1;
    test=test<<1;
/*      printf("i : %d, test : %ld\n", i, test);  */
  }
  if( (test & s) != 0 )
    printf("%d ", i);
}

partition *create_partition(int size){
  partition *p;
  
  p=malloc(sizeof(partition));
  p->size=size;
  p->current=0;
  p->sets=calloc(size,sizeof(unsigned long));
  return(p);
}

partition *add_set(unsigned long s, partition *p){
  assert( p->current < p->size );
  p->sets[p->current]=s;
  p->current++;
  return(p);
}

partition *remove_last_set(partition *p){
  assert( p->current > 0 );
  p->sets[p->current]=empty_set();
  p->current--;
  return(p);
}

void print_partition_old(partition *p){
  int i;

  for(i=0;i<p->current;i++){
    print_set_no_newline(p->sets[i],SET_MAX_SIZE);
    //num_partitions++;
    printf("| ");
  }
  printf("\n");
}

void print_partition(partition *p){
  int i;
  num_partitions++;
  //for(i=0;i<p->current;i++)
  //  printf("%ld ", p->sets[i]);
  //printf("\n");
}

/* generate connected partitions */

void gen_partitions_prune_aux(unsigned long comp_tmp,unsigned long neighbours,unsigned long forbidden){
  int v, vn,i;
  unsigned long new_neighbours, new_forbidden, new_comp_tmp;

  if( neighbours != 0 ){
    /* possibility to extend comp_tmp with a vertex from neighbours */

    vn = set_get_first(neighbours);
    comp_tmp=add_to_set(vn,comp_tmp);

    new_neighbours=remove_from_set(vn, neighbours);
    for(i=0;i<g->degrees[vn];i++)
      new_neighbours=add_to_set_no_fail(g->links[vn][i],new_neighbours);
    new_neighbours=set_difference(new_neighbours,forbidden);
    
    forbidden=add_to_set(vn,forbidden);

    nodes_to_place=remove_from_set(vn,nodes_to_place);

    /*  Recursive call, extend comp_tmp with first neighbour" */
    gen_partitions_prune_aux(comp_tmp, new_neighbours, forbidden);

    nodes_to_place=add_to_set(vn,nodes_to_place);
    comp_tmp=remove_from_set(vn,comp_tmp);
        
    /* exploring the possibility to extend comp_tmp without using vn */
    new_neighbours=remove_from_set(vn,neighbours);
        
    /*  rec call, extend comp_tmp with other neighour */
    gen_partitions_prune_aux(comp_tmp, new_neighbours, forbidden);

    forbidden=remove_from_set(vn,forbidden);
  }
  else{
    /* all possibilities to extend comp_tmp have been considered, */
    /* new consider comp_tmp as a whole set */

    /* if comp_tmp has only one element we do not want it as a whole set */
    //if( has_length_one(comp_tmp) )
    //  return;

    if( nodes_to_place == 0 ){
      /* we have a partition */
      part=add_set(comp_tmp,part);             
      print_partition(part);
      remove_last_set(part);
      return;
    }
    /* start a new set with first vertex from nodes_to_place */

    v=set_get_first(nodes_to_place);
    nodes_to_place=remove_from_set(v,nodes_to_place);
                
    add_set(comp_tmp,part);
    part_union=set_union(comp_tmp,part_union);

    new_forbidden=add_to_set(v,part_union);
    new_comp_tmp=empty_set();
    new_comp_tmp=add_to_set(v,new_comp_tmp);

    new_neighbours=empty_set();
    for(i=0;i<g->degrees[v];i++)
      new_neighbours=add_to_set_no_fail(g->links[v][i],new_neighbours);
    new_neighbours=set_difference(new_neighbours,new_forbidden);
        
        
    /* rec call, comp_tmp is a new set */
    gen_partitions_prune_aux(new_comp_tmp, new_neighbours, new_forbidden);

    nodes_to_place=add_to_set(v,nodes_to_place);

    part=remove_last_set(part);
    part_union=set_difference(part_union,comp_tmp);
  }
}

int main(){
/*   graph *g; */
  unsigned long s,s2,comp_tmp,neighbours,forbidden;
  int i;
  partition *p;

  g=graph_from_file(stdin);
  nodes_to_place=empty_set();
  for(i=1;i<g->n;i++)
    nodes_to_place=add_to_set(i,nodes_to_place);
  part=create_partition(g->n);
  part_union=empty_set();

  comp_tmp=empty_set();
  comp_tmp=add_to_set(0,comp_tmp);
  neighbours=empty_set();
  for(i=0;i<g->degrees[0];i++)
    neighbours=add_to_set(g->links[0][i],neighbours);
  forbidden=empty_set();
  forbidden=add_to_set(0,forbidden);
  /*  printf("begin\n"); */
  gen_partitions_prune_aux(comp_tmp, neighbours, forbidden);
/*   s=empty_set(); */
/*   s=add_to_set(0,s); */
/*   s=add_to_set(2,s); */
/*   s=remove_from_set(0,s); */
/*   s=add_to_set(35,s); */
/*   s2=empty_set(); */
/*   s2=add_to_set(2,s2); */
/*   s2=add_to_set(3,s2); */
/*   s2=add_to_set(4,s2); */
/*   part=add_set(s,part); */
/*   part=add_set(s2,part); */
/*   print_partition(part); */
/*   print_set(set_difference(s,s2),35); */
/*   printf("smallest elt : %d\n", set_get_first(s)); */
/*   p=create_partition(2); */
/*   p->sets[0]=s; */
/*   print_partition(p); */
  printf("%d",num_partitions);
  free_graph(g);
  return(0);
}
