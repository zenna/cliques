def evaluate_all_partitions(graph,partitions,t):
    scores = []
    for partition in partitions:
        scores.append(evaluate_partition(graph,partition,t))
    
    return scores

#Evaluate quality measure
def evaluate_partition(graph,partition,t):
    m = graph.number_of_edges()
    first_term = 0.0
    second_term = 0.0
    
    for part in partition:
        for i in part:
            for j in part:
                first_term = first_term + (float(graph.k(i))*float(graph.k(j)))
                second_term = second_term + float(graph.A(i,j))

    two_m = 2.0*float(m)
    first_term = first_term/(two_m*two_m)
    second_term = second_term /two_m

    t = float(t)
    R = (1-t) - first_term + t * second_term
    return R

import math
def combinatorial_stability(graph,partition,time):
    N = float(len(graph))
    
    mean_k = 0.0
    for k in graph.degree().values():
        mean_k = mean_k = k
    mean_k = mean_k / len(graph)
    
    R = 0.0
    for part in partition:
        for i in part:
            for j in part:
                if i == j:
                    laplace = float(graph.degree(i))
                else:
                    laplace = 0.0              
                R = R + math.exp(time*laplace/mean_k) * 1.0/N - 1.0/(N*N)

    return R
