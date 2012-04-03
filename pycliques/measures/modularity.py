def evaluate_all_partitions(graph,partitions):
    scores = []
    for partition in partitions:
        scores.append(evaluate_partition(graph,partition))
    
    return scores

#Evaluate quality measure
def evaluate_partition(graph,partition):
    m = graph.number_of_edges()
    Q = float(0)
    for part in partition:
        for i in part:
            for j in part:
                Q = Q + float(graph.A(i,j)) - (float(graph.k(i))*float(graph.k(j)))/(float(2*m))

    Q = Q/(2*m)
    return Q
