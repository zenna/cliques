import math
def variation_of_information(partitions,n):
    #TODO calculate n
    all_var_infs = []
    for partition_1 in partitions:
        row_var_infs = []
        for partition_2 in partitions:
            var_inf = 0.0
            for i in partition_1:
                for j in partition_2:
                    a = set(i)
                    b = set(j)
                    n_ij = len(a.intersection(b))
                    if n_ij !=0: #Avoid log(0)
                        log_term = math.log(float(n_ij*n_ij)/float(len(i)*len(b)))
                        var_inf = var_inf + n_ij*log_term
            
            var_inf = -var_inf/float(n)
            row_var_infs.append(var_inf)
        all_var_infs.append(row_var_infs)
    
    return all_var_infs

import math
def vi_two_parts(partitions,n):
    var_inf = 0.0           
    for i in partitions[0]:
        for j in partitions[1]:
            a = set(i)
            b = set(j)
            n_ij = len(a.intersection(b))
            if n_ij !=0: #Avoid log(0)
                log_term = math.log(float(n_ij*n_ij)/float(len(i)*len(b)))
                var_inf = var_inf + n_ij*log_term
    
    var_inf = -var_inf/float(n)
    
    return var_inf