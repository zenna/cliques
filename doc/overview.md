# Robustness

The question we should answer is: *"When is a partition stable?"*.
This question is important to 1. provide confidence (or otherwise) when partitions are found through an optimisation process 2. understand real networks which may exploit, or depend upon robustness.  Concretely, we would like to be able to say:

1. The robustness of a partition is best described as property $X$ with respect to the landscape.
2. Our justification for this definition is $Y$, which we can show this directly with small artificial and natural graphs.
4. For larger graphs we can approximate this with measure $Z$

## General Approach
Our general approach is to understand robustness of a partition as dependent upon its relevant position in the space of *all* partitions.
With a small enough graph, we can find these by explicit enumeration.
The overall strategy is then.

1. Enumerate all partitions
2. Weigh each partition by a measure (e.g. stability)
3. Connect partitions by some notion of locality to form a partition graph
4. Analyse the graph

## Robustness measures
By stable, or robust, we informally mean is robust to perturbation.
Formalising this forces us to address:

__Key question 1:__ What are our candidate measures which we can compute directly with the full landscape

- Do they capture our intuition?
- How can we justify this as a ground truth?
- Is it reasonable to consider a single value, or should we consider distributions?
- How does this connect to existing literature on robustness, combinatorial or otherwise?

__Key question 2:__ To what extent is the stability of a partition dependent on a) the measure b) the moveset

__Key question 3:__ Can we find natural graphs which exhibit dynamics which depend upon robustness?

__Key question 4:__ What's the connection between robustness and timescales?

## TODO

### Code

- Codereview: cliques library, file by file
- Make cliques compile.
- Document file formats