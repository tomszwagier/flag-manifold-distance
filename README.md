# flag-manifold-distance
Authors' implementation of the paper ["Rethinking the Riemannian Logarithm on Flag Manifolds as an Orthogonal Alignment Problem"](https://link.springer.com/chapter/10.1007/978-3-031-38271-0_37), accepted with oral to the [GSI'23 6th International Conference on Geometric Science of Information](https://conference-gsi.org/).

![Illustrative summary of the paper](illustrative_summary.png)

### Installation
You can create a conda environment and then install the required packages by running the following commands on the Anaconda prompt.
```python
conda create -n flag-manifold-distance python=3.9
pip install -r requirements.txt
```

### Notes
- The code has been open-sourced for reproducibility purposes.
- An improvement of its structure, docstring and readability, as well as the design of an explanatory notebook with more details and illustrations will be considered between now and the conference.
- An integration of the code into [`geomstats`](https://github.com/geomstats/geomstats), an open-source Python package for geometric statistics, will be considered in the future. 


### Citation
```bibtex
@InProceedings{szwagier_flaglog_2023,
author="Szwagier, Tom
and Pennec, Xavier",
editor="Nielsen, Frank
and Barbaresco, Fr{\'e}d{\'e}ric",
title="Rethinking the Riemannian Logarithm on Flag Manifolds as an Orthogonal Alignment Problem",
booktitle="Geometric Science of Information",
year="2023",
publisher="Springer Nature Switzerland",
address="Cham",
pages="375--383",
abstract="Flags are sequences of nested linear subspaces of increasing dimension. They belong to smooth manifolds generalizing Grassmannians and bring a richer multi-scale point of view to the traditional subspace methods in statistical analysis. Hence, there is an increasing interest in generalizing the formulae and statistical methods already developed for Grassmannians to flag manifolds. In particular, it is critical to compute accurately and efficiently the geodesic distance and the logarithm due to their fundamental importance in geometric statistics. However, there is no explicit expression known in the case of flags. In this work, we exploit the homogeneous quotient space structure of flag manifolds and rethink the geodesic endpoint problem as an alignment of orthogonal matrices on their equivalence classes. The relaxed problem with the Frobenius metric surprisingly enjoys an explicit solution. This is the key to modify a previously proposed algorithm. We show that our explicit alignment step brings drastic improvements in accuracy, speed and radius of convergence, in addition to overcoming the combinatorial issues raised by the non-connectedness of the equivalence classes.",
isbn="978-3-031-38271-0"
}
