# Algorithms for Finding Orientations of Supersingular Elliptic Curves and Quaternion Orders.

This code accompanies the paper "Finding Orientations of Supersingular Elliptic Curves and Quaternion Orders", available on [arXiv](https://arxiv.org/abs/2308.11539) and [eprint](https://eprint.iacr.org/2023/1268). See the paper for the list of authors.  

It includes the following:

1. A brute force algorithm for solving the $\mathfrak{O}$-Orienting problem - Given a curve $E$ which can be oriented by a quadratic order $\mathfrak{O}$, find the orientation.  
2. Algorithms 5.1 and 5.2 for the *Quaternion Embedding Problem* - Finding embeddings of a quadratic order within a maximal quaternion order in $B_{p,\infty}$. 
    + Includes rerandomization of the basis of the order, making it very fast for any quadratic order with discriminant less than $O(p)$.
    + Algorithm 5.3 for checking if a solution gives a primitive embedding.
    + Experimental results on finding primitive embeddings (i.e. orientations) as described in Section 5.5.
3. Examples for Lemma 5.3 - Checking properties of the Hermite normal form of basis of quaternion orders.

The code is explained with examples in jupyter notebook files which can be found in the `Notebooks` directory. However, if you just need the code in `.py` files, it is in the repository root. Tested in Sagemath 10.0.