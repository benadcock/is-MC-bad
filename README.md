# is-MC-bad
Code used to generate the figures in the article "Is Monte Carlo a bad sampling strategy for learning smooth functions in high dimensions?"

## Article's abstract

This paper concerns the approximation of smooth, high-dimensional functions on bounded hypercubes from limited samples using polynomials. This task lies at the heart of many applications in computational science and engineering -- notably, those arising from parametric modelling and computational uncertainty quantification. It is common to use Monte Carlo sampling in such applications, so as not to succumb to the curse of dimensionality. However, it is well known that such a strategy is theoretically suboptimal. Specifically, there are many polynomial spaces of dimension $n$ for which the sample complexity scales log-quadratically, i.e., like $c \cdot n^2 \cdot \log(n)$ as $n \rightarrow \infty$. This well-documented phenomenon has led to a concerted effort over the last decade to design improved, in fact, near-optimal strategies, whose sample complexities scale log-linearly, or even linearly in $n$.


Paradoxically, in this work we demonstrate that Monte Carlo is actually a perfectly good strategy in high dimensions, despite this apparent suboptimality. We first document this phenomenon empirically via several numerical examples. Next, we present a theoretical analysis that resolves this seeming contradiction in the case of holomorphic functions of infinitely-many variables. We show that there is a least-squares approximation based on $m$ Monte Carlo samples whose error decays algebraically fast in $m/\log(m)$, with a rate that is the same as that of the best $n$-term polynomial approximation. This result is non-constructive, since it assumes knowledge of a suitable polynomial subspace in which to perform the approximation. We next present a compressed sensing-based scheme that achieves the same rate, except for a larger polylogarithmic factor.  This scheme is practical, and numerically it performs as well as or better than well-known adaptive least-squares schemes. 

Overall,  our findings in this paper demonstrate that Monte Carlo sampling is eminently suitable for smooth function approximation tasks when the dimension is sufficiently high. Hence the benefits of improved sampling strategy are generically limited to lower-dimensional settings.

## Code organization

### Disclaimer

Since the experiments are random in nature, the results obtained from running the scripts may differ slighlty from those in the paper.
