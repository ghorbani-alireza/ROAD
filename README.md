A ROAD to Classification in High-Dimensional Space
================

This repository contains all the R code for the high dimensional classification method introduced in the paper “A ROAD to Classification in High Dimensional Space (2010)” which can be found [here](http://yangfeng.hosting.nyu.edu/publication/fan-2012-road/fan-2012-road.pdf).

**Abstract:** For high dimensional classification, it is well known that naively performing the Fisher discriminant rule leads to poor results due to diverging spectra and accumulation of noise. Therefore, researchers proposed independence rules to circumvent the diverging spectra, and sparse independence rules to mitigate the issue of accumulation of noise. However, in
biological applications, often a group of correlated genes are responsible for clinical outcomes, and the use of the covariance information can significantly reduce misclassification rates. In theory, the extent of such error rate reductions is unveiled by comparing the misclassification rates of the Fisher discriminant rule and the independence rule. To materialize the gain on the basis of finite samples, a regularized optimal affine discriminant (ROAD) is proposed. The ROAD
selects an increasing number of features as the regularization relaxes. Further benefits can be achieved when a screening method is employed to narrow the feature pool before applying the ROAD method. An efficient constrained co-ordinate descent algorithm is also developed to solve the associated optimization problems. Sampling properties of oracle type are established. Simulation studies and real data analysis support our theoretical results and demonstrate the
advantages of the new classification procedure under a variety of correlation structures. A delicate result on continuous piecewise linear solution paths for the ROAD optimization problem at the population level justifies the linear interpolation of the constrained co-ordinate descent algorithm.

**Keywords:** Fisher discriminant; High dimensional classification; Independence rule; Linear discriminant analysis; Regularized optimal affine discriminant
