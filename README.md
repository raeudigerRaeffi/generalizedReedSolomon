# generalizedReedSolomon
This is an implementation of a special kind of Reed Solomon Code called Generalized Reed Solomon code according to the paper "A Class of Generalized RS-Codes with Faster
Encoding and Decoding Algorithms" by Amin Shokrollahi. A version of the paper can be found at https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.306.4767&rep=rep1&type=pdf.

This special code can be sped up by a factor of p, since the code can be parallelized over p components, as well as allows to construct codes of arbitrary length as long as n%p == 0 is fulfilled. Whereas the standard Reed Solomon code imposes n= q^m -1.
