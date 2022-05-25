# Generalized Reed Solomon Code
This is an implementation of a special kind of Reed Solomon Code called Generalized Reed Solomon code according to the paper "A Class of Generalized RS-Codes with Faster
Encoding and Decoding Algorithms" by Amin Shokrollahi. A version of the paper can be found at https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.306.4767&rep=rep1&type=pdf.

This special code can be sped up by a factor of p, since the code can be parallelized over p components. Additionally it allows to construct codes of arbitrary length as long as n%p == 0 is fulfilled. Whereas the standard Reed Solomon code imposes n= q^m -1.

## Documentation
n = k + r <br/>
k = message symbols sent <br/>
r = number of parity symbols <br/>
p = speed up factor <br/>

Currently the code can also handle codes where n %p != 0 however only if r is still divisble by p.
```
class Generalized_Reed_Solomon(
field_size:int, #q characterisitc of the field
message_length:int, # n
payload_length:int, # k
symbol_size:int, #m exponent of field. Order of field is q^m
p_factor:int, #p speedup factor
irr_poly=None, # optional in case that the galois libary has not a precomputed irreducible polynominal
debug=False # can be set to true in order to have the steps printed to the console
)
```
## Usage
```
import generalizedReedSolomon
import numpy as np

org_info = [30, 20, 20, 9, 25, 1, 7, 2, 0, 18, 17, 22, 11,0,0]
org_compare = np.copy(org_info)
test_normal_rs = generalizedReedSolomon.Generalized_Reed_Solomon(field_size=31,
                                                                 message_length=27,
                                                                 payload_length=15,
                                                                 symbol_size=1,
                                                                 p_factor=3)
normal_msg = test_normal_rs.encode(org_info)
normal_msg[6] = 10
normal_msg[9] = 12
normal_msg[10] = 0
normal_msg[16] = 4
normal_msg[23] = 10

dec_test_normal = test_normal_rs.decode(normal_msg)
print("dec",dec_test_normal)
print("org",org_info)
print("equal",org_info == dec_test_normal)

#dec [30, 20, 20, 9, 25, 1, 7, 2, 0, 18, 17, 22, 11, 0, 0]
#org [30, 20, 20, 9, 25, 1, 7, 2, 0, 18, 17, 22, 11, 0, 0]
#equal True
```
## Dependencies
galois == 0.0.27 <br/>
numpy == 1.19.5

## TODO
1. Add case for r % p != 0 <br/>
2. Add parallelization
3. Add functionality to directly input blocks into encode and output blocks for decode
