from xmlrpc.client import Boolean
import galois as gl

from . import helper

class Base_Reed_Solomon():
    def __init__(self, field_size:int, message_length:int, payload_length:int,symbol_size:int,irr_poly=None,debug=True,p=1):
        self.field_size = field_size
        if irr_poly ==None:
            try:
                self.galois_field = gl.GF(field_size**symbol_size)#TODO remove power asap
            except LookupError:
                print("need to find irreducible polynominal manually")
                irreducible_polynominal = gl.irreducible_poly(field_size,symbol_size)
                self.galois_field = gl.GF(field_size**symbol_size,irreducible_poly=irreducible_polynominal)   
        else:
            self.galois_field = gl.GF(field_size**symbol_size,irreducible_poly=irr_poly)   
        self.message_length = message_length
        self.payload_length = payload_length
        self.two_s = message_length - payload_length
        self.helper = helper.Galois_Helper(self.galois_field,p,debug)
    
    

    def convert_to_symbol_array(self,matrix):
        output = []
        print(matrix.shape)
        for i in range(0,matrix.shape[1]):
            output.append(int(self.galois_field.Vector(matrix[:,i])))
        return output