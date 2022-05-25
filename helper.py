import galois as gl
import numpy as np
from multiprocessing import Pool

class Galois_Helper():

    def __init__(self,_field,p,debug_mode) -> None:
        self.field = _field
        self.base_field = gl.Field(_field.characteristic)
        self.debug_active = debug_mode
        self.pool = Pool(p)

    def get_power_of_field_element(self,el): 
        divider = self.field.primitive_element
        i = 1
        while int(el/divider) != 1:
            i+= 1
            divider *= self.field.primitive_element
        
        return i


    def get_nth_unity_root_of_field(self,n):
        power_of_one = (self.field.characteristic **self.field.degree)-1 #self.get_power_of_field_element(self.field(1))
        if power_of_one % n != 0:
            raise ValueError("failure in getting nth root") 
        return self.field.primitive_element**(power_of_one//n)
    

    def fft_on_matrix(self,matrix):
        n = len(matrix)
        output = []
        for i in range(0,n):
            output_vec = self.field(matrix[0])
            for j in range(1,n):
                p_factor = i * j
                p = self.get_nth_unity_root_of_field(n)
                vec = self.field(matrix[j])
                output_vec = output_vec + (p**p_factor * vec)
            output.append([int(x) for x in output_vec])
        return output
    
    def ifft_on_matrix(self,matrix): #aus irgend nem grund sind hier die sachen vertauscht dafq
        n = len(matrix)
        output = []
        for i in range(0,n):
            output_vec = self.field(matrix[0])
            for j in range(1,n):
                p_factor = -i * j

                p = self.get_nth_unity_root_of_field(n)
                vec = self.field(matrix[j])
                output_vec = output_vec + (p**p_factor * vec)
            output_vec= int(-self.base_field(1)/ self.base_field(n % self.base_field.characteristic)) * output_vec
            output.append([int(x) for x in output_vec])
        return output

    
    def debug_print(self,*args):
        if self.debug_active:
            print(*args)




    
