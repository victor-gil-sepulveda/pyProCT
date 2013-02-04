'''
Created on 04/06/2012

@author: victor
'''
import json
from pyproclust.tools.commonTools import convert

class ProtocolParameters():
    """
    Stores the values of the parameters used to execute the protocol. 
    """
    def __init__(self, params_dic):
        self.params_dic = params_dic

    def __getitem__(self,key):
        return self.params_dic[key]
    
    def __str__(self):
        return json.dumps(self.params_dic, sort_keys = False, indent=4, separators=(',', ': '))        
    
    @classmethod
    def get_params_from_json(cls, json_string):
        return convert(json.loads(json_string));
    
    @classmethod
    def get_default_params(cls, source):
        source_json_string = open(source,"r")
        return ProtocolParameters.get_params_from_json(source_json_string)
    