'''
Created on 04/06/2012

@author: victor
'''
import json
from pyproct.tools.commonTools import convert_to_utf8, get_parameter_value

class ProtocolParameters():
    """
    Stores the values of the parameters used to execute the protocol.
    """
    def __init__(self, params_dic):
        self.params_dic = params_dic

    def __getitem__(self,key):
        return self.params_dic[key]

    def __contains__(self, key):
        return key in self.params_dic

    def __str__(self):
        return json.dumps(self.params_dic, sort_keys = False, indent = 4, separators = (',', ': '))

    def get_value(self, key_description, default_value):
        get_parameter_value(key_description, self.params_dic, default_value)

    @classmethod
    def get_params_from_json(cls, json_string):
        return ProtocolParameters(convert_to_utf8(json.loads(json_string)));

    @classmethod
    def get_default_params(cls, source):
        return ProtocolParameters.get_params_from_json(open(source,"r").read())



