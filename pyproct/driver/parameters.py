"""
Created on 04/06/2012

@author: victor
"""
import json
from pyproct.tools.commonTools import convert_to_utf8, get_parameter_value


class ProtocolParameters():
    """
    Stores the values of the parameters used to execute the protocol.
    """

    def __init__(self, params_dic):
        self.params_dic = params_dic

    def __getitem__(self,key):
        if isinstance(self.params_dic[key], dict):
            return ProtocolParameters(self.params_dic[key])
        else:
            return self.params_dic[key]

    def __setitem__(self,key,value):
        self.params_dic[key] = value

    def __contains__(self, key):
        return key in self.params_dic

    def __iter__(self):
        for k in self.keys():
            yield k

    def __str__(self):
        return json.dumps(ProtocolParameters.to_dict(self.params_dic), sort_keys = False, indent = 4, separators = (',', ': '))

    def get_value(self, key_description, default_value = None):
        return get_parameter_value(key_description, self.params_dic, default_value)

    def keys(self):
        return [k for k in self.params_dic.keys()]

    @classmethod
    def get_params_from_json(cls, json_string):
        return ProtocolParameters(convert_to_utf8(json.loads(json_string)));

    @classmethod
    def get_default_params(cls, source):
        return ProtocolParameters.get_params_from_json(open(source,"r").read())

    @classmethod
    def empty(cls):
        return ProtocolParameters({})

    @classmethod
    def to_dict(cls, proParam):
        new_dic = {}
        for k in proParam.keys():
            item  = proParam[k]
            if isinstance(item, ProtocolParameters):
                new_dic[k] = ProtocolParameters.to_dict(item)

            elif isinstance(item, dict):
                new_dic[k] = ProtocolParameters.to_dict(ProtocolParameters(proParam[k]))

            else:
                new_dic[k] = proParam[k]
        return new_dic







