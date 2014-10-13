"""
Created on Sep 3, 2014

@author: victor
"""
import glob
import os
from pyproct.data.handler.dataSource import DataSource
import json
from pyproct.tools.commonTools import convert_to_utf8
import copy

class SourceGenerator(object):
    """
    Performs source expansion and yields data sources.
    """
    def __init__(self, source_list):
        self.source_list = [DataSource(source) for source in SourceGenerator.inflate_source_list(source_list)]
        
    @classmethod
    def inflate_source_list(cls, source_list):    
        inflated_list = []
        for source in source_list:
            
            if isinstance(source, basestring):
                # is string-like
                _, ext = os.path.splitext(source)
                if ext == ".lst": #pyProCT file list
                    try:
                        # If it is a list, we keep inflating.
                        # Infinite recursivity (circular references for instance) is not
                        # checked.
                        inflated_list.extend(cls.inflate_source_list(cls.get_sources_from_file_list(source)))
                    except Exception, e:
                        print "[ERROR SourceGenerator::init] Impossible to read list file %s"%source
                        print e.message
                        exit()
                else:
                    # Is a normal string
                    inflated_list.extend(cls.do_glob(source))
            else:
                # is dict-like
                inflated_list.extend(cls.get_sources_from_dictionary(source))
        return inflated_list
    
    @classmethod
    def do_glob(cls, path):
        paths = glob.glob(path)
        if len(paths) > 0:
            return paths
        else:
            print "[ERROR SourceGenerator::init] Impossible to find one or all of this files: %s"%path
            exit()
    
    @classmethod                
    def get_sources_from_file_list(cls, list_file):
        """
        Gets the sources from a ".lst" file. Each line of this file describes a source using:
            - A single string: the source path.
            - A json object with at least a "source" keyword with the source path.
        
        :param list_file: The path of the list file.
        
        :return: An array with the contents of the file (strings and dictionaries) 
        """
        return [convert_to_utf8(json.loads(line)) for line in open(list_file,"r")]
    
    @classmethod
    def get_sources_from_dictionary(cls, info_dict):
        """
        If a dictionary source string is a glob and can be inflated, creates a dictionary copy
        with a different "source" key contents.
        
        :param info_dict: A dictionary with at least "source" keyword. 
        
        :return: An array with the same array or more in case it was inflated.
        """
        inflated_dics = []
        for path in cls.do_glob(info_dict["source"]):
            clone = copy.deepcopy(info_dict)
            clone["source"] = path
            inflated_dics.append(clone)
        return inflated_dics
