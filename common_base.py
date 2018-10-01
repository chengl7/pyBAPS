#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:48:16 2018

common list for functions

@author: lcheng
"""

class funcList:
    FUNC_NAME_LIST = ('f')
    FUNC_NAME_DICT = {}
    
    @classmethod
    def init(cls):
        cls.FUNC_NAME_DICT = {cls.FUNC_NAME_LIST[i]:i \
                              for i in range(len(cls.FUNC_NAME_LIST))}
    
    @classmethod
    def get_func_ind(cls,funcstr):
        if str in cls.FUNC_NAME_DICT:
            return cls.FUNC_NAME_DICT[funcstr]
        else:
            raise KeyError(f'Given function name {funcstr} is not supported.')
    
    @classmethod
    def get_func_name(cls, index):
        return cls.FUNC_NAME_LIST[index]
    
    # to be called by the worker to get the actual function object
    @classmethod
    def get_func_dict(cls):
        import sys
        return { cls.FUNC_NAME_LIST[i]:getattr(sys.modules['__main__'],cls.FUNC_NAME_LIST[i]) \
                for i in range(len(cls.FUNC_NAME_LIST)) }