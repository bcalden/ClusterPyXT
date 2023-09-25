"""
File: cpxt/core/stages.py
Author: Brian Alden
Date: 29 Mar 2022
Description: This file contains the Stage enumeration. It is used to make 
             references to a clusters progress through `cpxt` without using
             integers and/or strings.

"""
from enum import IntEnum, Enum

class Stage(IntEnum):
    """
    This class is used so we can reference what stage the of pipeline a 
    particular cluster is at. It is an Integer Enumeration class and as 
    such there are no methods or attributes.
    """
    none = -1
    init = 0
    one = 1
    two = 2
    three = 3
    four = 4
    five = 5
    spec = 6

class MessageType(Enum):
    """
    This class is used to reference the type of message being sent to the
    user and logged to the log file.
    """
    preparation = "preparation"
    completion = "completion"