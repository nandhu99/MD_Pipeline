#!/usr/bin/env python

from Example_services import *
from ZSI import dispatch

def echo(message):
    response = EchoResponseWrapper()
    response._Message = message
    return response

if __name__ == '__main__':
    dispatch.AsServer(port=8080)
    
