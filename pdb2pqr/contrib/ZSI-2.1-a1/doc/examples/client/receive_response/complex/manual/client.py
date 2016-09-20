#!/usr/bin/env python

import sys

from ZSI import ServiceProxy

import ComplexTypes as MyComplexTypes

def main():
    server = ServiceProxy('../binding.wsdl',
                          typesmodule=MyComplexTypes,
                          tracefile=sys.stdout)
    
    user = server.GetUser('john_doe')
    print '   Age: %d' % user.Age
    print '  Name: %s' % user.Name
    print 'UserId: %s' % user.UserId
        
if __name__ == '__main__':
    main()
