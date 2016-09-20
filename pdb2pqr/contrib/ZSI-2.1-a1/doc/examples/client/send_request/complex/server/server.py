#!/usr/bin/env python

from Registration_services import RegisterUserResponseWrapper
from ZSI import dispatch

import ComplexTypes as MyComplexTypes


def RegisterUser(user):
    response = RegisterUserResponseWrapper()
    response._Message = "OK"
    return response


if __name__ == '__main__':
    dispatch.AsServer(port=8080, typesmodule=(MyComplexTypes,))
