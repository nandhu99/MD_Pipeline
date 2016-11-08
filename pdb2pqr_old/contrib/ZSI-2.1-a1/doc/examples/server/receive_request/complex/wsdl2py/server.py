#!/usr/bin/env python

import Registration_services_types
from Registration_services import RegisterUserResponseWrapper
from ZSI import dispatch


def RegisterUser(user):
    response = RegisterUserResponseWrapper()
    response._Message = "OK"
    return response

if __name__ == '__main__':
    dispatch.AsServer(port=8080, typesmodule=(Registration_services_types,))
