import os
import sys

path = '/'
if path not in sys.path:
    sys.path.append(path)

os.environ['DJANGO_SETTINGS_MODULE'] = 'igem.settings'

import django.core.handlers.wsgi

application = django.core.handlers.wsgi.WSGIHandler()



#TODO: add robots.txt and favicon.ico
