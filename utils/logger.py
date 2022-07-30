# -*- coding: utf-8 -*-

"""
Logger.
"""


# Import modules
import logging.config
from warnings import filterwarnings


def load_logging_cfg():
    # Define three log output formats
    standard_format = '%(asctime)-15s  [%(threadName)s:%(thread)d, task_id:%(name)s, %(filename)s:%(lineno)d]' \
                  '\n[%(levelname)s]  %(message)s\n' 
    simple_format = '[%(levelname)s]  %(asctime)-15s\n%(message)s\n'
    #id_simple_format = '[%(levelname)s][%(asctime)s]%(message)s'
    # Full path of log file
    logfile = DEFAULT['logfile']
    # Configuration of the logger
    LOGGING_DIC = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'standard': {
                'format': standard_format
            },
            'simple': {
                'format': simple_format
            },
        },
        'filters': {},
        'handlers': {
            # Print logs to terminal
            'console': {
                'level': 'DEBUG',
                'class': 'logging.StreamHandler',
                'formatter': 'simple',
            },
            # Print logs into file, collect logs above 'INFO' 
            'default': {
                'level': 'DEBUG',
                'class': 'logging.handlers.RotatingFileHandler', 
                'formatter': 'standard',
                'filename': logfile, 
                'maxBytes': 1024*1024*5, 
                'backupCount': 5,
                'encoding': 'utf-8',
            },
        },
        'loggers': {
            # Using logging.getLogger(__name__) to get the configuration
            '': {
                'handlers': ['default', 'console'],
                'level': 'DEBUG',
                'propagate': True,
            },
        },
    }    
    logging.config.dictConfig(LOGGING_DIC)
    global logger
    logger = logging.getLogger(__name__)