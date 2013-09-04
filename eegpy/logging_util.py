# -*- coding: utf-8 -*-
"""
Abstraction layer for logging
"""
import logging
import sys

def set_handler_for_loggers(loggers, level = logging.DEBUG, filename = None):
    if filename == None:
        script_name = sys.argv[0]
        if len(script_name) > 0:
            filename = script_name + ".log"
        else:
            raise ValueError("Please give a filename, couldn't determine good" +\
                             "name automatically.")
    log_handler = logging.handlers.RotatingFileHandler(filename, 
                                maxBytes = 1024*200, backupCount=5, 
                                encoding = "UTF-8")
    log_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    for logger in loggers:
        logger.addHandler(log_handler)
        logger.setLevel(level)

def get_loggers(logger_names, *args, **kwargs):
    loggers = [logging.getLogger(logger_name) for logger_name in logger_names]
    set_handler_for_loggers(loggers, *args, **kwargs)
    return loggers
    
def prepare_logger(logger_name, all_loggers=None, *args, **kwargs):
    """Prepares all desired loggers and returns the logger defined by logger_name"""
    all_loggers = list([] if all_loggers is None else all_loggers)
    try:
        script_logger_index = all_loggers.index(".")
    except ValueError:
        script_logger_index = 0
        all_loggers.insert(0,".")        
    all_loggers[script_logger_index] = logger_name
    loggers = [logging.getLogger(ln) for ln in all_loggers]
    set_handler_for_loggers(loggers, *args, **kwargs)
    return logging.getLogger(logger_name)

