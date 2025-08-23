
# mp_logging.py
import sys
import logging
import multiprocessing as mp
from logging.handlers import QueueHandler, QueueListener

_LEVEL_MAP = {
    "Norm": logging.INFO,
    "Pass": logging.INFO,
    "Header": logging.INFO,
    "Warn": logging.WARNING,
    "Fail": logging.ERROR,
}

def setup_mp_logging(logfile: str, level=logging.INFO, console: bool = True):
    """
    Main-process setup: start a QueueListener that writes to logfile (and console).
    Also put a QueueHandler on the root logger so main process logs go to the queue too.
    Returns (queue, listener).
    """
    log_queue = mp.Queue()

    fmt = logging.Formatter("%(asctime)s [%(processName)s] %(levelname)s: %(message)s")

    file_h = logging.FileHandler(logfile, mode="a", encoding="utf-8")
    file_h.setFormatter(fmt)
    handlers = [file_h]

    if console:
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(fmt)
        handlers.append(sh)

    listener = QueueListener(log_queue, *handlers, respect_handler_level=True)
    listener.start()

    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(level)
    root.addHandler(QueueHandler(log_queue))

    return log_queue, listener

def init_worker_logging(log_queue, level=logging.INFO):
    """Worker-process initializer: attach a QueueHandler so worker logs go to the queue."""
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(level)
    root.addHandler(QueueHandler(log_queue))

def compat_logger(_logfile_unused, message: str = "", message_type: str = "Norm",
                  mode: str = "timestamp", s_type: str = "plain"):
    """Drop-in replacement for your utils.logger()."""
    level = _LEVEL_MAP.get(message_type, logging.INFO)
    logging.log(level, message)

def time_print_to_logger(message: str, message_type: str = "Norm"):
    level = _LEVEL_MAP.get(message_type, logging.INFO)
    logging.log(level, message)
