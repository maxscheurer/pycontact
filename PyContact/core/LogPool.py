# from:
# Shortcut to multiprocessing's logger

import multiprocessing
from multiprocessing.pool import Pool
import traceback


def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class LogExceptions(object):
    def __init__(self, argcallable):
        self.__callable = argcallable

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception:
            error(traceback.format_exc())
            raise

        return result


class LoggingPool(Pool):
    """Used for multiprocessing logging."""
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)
