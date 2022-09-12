import logging

logger_sctriangulate = logging.getLogger(__name__)  # It is a singleton, so multiple call with same name will always reference to same object
                                                    # no matter it is imported to main_class, or submodule, it's all the same object