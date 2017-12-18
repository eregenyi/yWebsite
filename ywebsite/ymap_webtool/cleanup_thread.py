from threading import Thread
import time

from .tasks import *

class CleaningLady(Thread):

    HOUR_IN_SECONDS = 3600

    def __init__(self, input_path, output_path, archive_path, seconds_to_wait: int=HOUR_IN_SECONDS, logger=None):
        super().__init__()
        self.input_path = input_path
        self.output_path = output_path
        self.archive_path = archive_path
        self.seconds_to_wait = seconds_to_wait
        self.logger = logger

    def run(self) -> None:

        if self.logger is not None:
            self.logger.debug("Cleaning lady is sleeping " + str(self.seconds_to_wait) + " seconds..")

        time.sleep(self.seconds_to_wait)

        if self.logger is not None:
            self.logger.debug("Cleaning lady is woke af..")

        clean_up(self.input_path, self.output_path, self.archive_path)

        if self.logger is not None:
            self.logger.debug("Cleaning lady cleaned up.")

