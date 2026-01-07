from time import sleep
import elliptec
import sys

controller = elliptec.Controller("/dev/ttyUSB0")
rotator_1 = elliptec.Rotator(controller, address="0")
rotator_2 = elliptec.Rotator(controller, address="1")


class Polarization_Waveplates:
    def __init__(self, port: str = "/dev/ttyUSB0"):
        try:
            self.controller = elliptec.Controller(port)
            self.l_2 = elliptec.Rotator(self.controller, address="0")
            self.l_4 = elliptec.Rotator(self.controller, address="1")
            self.l_2.home()
            sleep(2)
            self.l_4.home()
            sleep(2)
        except Exception as e:
            print("Wasn't able to connect to waveplates, exiting now")
            print(f"Whole error: \n {e}")
            self.controller.close_connection()
            sys.exit()

    def close(self):
        self.controller.close_connection()
