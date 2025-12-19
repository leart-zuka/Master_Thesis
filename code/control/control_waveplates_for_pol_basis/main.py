import sys
from gui import WaveplateGUI
from waveplates import Polarization_Waveplates
from PyQt5.QtWidgets import QApplication
# Settings

waveplate_port = "COM3"  # Windows
waveplate_port = "/dev/ttyUSB0"  # Linux


if __name__ == "__main__":
    polarization_waveplates = Polarization_Waveplates(port=waveplate_port)
    app = QApplication(sys.argv)
    win = WaveplateGUI(polarization_waveplates)
    win.show()
    sys.exit(app.exec_())
