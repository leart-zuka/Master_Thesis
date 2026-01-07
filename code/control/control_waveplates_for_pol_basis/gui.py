import sys
import json
from pathlib import Path
from time import sleep

from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QPushButton,
    QLabel,
    QDoubleSpinBox,
    QFileDialog,
    QMessageBox,
    QGridLayout,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QFrame,
)


from PyQt5.QtCore import Qt
from PyQt5 import QtGui

from waveplates import Polarization_Waveplates
import typing


SETTINGS_FILE = Path("waveplate_settings.json")


class WaveplateGUI(QMainWindow):
    def __init__(self, polarization_waveplates: Polarization_Waveplates):
        super().__init__()
        self.setWindowTitle("Waveplate Control")
        self.resize(1000, 600)

        self.polarizations = ["H", "D", "R", "V", "A", "L"]

        # internal state
        self.settings = {pol: {"QWP": 0.0, "HWP": 0.0} for pol in self.polarizations}
        self.current_pol = None
        self.polarization_waveplates = polarization_waveplates

        self._build_ui()
        self._apply_style()

    # ---------------- UI ---------------- #

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)

        main_layout = QHBoxLayout(central)

        # Left: polarization buttons
        left_frame = QFrame()
        left_frame.setFrameShape(QFrame.StyledPanel)
        left_layout = QGridLayout(left_frame)

        self.pol_buttons = {}
        for i, pol in enumerate(self.polarizations):
            btn = QPushButton(pol)
            btn.setCheckable(True)
            btn.setMinimumSize(100, 80)
            btn.clicked.connect(lambda _, p=pol: self.select_polarization(p))
            self.pol_buttons[pol] = btn
            left_layout.addWidget(btn, i // 3, i % 3)

        # Right: angle editor
        right_frame = QFrame()
        right_frame.setFrameShape(QFrame.StyledPanel)
        right_layout = QVBoxLayout(right_frame)

        title = QLabel("Waveplate Angles")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size: 18px; font-weight: bold;")

        form = QFormLayout()

        self.qwp_spin = QDoubleSpinBox()
        self.qwp_spin.setRange(-360.0, 360.0)
        self.qwp_spin.setDecimals(3)
        self.qwp_spin.setSuffix(" °")

        self.hwp_spin = QDoubleSpinBox()
        self.hwp_spin.setRange(-360.0, 360.0)
        self.hwp_spin.setDecimals(3)
        self.hwp_spin.setSuffix(" °")

        form.addRow("QWP:", self.qwp_spin)
        form.addRow("HWP:", self.hwp_spin)

        save_btn = QPushButton("Save to Polarization")
        save_btn.clicked.connect(self.save_current_angles)

        io_layout = QHBoxLayout()
        load_btn = QPushButton("Load JSON")
        save_json_btn = QPushButton("Save JSON")
        load_btn.clicked.connect(self.load_json)
        save_json_btn.clicked.connect(self.save_json)

        exit_btn = QPushButton("Exit")
        exit_btn.setMinimumSize(120, 40)
        exit_btn.clicked.connect(self.close)
        exit_btn.setShortcut("Ctrl+Q")

        io_layout.addWidget(load_btn)
        io_layout.addWidget(save_json_btn)

        right_layout.addWidget(title)
        right_layout.addLayout(form)
        right_layout.addWidget(save_btn)
        right_layout.addStretch()
        right_layout.addLayout(io_layout)
        right_layout.addWidget(exit_btn)

        # proportions
        main_layout.addWidget(left_frame, 3)
        main_layout.addWidget(right_frame, 2)

    # ---------------- Logic ---------------- #

    def select_polarization(self, pol):
        self.current_pol = pol
        angles = self.settings[pol]

        self.qwp_spin.setValue(angles["QWP"])
        self.hwp_spin.setValue(angles["HWP"])

        for p, btn in self.pol_buttons.items():
            btn.setChecked(p == pol)

        # 🔧 APPLY POLARIZATION IMMEDIATELY
        self.apply_polarization(pol)

    def apply_polarization(self, pol):
        angles = self.settings[pol]
        self.polarization_waveplates.l_2.set_angle(angles["HWP"])
        sleep(2)
        self.polarization_waveplates.l_4.set_angle(angles["QWP"])

    def save_current_angles(self):
        if self.current_pol is None:
            QMessageBox.warning(self, "No Selection", "Select a polarization first.")
            return

        self.settings[self.current_pol] = {
            "QWP": self.qwp_spin.value(),
            "HWP": self.hwp_spin.value(),
        }

        self.flash_button(self.pol_buttons[self.current_pol])

    def save_json(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Settings", str(SETTINGS_FILE), "JSON (*.json)"
        )
        if not path:
            return

        with open(path, "w") as f:
            json.dump(self.settings, f, indent=4)

    def load_json(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Settings", "", "JSON (*.json)"
        )
        if not path:
            return

        with open(path, "r") as f:
            self.settings = json.load(f)

        QMessageBox.information(self, "Loaded", "Settings loaded successfully.")

    # ---------------- Styling ---------------- #

    def _apply_style(self):
        self.setStyleSheet("""
                QWidget {
                    background-color: #1e1e1e;
                    color: #dddddd;
                    font-size: 14px;
                }

                QPushButton {
                    background-color: #2d2d2d;
                    border: 1px solid #444;
                    border-radius: 6px;
                    padding: 10px;
                }

                QPushButton:hover {
                    background-color: #3a3a3a;
                }

                QPushButton:checked {
                    background-color: #005f87;
                    border: 1px solid #00bfff;
                }

                QPushButton[saved="true"] {
                    background-color: #2ecc71;
                    color: #000000;
                    border: 1px solid #27ae60;
                }

                QFrame {
                    border-radius: 8px;
                }

                QDoubleSpinBox {
                    background-color: #2a2a2a;
                    border: 1px solid #444;
                    border-radius: 4px;
                    padding: 4px;
                }

                QPushButton#exitButton {
                    background-color: #8b0000;
                    border: 1px solid #cc4444;
                    color: #ffffff;
                }

                QPushButton#exitButton:hover {
                    background-color: #b22222;
                }
            """)

    def flash_button(self, button, duration_ms=300):
        button.setProperty("saved", True)
        button.style().unpolish(button)
        button.style().polish(button)

        QTimer.singleShot(duration_ms, lambda: self._clear_flash(button))

    def _clear_flash(self, button):
        button.setProperty("saved", False)
        button.style().unpolish(button)
        button.style().polish(button)

    def closeEvent(self, event):
        try:
            self.polarization_waveplates.close()
        except Exception as e:
            print("Shutdown error:", e)
        event.accept()


# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     win = WaveplateGUI()
#     win.show()
#     sys.exit(app.exec_())
