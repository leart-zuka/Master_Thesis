from flask import Flask, request, jsonify, render_template
import json
from pathlib import Path
import threading
from collections import deque
from waveplates import Polarization_Waveplates
import time

app = Flask(__name__)

SETTINGS_FILE = Path("waveplate_settings.json")
POLARIZATIONS = ["H", "V", "D", "A", "R", "L"]

# ---------------- Hardware ---------------- #

waveplate_port = "/dev/ttyUSB0"
polarization_waveplates = Polarization_Waveplates(port=waveplate_port)

is_moving = False
lock = threading.Lock()

# ---------------- Settings ---------------- #

# store last 100 messages
log = deque(maxlen=100)


def log_msg(msg):
    log.append(msg)
    print(msg)  # still print to terminal


def load_settings():
    if SETTINGS_FILE.exists():
        with open(SETTINGS_FILE) as f:
            return json.load(f)
    return {p: {"QWP": 0.0, "HWP": 0.0} for p in POLARIZATIONS}


def save_settings(settings):
    with open(SETTINGS_FILE, "w") as f:
        json.dump(settings, f, indent=4)


settings = load_settings()

# ---------------- Logic ---------------- #


def apply_polarization(pol):
    angles = settings[pol]
    log_msg(f"Applying {pol}: HWP={angles['HWP']} QWP={angles['QWP']}")

    polarization_waveplates.l_2.set_angle(angles["HWP"])

    time.sleep(2)
    _finish_qwp_move(angles["QWP"])


def _finish_qwp_move(qwp_angle):
    # global is_moving
    polarization_waveplates.l_4.set_angle(qwp_angle)
    log_msg(f"QWP move complete to {qwp_angle}")


# ---------------- Routes ---------------- #


@app.route("/")
def index():
    return render_template("index.html", settings=settings)


@app.route("/logs")
def get_logs():
    # return all current log messages
    return jsonify(list(log))


@app.route("/select", methods=["POST"])
def select():
    pol = request.json["pol"]
    apply_polarization(pol)
    return jsonify(success=True)


@app.route("/save", methods=["POST"])
def save():
    pol = request.json["pol"]
    settings[pol]["QWP"] = float(request.json["QWP"])
    settings[pol]["HWP"] = float(request.json["HWP"])
    save_settings(settings)
    return jsonify(success=True)


@app.route("/shutdown", methods=["POST"])
def shutdown():
    print("Shutdown requested")
    return jsonify(success=True)


# ---------------- Run ---------------- #

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
