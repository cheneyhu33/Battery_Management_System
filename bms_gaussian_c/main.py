"""bms — web view of the C thermal model.

The model is `bms_gaussian.c`, compiled for this board and run as a child
process. This file contains no physics at all: it reads the binary frame
stream the child writes to stdout and relays each frame to the browser. There
is therefore exactly one implementation of the model, the C one, and the page
shows what that implementation computed.

    stdout : one binary frame every STREAM_EVERY steps  (see bms_gaussian.c)
    stderr : the model's own log, which lands in `arduino-app-cli app logs`

Frame layout, little-endian, no padding:

    char   magic[4]   "BMSF"
    int32  nodes, step, final
    double t, t_max, t_center, t_surf, energy, t_amb, t_limit
    float  xy[P*P], xz[P*P], yz[P*P]      temperatures as rises above ambient

The three planes are the slices through the cell centre, flattened row-major
as [first index][second index], which is what the page expects.
"""

import base64
import os
import struct
import subprocess
import threading
import time

import numpy as np

from arduino.app_bricks.web_ui import WebUI
from arduino.app_utils import App

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "bin", "bms_gaussian")
STREAM_EVERY = 10           # model steps per frame

HEADER = struct.Struct("<4siii7d")
MAGIC = b"BMSF"
PLANE_KEYS = ("xy", "xz", "yz")

ui = WebUI()

_lock = threading.Lock()
_latest = None              # most recent complete frame, or None
_published_step = -1
_finished = False


def _read_exact(stream, n):
    buf = bytearray()
    while len(buf) < n:
        chunk = stream.read(n - len(buf))
        if not chunk:
            raise EOFError
        buf.extend(chunk)
    return bytes(buf)


def _encode_plane(plane, scale):
    """Quantise one 2D slice to 0..255 against this frame's peak, then base64.

    rint rather than a bare astype: astype truncates, which biases the whole
    field low by up to 1 LSB and drops the full-scale cell to 254.
    """
    q = np.clip(np.rint(plane * scale), 0.0, 255.0).astype(np.uint8)
    return base64.b64encode(np.ascontiguousarray(q).tobytes()).decode("ascii")


def _read_frame(stream):
    """Read one frame and turn it into the dict the page consumes."""
    magic, nodes, step, final, t, t_max, t_center, t_surf, energy, t_amb, t_limit = \
        HEADER.unpack(_read_exact(stream, HEADER.size))
    if magic != MAGIC:
        raise ValueError("bad frame magic %r" % (magic,))

    n = nodes * nodes
    planes = {}
    for key in PLANE_KEYS:
        raw = _read_exact(stream, 4 * n)
        planes[key] = np.frombuffer(raw, dtype="<f4").reshape(nodes, nodes)

    scale = 255.0 / t_max if t_max > 0.0 else 0.0
    return {
        "step": step,
        "t": t,
        "state": "stopped" if final else "running",
        "t_limit": t_limit,
        "t_amb": t_amb,
        "t_max": t_amb + t_max,
        "t_center": t_amb + t_center,
        "t_surf": t_amb + t_surf,
        "energy": energy,
        "nodes": nodes,
        "planes": {k: _encode_plane(v, scale) for k, v in planes.items()},
    }


def _reader(proc):
    global _latest, _finished
    try:
        while True:
            frame = _read_frame(proc.stdout)
            with _lock:
                _latest = frame
            if frame["state"] == "stopped":
                break
    except EOFError:
        print("frame stream ended without a final frame")
    except Exception as exc:                      # noqa: BLE001
        print("frame stream error: %r" % (exc,))
    finally:
        with _lock:
            _finished = True
    print("model finished; serving the final frame")


def publish_frame():
    """Push the newest frame to every connected browser, at most one per loop."""
    global _published_step
    with _lock:
        frame = _latest
    if frame is None or frame["step"] == _published_step:
        return
    _published_step = frame["step"]
    ui.send_message("heatmap", frame)


def on_request_state(sid, data):
    """A browser that connects mid-run asks for the current frame."""
    with _lock:
        return _latest


ui.on_message("request_state", on_request_state)


# =========================================================
# START THE MODEL
# =========================================================

print("=" * 72)
print("bms - web view of the C model")
print("=" * 72)
print("model    : %s" % BINARY)
print("stream   : one frame every %d steps" % STREAM_EVERY)
print("physics  : none in this process; see bms_gaussian.c")
print("-" * 72)

if not os.path.isfile(BINARY):
    raise SystemExit("model binary not found at %s" % BINARY)
if not os.access(BINARY, os.X_OK):
    raise SystemExit("model binary is not executable: %s" % BINARY)

_proc = None


def _start_model():
    """Launch the model on the first loop() call, not at import.

    App.run() takes about twenty seconds to bring the web server up. The whole
    model takes about as long, so starting it at import would mean it finished
    before the page was reachable, and the browser would only ever see the
    final frame. loop() does not run until the server is ready.
    """
    global _proc
    if _proc is not None:
        return
    # stderr is inherited, so the model's banner and per-500-step lines land in
    # the app log next to this process's own output.
    _proc = subprocess.Popen(
        [BINARY, "--stream", str(STREAM_EVERY)],
        stdout=subprocess.PIPE,
        stderr=None,
        bufsize=0,
    )
    threading.Thread(target=_reader, args=(_proc,), daemon=True).start()
    print("model started, pid %d" % _proc.pid)


def loop():
    """Called repeatedly by the App framework."""
    _start_model()
    publish_frame()
    if _finished:
        # Keep the final frame on screen; the reader has already stored it.
        time.sleep(0.25)
        return
    time.sleep(0.05)


App.run(user_loop=loop)
