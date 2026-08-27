"""Windows input tap for magneto_session — Raw Input (WM_INPUT) listener.

Counterpart of the macOS Quartz event tap: a message-only window registered
with RIDEV_INPUTSINK receives raw mouse deltas and key events system-wide,
including while CS2 holds raw-input focus (a low-level hook would not see
the deltas the game sees; WM_INPUT does, pre-acceleration).

Writes the same event vocabulary as the mac tap into the capture file:
  m dx dy | L / Lu / R | w a s d j downs, Wu Au Su Du Ju ups | W2/W2u C2/C2u
Timestamps are time.perf_counter_ns() (QPC), written as integer nanoseconds —
capture files started on Windows carry a `t_ns` header instead of `t_ticks`.

API: start(get_file) once, then set_enabled(bool) around live play.
Events are dropped (zero work beyond the callback) while disabled.
"""
import ctypes
import ctypes.wintypes as wt
import threading
import time

user32 = ctypes.windll.user32
kernel32 = ctypes.windll.kernel32

WM_INPUT = 0x00FF
RIDEV_INPUTSINK = 0x00000100
RID_INPUT = 0x10000003
RIM_TYPEMOUSE = 0
RIM_TYPEKEYBOARD = 1
MOUSE_MOVE_ABSOLUTE = 0x0001
RI_MOUSE_LEFT_BUTTON_DOWN = 0x0001
RI_MOUSE_LEFT_BUTTON_UP = 0x0002
RI_MOUSE_RIGHT_BUTTON_DOWN = 0x0004
RI_KEY_BREAK = 0x0001
HWND_MESSAGE = -3

WASD = {0x41: "a", 0x53: "s", 0x44: "d", 0x57: "w", 0x20: "j",  # j = jump (space)
        0x43: "c", 0x51: "q", 0x45: "e", 0x58: "x"}  # ability keys (Val util)
MODS = {0x10: "W2", 0x11: "C2"}  # shift = walkmod, ctrl = crouchmod

LRESULT = ctypes.c_ssize_t
WNDPROC = ctypes.WINFUNCTYPE(LRESULT, wt.HWND, wt.UINT, wt.WPARAM, wt.LPARAM)
user32.DefWindowProcW.restype = LRESULT
user32.DefWindowProcW.argtypes = [wt.HWND, wt.UINT, wt.WPARAM, wt.LPARAM]


class RAWINPUTDEVICE(ctypes.Structure):
    _fields_ = [("usUsagePage", wt.USHORT), ("usUsage", wt.USHORT),
                ("dwFlags", wt.DWORD), ("hwndTarget", wt.HWND)]


class RAWINPUTHEADER(ctypes.Structure):
    _fields_ = [("dwType", wt.DWORD), ("dwSize", wt.DWORD),
                ("hDevice", wt.HANDLE), ("wParam", wt.WPARAM)]


class _BUTTONS(ctypes.Structure):
    _fields_ = [("usButtonFlags", wt.USHORT), ("usButtonData", wt.USHORT)]


class _BU(ctypes.Union):
    _anonymous_ = ("s",)
    _fields_ = [("ulButtons", wt.ULONG), ("s", _BUTTONS)]


class RAWMOUSE(ctypes.Structure):
    _anonymous_ = ("u",)
    _fields_ = [("usFlags", wt.USHORT), ("u", _BU), ("ulRawButtons", wt.ULONG),
                ("lLastX", wt.LONG), ("lLastY", wt.LONG),
                ("ulExtraInformation", wt.ULONG)]


class RAWKEYBOARD(ctypes.Structure):
    _fields_ = [("MakeCode", wt.USHORT), ("Flags", wt.USHORT),
                ("Reserved", wt.USHORT), ("VKey", wt.USHORT),
                ("Message", wt.UINT), ("ExtraInformation", wt.ULONG)]


class _RIDATA(ctypes.Union):
    _fields_ = [("mouse", RAWMOUSE), ("keyboard", RAWKEYBOARD)]


class RAWINPUT(ctypes.Structure):
    _fields_ = [("header", RAWINPUTHEADER), ("data", _RIDATA)]


_enabled = threading.Event()
_get_file = None
_key_down = {}  # VKey -> bool, to drop typematic auto-repeat
_wndproc_ref = None  # keep the callback alive for the window's lifetime


def set_enabled(on):
    if on:
        _enabled.set()
    else:
        _enabled.clear()
        _key_down.clear()


def _handle_input(lparam):
    f = _get_file()
    if f is None:
        return
    size = wt.UINT(0)
    user32.GetRawInputData(ctypes.c_void_p(lparam), RID_INPUT, None,
                           ctypes.byref(size), ctypes.sizeof(RAWINPUTHEADER))
    buf = ctypes.create_string_buffer(size.value)
    if user32.GetRawInputData(ctypes.c_void_p(lparam), RID_INPUT, buf,
                              ctypes.byref(size),
                              ctypes.sizeof(RAWINPUTHEADER)) != size.value:
        return
    ri = ctypes.cast(buf, ctypes.POINTER(RAWINPUT)).contents
    t = time.perf_counter_ns()
    try:
        if ri.header.dwType == RIM_TYPEMOUSE:
            mo = ri.data.mouse
            bf = mo.usButtonFlags
            if bf & RI_MOUSE_LEFT_BUTTON_DOWN:
                f.write(f"{t},L,0,0\n")
            if bf & RI_MOUSE_LEFT_BUTTON_UP:
                f.write(f"{t},Lu,0,0\n")
            if bf & RI_MOUSE_RIGHT_BUTTON_DOWN:
                f.write(f"{t},R,0,0\n")
            if not (mo.usFlags & MOUSE_MOVE_ABSOLUTE) and (mo.lLastX or mo.lLastY):
                f.write(f"{t},m,{mo.lLastX},{mo.lLastY}\n")
        elif ri.header.dwType == RIM_TYPEKEYBOARD:
            kb = ri.data.keyboard
            vk = kb.VKey
            down = not (kb.Flags & RI_KEY_BREAK)
            k = WASD.get(vk)
            mod = MODS.get(vk)
            if k is None and mod is None:
                return
            if down == _key_down.get(vk, False):
                return  # typematic repeat
            _key_down[vk] = down
            if k:
                f.write(f"{t},{k if down else k.upper() + 'u'},0,0\n")
            else:
                f.write(f"{t},{mod if down else mod + 'u'},0,0\n")
    except Exception:
        pass


def _thread(on_fail):
    global _wndproc_ref

    def wndproc(hwnd, msg, wparam, lparam):
        if msg == WM_INPUT and _enabled.is_set():
            _handle_input(lparam)
        return user32.DefWindowProcW(hwnd, msg, wparam, lparam)

    _wndproc_ref = WNDPROC(wndproc)

    class WNDCLASS(ctypes.Structure):
        _fields_ = [("style", wt.UINT), ("lpfnWndProc", WNDPROC),
                    ("cbClsExtra", ctypes.c_int), ("cbWndExtra", ctypes.c_int),
                    ("hInstance", wt.HINSTANCE), ("hIcon", wt.HANDLE),
                    ("hCursor", wt.HANDLE), ("hbrBackground", wt.HANDLE),
                    ("lpszMenuName", wt.LPCWSTR), ("lpszClassName", wt.LPCWSTR)]

    wc = WNDCLASS()
    wc.lpfnWndProc = _wndproc_ref
    wc.lpszClassName = "MagnetoRawInput"
    wc.hInstance = kernel32.GetModuleHandleW(None)
    if not user32.RegisterClassW(ctypes.byref(wc)):
        on_fail("RegisterClassW failed")
        return
    hwnd = user32.CreateWindowExW(0, wc.lpszClassName, None, 0, 0, 0, 0, 0,
                                  wt.HWND(HWND_MESSAGE), None, wc.hInstance, None)
    if not hwnd:
        on_fail("CreateWindowExW failed")
        return
    devices = (RAWINPUTDEVICE * 2)(
        RAWINPUTDEVICE(0x01, 0x02, RIDEV_INPUTSINK, hwnd),  # mouse
        RAWINPUTDEVICE(0x01, 0x06, RIDEV_INPUTSINK, hwnd),  # keyboard
    )
    if not user32.RegisterRawInputDevices(devices, 2,
                                          ctypes.sizeof(RAWINPUTDEVICE)):
        on_fail("RegisterRawInputDevices failed")
        return
    msg = wt.MSG()
    while user32.GetMessageW(ctypes.byref(msg), None, 0, 0) > 0:
        user32.TranslateMessage(ctypes.byref(msg))
        user32.DispatchMessageW(ctypes.byref(msg))


def start(get_file, on_fail=None):
    """get_file: callable returning the open capture file or None."""
    global _get_file
    _get_file = get_file
    fail = on_fail or (lambda why: print(f"TAP_FAILED:{why}", flush=True))
    threading.Thread(target=_thread, args=(fail,), daemon=True).start()
