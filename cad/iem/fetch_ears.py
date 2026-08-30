#!/usr/bin/env python3
"""
fetch_ears.py -- download public 3D ear/head meshes for the IEM try-on pipeline.

Two sources, both scriptable, both fetched into cad/iem/ears/ (gitignored):

  hutubs   TU Berlin HUTUBS database, 96 subjects, CC-BY 4.0.
           DSpace 9 REST API on api-depositonce.tu-berlin.de.  The whole
           "3D head meshes.zip" bitstream (149 MB) holds every subject, so
           there is nothing per-subject to iterate: one GET, one unzip.
           Meshes are PLY, head + separately-scanned pinna (0.05 mm Artec
           Space Spider), in the HUTUBS head-centred frame.

  sonicom  SONICOM HRTF dataset, 200+ subjects, hosted on a CrushFTP
           WebInterface at transfer.ic.ac.uk:9090.  The browser UI is
           JavaScript, but the backend exposes an unauthenticated JSON
           listing endpoint:

               POST /WebInterface/function/
                    command=getXMLListing&format=JSONOBJ&path=/some/dir/

           and plain GET on the file path for content.  We take
           SYNTHETIC_HRTF/PXXXX_preprocessed.stl rather than the raw
           3DSCAN/PXXXX.stl: the preprocessed mesh is already aligned to the
           Frankfurt plane, beheaded, smoothed and de-haired, and -- unlike
           PXXXX_plugged.stl -- still has the ear canal open, which is what
           our canal-aperture detector keys on.  ~18 MB per subject.

Usage:
    python fetch_ears.py --source hutubs
    python fetch_ears.py --source sonicom --n 40
    python fetch_ears.py --source all --n 40
"""

from __future__ import annotations

import argparse
import json
import os
import ssl
import sys
import urllib.parse
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
EARS = os.path.join(HERE, "ears")

HUTUBS_ITEM = "dc2a3076-a291-417e-97f0-7697e332c960"
HUTUBS_API = "https://api-depositonce.tu-berlin.de/server/api"
HUTUBS_MESH_BITSTREAM = "3D head meshes.zip"

SONICOM_HOST = "https://transfer.ic.ac.uk:9090"
SONICOM_ROOT = "/2022_SONICOM-HRTF-DATASET"

# transfer.ic.ac.uk:9090 serves a cert chain that the stdlib store rejects on
# macOS; the payload is public research data with no credentials in flight, so
# an unverified context is acceptable here and nowhere else in this repo.
_NOVERIFY = ssl.create_default_context()
_NOVERIFY.check_hostname = False
_NOVERIFY.verify_mode = ssl.CERT_NONE


def _get(url, ctx=None, timeout=900):
    req = urllib.request.Request(url, headers={"User-Agent": "magneto-fit/1.0"})
    return urllib.request.urlopen(req, timeout=timeout, context=ctx)


def _download(url, dest, ctx=None):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return False
    tmp = dest + ".part"
    with _get(url, ctx) as r, open(tmp, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    os.replace(tmp, dest)
    return True


# --------------------------------------------------------------------------- #
# HUTUBS
# --------------------------------------------------------------------------- #

def fetch_hutubs():
    out = os.path.join(EARS, "hutubs")
    os.makedirs(out, exist_ok=True)

    with _get(f"{HUTUBS_API}/core/items/{HUTUBS_ITEM}/bundles") as r:
        bundles = json.load(r)["_embedded"]["bundles"]
    orig = next(b for b in bundles if b["name"] == "ORIGINAL")

    with _get(orig["_links"]["bitstreams"]["href"] + "?size=100") as r:
        bits = json.load(r)["_embedded"]["bitstreams"]
    mesh = next(b for b in bits if b["name"] == HUTUBS_MESH_BITSTREAM)

    zpath = os.path.join(out, "3D_head_meshes.zip")
    print(f"hutubs: {mesh['name']} ({mesh['sizeBytes']/1e6:.0f} MB)")
    if _download(mesh["_links"]["content"]["href"], zpath):
        print("hutubs: downloaded")
    else:
        print("hutubs: already present")

    with zipfile.ZipFile(zpath) as z:
        z.extractall(out)
    n = sum(len(f) for _, _, f in os.walk(out))
    print(f"hutubs: {n} files extracted under {out}")


# --------------------------------------------------------------------------- #
# SONICOM
# --------------------------------------------------------------------------- #

def sonicom_list(path):
    """CrushFTP JSON directory listing.  path must end in '/'."""
    body = urllib.parse.urlencode({
        "command": "getXMLListing",
        "format": "JSONOBJ",
        "path": path,
        "random": "0.1",
    }).encode()
    req = urllib.request.Request(
        f"{SONICOM_HOST}/WebInterface/function/", data=body,
        headers={"User-Agent": "magneto-fit/1.0"})
    with urllib.request.urlopen(req, timeout=120, context=_NOVERIFY) as r:
        return json.load(r)["listing"]


def fetch_sonicom(n_subjects):
    out = os.path.join(EARS, "sonicom")
    os.makedirs(out, exist_ok=True)

    subs = [e["name"] for e in sonicom_list(SONICOM_ROOT + "/")
            if e["type"] == "DIR" and e["name"].startswith("P0")]
    subs.sort()
    subs = subs[:n_subjects]
    print(f"sonicom: {len(subs)} subjects targeted")

    got = 0
    for s in subs:
        dest = os.path.join(out, f"{s}_preprocessed.stl")
        url = f"{SONICOM_HOST}{SONICOM_ROOT}/{s}/SYNTHETIC_HRTF/{s}_preprocessed.stl"
        try:
            fresh = _download(url, dest, _NOVERIFY)
            got += 1
            print(f"  {s} {'ok' if fresh else 'cached'} "
                  f"{os.path.getsize(dest)/1e6:.1f} MB", flush=True)
        except Exception as e:                      # noqa: BLE001
            print(f"  {s} FAILED {e}", flush=True)
    print(f"sonicom: {got}/{len(subs)} meshes in {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=("hutubs", "sonicom", "all"), default="all")
    ap.add_argument("--n", type=int, default=40,
                    help="number of SONICOM subjects (HUTUBS is all-or-nothing)")
    a = ap.parse_args()
    os.makedirs(EARS, exist_ok=True)
    if a.source in ("hutubs", "all"):
        fetch_hutubs()
    if a.source in ("sonicom", "all"):
        fetch_sonicom(a.n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
