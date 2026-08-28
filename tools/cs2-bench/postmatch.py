#!/usr/bin/env python3
"""Post-match autopilot: fetch Riot match details for a just-finished capture,
then run val_enrich. Called by val_session after each report (detached), or
manually: postmatch.py <sessions/...base>

Auto-detects the player from the local Riot client (lockfile -> entitlements).
Uses curl.exe for pd.* requests (python TLS gets Cloudflare-blocked).
Retries while Riot's history indexer lags behind the match end.
"""
import base64, json, os, ssl, subprocess, sys, time, urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
PLAT = base64.b64encode(json.dumps({
    "platformType": "PC", "platformOS": "Windows",
    "platformOSVersion": "10.0.19042.1.256.64bit",
    "platformChipset": "Unknown"}).encode()).decode()
UA = "ShooterGame/18 Windows/10.0.19043.1.256.64bit"


def client_version():
    try:
        return json.loads(urllib.request.urlopen(
            "https://valorant-api.com/v1/version", timeout=10
        ).read())["data"]["riotClientVersion"]
    except Exception:
        return "release-13.04-shipping-18-5304478"


def tokens():
    lock = open(os.path.join(os.environ["LOCALAPPDATA"], "Riot Games",
                             "Riot Client", "Config", "lockfile")).read().strip().split(":")
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(f"https://127.0.0.1:{lock[2]}/entitlements/v1/token")
    req.add_header("Authorization", "Basic " +
                   base64.b64encode(f"riot:{lock[3]}".encode()).decode())
    ent = json.loads(urllib.request.urlopen(req, timeout=5, context=ctx).read())
    return ent["subject"], ent["accessToken"], ent["token"]


def pd_get(path, out, access, jwt, ver):
    r = subprocess.run(
        ["curl.exe", "-s", "-o", out, "-w", "%{http_code}",
         f"https://pd.na.a.pvp.net{path}",
         "-H", f"Authorization: Bearer {access}",
         "-H", f"X-Riot-Entitlements-JWT: {jwt}",
         "-H", f"X-Riot-ClientPlatform: {PLAT}",
         "-H", f"X-Riot-ClientVersion: {ver}",
         "-H", f"User-Agent: {UA}"],
        capture_output=True, text=True, timeout=30)
    return r.stdout.strip() == "200"


def main(base):
    base = base[:-len(".cap.csv")] if base.endswith(".cap.csv") else base
    anchor = float(open(base + ".anchor").read().strip())
    puuid, access, jwt = tokens()
    ver = client_version()
    tmp = base + ".hist.tmp"
    mid = None
    for attempt in range(8):  # indexer can lag minutes behind match end
        if pd_get(f"/match-history/v1/history/{puuid}?startIndex=0&endIndex=5",
                  tmp, access, jwt, ver):
            try:
                hist = json.load(open(tmp)).get("History", [])
            except Exception:
                hist = []
            for h in hist:
                start = (h.get("GameStartTime") or 0) / 1000.0
                if abs(start - anchor) < 1800:  # started within 30 min of capture
                    mid = h["MatchID"]
                    break
        if mid:
            break
        time.sleep(60)
    try:
        os.remove(tmp)
    except Exception:
        pass
    if not mid:
        print("POSTMATCH: match never appeared in history (remake or lag)")
        return
    out = base + ".match.json"
    if not pd_get(f"/match-details/v1/matches/{mid}", out, access, jwt, ver):
        print("POSTMATCH: details fetch failed")
        return
    print(f"POSTMATCH: details saved for {mid}")
    r = subprocess.run([sys.executable, os.path.join(BASE, "val_enrich.py"),
                        base, puuid], capture_output=True, text=True, cwd=BASE)
    print(r.stderr.strip() or "POSTMATCH: enrich done")


if __name__ == "__main__":
    main(sys.argv[1])
